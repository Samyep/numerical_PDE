"""Algebraic regression tests for the 2-D Euler HCFL core."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import euler_2d_common as C  # noqa: E402


def conserved(rho: torch.Tensor, u: torch.Tensor, v: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    energy = p / (C.GAMMA - 1.0) + 0.5 * rho * (u.square() + v.square())
    return torch.stack((rho, rho * u, rho * v, energy), dim=-1)


def random_state(shape: tuple[int, ...], seed: int = 1) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    rho = 0.5 + 1.5 * torch.rand(shape, generator=generator)
    u = -0.6 + 1.2 * torch.rand(shape, generator=generator)
    v = -0.6 + 1.2 * torch.rand(shape, generator=generator)
    p = 0.4 + 1.6 * torch.rand(shape, generator=generator)
    return conserved(rho, u, v, p)


def test_orientation_and_physical_flux() -> None:
    state = random_state((2, 5, 7), seed=11)
    assert torch.equal(C.orient_state(C.orient_state(state, "y"), "y"), state)
    rho = state[..., 0]
    mx = state[..., 1]
    my = state[..., 2]
    energy = state[..., 3]
    velocity_y = my / rho
    pressure = C.pressure_raw(state)
    expected = torch.stack(
        (my, mx * velocity_y, my * velocity_y + pressure, velocity_y * (energy + pressure)),
        dim=-1,
    )
    actual = C.deorient_flux(
        C.physical_flux_oriented(C.orient_state(state, "y")), "y"
    )
    assert torch.allclose(actual, expected, rtol=2.0e-6, atol=2.0e-6)


def test_roe_reconstruction() -> None:
    left = random_state((2, 4, 6), seed=21)
    right = random_state((2, 4, 6), seed=22)
    matrix, strengths, speeds = C.roe_waves_oriented(left, right)
    reconstruction = torch.einsum("...ij,...j->...i", matrix, strengths)
    assert torch.allclose(reconstruction, right - left, rtol=2.0e-5, atol=2.0e-5)
    assert bool((speeds >= 0.0).all())


def test_uniform_consistency() -> None:
    rho = torch.full((2, 6, 8), 1.2)
    u = torch.full_like(rho, 0.3)
    v = torch.full_like(rho, -0.2)
    p = torch.full_like(rho, 0.9)
    state = conserved(rho, u, v, p)
    model = C.HCFL2DEuler(
        np.asarray([1.0, 0.0, 0.0, 1.0], dtype=np.float32),
        np.ones(4, dtype=np.float32),
        stencil_cells=4,
    )
    raw_x, raw_y = model.raw_fluxes(state)
    projected_x, projected_y = model.project_raw_fluxes(state, (raw_x, raw_y))
    expected_x = C.physical_flux_oriented(C.orient_state(state, "x"))
    expected_y = C.deorient_flux(
        C.physical_flux_oriented(C.orient_state(state, "y")), "y"
    )
    assert torch.allclose(raw_x[..., 1:-1, :], expected_x[..., :-1, :], atol=2.0e-6)
    assert torch.allclose(raw_y[..., 1:-1, :, :], expected_y[..., :-1, :, :], atol=2.0e-6)
    assert torch.allclose(projected_x, raw_x, atol=2.0e-6)
    assert torch.allclose(projected_y, raw_y, atol=2.0e-6)


def test_hard_projection_and_conservation() -> None:
    state = random_state((2, 9, 11), seed=31)
    model = C.HCFL2DEuler(
        np.asarray([1.0, 0.0, 0.0, 1.0], dtype=np.float32),
        np.ones(4, dtype=np.float32),
        stencil_cells=6,
    )
    with torch.no_grad():
        model.flux_net.net[-1].bias.copy_(torch.tensor([0.7, -0.4, 0.5, -0.8]))
    flux_x, flux_y = model.projected_fluxes(state)
    residual_x = C.interface_entropy_residual_oriented(
        flux_x.double(), C.orient_state(state.double(), "x")
    )
    residual_y = C.interface_entropy_residual_oriented(
        C.orient_state(flux_y.double(), "y"), C.orient_state(state.double(), "y")
    )
    assert float(residual_x.detach().max()) <= 2.0e-6
    assert float(residual_y.detach().max()) <= 2.0e-6

    dt = 1.0e-5
    updated = C.flux_divergence(state, flux_x, flux_y, dt)
    closure = (
        updated.double().sum(dim=(-3, -2))
        - state.double().sum(dim=(-3, -2))
        + C.boundary_conservation_increment(state.double(), dt)
    )
    assert float(closure.detach().abs().max()) <= 2.0e-5


def test_low_order_short_step_is_admissible_and_entropy_stable() -> None:
    state = random_state((2, 12, 12), seed=41)
    dt = 0.05 / C.maximum_2d_rate(state)
    updated = C.flux_divergence(state, *C.hll_fluxes(state), dt)
    assert bool(C.admissible(updated).all())
    balance = (
        C.total_entropy(updated)
        - C.total_entropy(state)
        + C.boundary_entropy_increment(state, dt)
    )
    assert float(balance.max()) <= C.ENTROPY_TOL


if __name__ == "__main__":
    tests = [
        test_orientation_and_physical_flux,
        test_roe_reconstruction,
        test_uniform_consistency,
        test_hard_projection_and_conservation,
        test_low_order_short_step_is_admissible_and_entropy_stable,
    ]
    for test in tests:
        test()
        print({"passed": test.__name__}, flush=True)
