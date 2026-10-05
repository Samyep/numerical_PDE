"""Regression tests for the 18-cell transverse HCFL variants."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "euler_2d_clawpack_quadrants"
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(HERE))
import euler_2d_common as C  # noqa: E402
import transverse_models as M  # noqa: E402


MEAN = np.asarray([1.0, 0.0, 0.0, 1.0], dtype=np.float32)
STD = np.ones(4, dtype=np.float32)


def conserved(rho: torch.Tensor, u: torch.Tensor, v: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    energy = p / (C.GAMMA - 1.0) + 0.5 * rho * (u.square() + v.square())
    return torch.stack((rho, rho * u, rho * v, energy), dim=-1)


def random_state(shape: tuple[int, ...], seed: int) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    rho = 0.5 + 1.5 * torch.rand(shape, generator=generator)
    u = -0.6 + 1.2 * torch.rand(shape, generator=generator)
    v = -0.6 + 1.2 * torch.rand(shape, generator=generator)
    p = 0.4 + 1.6 * torch.rand(shape, generator=generator)
    return conserved(rho, u, v, p)


def test_patch_shape_and_boundary_replication() -> None:
    state = random_state((2, 5, 8), seed=1)
    patch = M.oriented_primitive_patch(state)
    assert patch.shape == (2, 5, 7, 3, 6, 4)
    assert torch.equal(patch[:, 0, :, 0], patch[:, 0, :, 1])
    assert torch.equal(patch[:, -1, :, 1], patch[:, -1, :, 2])


def test_parameter_budgets_are_matched() -> None:
    normal = M.make_model(MEAN, STD, "normal6wide")
    flat = M.make_model(MEAN, STD, "flat18")
    gated = M.make_model(MEAN, STD, "gated18")
    normal_count = C.parameter_count(normal)
    flat_count = C.parameter_count(flat)
    gated_count = C.parameter_count(gated)
    assert normal_count == 10804
    assert flat_count == 10804
    assert gated_count == 10872
    assert abs(flat_count - gated_count) / flat_count < 0.01


def test_gated_model_has_exact_one_dimensional_reduction() -> None:
    line = random_state((2, 1, 11), seed=2)
    state = line.expand(-1, 7, -1, -1).clone()
    model = M.make_model(MEAN, STD, "gated18")
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.uniform_(-0.2, 0.2)
    coefficients, centre_only, gate = model.flux_net.coefficient_details(state)
    assert torch.equal(gate, torch.zeros_like(gate))
    assert torch.equal(coefficients, centre_only)


def test_flat_model_receives_real_transverse_information() -> None:
    state = random_state((1, 7, 10), seed=3)
    model = M.make_model(MEAN, STD, "flat18")
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.uniform_(-0.15, 0.15)
    coefficients, centre_only, gate = model.flux_net.coefficient_details(state)
    assert float(gate.max()) > 0.0
    assert not torch.allclose(coefficients, centre_only)


def test_hard_projection_and_conservation_for_both_variants() -> None:
    state = random_state((2, 9, 11), seed=4)
    for variant in M.MODEL_TYPES:
        model = M.make_model(MEAN, STD, variant)
        with torch.no_grad():
            for parameter in model.parameters():
                parameter.uniform_(-0.08, 0.08)
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


if __name__ == "__main__":
    tests = [
        test_patch_shape_and_boundary_replication,
        test_parameter_budgets_are_matched,
        test_gated_model_has_exact_one_dimensional_reduction,
        test_flat_model_receives_real_transverse_information,
        test_hard_projection_and_conservation_for_both_variants,
    ]
    for test in tests:
        test()
        print({"passed": test.__name__}, flush=True)
