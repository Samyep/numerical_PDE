"""Scientific invariants for six-cell HCFL plus fixed transverse transport."""

from __future__ import annotations

import numpy as np
import torch

import euler_2d_common as C
import fixed_transverse as T


MEAN = np.array([1.0, 0.0, 0.0, 1.0], dtype=np.float32)
STD = np.ones(4, dtype=np.float32)


def conserved(rho: torch.Tensor, u: torch.Tensor, v: torch.Tensor, p: torch.Tensor):
    energy = p / (C.GAMMA - 1.0) + 0.5 * rho * (u.square() + v.square())
    return torch.stack((rho, rho * u, rho * v, energy), dim=-1)


def smooth_state(batch: int = 2, cells: int = 12) -> torch.Tensor:
    y, x = torch.meshgrid(
        torch.linspace(0.0, 1.0, cells),
        torch.linspace(0.0, 1.0, cells),
        indexing="ij",
    )
    rho = 1.0 + 0.08 * torch.sin(2.0 * torch.pi * x) * torch.cos(torch.pi * y)
    u = 0.18 + 0.04 * torch.cos(2.0 * torch.pi * y)
    v = -0.11 + 0.03 * torch.sin(2.0 * torch.pi * x)
    p = 1.0 + 0.06 * torch.cos(torch.pi * (x + y))
    return conserved(rho, u, v, p).unsqueeze(0).repeat(batch, 1, 1, 1)


def test_transverse_split_reconstructs_roe_action() -> None:
    state = smooth_state(batch=1, cells=8)
    oriented = C.orient_state(state, "x")
    left = oriented[..., :-1, :]
    right = oriented[..., 1:, :]
    fluctuation = 0.13 * right - 0.07 * left
    down, up = T.transverse_roe_split_oriented(left, right, fluctuation)

    # B^-q + B^+q = B_roe q.  A centred finite-difference Jacobian action is
    # an independent check of the analytic characteristic decomposition.
    epsilon = 2.0e-4
    average = 0.5 * (left + right)
    numerical = (
        C.physical_flux_oriented(
            C.orient_state(average + epsilon * fluctuation, "y")
        )
        - C.physical_flux_oriented(
            C.orient_state(average - epsilon * fluctuation, "y")
        )
    ) / (2.0 * epsilon)
    numerical = C.orient_state(numerical, "y")
    assert torch.isfinite(down).all()
    assert torch.isfinite(up).all()
    # Roe linearization is evaluated at a Roe average, not the arithmetic
    # state used by the finite-difference check, so use a modest tolerance.
    assert torch.allclose(down + up, numerical, rtol=7.0e-2, atol=2.0e-2)


def test_zero_logits_are_standard_roe_and_multipliers_are_nonnegative() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=1, cells=10)
    oriented = C.orient_state(state, "x")
    multipliers = model.flux_net.multipliers(oriented)
    learned = model.flux_net.forward_oriented(oriented)
    standard = model.flux_net.standard_roe_faces_oriented(oriented)
    assert torch.equal(multipliers, torch.ones_like(multipliers))
    assert torch.allclose(learned, standard, rtol=0.0, atol=0.0)
    assert bool((multipliers >= 0.0).all())
    assert bool((multipliers <= 2.0).all())


def test_exact_one_dimensional_reduction() -> None:
    model = T.make_model(MEAN, STD)
    cells = 14
    x = torch.linspace(0.0, 1.0, cells)
    rho = 1.0 + 0.1 * torch.sin(2.0 * torch.pi * x)
    u = 0.2 + 0.03 * torch.cos(2.0 * torch.pi * x)
    v = torch.zeros_like(x)
    p = 1.0 + 0.05 * torch.sin(torch.pi * x)
    row = conserved(rho, u, v, p)
    state = row[None, None].repeat(1, 9, 1, 1)
    dt = 7.0e-4

    normal_x, normal_y = model.normal_raw_fluxes(state)
    full_x, full_y = model.raw_fluxes(state, dt=dt)
    base_update = C.flux_divergence(state, normal_x, normal_y, dt)
    full_update = C.flux_divergence(state, full_x, full_y, dt)
    assert torch.allclose(full_x, normal_x, rtol=0.0, atol=2.0e-7)
    assert torch.allclose(full_update, base_update, rtol=0.0, atol=2.0e-7)


def test_vectorized_corner_assembly_matches_direct_scatter() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=1, cells=7)
    oriented = C.orient_state(state, "x")
    normal_flux = model.flux_net.forward_oriented(oriented)
    dt = 6.0e-4
    vectorized = T.fixed_transverse_correction_oriented(
        oriented, normal_flux, dt
    )

    padded = torch.cat(
        (oriented[..., :1, :], oriented, oriented[..., -1:, :]), dim=-2
    )
    left = padded[..., :-1, :]
    right = padded[..., 1:, :]
    minus = normal_flux - C.physical_flux_oriented(left)
    plus = C.physical_flux_oriented(right) - normal_flux
    down_minus, up_minus = T.transverse_roe_split_oriented(left, right, minus)
    down_plus, up_plus = T.transverse_roe_split_oriented(left, right, plus)
    down_cell = down_plus[..., :-1, :] + down_minus[..., 1:, :]
    up_cell = up_plus[..., :-1, :] + up_minus[..., 1:, :]
    factor = -0.5 * dt / (C.DOMAIN_LENGTH / oriented.shape[-2])
    direct = torch.zeros_like(vectorized)
    transverse_cells = oriented.shape[-3]
    for row in range(transverse_cells):
        direct[..., row, :, :] += factor * down_cell[..., row, :, :]
        direct[..., row + 1, :, :] += factor * up_cell[..., row, :, :]
    # Contributions from the constant-extrapolation ghost rows.
    direct[..., 0, :, :] += factor * up_cell[..., 0, :, :]
    direct[..., -1, :, :] += factor * down_cell[..., -1, :, :]
    assert torch.allclose(vectorized, direct, rtol=2.0e-6, atol=2.0e-7)


def test_conservation_and_final_hard_entropy_projection() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=2, cells=11)
    dt = 4.0e-4
    raw = model.raw_fluxes(state, dt=dt)
    projected = model.project_raw_fluxes(state, raw)
    candidate = C.flux_divergence(state, *projected, dt)

    boundary = dt / (C.DOMAIN_LENGTH / state.shape[-2]) * (
        (projected[0][..., -1, :] - projected[0][..., 0, :]).sum(dim=-2)
        + (projected[1][..., -1, :, :] - projected[1][..., 0, :, :]).sum(dim=-2)
    )
    closure = candidate.sum(dim=(-3, -2)) - state.sum(dim=(-3, -2)) + boundary
    assert float(closure.detach().abs().max()) < 2.0e-5

    residual_x = C.interface_entropy_residual_oriented(
        projected[0].double(), C.orient_state(state.double(), "x")
    )
    residual_y = C.interface_entropy_residual_oriented(
        C.orient_state(projected[1].double(), "y"),
        C.orient_state(state.double(), "y"),
    )
    assert float(residual_x.detach().max()) <= 2.0e-6
    assert float(residual_y.detach().max()) <= 2.0e-6


def test_transverse_transport_has_no_trainable_parameters_and_backpropagates() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=1, cells=9)
    dt = 5.0e-4
    output = C.flux_divergence(
        state,
        *model.project_raw_fluxes_training(state, model.raw_fluxes(state, dt=dt)),
        dt,
    )
    output.square().mean().backward()
    assert model.transverse_trainable_parameters == 0
    assert all(parameter.grad is not None for parameter in model.parameters())
    assert all(torch.isfinite(parameter.grad).all() for parameter in model.parameters())


def test_inference_ablation_disables_only_fixed_transverse_term() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=1, cells=9)
    normal = model.normal_raw_fluxes(state)
    model.transverse_enabled = False
    ablated = model.raw_fluxes(state, dt=5.0e-4)
    assert all(
        torch.equal(candidate, expected)
        for candidate, expected in zip(ablated, normal)
    )
    assert sum(parameter.numel() for parameter in model.parameters()) == 7348


def test_shared_model_and_fixed_transport_are_xy_rotation_equivariant() -> None:
    model = T.make_model(MEAN, STD)
    state = smooth_state(batch=1, cells=9)
    dt = 5.0e-4
    flux_x, flux_y = model.raw_fluxes(state, dt=dt)
    rotated_state = C.orient_state(state, "y")
    rotated_flux_x, _ = model.raw_fluxes(rotated_state, dt=dt)
    expected = C.orient_state(flux_y, "y")
    assert torch.allclose(rotated_flux_x, expected, rtol=2.0e-6, atol=2.0e-7)
