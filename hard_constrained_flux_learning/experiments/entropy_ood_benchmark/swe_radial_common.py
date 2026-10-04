"""Shared 2-D radial-dam utilities for the entropy/OOD benchmark.

The task follows PDEBench's radial dam-break geometry and boundary condition,
but the trajectories are generated locally with an independently implemented
MC-HLL/SSP-RK2 reference solver.  State ordering in this file is always
``(h, hu, hv)``.  Conversion to clawNO's ``(hu, hv, h)`` convention is kept at
the operator-model boundary.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import torch
from torch import nn


G = 1.0
DOMAIN_LOWER = -2.5
DOMAIN_UPPER = 2.5
DOMAIN_LENGTH = DOMAIN_UPPER - DOMAIN_LOWER
N_COARSE = 32
N_REFERENCE = 128
RESTRICTION = N_REFERENCE // N_COARSE
SNAPSHOT_DT = 0.04
ORDINARY_STATES = 25
LONG_STATES = 49
H_FLOOR = 1.0e-5
ENTROPY_TOL = 5.0e-7


SPLIT_SPECS: dict[str, dict[str, Any]] = {
    "train": {"count": 100, "seed": 61001, "kind": "id", "states": ORDINARY_STATES},
    "validation": {"count": 20, "seed": 62001, "kind": "id", "states": ORDINARY_STATES},
    "test_id": {"count": 20, "seed": 63001, "kind": "id", "states": ORDINARY_STATES},
    "test_radius_ood": {
        "count": 20,
        "seed": 64001,
        "kind": "radius_ood",
        "states": ORDINARY_STATES,
    },
    "test_height_ood": {
        "count": 20,
        "seed": 65001,
        "kind": "height_ood",
        "states": ORDINARY_STATES,
    },
    "test_long": {"count": 10, "seed": 66001, "kind": "id", "states": LONG_STATES},
}


def primitive(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0].clamp_min(H_FLOOR)
    return torch.stack((h, state[..., 1] / h, state[..., 2] / h), dim=-1)


def physical_flux_oriented(state: torch.Tensor) -> torch.Tensor:
    """Physical flux when the last two components are normal/tangent momentum."""
    h = state[..., 0].clamp_min(H_FLOOR)
    mn = state[..., 1]
    mt = state[..., 2]
    un = mn / h
    ut = mt / h
    return torch.stack(
        (mn, mn * un + 0.5 * G * h * h, mn * ut), dim=-1
    )


def hll_pair(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    f_left = physical_flux_oriented(left)
    f_right = physical_flux_oriented(right)
    h_left = left[..., 0].clamp_min(H_FLOOR)
    h_right = right[..., 0].clamp_min(H_FLOOR)
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    c_left = torch.sqrt(G * h_left)
    c_right = torch.sqrt(G * h_right)
    speed_left = torch.minimum(u_left - c_left, u_right - c_right)
    speed_right = torch.maximum(u_left + c_left, u_right + c_right)
    denominator = (speed_right - speed_left).clamp_min(1.0e-12)
    middle = (
        speed_right[..., None] * f_left
        - speed_left[..., None] * f_right
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / denominator[..., None]
    return torch.where(
        (speed_left >= 0.0)[..., None],
        f_left,
        torch.where((speed_right <= 0.0)[..., None], f_right, middle),
    )


def orient_state(state: torch.Tensor, direction: str) -> torch.Tensor:
    if direction == "x":
        return state
    if direction == "y":
        return state.transpose(-3, -2)[..., [0, 2, 1]]
    raise ValueError(direction)


def deorient_flux(flux: torch.Tensor, direction: str) -> torch.Tensor:
    if direction == "x":
        return flux
    if direction == "y":
        return flux[..., [0, 2, 1]].transpose(-3, -2)
    raise ValueError(direction)


def hll_faces_oriented(state: torch.Tensor) -> torch.Tensor:
    padded = torch.cat((state[..., :1, :], state, state[..., -1:, :]), dim=-2)
    return hll_pair(padded[..., :-1, :], padded[..., 1:, :])


def hll_fluxes(state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    fx = hll_faces_oriented(orient_state(state, "x"))
    gy = deorient_flux(hll_faces_oriented(orient_state(state, "y")), "y")
    return fx, gy


def flux_divergence(
    state: torch.Tensor,
    flux_x: torch.Tensor,
    flux_y: torch.Tensor,
    dt: float,
) -> torch.Tensor:
    dx = DOMAIN_LENGTH / state.shape[-2]
    return state - (dt / dx) * (
        flux_x[..., 1:, :] - flux_x[..., :-1, :]
        + flux_y[..., 1:, :, :] - flux_y[..., :-1, :, :]
    )


def _minmod3(a: torch.Tensor, b: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
    same = (torch.sign(a) == torch.sign(b)) & (torch.sign(b) == torch.sign(c))
    magnitude = torch.minimum(torch.minimum(a.abs(), b.abs()), c.abs())
    return torch.where(same, torch.sign(a) * magnitude, torch.zeros_like(a))


def mc_hll_faces_oriented(state: torch.Tensor) -> torch.Tensor:
    previous = torch.cat((state[..., :1, :], state[..., :-1, :]), dim=-2)
    following = torch.cat((state[..., 1:, :], state[..., -1:, :]), dim=-2)
    backward = state - previous
    forward = following - state
    slope = _minmod3(2.0 * backward, 0.5 * (backward + forward), 2.0 * forward)

    left = state[..., :-1, :] + 0.5 * slope[..., :-1, :]
    right = state[..., 1:, :] - 0.5 * slope[..., 1:, :]
    left_cell = state[..., :-1, :]
    right_cell = state[..., 1:, :]
    left = torch.where((left[..., 0] > H_FLOOR)[..., None], left, left_cell)
    right = torch.where((right[..., 0] > H_FLOOR)[..., None], right, right_cell)
    interior = hll_pair(left, right)
    return torch.cat(
        (
            physical_flux_oriented(state[..., :1, :]),
            interior,
            physical_flux_oriented(state[..., -1:, :]),
        ),
        dim=-2,
    )


def reference_rhs(state: torch.Tensor) -> torch.Tensor:
    fx = mc_hll_faces_oriented(orient_state(state, "x"))
    gy = deorient_flux(mc_hll_faces_oriented(orient_state(state, "y")), "y")
    dx = DOMAIN_LENGTH / state.shape[-2]
    return -(
        fx[..., 1:, :] - fx[..., :-1, :]
        + gy[..., 1:, :, :] - gy[..., :-1, :, :]
    ) / dx


def reference_ssprk2(state: torch.Tensor, dt: float) -> torch.Tensor:
    stage = state + dt * reference_rhs(state)
    return 0.5 * state + 0.5 * (stage + dt * reference_rhs(stage))


def maximum_2d_rate(state: torch.Tensor) -> float:
    values = primitive(state)
    c = torch.sqrt(G * values[..., 0])
    dx = DOMAIN_LENGTH / state.shape[-2]
    rate = ((values[..., 1].abs() + c) / dx + (values[..., 2].abs() + c) / dx).max()
    return float(rate)


def restrict_reference(state: torch.Tensor) -> torch.Tensor:
    batch = state.shape[0]
    return state.reshape(
        batch,
        N_COARSE,
        RESTRICTION,
        N_COARSE,
        RESTRICTION,
        3,
    ).mean(dim=(2, 4))


def sample_parameters(count: int, seed: int, kind: str) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if kind == "id":
        radii = rng.uniform(0.3, 0.7, size=count)
        heights = np.full(count, 2.0)
    elif kind == "radius_ood":
        low_count = count // 2
        radii = np.concatenate(
            (rng.uniform(0.15, 0.25, size=low_count), rng.uniform(0.75, 0.85, size=count - low_count))
        )
        rng.shuffle(radii)
        heights = np.full(count, 2.0)
    elif kind == "height_ood":
        radii = rng.uniform(0.3, 0.7, size=count)
        heights = rng.uniform(2.5, 3.0, size=count)
    else:
        raise ValueError(kind)
    return radii.astype(np.float32), heights.astype(np.float32)


def radial_initial_state(
    radii: np.ndarray,
    heights: np.ndarray,
    cells: int,
    device: torch.device,
) -> torch.Tensor:
    coordinate = torch.linspace(
        DOMAIN_LOWER + 0.5 * DOMAIN_LENGTH / cells,
        DOMAIN_UPPER - 0.5 * DOMAIN_LENGTH / cells,
        cells,
        device=device,
        dtype=torch.float32,
    )
    yy, xx = torch.meshgrid(coordinate, coordinate, indexing="ij")
    radius = torch.sqrt(xx * xx + yy * yy)[None]
    thresholds = torch.as_tensor(radii, device=device)[:, None, None]
    inner = torch.as_tensor(heights, device=device)[:, None, None]
    h = torch.where(radius <= thresholds, inner, torch.ones_like(inner + radius))
    zeros = torch.zeros_like(h)
    return torch.stack((h, zeros, zeros), dim=-1)


@torch.no_grad()
def generate_reference_split(
    count: int,
    seed: int,
    kind: str,
    states: int,
    device: torch.device,
    batch_size: int = 10,
) -> dict[str, torch.Tensor]:
    radii, heights = sample_parameters(count, seed, kind)
    trajectories: list[torch.Tensor] = []
    for start in range(0, count, batch_size):
        stop = min(start + batch_size, count)
        current = radial_initial_state(
            radii[start:stop], heights[start:stop], N_REFERENCE, device
        )
        snapshots = [restrict_reference(current).cpu()]
        for _ in range(1, states):
            remaining = SNAPSHOT_DT
            while remaining > 1.0e-10:
                dt = min(remaining, 0.42 / max(maximum_2d_rate(current), 1.0e-12))
                accepted = False
                for _ in range(18):
                    candidate = reference_ssprk2(current, dt)
                    if bool(torch.isfinite(candidate).all()) and float(candidate[..., 0].min()) > H_FLOOR:
                        accepted = True
                        break
                    dt *= 0.5
                if not accepted:
                    raise RuntimeError("Reference solver could not retain positive depth")
                current = candidate
                remaining -= dt
            snapshots.append(restrict_reference(current).cpu())
        trajectories.append(torch.stack(snapshots, dim=1))
    return {
        "trajectory": torch.cat(trajectories, dim=0).contiguous(),
        "radius": torch.from_numpy(radii),
        "inner_height": torch.from_numpy(heights),
    }


def load_or_generate_dataset(cache_path: Path, device: torch.device) -> dict[str, Any]:
    if cache_path.exists():
        return torch.load(cache_path, map_location="cpu", weights_only=False)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {
        "metadata": {
            "generator": "locally implemented MC-HLL SSP-RK2",
            "reference_cells": N_REFERENCE,
            "coarse_cells": N_COARSE,
            "gravity": G,
            "domain": [DOMAIN_LOWER, DOMAIN_UPPER],
            "snapshot_dt": SNAPSHOT_DT,
            "boundary": "constant extrapolation",
        }
    }
    for name, spec in SPLIT_SPECS.items():
        print({"stage": "generate_swe_reference", "split": name, **spec}, flush=True)
        result[name] = generate_reference_split(device=device, **spec)
    torch.save(result, cache_path)
    return result


def entropy(state: torch.Tensor) -> torch.Tensor:
    values = primitive(state)
    h, u, v = values.unbind(dim=-1)
    return 0.5 * h * (u * u + v * v) + 0.5 * G * h * h


def entropy_variables_oriented(state: torch.Tensor) -> torch.Tensor:
    values = primitive(state)
    h, un, ut = values.unbind(dim=-1)
    return torch.stack((G * h - 0.5 * (un * un + ut * ut), un, ut), dim=-1)


def entropy_potential_oriented(state: torch.Tensor) -> torch.Tensor:
    values = primitive(state)
    h, un = values[..., 0], values[..., 1]
    return 0.5 * G * h * h * un


def physical_entropy_flux_oriented(state: torch.Tensor) -> torch.Tensor:
    values = primitive(state)
    h, un, ut = values.unbind(dim=-1)
    return un * (0.5 * h * (un * un + ut * ut) + G * h * h)


def roe_waves_oriented(
    left: torch.Tensor, right: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    h_left = left[..., 0].clamp_min(H_FLOOR)
    h_right = right[..., 0].clamp_min(H_FLOOR)
    root_left = torch.sqrt(h_left)
    root_right = torch.sqrt(h_right)
    denominator = (root_left + root_right).clamp_min(1.0e-8)
    un = (left[..., 1] / root_left + right[..., 1] / root_right) / denominator
    ut = (left[..., 2] / root_left + right[..., 2] / root_right) / denominator
    c = torch.sqrt(0.5 * G * (h_left + h_right)).clamp_min(1.0e-6)
    matrix = torch.stack(
        (
            torch.stack((torch.ones_like(un), un - c, ut), dim=-1),
            torch.stack((torch.zeros_like(un), torch.zeros_like(un), torch.ones_like(un)), dim=-1),
            torch.stack((torch.ones_like(un), un + c, ut), dim=-1),
        ),
        dim=-1,
    )
    # The stack above directly builds columns [r_minus, r_shear, r_plus].
    jump = right - left
    alpha = torch.linalg.solve(matrix, jump.unsqueeze(-1)).squeeze(-1)
    eigenvalues = torch.stack((un - c, un, un + c), dim=-1)
    absolute = eigenvalues.abs()
    entropy_fix = (0.10 * c).clamp_min(1.0e-6)[..., None]
    fixed = torch.where(
        absolute < entropy_fix,
        0.5 * (absolute * absolute / entropy_fix + entropy_fix),
        absolute,
    )
    return matrix, alpha, fixed


def interface_entropy_residual_oriented(
    flux: torch.Tensor, state: torch.Tensor
) -> torch.Tensor:
    left = state[..., :-1, :]
    right = state[..., 1:, :]
    normal = entropy_variables_oriented(right) - entropy_variables_oriented(left)
    bound = entropy_potential_oriented(right) - entropy_potential_oriented(left)
    return (normal * flux[..., 1:-1, :]).sum(dim=-1) - bound


def strict_entropy_projection_oriented(
    proposal: torch.Tensor, state: torch.Tensor
) -> torch.Tensor:
    """Project only interior learned faces; physical boundary fluxes are fixed."""
    interior = proposal[..., 1:-1, :]
    state64 = state.double()
    left = state64[..., :-1, :]
    right = state64[..., 1:, :]
    normal = entropy_variables_oriented(right) - entropy_variables_oriented(left)
    bound = entropy_potential_oriented(right) - entropy_potential_oriented(left)
    norm_squared = (normal * normal).sum(dim=-1)
    epsilon = torch.finfo(proposal.dtype).eps
    output = interior
    for _ in range(4):
        flux64 = output.double()
        residual = (normal * flux64).sum(dim=-1) - bound
        scale = bound.abs() + (normal * flux64).abs().sum(dim=-1)
        target = -8.0 * epsilon * scale
        active = (residual > target) & (norm_squared > 1.0e-14)
        if not bool(active.any()):
            break
        coefficient = torch.zeros_like(residual)
        coefficient[active] = (residual[active] - target[active]) / norm_squared[active]
        output = (flux64 - coefficient[..., None] * normal).to(proposal.dtype)
    return torch.cat((proposal[..., :1, :], output, proposal[..., -1:, :]), dim=-2)


def training_entropy_projection_oriented(
    proposal: torch.Tensor, state: torch.Tensor
) -> torch.Tensor:
    """Differentiable same-dtype projection used only inside optimizer steps.

    Deployment and every reported safety metric use the stricter float64
    iterative projection above.  Avoiding consumer-GPU float64 here changes
    optimization cost, not the accepted forward solver.
    """
    interior = proposal[..., 1:-1, :]
    left = state[..., :-1, :]
    right = state[..., 1:, :]
    normal = entropy_variables_oriented(right) - entropy_variables_oriented(left)
    bound = entropy_potential_oriented(right) - entropy_potential_oriented(left)
    residual = (normal * interior).sum(dim=-1) - bound
    norm_squared = (normal * normal).sum(dim=-1)
    scale = bound.abs() + (normal * interior).abs().sum(dim=-1)
    target = -16.0 * torch.finfo(proposal.dtype).eps * scale
    coefficient = torch.where(
        (residual > target) & (norm_squared > 1.0e-12),
        (residual - target) / norm_squared.clamp_min(1.0e-12),
        torch.zeros_like(residual),
    )
    projected = interior - coefficient[..., None] * normal
    return torch.cat((proposal[..., :1, :], projected, proposal[..., -1:, :]), dim=-2)


def stencil_shifts(stencil_cells: int) -> tuple[int, ...]:
    if stencil_cells == 4:
        return (-2, -1, 0, 1)
    if stencil_cells == 6:
        return (-3, -2, -1, 0, 1, 2)
    raise ValueError("Only the recorded four- and six-cell stencils are supported")


class DirectionalHLLRoeCorrection(nn.Module):
    """One shared, orientation-aware flux network for x and y faces."""

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ) -> None:
        super().__init__()
        self.shifts = stencil_shifts(stencil_cells)
        self.net = nn.Sequential(
            nn.Linear(3 * stencil_cells, width),
            nn.Tanh(),
            nn.Linear(width, width),
            nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def coefficients(self, state: torch.Tensor) -> torch.Tensor:
        # state is oriented [batch, transverse, normal, channels].
        values = primitive(state)
        cells = values.shape[-2]
        interfaces = torch.arange(1, cells, device=state.device)
        offsets = torch.tensor(self.shifts, device=state.device)
        indices = (interfaces[:, None] + offsets[None, :]).clamp(0, cells - 1)
        gathered = values.index_select(-2, indices.reshape(-1))
        gathered = gathered.reshape(*values.shape[:-2], cells - 1, len(self.shifts), 3)
        features = ((gathered - self.mean) / self.std).flatten(start_dim=-2)
        return torch.tanh(self.net(features))

    def forward_oriented(self, state: torch.Tensor) -> torch.Tensor:
        base = hll_faces_oriented(state)
        left = state[..., :-1, :]
        right = state[..., 1:, :]
        matrix, alpha, speeds = roe_waves_oriented(left, right)
        correction = torch.einsum(
            "...ij,...j->...i", matrix, self.coefficients(state) * speeds * alpha
        )
        interior = base[..., 1:-1, :] - 0.5 * correction
        return torch.cat((base[..., :1, :], interior, base[..., -1:, :]), dim=-2)


class HCFL2D(nn.Module):
    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ) -> None:
        super().__init__()
        self.flux_net = DirectionalHLLRoeCorrection(mean, std, width, stencil_cells)
        self.stencil_cells = stencil_cells

    def raw_fluxes(self, state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x_state = orient_state(state, "x")
        y_state = orient_state(state, "y")
        fx = self.flux_net.forward_oriented(x_state)
        gy = deorient_flux(self.flux_net.forward_oriented(y_state), "y")
        return fx, gy

    def project_raw_fluxes(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        raw_x, raw_y = raw_fluxes
        x_state = orient_state(state, "x")
        y_state = orient_state(state, "y")
        fx = strict_entropy_projection_oriented(raw_x, x_state)
        gy_oriented = strict_entropy_projection_oriented(orient_state(raw_y, "y"), y_state)
        return fx, deorient_flux(gy_oriented, "y")

    def project_raw_fluxes_training(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        raw_x, raw_y = raw_fluxes
        x_state = orient_state(state, "x")
        y_state = orient_state(state, "y")
        fx = training_entropy_projection_oriented(raw_x, x_state)
        gy_oriented = training_entropy_projection_oriented(orient_state(raw_y, "y"), y_state)
        return fx, deorient_flux(gy_oriented, "y")

    def projected_fluxes(self, state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.project_raw_fluxes(state, self.raw_fluxes(state))

    def projected_step(self, state: torch.Tensor, dt: float) -> torch.Tensor:
        return flux_divergence(state, *self.projected_fluxes(state), dt)

    def feasibility_loss(self, state: torch.Tensor) -> torch.Tensor:
        return self.feasibility_loss_from_raw(state, self.raw_fluxes(state))

    def feasibility_loss_from_raw(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> torch.Tensor:
        raw_x, raw_y = raw_fluxes
        residual_x = interface_entropy_residual_oriented(raw_x, orient_state(state, "x"))
        residual_y = interface_entropy_residual_oriented(
            orient_state(raw_y, "y"), orient_state(state, "y")
        )
        return torch.relu(residual_x).square().mean() + torch.relu(residual_y).square().mean()


def total_entropy(state: torch.Tensor) -> torch.Tensor:
    return entropy(state.double()).sum(dim=(-2, -1))


def boundary_entropy_increment(state: torch.Tensor, dt: float) -> torch.Tensor:
    lam = dt / (DOMAIN_LENGTH / state.shape[-2])
    x_state = orient_state(state, "x")
    y_state = orient_state(state, "y")
    qx = physical_entropy_flux_oriented(x_state)
    qy = physical_entropy_flux_oriented(y_state)
    return lam * (
        (qx[..., -1] - qx[..., 0]).sum(dim=-1)
        + (qy[..., -1] - qy[..., 0]).sum(dim=-1)
    )


def admissible(state: torch.Tensor) -> torch.Tensor:
    return torch.isfinite(state).all(dim=(-3, -2, -1)) & (
        state[..., 0].amin(dim=(-2, -1)) >= H_FLOOR
    )


def _global_depth_blend(
    state: torch.Tensor,
    high_fluxes: tuple[torch.Tensor, torch.Tensor],
    low_fluxes: tuple[torch.Tensor, torch.Tensor],
    dt: float,
    bisections: int = 36,
) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
    low_state = flux_divergence(state, *low_fluxes, dt)
    high_state = flux_divergence(state, *high_fluxes, dt)
    if not bool(admissible(low_state).all()):
        raise RuntimeError("Low-order 2-D SWE step is not depth-admissible")
    theta = torch.ones(state.shape[0], device=state.device, dtype=state.dtype)
    need = ~admissible(high_state)
    if bool(need.any()):
        lo = torch.zeros(int(need.sum()), device=state.device, dtype=torch.float64)
        hi = torch.ones_like(lo)
        base = low_state[need].double()
        delta = (high_state - low_state)[need].double()
        for _ in range(bisections):
            mid = 0.5 * (lo + hi)
            trial = base + mid[:, None, None, None] * delta
            okay = admissible(trial)
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        theta[need] = (lo.to(state.dtype) - 1.0e-6).clamp_min(0.0)
    shaped = theta[:, None, None, None]
    result = tuple(low + shaped * (high - low) for high, low in zip(high_fluxes, low_fluxes))
    rounded_state = flux_divergence(state, result[0], result[1], dt)
    rounded_bad = ~admissible(rounded_state)
    if bool(rounded_bad.any()):
        theta[rounded_bad] = 0.0
        shaped = theta[:, None, None, None]
        result = tuple(low + shaped * (high - low) for high, low in zip(high_fluxes, low_fluxes))
    return (result[0], result[1]), theta


def _global_entropy_blend(
    state: torch.Tensor,
    high_fluxes: tuple[torch.Tensor, torch.Tensor],
    low_fluxes: tuple[torch.Tensor, torch.Tensor],
    dt: float,
    bisections: int = 42,
) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
    low_state = flux_divergence(state, *low_fluxes, dt)
    high_state = flux_divergence(state, *high_fluxes, dt)
    target = total_entropy(state) - boundary_entropy_increment(state, dt)
    beta = torch.ones(state.shape[0], device=state.device, dtype=state.dtype)
    need = total_entropy(high_state) > target + ENTROPY_TOL
    if bool(need.any()):
        lo = torch.zeros(int(need.sum()), device=state.device, dtype=torch.float64)
        hi = torch.ones_like(lo)
        base = low_state[need].double()
        delta = (high_state - low_state)[need].double()
        target_selected = target[need]
        for _ in range(bisections):
            mid = 0.5 * (lo + hi)
            trial = base + mid[:, None, None, None] * delta
            okay = total_entropy(trial) <= target_selected
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        beta[need] = (lo.to(state.dtype) - 1.0e-6).clamp_min(0.0)
    shaped = beta[:, None, None, None]
    result = tuple(low + shaped * (high - low) for high, low in zip(high_fluxes, low_fluxes))
    rounded_state = flux_divergence(state, result[0], result[1], dt)
    rounded_balance = (
        total_entropy(rounded_state)
        - total_entropy(state)
        + boundary_entropy_increment(state, dt)
    )
    rounded_bad = rounded_balance > ENTROPY_TOL
    if bool(rounded_bad.any()):
        beta[rounded_bad] = 0.0
        shaped = beta[:, None, None, None]
        result = tuple(low + shaped * (high - low) for high, low in zip(high_fluxes, low_fluxes))
    return (result[0], result[1]), beta


@torch.no_grad()
def advance_hll_interval(
    state: torch.Tensor, interval: float = SNAPSHOT_DT, cfl: float = 0.42
) -> tuple[torch.Tensor, dict[str, float]]:
    remaining = interval
    stats: defaultdict[str, float] = defaultdict(float)
    stats["maximum_entropy_balance"] = -float("inf")
    stats["maximum_conservation_closure"] = 0.0
    while remaining > 1.0e-10:
        dt = min(remaining, cfl / max(maximum_2d_rate(state), 1.0e-12))
        for _ in range(20):
            fluxes = hll_fluxes(state)
            candidate = flux_divergence(state, *fluxes, dt)
            balance = total_entropy(candidate) - total_entropy(state) + boundary_entropy_increment(state, dt)
            if bool(admissible(candidate).all()) and bool((balance <= ENTROPY_TOL).all()):
                break
            dt *= 0.5
            stats["dt_halvings"] += 1.0
        else:
            raise RuntimeError("Could not establish the low-order SWE premise")
        stats["maximum_entropy_balance"] = max(
            stats["maximum_entropy_balance"], float(balance.max())
        )
        conservation = (
            candidate.double().sum(dim=(-3, -2))
            - state.double().sum(dim=(-3, -2))
            + _boundary_conservation_increment(state, dt).double()
        )
        conservation_scale = state.double().abs().sum(dim=(-3, -2)).clamp_min(1.0)
        stats["maximum_conservation_closure"] = max(
            stats["maximum_conservation_closure"],
            float((conservation.abs() / conservation_scale).max()),
        )
        stats["substeps"] += 1.0
        state = candidate
        remaining -= dt
    return state, dict(stats)


@torch.no_grad()
def advance_hcfl_interval(
    model: HCFL2D,
    state: torch.Tensor,
    interval: float = SNAPSHOT_DT,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, float]]:
    remaining = interval
    stats: defaultdict[str, float] = defaultdict(float)
    stats["minimum_depth_theta"] = 1.0
    stats["minimum_entropy_beta"] = 1.0
    stats["maximum_entropy_balance"] = -float("inf")
    stats["maximum_interface_residual"] = -float("inf")
    stats["maximum_conservation_closure"] = 0.0
    while remaining > 1.0e-10:
        dt = min(remaining, cfl / max(maximum_2d_rate(state), 1.0e-12))
        for _ in range(20):
            low_fluxes = hll_fluxes(state)
            low_state = flux_divergence(state, *low_fluxes, dt)
            low_balance = total_entropy(low_state) - total_entropy(state) + boundary_entropy_increment(state, dt)
            if bool(admissible(low_state).all()) and bool((low_balance <= ENTROPY_TOL).all()):
                break
            dt *= 0.5
            stats["low_order_dt_halvings"] += 1.0
        else:
            raise RuntimeError("Could not establish the low-order SWE premise")

        high_fluxes = model.projected_fluxes(state)
        depth_fluxes, theta = _global_depth_blend(state, high_fluxes, low_fluxes, dt)
        final_fluxes, beta = _global_entropy_blend(state, depth_fluxes, low_fluxes, dt)
        candidate = flux_divergence(state, *final_fluxes, dt)
        if not bool(admissible(candidate).all()):
            raise RuntimeError("HCFL 2-D deployment lost positive depth")
        balance = total_entropy(candidate) - total_entropy(state) + boundary_entropy_increment(state, dt)

        oriented_x = final_fluxes[0]
        oriented_y = orient_state(final_fluxes[1], "y")
        residual_x = interface_entropy_residual_oriented(oriented_x.double(), orient_state(state.double(), "x"))
        residual_y = interface_entropy_residual_oriented(oriented_y.double(), orient_state(state.double(), "y"))
        stats["minimum_depth_theta"] = min(stats["minimum_depth_theta"], float(theta.min()))
        stats["minimum_entropy_beta"] = min(stats["minimum_entropy_beta"], float(beta.min()))
        stats["depth_limiter_active"] += float((theta < 1.0 - 1.0e-7).sum())
        stats["entropy_limiter_active"] += float((beta < 1.0 - 1.0e-7).sum())
        stats["batch_substeps"] += float(beta.numel())
        stats["maximum_entropy_balance"] = max(stats["maximum_entropy_balance"], float(balance.max()))
        stats["maximum_interface_residual"] = max(
            stats["maximum_interface_residual"], float(residual_x.max()), float(residual_y.max())
        )
        conservation = (
            candidate.double().sum(dim=(-3, -2))
            - state.double().sum(dim=(-3, -2))
            + _boundary_conservation_increment(state, dt).double()
        )
        conservation_scale = state.double().abs().sum(dim=(-3, -2)).clamp_min(1.0)
        stats["maximum_conservation_closure"] = max(
            stats["maximum_conservation_closure"],
            float((conservation.abs() / conservation_scale).max()),
        )
        stats["substeps"] += 1.0
        state = candidate
        remaining -= dt
    return state, dict(stats)


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}


def merge_stats(total: defaultdict[str, float], step: dict[str, float]) -> None:
    for key, value in step.items():
        if key.startswith("minimum_"):
            total[key] = min(total.get(key, 1.0), value)
        elif key.startswith("maximum_"):
            total[key] = max(total.get(key, -float("inf")), value)
        else:
            total[key] += value


def conserved_to_claw(state: torch.Tensor) -> torch.Tensor:
    return state[..., [1, 2, 0]]


def claw_to_conserved(state: torch.Tensor) -> torch.Tensor:
    return state[..., [2, 0, 1]]


def primitive_channel_scale(training: torch.Tensor) -> torch.Tensor:
    values = primitive(training)
    return values.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-5)


def primitive_nrmse(
    prediction: torch.Tensor, reference: torch.Tensor, scale: torch.Tensor
) -> float:
    normalized = (primitive(prediction) - primitive(reference)) / scale.to(prediction.device)
    return float(torch.sqrt((normalized.double().square()).mean()))


@dataclass
class SequenceAudit:
    completed: bool
    minimum_depth: float
    rollout_nrmse: float | None
    final_nrmse: float | None
    maximum_conservation_residual: float
    maximum_entropy_balance: float | None
    entropy_violation_rate: float | None
    height_tv_excess: float | None
    height_range_overshoot: float | None


def _boundary_conservation_increment(state: torch.Tensor, dt: float) -> torch.Tensor:
    lam = dt / (DOMAIN_LENGTH / state.shape[-2])
    x_state = orient_state(state, "x")
    y_state = orient_state(state, "y")
    fx = physical_flux_oriented(x_state)
    fy = deorient_flux(physical_flux_oriented(y_state), "y")
    return lam * (
        (fx[..., -1, :] - fx[..., 0, :]).sum(dim=-2)
        + (fy[..., -1, :, :] - fy[..., 0, :, :]).sum(dim=-2)
    )


def height_total_variation(state: torch.Tensor) -> torch.Tensor:
    h = state[..., 0]
    return (h[..., 1:] - h[..., :-1]).abs().sum(dim=(-2, -1)) + (
        h[..., 1:, :] - h[..., :-1, :]
    ).abs().sum(dim=(-2, -1))


def audit_sequence(
    prediction: torch.Tensor,
    reference: torch.Tensor,
    scale: torch.Tensor,
    dt: float = SNAPSHOT_DT,
) -> SequenceAudit:
    """Black-box saved-frame audit; prediction/reference have [T,Y,X,3]."""
    finite = bool(torch.isfinite(prediction).all())
    minimum_depth = float(prediction[..., 0].min()) if finite else -float("inf")
    completed = finite and minimum_depth >= H_FLOOR
    if completed:
        rollout = primitive_nrmse(prediction[1:], reference[1:], scale)
        final = primitive_nrmse(prediction[-1:], reference[-1:], scale)
    else:
        rollout = None
        final = None

    initial_sum = prediction[0].double().sum(dim=(-3, -2)) if finite else torch.zeros(3, dtype=torch.float64)
    cumulative = torch.zeros(3, dtype=torch.float64, device=prediction.device)
    conservation_max = 0.0
    entropy_max = -float("inf")
    violations = 0
    for index in range(prediction.shape[0] - 1):
        left = prediction[index]
        right = prediction[index + 1]
        if not bool(torch.isfinite(left).all() and torch.isfinite(right).all()):
            conservation_max = float("inf")
            entropy_max = float("inf")
            violations += 1
            continue
        midpoint = 0.5 * (left + right)
        cumulative = cumulative + _boundary_conservation_increment(midpoint[None], dt)[0].double()
        residual = right.double().sum(dim=(-3, -2)) - initial_sum + cumulative
        normalization = torch.maximum(
            prediction[0].double().abs().sum(dim=(-3, -2)), torch.ones(3, device=prediction.device, dtype=torch.float64)
        )
        conservation_max = max(conservation_max, float((residual.abs() / normalization).max()))
        balance = total_entropy(right[None])[0] - total_entropy(left[None])[0] + boundary_entropy_increment(midpoint[None], dt)[0]
        entropy_max = max(entropy_max, float(balance))
        violations += int(float(balance) > ENTROPY_TOL)

    if completed:
        tv_prediction = float(height_total_variation(prediction[-1:])[0])
        tv_reference = float(height_total_variation(reference[-1:])[0])
        tv_excess = max(0.0, tv_prediction - tv_reference) / max(tv_reference, 1.0e-12)
        ref_min = float(reference[-1, ..., 0].min())
        ref_max = float(reference[-1, ..., 0].max())
        pred_min = float(prediction[-1, ..., 0].min())
        pred_max = float(prediction[-1, ..., 0].max())
        height_overshoot = max(ref_min - pred_min, pred_max - ref_max, 0.0)
    else:
        tv_excess = None
        height_overshoot = None
    return SequenceAudit(
        completed=completed,
        minimum_depth=minimum_depth,
        rollout_nrmse=rollout,
        final_nrmse=final,
        maximum_conservation_residual=conservation_max,
        # Entropy is not defined once depth is nonpositive.  The failure and
        # minimum depth remain reported, but post-failure entropy is NA.
        maximum_entropy_balance=entropy_max if completed else None,
        entropy_violation_rate=(
            violations / max(1, prediction.shape[0] - 1) if completed else None
        ),
        height_tv_excess=tv_excess,
        height_range_overshoot=height_overshoot,
    )
