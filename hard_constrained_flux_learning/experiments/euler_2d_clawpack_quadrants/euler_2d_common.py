"""Shared numerical core for the audited 2-D Euler quadrant experiment.

State order is always ``(rho, rho*u, rho*v, E)``.  Directional routines see
``(rho, m_normal, m_tangent, E)`` so one flux network can be shared by x and y.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import torch
from torch import nn


GAMMA = 1.4
DOMAIN_LENGTH = 1.0
N_COARSE = 64
RHO_FLOOR = 1.0e-7
PRESSURE_FLOOR = 1.0e-7
ENTROPY_TOL = 5.0e-7


def pressure_raw(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0]
    safe_rho = rho.clamp_min(RHO_FLOOR)
    kinetic = 0.5 * (
        state[..., 1].square() + state[..., 2].square()
    ) / safe_rho
    return (GAMMA - 1.0) * (state[..., 3] - kinetic)


def primitive(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    return torch.stack(
        (
            rho,
            state[..., 1] / rho,
            state[..., 2] / rho,
            pressure_raw(state),
        ),
        dim=-1,
    )


def physical_flux_oriented(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    normal_momentum = state[..., 1]
    tangent_momentum = state[..., 2]
    energy = state[..., 3]
    normal_velocity = normal_momentum / rho
    tangent_velocity = tangent_momentum / rho
    pressure = pressure_raw(state)
    return torch.stack(
        (
            normal_momentum,
            normal_momentum * normal_velocity + pressure,
            normal_momentum * tangent_velocity,
            normal_velocity * (energy + pressure),
        ),
        dim=-1,
    )


def orient_state(state: torch.Tensor, direction: str) -> torch.Tensor:
    if direction == "x":
        return state
    if direction == "y":
        return state.transpose(-3, -2)[..., [0, 2, 1, 3]]
    raise ValueError(direction)


def deorient_flux(flux: torch.Tensor, direction: str) -> torch.Tensor:
    if direction == "x":
        return flux
    if direction == "y":
        return flux[..., [0, 2, 1, 3]].transpose(-3, -2)
    raise ValueError(direction)


def _sound_speed(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    pressure = pressure_raw(state).clamp_min(PRESSURE_FLOOR)
    return torch.sqrt(GAMMA * pressure / rho)


def hll_pair(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    left_flux = physical_flux_oriented(left)
    right_flux = physical_flux_oriented(right)
    left_rho = left[..., 0].clamp_min(RHO_FLOOR)
    right_rho = right[..., 0].clamp_min(RHO_FLOOR)
    left_velocity = left[..., 1] / left_rho
    right_velocity = right[..., 1] / right_rho
    left_sound = _sound_speed(left)
    right_sound = _sound_speed(right)
    speed_left = torch.minimum(
        left_velocity - left_sound, right_velocity - right_sound
    )
    speed_right = torch.maximum(
        left_velocity + left_sound, right_velocity + right_sound
    )
    denominator = (speed_right - speed_left).clamp_min(1.0e-12)
    middle = (
        speed_right[..., None] * left_flux
        - speed_left[..., None] * right_flux
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / denominator[..., None]
    return torch.where(
        (speed_left >= 0.0)[..., None],
        left_flux,
        torch.where((speed_right <= 0.0)[..., None], right_flux, middle),
    )


def _signed_safe(value: torch.Tensor, epsilon: float = 1.0e-12) -> torch.Tensor:
    sign = torch.where(value >= 0.0, torch.ones_like(value), -torch.ones_like(value))
    return torch.where(value.abs() >= epsilon, value, sign * epsilon)


def hllc_pair(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    """Toro HLLC flux with a passively advected tangential velocity."""
    rho_left = left[..., 0].clamp_min(RHO_FLOOR)
    rho_right = right[..., 0].clamp_min(RHO_FLOOR)
    un_left = left[..., 1] / rho_left
    un_right = right[..., 1] / rho_right
    ut_left = left[..., 2] / rho_left
    ut_right = right[..., 2] / rho_right
    p_left = pressure_raw(left).clamp_min(PRESSURE_FLOOR)
    p_right = pressure_raw(right).clamp_min(PRESSURE_FLOOR)
    c_left = torch.sqrt(GAMMA * p_left / rho_left)
    c_right = torch.sqrt(GAMMA * p_right / rho_right)
    speed_left = torch.minimum(un_left - c_left, un_right - c_right)
    speed_right = torch.maximum(un_left + c_left, un_right + c_right)

    contact_denominator = _signed_safe(
        rho_left * (speed_left - un_left)
        - rho_right * (speed_right - un_right)
    )
    contact_speed = (
        p_right
        - p_left
        + rho_left * un_left * (speed_left - un_left)
        - rho_right * un_right * (speed_right - un_right)
    ) / contact_denominator

    def star_state(
        state: torch.Tensor,
        rho: torch.Tensor,
        un: torch.Tensor,
        ut: torch.Tensor,
        pressure: torch.Tensor,
        wave_speed: torch.Tensor,
    ) -> torch.Tensor:
        wave_jump = _signed_safe(wave_speed - contact_speed)
        density = rho * (wave_speed - un) / wave_jump
        specific_energy = state[..., 3] / rho
        energy = density * (
            specific_energy
            + (contact_speed - un)
            * (
                contact_speed
                + pressure / _signed_safe(rho * (wave_speed - un))
            )
        )
        return torch.stack(
            (density, density * contact_speed, density * ut, energy), dim=-1
        )

    star_left = star_state(
        left, rho_left, un_left, ut_left, p_left, speed_left
    )
    star_right = star_state(
        right, rho_right, un_right, ut_right, p_right, speed_right
    )
    flux_left = physical_flux_oriented(left)
    flux_right = physical_flux_oriented(right)
    star_flux_left = flux_left + speed_left[..., None] * (star_left - left)
    star_flux_right = flux_right + speed_right[..., None] * (star_right - right)
    return torch.where(
        (speed_left >= 0.0)[..., None],
        flux_left,
        torch.where(
            (contact_speed >= 0.0)[..., None],
            star_flux_left,
            torch.where(
                (speed_right > 0.0)[..., None], star_flux_right, flux_right
            ),
        ),
    )


def _faces_oriented(state: torch.Tensor, pair_flux) -> torch.Tensor:
    padded = torch.cat((state[..., :1, :], state, state[..., -1:, :]), dim=-2)
    return pair_flux(padded[..., :-1, :], padded[..., 1:, :])


def hll_faces_oriented(state: torch.Tensor) -> torch.Tensor:
    return _faces_oriented(state, hll_pair)


def hllc_faces_oriented(state: torch.Tensor) -> torch.Tensor:
    return _faces_oriented(state, hllc_pair)


def _directional_fluxes(state: torch.Tensor, face_flux) -> tuple[torch.Tensor, torch.Tensor]:
    flux_x = face_flux(orient_state(state, "x"))
    flux_y = deorient_flux(face_flux(orient_state(state, "y")), "y")
    return flux_x, flux_y


def hll_fluxes(state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    return _directional_fluxes(state, hll_faces_oriented)


def hllc_fluxes(state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    return _directional_fluxes(state, hllc_faces_oriented)


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


def maximum_2d_rate(state: torch.Tensor) -> float:
    values = primitive(state)
    sound = torch.sqrt(
        GAMMA
        * values[..., 3].clamp_min(PRESSURE_FLOOR)
        / values[..., 0].clamp_min(RHO_FLOOR)
    )
    dx = DOMAIN_LENGTH / state.shape[-2]
    rate = (
        (values[..., 1].abs() + sound) / dx
        + (values[..., 2].abs() + sound) / dx
    ).max()
    return float(rate)


def entropy(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    pressure = pressure_raw(state).clamp_min(PRESSURE_FLOOR)
    specific_entropy = torch.log(pressure) - GAMMA * torch.log(rho)
    return -rho * specific_entropy / (GAMMA - 1.0)


def entropy_variables_oriented(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    un = state[..., 1] / rho
    ut = state[..., 2] / rho
    pressure = pressure_raw(state).clamp_min(PRESSURE_FLOOR)
    specific_entropy = torch.log(pressure) - GAMMA * torch.log(rho)
    first = (
        (GAMMA - specific_entropy) / (GAMMA - 1.0)
        - rho * (un.square() + ut.square()) / (2.0 * pressure)
    )
    return torch.stack(
        (first, rho * un / pressure, rho * ut / pressure, -rho / pressure),
        dim=-1,
    )


def entropy_potential_oriented(state: torch.Tensor) -> torch.Tensor:
    return state[..., 1]


def physical_entropy_flux_oriented(state: torch.Tensor) -> torch.Tensor:
    rho = state[..., 0].clamp_min(RHO_FLOOR)
    return state[..., 1] / rho * entropy(state)


def roe_waves_oriented(
    left: torch.Tensor, right: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    rho_left = left[..., 0].clamp_min(RHO_FLOOR)
    rho_right = right[..., 0].clamp_min(RHO_FLOOR)
    root_left = torch.sqrt(rho_left)
    root_right = torch.sqrt(rho_right)
    denominator = (root_left + root_right).clamp_min(1.0e-12)
    un_left = left[..., 1] / rho_left
    un_right = right[..., 1] / rho_right
    ut_left = left[..., 2] / rho_left
    ut_right = right[..., 2] / rho_right
    p_left = pressure_raw(left).clamp_min(PRESSURE_FLOOR)
    p_right = pressure_raw(right).clamp_min(PRESSURE_FLOOR)
    h_left = (left[..., 3] + p_left) / rho_left
    h_right = (right[..., 3] + p_right) / rho_right
    un = (root_left * un_left + root_right * un_right) / denominator
    ut = (root_left * ut_left + root_right * ut_right) / denominator
    enthalpy = (root_left * h_left + root_right * h_right) / denominator
    speed = torch.sqrt(
        ((GAMMA - 1.0) * (enthalpy - 0.5 * (un.square() + ut.square())))
        .clamp_min(1.0e-12)
    )
    ones = torch.ones_like(un)
    zeros = torch.zeros_like(un)
    acoustic_minus = torch.stack(
        (ones, un - speed, ut, enthalpy - un * speed), dim=-1
    )
    contact = torch.stack(
        (ones, un, ut, 0.5 * (un.square() + ut.square())), dim=-1
    )
    shear = torch.stack((zeros, zeros, ones, ut), dim=-1)
    acoustic_plus = torch.stack(
        (ones, un + speed, ut, enthalpy + un * speed), dim=-1
    )
    matrix = torch.stack(
        (acoustic_minus, contact, shear, acoustic_plus), dim=-1
    )
    jump = right - left
    strengths = torch.linalg.solve(matrix, jump.unsqueeze(-1)).squeeze(-1)
    eigenvalues = torch.stack(
        (un - speed, un, un, un + speed), dim=-1
    )
    magnitude = eigenvalues.abs()
    entropy_fix = (0.1 * speed).clamp_min(1.0e-6)[..., None]
    fixed = torch.where(
        magnitude < entropy_fix,
        0.5 * (magnitude.square() / entropy_fix + entropy_fix),
        magnitude,
    )
    return matrix, strengths, fixed


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
    interior = proposal[..., 1:-1, :]
    state64 = state.double()
    left = state64[..., :-1, :]
    right = state64[..., 1:, :]
    normal = entropy_variables_oriented(right) - entropy_variables_oriented(left)
    bound = entropy_potential_oriented(right) - entropy_potential_oriented(left)
    norm_squared = normal.square().sum(dim=-1)
    epsilon = torch.finfo(proposal.dtype).eps
    output = interior
    for _ in range(5):
        flux64 = output.double()
        residual = (normal * flux64).sum(dim=-1) - bound
        scale = bound.abs() + (normal * flux64).abs().sum(dim=-1)
        target = -8.0 * epsilon * scale
        active = (residual > target) & (norm_squared > 1.0e-14)
        if not bool(active.any()):
            break
        coefficient = torch.zeros_like(residual)
        coefficient[active] = (
            residual[active] - target[active]
        ) / norm_squared[active]
        output = (flux64 - coefficient[..., None] * normal).to(proposal.dtype)
    return torch.cat((proposal[..., :1, :], output, proposal[..., -1:, :]), dim=-2)


def training_entropy_projection_oriented(
    proposal: torch.Tensor, state: torch.Tensor
) -> torch.Tensor:
    interior = proposal[..., 1:-1, :]
    left = state[..., :-1, :]
    right = state[..., 1:, :]
    normal = entropy_variables_oriented(right) - entropy_variables_oriented(left)
    bound = entropy_potential_oriented(right) - entropy_potential_oriented(left)
    residual = (normal * interior).sum(dim=-1) - bound
    norm_squared = normal.square().sum(dim=-1)
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
    raise ValueError("Only interface-centred four- and six-cell stencils are locked")


class DirectionalHLLCRoeCorrection(nn.Module):
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
            nn.Linear(4 * stencil_cells, width),
            nn.Tanh(),
            nn.Linear(width, width),
            nn.Tanh(),
            nn.Linear(width, 4),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def coefficients(self, state: torch.Tensor) -> torch.Tensor:
        values = primitive(state)
        cells = values.shape[-2]
        interfaces = torch.arange(1, cells, device=state.device)
        offsets = torch.tensor(self.shifts, device=state.device)
        indices = (interfaces[:, None] + offsets[None, :]).clamp(0, cells - 1)
        gathered = values.index_select(-2, indices.reshape(-1))
        gathered = gathered.reshape(
            *values.shape[:-2], cells - 1, len(self.shifts), 4
        )
        features = ((gathered - self.mean) / self.std).flatten(start_dim=-2)
        return torch.tanh(self.net(features))

    def forward_oriented(self, state: torch.Tensor) -> torch.Tensor:
        base = hllc_faces_oriented(state)
        left = state[..., :-1, :]
        right = state[..., 1:, :]
        matrix, strengths, speeds = roe_waves_oriented(left, right)
        correction = torch.einsum(
            "...ij,...j->...i",
            matrix,
            self.coefficients(state) * speeds * strengths,
        )
        interior = base[..., 1:-1, :] - 0.5 * correction
        return torch.cat((base[..., :1, :], interior, base[..., -1:, :]), dim=-2)


class HCFL2DEuler(nn.Module):
    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
        stencil_cells: int = 4,
    ) -> None:
        super().__init__()
        self.flux_net = DirectionalHLLCRoeCorrection(
            mean, std, width, stencil_cells
        )
        self.stencil_cells = stencil_cells

    def raw_fluxes(self, state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        raw_x = self.flux_net.forward_oriented(orient_state(state, "x"))
        raw_y = deorient_flux(
            self.flux_net.forward_oriented(orient_state(state, "y")), "y"
        )
        return raw_x, raw_y

    def project_raw_fluxes(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        raw_x, raw_y = raw_fluxes
        flux_x = strict_entropy_projection_oriented(
            raw_x, orient_state(state, "x")
        )
        flux_y_oriented = strict_entropy_projection_oriented(
            orient_state(raw_y, "y"), orient_state(state, "y")
        )
        return flux_x, deorient_flux(flux_y_oriented, "y")

    def project_raw_fluxes_training(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, torch.Tensor]:
        raw_x, raw_y = raw_fluxes
        flux_x = training_entropy_projection_oriented(
            raw_x, orient_state(state, "x")
        )
        flux_y_oriented = training_entropy_projection_oriented(
            orient_state(raw_y, "y"), orient_state(state, "y")
        )
        return flux_x, deorient_flux(flux_y_oriented, "y")

    def projected_fluxes(self, state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.project_raw_fluxes(state, self.raw_fluxes(state))

    def feasibility_loss_from_raw(
        self,
        state: torch.Tensor,
        raw_fluxes: tuple[torch.Tensor, torch.Tensor],
    ) -> torch.Tensor:
        raw_x, raw_y = raw_fluxes
        residual_x = interface_entropy_residual_oriented(
            raw_x, orient_state(state, "x")
        )
        residual_y = interface_entropy_residual_oriented(
            orient_state(raw_y, "y"), orient_state(state, "y")
        )
        return torch.relu(residual_x).square().mean() + torch.relu(
            residual_y
        ).square().mean()


def admissible(state: torch.Tensor) -> torch.Tensor:
    return (
        torch.isfinite(state).all(dim=(-3, -2, -1))
        & (state[..., 0].amin(dim=(-2, -1)) >= RHO_FLOOR)
        & (pressure_raw(state).amin(dim=(-2, -1)) >= PRESSURE_FLOOR)
    )


def total_entropy(state: torch.Tensor) -> torch.Tensor:
    return entropy(state.double()).sum(dim=(-2, -1))


def boundary_entropy_increment(state: torch.Tensor, dt: float) -> torch.Tensor:
    factor = dt / (DOMAIN_LENGTH / state.shape[-2])
    x_state = orient_state(state, "x")
    y_state = orient_state(state, "y")
    entropy_x = physical_entropy_flux_oriented(x_state)
    entropy_y = physical_entropy_flux_oriented(y_state)
    return factor * (
        (entropy_x[..., -1] - entropy_x[..., 0]).sum(dim=-1)
        + (entropy_y[..., -1] - entropy_y[..., 0]).sum(dim=-1)
    )


def boundary_conservation_increment(state: torch.Tensor, dt: float) -> torch.Tensor:
    factor = dt / (DOMAIN_LENGTH / state.shape[-2])
    flux_x = physical_flux_oriented(orient_state(state, "x"))
    flux_y = physical_flux_oriented(orient_state(state, "y"))
    x_increment = (flux_x[..., -1, :] - flux_x[..., 0, :]).sum(dim=-2)
    y_oriented = (flux_y[..., -1, :] - flux_y[..., 0, :]).sum(dim=-2)
    y_increment = y_oriented[..., [0, 2, 1, 3]]
    return factor * (x_increment + y_increment)


def _global_admissibility_blend(
    state: torch.Tensor,
    high_fluxes: tuple[torch.Tensor, torch.Tensor],
    low_fluxes: tuple[torch.Tensor, torch.Tensor],
    dt: float,
    bisections: int = 40,
) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
    low_state = flux_divergence(state, *low_fluxes, dt)
    high_state = flux_divergence(state, *high_fluxes, dt)
    if not bool(admissible(low_state).all()):
        raise RuntimeError("Low-order Euler step is not admissible")
    theta = torch.ones(state.shape[0], device=state.device, dtype=state.dtype)
    need = ~admissible(high_state)
    if bool(need.any()):
        lo = torch.zeros(int(need.sum()), device=state.device, dtype=torch.float64)
        hi = torch.ones_like(lo)
        base = low_state[need].double()
        delta = (high_state - low_state)[need].double()
        for _ in range(bisections):
            middle = 0.5 * (lo + hi)
            trial = base + middle[:, None, None, None] * delta
            okay = admissible(trial)
            lo = torch.where(okay, middle, lo)
            hi = torch.where(okay, hi, middle)
        theta[need] = (lo.to(state.dtype) - 1.0e-6).clamp_min(0.0)
    shaped = theta[:, None, None, None]
    result = tuple(
        low + shaped * (high - low)
        for high, low in zip(high_fluxes, low_fluxes)
    )
    rounded = flux_divergence(state, result[0], result[1], dt)
    bad = ~admissible(rounded)
    if bool(bad.any()):
        theta[bad] = 0.0
        shaped = theta[:, None, None, None]
        result = tuple(
            low + shaped * (high - low)
            for high, low in zip(high_fluxes, low_fluxes)
        )
    return (result[0], result[1]), theta


def _global_entropy_blend(
    state: torch.Tensor,
    high_fluxes: tuple[torch.Tensor, torch.Tensor],
    low_fluxes: tuple[torch.Tensor, torch.Tensor],
    dt: float,
    bisections: int = 44,
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
            middle = 0.5 * (lo + hi)
            trial = base + middle[:, None, None, None] * delta
            okay = total_entropy(trial) <= target_selected
            lo = torch.where(okay, middle, lo)
            hi = torch.where(okay, hi, middle)
        beta[need] = (lo.to(state.dtype) - 1.0e-6).clamp_min(0.0)
    shaped = beta[:, None, None, None]
    result = tuple(
        low + shaped * (high - low)
        for high, low in zip(high_fluxes, low_fluxes)
    )
    candidate = flux_divergence(state, result[0], result[1], dt)
    balance = (
        total_entropy(candidate)
        - total_entropy(state)
        + boundary_entropy_increment(state, dt)
    )
    bad = balance > ENTROPY_TOL
    if bool(bad.any()):
        beta[bad] = 0.0
        shaped = beta[:, None, None, None]
        result = tuple(
            low + shaped * (high - low)
            for high, low in zip(high_fluxes, low_fluxes)
        )
    return (result[0], result[1]), beta


def _conservation_closure(
    before: torch.Tensor, after: torch.Tensor, dt: float
) -> float:
    closure = (
        after.double().sum(dim=(-3, -2))
        - before.double().sum(dim=(-3, -2))
        + boundary_conservation_increment(before.double(), dt)
    )
    scale = before.double().abs().sum(dim=(-3, -2)).clamp_min(1.0)
    return float((closure.abs() / scale).max())


@torch.no_grad()
def advance_hllc_interval(
    state: torch.Tensor, interval: float, cfl: float = 0.42
) -> tuple[torch.Tensor, dict[str, float]]:
    remaining = interval
    statistics: defaultdict[str, float] = defaultdict(float)
    statistics["maximum_conservation_closure"] = 0.0
    while remaining > 1.0e-11:
        dt = min(remaining, cfl / max(maximum_2d_rate(state), 1.0e-12))
        for _ in range(24):
            fluxes = hllc_fluxes(state)
            candidate = flux_divergence(state, *fluxes, dt)
            if bool(admissible(candidate).all()):
                break
            dt *= 0.5
            statistics["dt_halvings"] += 1.0
        else:
            raise RuntimeError("HLLC could not retain an admissible Euler state")
        statistics["maximum_conservation_closure"] = max(
            statistics["maximum_conservation_closure"],
            _conservation_closure(state, candidate, dt),
        )
        statistics["substeps"] += 1.0
        state = candidate
        remaining -= dt
    return state, dict(statistics)


@torch.no_grad()
def advance_hcfl_interval(
    model: HCFL2DEuler,
    state: torch.Tensor,
    interval: float,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, float]]:
    remaining = interval
    statistics: defaultdict[str, float] = defaultdict(float)
    statistics["minimum_positivity_theta"] = 1.0
    statistics["minimum_entropy_beta"] = 1.0
    statistics["maximum_entropy_balance"] = -float("inf")
    statistics["maximum_interface_residual"] = -float("inf")
    statistics["maximum_conservation_closure"] = 0.0
    while remaining > 1.0e-11:
        dt = min(remaining, cfl / max(maximum_2d_rate(state), 1.0e-12))
        for _ in range(24):
            low_fluxes = hll_fluxes(state)
            low_state = flux_divergence(state, *low_fluxes, dt)
            low_balance = (
                total_entropy(low_state)
                - total_entropy(state)
                + boundary_entropy_increment(state, dt)
            )
            if bool(admissible(low_state).all()) and bool(
                (low_balance <= ENTROPY_TOL).all()
            ):
                break
            dt *= 0.5
            statistics["low_order_dt_halvings"] += 1.0
        else:
            raise RuntimeError("Could not establish the low-order Euler premise")

        high_fluxes = model.projected_fluxes(state)
        positive_fluxes, theta = _global_admissibility_blend(
            state, high_fluxes, low_fluxes, dt
        )
        final_fluxes, beta = _global_entropy_blend(
            state, positive_fluxes, low_fluxes, dt
        )
        candidate = flux_divergence(state, *final_fluxes, dt)
        if not bool(admissible(candidate).all()):
            raise RuntimeError("HCFL deployment lost Euler admissibility")
        balance = (
            total_entropy(candidate)
            - total_entropy(state)
            + boundary_entropy_increment(state, dt)
        )
        residual_x = interface_entropy_residual_oriented(
            final_fluxes[0].double(), orient_state(state.double(), "x")
        )
        residual_y = interface_entropy_residual_oriented(
            orient_state(final_fluxes[1].double(), "y"),
            orient_state(state.double(), "y"),
        )
        statistics["minimum_positivity_theta"] = min(
            statistics["minimum_positivity_theta"], float(theta.min())
        )
        statistics["minimum_entropy_beta"] = min(
            statistics["minimum_entropy_beta"], float(beta.min())
        )
        statistics["positivity_limiter_active"] += float(
            (theta < 1.0 - 1.0e-7).sum()
        )
        statistics["entropy_limiter_active"] += float(
            (beta < 1.0 - 1.0e-7).sum()
        )
        statistics["batch_substeps"] += float(beta.numel())
        statistics["maximum_entropy_balance"] = max(
            statistics["maximum_entropy_balance"], float(balance.max())
        )
        statistics["maximum_interface_residual"] = max(
            statistics["maximum_interface_residual"],
            float(residual_x.max()),
            float(residual_y.max()),
        )
        statistics["maximum_conservation_closure"] = max(
            statistics["maximum_conservation_closure"],
            _conservation_closure(state, candidate, dt),
        )
        statistics["substeps"] += 1.0
        state = candidate
        remaining -= dt
    return state, dict(statistics)


def merge_statistics(total: defaultdict[str, float], step: dict[str, float]) -> None:
    for key, value in step.items():
        if key.startswith("minimum_"):
            total[key] = min(total.get(key, 1.0), value)
        elif key.startswith("maximum_"):
            total[key] = max(total.get(key, -float("inf")), value)
        else:
            total[key] += value


def primitive_scale(training: torch.Tensor) -> torch.Tensor:
    values = primitive(training)
    return values.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-5)


def primitive_nmae(
    prediction: torch.Tensor, reference: torch.Tensor, scale: torch.Tensor
) -> float:
    normalized = (primitive(prediction) - primitive(reference)) / scale.to(
        prediction.device
    )
    return float(normalized.double().abs().mean())


def primitive_nrmse(
    prediction: torch.Tensor, reference: torch.Tensor, scale: torch.Tensor
) -> float:
    normalized = (primitive(prediction) - primitive(reference)) / scale.to(
        prediction.device
    )
    return float(torch.sqrt(normalized.double().square().mean()))


def primitive_channel_mae(
    prediction: torch.Tensor, reference: torch.Tensor
) -> torch.Tensor:
    difference = (primitive(prediction) - primitive(reference)).double().abs()
    return difference.reshape(-1, 4).mean(dim=0)


def total_variation_density(state: torch.Tensor) -> torch.Tensor:
    density = state[..., 0]
    return (
        (density[..., 1:] - density[..., :-1]).abs().sum(dim=(-2, -1))
        + (density[..., 1:, :] - density[..., :-1, :]).abs().sum(dim=(-2, -1))
    )


@dataclass
class AccuracyAudit:
    completed: bool
    minimum_density: float
    minimum_pressure: float
    nmae: float | None
    nrmse: float | None
    density_mae: float | None
    x_velocity_mae: float | None
    y_velocity_mae: float | None
    pressure_mae: float | None
    density_tv_excess: float | None
    density_range_overshoot: float | None
    pressure_range_overshoot: float | None


def audit_accuracy(
    prediction: torch.Tensor,
    reference: torch.Tensor,
    scale: torch.Tensor,
) -> AccuracyAudit:
    finite = bool(torch.isfinite(prediction).all())
    minimum_density = float(prediction[..., 0].min()) if finite else -float("inf")
    minimum_pressure = (
        float(pressure_raw(prediction).min()) if finite else -float("inf")
    )
    completed = (
        finite
        and minimum_density >= RHO_FLOOR
        and minimum_pressure >= PRESSURE_FLOOR
    )
    if not completed:
        return AccuracyAudit(
            completed=False,
            minimum_density=minimum_density,
            minimum_pressure=minimum_pressure,
            nmae=None,
            nrmse=None,
            density_mae=None,
            x_velocity_mae=None,
            y_velocity_mae=None,
            pressure_mae=None,
            density_tv_excess=None,
            density_range_overshoot=None,
            pressure_range_overshoot=None,
        )
    channel_mae = primitive_channel_mae(prediction[1:], reference[1:])
    prediction_primitive = primitive(prediction[-1])
    reference_primitive = primitive(reference[-1])
    tv_prediction = float(total_variation_density(prediction[-1:])[0])
    tv_reference = float(total_variation_density(reference[-1:])[0])
    density_overshoot = max(
        float(reference_primitive[..., 0].min() - prediction_primitive[..., 0].min()),
        float(prediction_primitive[..., 0].max() - reference_primitive[..., 0].max()),
        0.0,
    )
    pressure_overshoot = max(
        float(reference_primitive[..., 3].min() - prediction_primitive[..., 3].min()),
        float(prediction_primitive[..., 3].max() - reference_primitive[..., 3].max()),
        0.0,
    )
    return AccuracyAudit(
        completed=True,
        minimum_density=minimum_density,
        minimum_pressure=minimum_pressure,
        nmae=primitive_nmae(prediction[1:], reference[1:], scale),
        nrmse=primitive_nrmse(prediction[1:], reference[1:], scale),
        density_mae=float(channel_mae[0]),
        x_velocity_mae=float(channel_mae[1]),
        y_velocity_mae=float(channel_mae[2]),
        pressure_mae=float(channel_mae[3]),
        density_tv_excess=max(0.0, tv_prediction - tv_reference)
        / max(tv_reference, 1.0e-12),
        density_range_overshoot=density_overshoot,
        pressure_range_overshoot=pressure_overshoot,
    )


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {
        key: value.detach().cpu().clone() for key, value in model.state_dict().items()
    }


def parameter_count(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters())


def stack_mean(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=np.float64)
    return float(array.mean())
