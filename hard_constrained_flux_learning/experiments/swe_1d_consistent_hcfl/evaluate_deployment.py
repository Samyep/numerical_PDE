"""Held-out periodic/nonperiodic deployment audit for consistent SWE HCFL.

The neural flux is used only at interior interfaces for transmissive runs.
The two physical boundary fluxes are classical SWE fluxes, so no boundary
condition is learned.  Fine HLL-2048 trajectories are conservatively
restricted to 512 cells only for scoring; the light-gray plot curves retain
all 2048 points.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

import run_swe_consistent as base


TARGET_CELLS = 512
REFERENCE_CELLS = 2048
EVAL_SNAPSHOTS = 51
PERIODIC_CASES = {
    "smooth_wave": "smooth wave",
    "dam_break": "periodic dam break",
    "counterflow_collision": "counterflow collision",
    "supercritical_right": "right-going supercritical",
}
NONPERIODIC_CASES = {
    "centered_dam": "centered dam break",
    "transcritical_right": "right-going transcritical",
    "left_exit": "left-going wave exits",
    "right_exit": "right-going wave exits",
}
COLORS = {4: "#0072B2", 6: "#D55E00"}
FAMILIES = {
    "hll_roe": "HLL + Roe correction",
    "central_nonnegative": "Central + nonnegative Roe + feasibility",
}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def initial_condition(name: str, cells: int) -> np.ndarray:
    x = (np.arange(cells, dtype=np.float64) + 0.5) / cells
    if name == "smooth_wave":
        h = 1.05 + 0.18 * np.sin(2.0 * np.pi * x) + 0.04 * np.sin(
            4.0 * np.pi * x + 0.3
        )
        u = 0.42 * np.sin(2.0 * np.pi * x + 0.8)
    elif name == "dam_break":
        inside = (x >= 0.24) & (x < 0.64)
        h = np.where(inside, 1.85, 0.55)
        u = np.where(inside, 0.12, -0.08)
    elif name == "counterflow_collision":
        h = 1.0 + 0.08 * np.cos(2.0 * np.pi * x)
        u = np.where(x < 0.5, 1.65, -1.65)
    elif name == "supercritical_right":
        h = np.where((x >= 0.20) & (x < 0.58), 1.35, 0.72)
        u = 1.38 * np.sqrt(base.G * h)
    elif name == "centered_dam":
        h = np.where(x < 0.48, 1.9, 0.58)
        u = np.zeros_like(x)
    elif name == "transcritical_right":
        h = np.where(x < 0.46, 1.25, 0.72)
        froude = np.where(x < 0.46, 0.82, 1.24)
        u = froude * np.sqrt(base.G * h)
    elif name == "left_exit":
        h = 0.82 + 0.52 * np.exp(-((x - 0.055) / 0.026) ** 2)
        u = -1.32 * np.sqrt(base.G * h)
    elif name == "right_exit":
        h = 0.82 + 0.52 * np.exp(-((x - 0.945) / 0.026) ** 2)
        u = 1.32 * np.sqrt(base.G * h)
    else:
        raise KeyError(name)
    return base.primitive_to_conservative(h, u)


def np_nonperiodic_rhs(state: np.ndarray, dx: float) -> np.ndarray:
    interior = base.np_hll_pair(state[:, :-1], state[:, 1:])
    flux = np.concatenate(
        [base.np_flux(state[:, :1]), interior, base.np_flux(state[:, -1:])],
        axis=1,
    )
    return -(flux[:, 1:] - flux[:, :-1]) / dx


def np_nonperiodic_ssprk2(
    state: np.ndarray,
    dt: float,
    dx: float,
) -> np.ndarray:
    first = state + dt * np_nonperiodic_rhs(state, dx)
    if not np.isfinite(first).all() or (first[..., 0] <= 0.0).any():
        raise RuntimeError("Nonperiodic HLL reference failed at RK stage one")
    second = first + dt * np_nonperiodic_rhs(first, dx)
    result = 0.5 * state + 0.5 * second
    if not np.isfinite(result).all() or (result[..., 0] <= 0.0).any():
        raise RuntimeError("Nonperiodic HLL reference failed at RK completion")
    return result


def strict_hll_rollout(
    initial: np.ndarray,
    periodic: bool,
    snapshots: int = EVAL_SNAPSHOTS,
) -> tuple[np.ndarray, dict[str, Any]]:
    state = np.asarray(initial, dtype=np.float64)[None].copy()
    cells = state.shape[1]
    dx = 1.0 / cells
    saved = [state[0].copy()]
    substeps = 0
    minimum_depth = float(state[..., 0].min())
    started = time.perf_counter()
    for _ in range(1, snapshots):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            h = state[..., 0]
            u = state[..., 1] / h
            speed = float(np.max(np.abs(u) + np.sqrt(base.G * h)))
            dt = min(remaining, 0.20 * dx / max(speed, 1.0e-12))
            if periodic:
                state = base.np_periodic_ssprk2(state, dt, dx)
            else:
                state = np_nonperiodic_ssprk2(state, dt, dx)
            minimum_depth = min(minimum_depth, float(state[..., 0].min()))
            remaining -= dt
            substeps += 1
        saved.append(state[0].copy())
    return np.stack(saved), {
        "minimum_depth": minimum_depth,
        "substeps": substeps,
        "seconds": time.perf_counter() - started,
        "integrator": "HLL + SSP-RK2",
    }


def conservative_restrict(values: np.ndarray, target_cells: int) -> np.ndarray:
    source_cells = values.shape[-2]
    if source_cells % target_cells:
        raise ValueError("Source cells must be divisible by target cells")
    factor = source_cells // target_cells
    return values.reshape(*values.shape[:-2], target_cells, factor, 2).mean(-2)


def t_hll_pair(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    flux_left = base.t_flux(left)
    flux_right = base.t_flux(right)
    h_left = left[..., 0].clamp_min(base.H_FLOOR)
    h_right = right[..., 0].clamp_min(base.H_FLOOR)
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    c_left = torch.sqrt(base.G * h_left)
    c_right = torch.sqrt(base.G * h_right)
    speed_left = torch.minimum(u_left - c_left, u_right - c_right)
    speed_right = torch.maximum(u_left + c_left, u_right + c_right)
    middle = (
        speed_right[..., None] * flux_left
        - speed_left[..., None] * flux_right
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / (speed_right - speed_left).clamp_min(1.0e-10)[..., None]
    return torch.where(
        (speed_left >= 0.0)[..., None],
        flux_left,
        torch.where((speed_right <= 0.0)[..., None], flux_right, middle),
    )


def project_pairs(
    proposal: torch.Tensor,
    left: torch.Tensor,
    right: torch.Tensor,
) -> torch.Tensor:
    output = proposal
    epsilon = torch.finfo(proposal.dtype).eps
    left64 = left.double()
    right64 = right.double()
    normal = base.entropy_variables(right64) - base.entropy_variables(left64)
    bound = base.entropy_potential(right64) - base.entropy_potential(left64)
    norm_squared = (normal * normal).sum(dim=-1)
    for _ in range(4):
        flux64 = output.double()
        residual = (normal * flux64).sum(dim=-1) - bound
        scale = bound.abs() + (normal * flux64).abs().sum(dim=-1)
        target = -8.0 * epsilon * scale
        mask = (residual > target) & (norm_squared > 1.0e-14)
        if not bool(mask.any()):
            break
        amount = torch.zeros_like(residual)
        amount[mask] = (residual[mask] - target[mask]) / norm_squared[mask]
        output = (flux64 - amount[..., None] * normal).to(proposal.dtype)
    return output


def boundary_physical_flux(state: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    return base.t_flux(state[:, :1]), base.t_flux(state[:, -1:])


def padded_model_state(
    model: base.Solver,
    state: torch.Tensor,
) -> tuple[torch.Tensor, int]:
    shifts = tuple(model.flux_net.stencil_shifts)
    left_pad = max(shifts)
    right_pad = max(-shift for shift in shifts)
    padded = torch.cat(
        [
            state[:, :1].expand(-1, left_pad, -1),
            state,
            state[:, -1:].expand(-1, right_pad, -1),
        ],
        dim=1,
    )
    return padded, left_pad


def interior_learned_flux(model: base.Solver, state: torch.Tensor) -> torch.Tensor:
    padded, left_pad = padded_model_state(model, state)
    all_flux = model.flux(padded)
    return all_flux[:, left_pad : left_pad + state.shape[1] - 1]


def interior_raw_flux(model: base.Solver, state: torch.Tensor) -> torch.Tensor:
    padded, left_pad = padded_model_state(model, state)
    all_flux = model.raw_flux(padded)
    return all_flux[:, left_pad : left_pad + state.shape[1] - 1]


def nonperiodic_model_flux(model: base.Solver, state: torch.Tensor) -> torch.Tensor:
    left_boundary, right_boundary = boundary_physical_flux(state)
    return torch.cat(
        [left_boundary, interior_learned_flux(model, state), right_boundary],
        dim=1,
    )


def nonperiodic_low_flux(state: torch.Tensor) -> torch.Tensor:
    interior = project_pairs(
        t_hll_pair(state[:, :-1], state[:, 1:]),
        state[:, :-1],
        state[:, 1:],
    )
    left_boundary, right_boundary = boundary_physical_flux(state)
    return torch.cat([left_boundary, interior, right_boundary], dim=1)


def nonperiodic_step(
    state: torch.Tensor,
    flux: torch.Tensor,
    lam: float,
) -> torch.Tensor:
    return state - lam * (flux[:, 1:] - flux[:, :-1])


def boundary_entropy_target(state: torch.Tensor, lam: float) -> torch.Tensor:
    boundary_change = base.entropy_flux(state[:, -1]) - base.entropy_flux(
        state[:, 0]
    )
    return base.total_entropy(state) - lam * boundary_change.double()


def nonperiodic_depth_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    delta = high_flux - low_flux
    low_state = nonperiodic_step(state, low_flux, lam)
    if not bool(base.admissible(low_state).all()):
        raise RuntimeError("Nonperiodic low-order depth premise failed")
    batch, cells, _ = state.shape
    alpha = torch.ones(batch, cells + 1, dtype=state.dtype, device=state.device)
    for _ in range(12):
        flux = low_flux + alpha[..., None] * delta
        current = nonperiodic_step(state, flux, lam)
        good = base.admissible(current)
        if bool(good.all()):
            return flux, alpha
        bad = ~good
        count = int(bad.sum())
        lo = torch.zeros(count, dtype=torch.float64, device=state.device)
        hi = torch.ones(count, dtype=torch.float64, device=state.device)
        origin = low_state[bad].double()
        direction = (current - low_state)[bad].double()
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            okay = base.admissible(origin + mid[:, None] * direction)
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        theta = torch.ones(batch, cells, dtype=state.dtype, device=state.device)
        theta[bad] = torch.clamp(lo.to(state.dtype) - 1.0e-5, min=0.0)
        factors = torch.cat(
            [
                theta[:, :1],
                torch.minimum(theta[:, :-1], theta[:, 1:]),
                theta[:, -1:],
            ],
            dim=1,
        )
        alpha = alpha * factors
    flux = low_flux + alpha[..., None] * delta
    failed = (~base.admissible(nonperiodic_step(state, flux, lam))).any(dim=-1)
    if bool(failed.any()):
        alpha[failed] = 0.0
        flux = low_flux + alpha[..., None] * delta
    if not bool(base.admissible(nonperiodic_step(state, flux, lam)).all()):
        raise RuntimeError("Nonperiodic depth limiter failed")
    return flux, alpha


def nonperiodic_entropy_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    low_state = nonperiodic_step(state, low_flux, lam)
    high_state = nonperiodic_step(state, high_flux, lam)
    target = boundary_entropy_target(state, lam)
    beta = torch.ones(state.shape[0], dtype=state.dtype, device=state.device)
    need = base.total_entropy(high_state) > target
    if bool(need.any()):
        count = int(need.sum())
        lo = torch.zeros(count, dtype=torch.float64, device=state.device)
        hi = torch.ones(count, dtype=torch.float64, device=state.device)
        low_selected = low_state[need].double()
        direction = (high_state - low_state)[need].double()
        target_selected = target[need]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            okay = base.total_entropy(
                low_selected + mid[:, None, None] * direction
            ) <= target_selected
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)
        beta[need] = torch.clamp(lo.to(state.dtype) - 1.0e-6, min=0.0)
    flux = low_flux + beta[:, None, None] * (high_flux - low_flux)
    rounded_bad = base.total_entropy(nonperiodic_step(state, flux, lam)) > (
        target + base.FD_ENTROPY_TOLERANCE
    )
    if bool(rounded_bad.any()):
        beta[rounded_bad] = 0.0
        flux = low_flux + beta[:, None, None] * (high_flux - low_flux)
    return flux, beta


@torch.no_grad()
def advance_nonperiodic_snapshot(
    model: base.Solver,
    state: torch.Tensor,
) -> tuple[torch.Tensor, dict[str, float]]:
    dx = 1.0 / state.shape[1]
    remaining = base.DT_SNAPSHOT
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["maximum_entropy_residual"] = -float("inf")
    totals["maximum_entropy_balance_violation"] = -float("inf")
    while remaining > 1.0e-14:
        values = base.primitive(state)
        speed = float(
            (values[..., 1].abs() + torch.sqrt(base.G * values[..., 0])).max()
        )
        dt = min(remaining, 0.42 * dx / max(speed, 1.0e-12))
        for _ in range(24):
            lam = dt / dx
            low_flux = nonperiodic_low_flux(state)
            low_state = nonperiodic_step(state, low_flux, lam)
            target = boundary_entropy_target(state, lam)
            entropy_ok = bool(
                (
                    base.total_entropy(low_state)
                    <= target + base.FD_ENTROPY_TOLERANCE
                ).all()
            )
            if bool(base.admissible(low_state).all()) and entropy_ok:
                break
            dt *= 0.5
            totals["low_order_dt_halvings"] += 1.0
        else:
            raise RuntimeError("Could not establish nonperiodic low-order premise")

        high_flux = nonperiodic_model_flux(model, state)
        local_flux, alpha = nonperiodic_depth_limiter(
            state, high_flux, low_flux, lam
        )
        final_flux, beta = nonperiodic_entropy_limiter(
            state, local_flux, low_flux, lam
        )
        next_state = nonperiodic_step(state, final_flux, lam)
        if not bool(base.admissible(next_state).all()):
            raise RuntimeError("Safe nonperiodic SWE update lost positive depth")

        interior_residual = (
            (
                base.entropy_variables(state[:, 1:].double())
                - base.entropy_variables(state[:, :-1].double())
            )
            * final_flux[:, 1:-1].double()
        ).sum(-1) - (
            base.entropy_potential(state[:, 1:].double())
            - base.entropy_potential(state[:, :-1].double())
        )
        balance = base.total_entropy(next_state) - target
        totals["local_active"] += float((alpha < 1.0 - 1.0e-7).sum())
        totals["local_total"] += float(alpha.numel())
        totals["local_alpha_sum"] += float(alpha.sum())
        totals["fd_active"] += float((beta < 1.0 - 1.0e-7).sum())
        totals["fd_total"] += float(beta.numel())
        totals["fd_beta_sum"] += float(beta.sum())
        totals["fd_beta_min"] = min(totals["fd_beta_min"], float(beta.min()))
        totals["entropy_violations"] += float(
            (interior_residual > base.ENTROPY_RESIDUAL_TOLERANCE).sum()
        )
        totals["entropy_total"] += float(interior_residual.numel())
        totals["maximum_entropy_residual"] = max(
            totals["maximum_entropy_residual"], float(interior_residual.max())
        )
        totals["maximum_entropy_balance_violation"] = max(
            totals["maximum_entropy_balance_violation"], float(balance.max())
        )
        totals["substeps"] += 1.0
        state = next_state
        remaining -= dt
    return state, dict(totals)


@torch.no_grad()
def learned_rollout(
    model: base.Solver,
    initial: np.ndarray,
    periodic: bool,
) -> tuple[np.ndarray, dict[str, Any]]:
    state = torch.from_numpy(initial.astype(np.float32))[None]
    saved = [state[0].numpy().copy()]
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["maximum_entropy_residual"] = -float("inf")
    if periodic:
        totals["maximum_total_entropy_change"] = -float("inf")
    else:
        totals["maximum_entropy_balance_violation"] = -float("inf")
    minimum_depth = float(state[..., 0].min())
    started = time.perf_counter()
    for _ in range(1, EVAL_SNAPSHOTS):
        if periodic:
            state, step = base.advance_safe_snapshot(model, state)
        else:
            state, step = advance_nonperiodic_snapshot(model, state)
        base.merge_stats(totals, step)
        minimum_depth = min(minimum_depth, float(state[..., 0].min()))
        saved.append(state[0].numpy().copy())
    result = {
        "minimum_depth": minimum_depth,
        "entropy_violation_rate": totals["entropy_violations"]
        / max(totals["entropy_total"], 1.0),
        "maximum_entropy_residual": totals["maximum_entropy_residual"],
        "local_limiter_intervention_rate": totals["local_active"]
        / max(totals["local_total"], 1.0),
        "mean_local_alpha": totals["local_alpha_sum"]
        / max(totals["local_total"], 1.0),
        "fd_entropy_intervention_rate": totals["fd_active"]
        / max(totals["fd_total"], 1.0),
        "mean_fd_beta": totals["fd_beta_sum"]
        / max(totals["fd_total"], 1.0),
        "minimum_fd_beta": totals["fd_beta_min"],
        "substeps": int(totals["substeps"]),
        "seconds": time.perf_counter() - started,
        "integrator": "safe forward Euler",
    }
    if periodic:
        result["maximum_total_entropy_change"] = totals[
            "maximum_total_entropy_change"
        ]
    else:
        result["maximum_entropy_balance_violation"] = totals[
            "maximum_entropy_balance_violation"
        ]
    return np.stack(saved), result


def primitive_np(state: np.ndarray) -> np.ndarray:
    return np.stack([state[..., 0], state[..., 1] / state[..., 0]], axis=-1)


def total_variation(values: np.ndarray, periodic: bool) -> float:
    variation = float(np.abs(np.diff(values)).sum())
    if periodic:
        variation += float(abs(values[0] - values[-1]))
    return variation


def extrema_count(values: np.ndarray, periodic: bool) -> int:
    scale = max(float(np.ptp(values)), 1.0e-12)
    threshold = 1.0e-3 * scale
    if periodic:
        left = values - np.roll(values, 1)
        right = np.roll(values, -1) - values
    else:
        left = values[1:-1] - values[:-2]
        right = values[2:] - values[1:-1]
    return int(((left * right < 0.0) & (np.abs(left) > threshold) & (np.abs(right) > threshold)).sum())


def diagnostics(
    reference: np.ndarray,
    prediction: np.ndarray,
    state_std: np.ndarray,
    periodic: bool,
) -> dict[str, Any]:
    normalized = (prediction[1:] - reference[1:]) / state_std
    final_normalized = (prediction[-1] - reference[-1]) / state_std
    ref_primitive = primitive_np(reference[-1])
    pred_primitive = primitive_np(prediction[-1])
    tv_relative: list[float] = []
    range_violation: list[float] = []
    excess_extrema = 0
    for component in range(2):
        ref_values = ref_primitive[:, component]
        pred_values = pred_primitive[:, component]
        ref_tv = total_variation(ref_values, periodic)
        pred_tv = total_variation(pred_values, periodic)
        tv_relative.append((pred_tv - ref_tv) / max(ref_tv, 1.0e-12))
        ref_min = float(ref_values.min())
        ref_max = float(ref_values.max())
        span = max(ref_max - ref_min, 1.0e-12)
        range_violation.append(
            (max(float(pred_values.max()) - ref_max, 0.0)
             + max(ref_min - float(pred_values.min()), 0.0))
            / span
        )
        excess_extrema += max(
            extrema_count(pred_values, periodic)
            - extrema_count(ref_values, periodic),
            0,
        )
    mass = prediction[..., 0].sum(axis=1)
    momentum = prediction[..., 1].sum(axis=1)
    return {
        "rollout_nrmse": float(np.sqrt(np.mean(normalized * normalized))),
        "final_snapshot_nrmse": float(
            np.sqrt(np.mean(final_normalized * final_normalized))
        ),
        "minimum_depth": float(prediction[..., 0].min()),
        "mean_relative_final_tv_error": float(np.mean(tv_relative)),
        "mean_positive_final_tv_excess": float(
            np.mean(np.maximum(tv_relative, 0.0))
        ),
        "mean_normalized_range_violation": float(np.mean(range_violation)),
        "total_excess_significant_extrema": excess_extrema,
        "relative_mass_change": float(
            (mass[-1] - mass[0]) / max(abs(float(mass[0])), 1.0e-12)
        ),
        "relative_momentum_change": float(
            (momentum[-1] - momentum[0])
            / max(float(np.abs(prediction[0, :, 1]).sum()), 1.0e-12)
        ),
    }


@torch.no_grad()
def proposal_diagnostics(
    model: base.Solver,
    reference: np.ndarray,
    periodic: bool,
) -> dict[str, float]:
    """Verify the network is active and quantify pre-projection feasibility."""
    state = torch.from_numpy(reference.astype(np.float32))
    if periodic:
        raw = model.raw_flux(state)
        projected = model.flux(state)
        residual = base.entropy_residual(raw.double(), state.double())
        coefficients = model.flux_net.coefficients(state)
    else:
        raw = interior_raw_flux(model, state)
        projected = interior_learned_flux(model, state)
        left = state[:, :-1].double()
        right = state[:, 1:].double()
        residual = (
            (base.entropy_variables(right) - base.entropy_variables(left))
            * raw.double()
        ).sum(-1) - (
            base.entropy_potential(right) - base.entropy_potential(left)
        )
        padded, left_pad = padded_model_state(model, state)
        all_coefficients = model.flux_net.coefficients(padded)
        coefficients = all_coefficients[
            :, left_pad : left_pad + state.shape[1] - 1
        ]
    positive = torch.relu(residual)
    intervention = (projected - raw).abs().amax(dim=-1) > 1.0e-7
    output = {
        "raw_proposal_positive_residual_rate": float(
            (residual > 0.0).float().mean()
        ),
        "raw_proposal_violation_rate_at_1e-8": float(
            (residual > 1.0e-8).float().mean()
        ),
        "raw_proposal_mean_positive_residual": float(positive.mean()),
        "raw_proposal_maximum_residual": float(residual.max()),
        "hard_projection_intervention_rate": float(intervention.float().mean()),
    }
    if isinstance(model.flux_net, base.HLLRoeCorrectionFlux):
        output.update(
            {
                "mean_absolute_learned_roe_correction_coefficient": float(
                    coefficients.abs().mean()
                ),
                "maximum_absolute_learned_roe_correction_coefficient": float(
                    coefficients.abs().max()
                ),
            }
        )
    else:
        output.update(
            {
                "mean_absolute_roe_multiplier_change_from_one": float(
                    (coefficients - 1.0).abs().mean()
                ),
                "minimum_roe_multiplier": float(coefficients.min()),
                "maximum_roe_multiplier": float(coefficients.max()),
            }
        )
    return output


def aggregate(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for method in dict.fromkeys(row["method"] for row in rows):
        selected = [row for row in rows if row["method"] == method]
        output[method] = {
            "mean_rollout_nrmse": float(
                np.mean([row["rollout_nrmse"] for row in selected])
            ),
            "mean_final_snapshot_nrmse": float(
                np.mean([row["final_snapshot_nrmse"] for row in selected])
            ),
            "mean_positive_final_tv_excess": float(
                np.mean([row["mean_positive_final_tv_excess"] for row in selected])
            ),
            "total_excess_significant_extrema": float(
                np.sum([row["total_excess_significant_extrema"] for row in selected])
            ),
            "minimum_depth": float(min(row["minimum_depth"] for row in selected)),
        }
        for key in (
            "raw_proposal_positive_residual_rate",
            "raw_proposal_violation_rate_at_1e-8",
            "hard_projection_intervention_rate",
            "run_local_limiter_intervention_rate",
            "run_mean_local_alpha",
            "run_fd_entropy_intervention_rate",
            "run_mean_fd_beta",
        ):
            if all(
                key in row and row[key] not in ("", None) for row in selected
            ):
                output[method][f"mean_{key}"] = float(
                    np.mean([row[key] for row in selected])
                )
        if all(
            "run_minimum_fd_beta" in row
            and row["run_minimum_fd_beta"] not in ("", None)
            for row in selected
        ):
            output[method]["minimum_run_fd_beta"] = float(
                min(row["run_minimum_fd_beta"] for row in selected)
            )
        if all(
            "run_maximum_entropy_residual" in row
            and row["run_maximum_entropy_residual"] not in ("", None)
            for row in selected
        ):
            output[method]["maximum_run_entropy_residual"] = float(
                max(row["run_maximum_entropy_residual"] for row in selected)
            )
    return output


def wave_coverage(data: torch.Tensor) -> dict[str, float]:
    values = base.primitive(data)
    velocity = values[..., 1]
    gravity_speed = torch.sqrt(base.G * values[..., 0])
    froude = velocity.abs() / gravity_speed
    right_supercritical = velocity - gravity_speed > 0.0
    left_supercritical = velocity + gravity_speed < 0.0
    subcritical = ~(right_supercritical | left_supercritical)
    total = float(froude.numel())
    return {
        "maximum_absolute_froude": float(froude.max()),
        "median_absolute_froude": float(froude.median()),
        "subcritical_fraction": float(subcritical.sum()) / total,
        "right_supercritical_fraction": float(right_supercritical.sum()) / total,
        "left_supercritical_fraction": float(left_supercritical.sum()) / total,
    }


def evaluate_cases(
    cases: dict[str, str],
    periodic: bool,
    models: dict[str, base.Solver],
    state_std: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    trajectories: dict[str, Any] = {}
    for case, display in cases.items():
        print(json.dumps({"stage": "evaluate", "case": case}), flush=True)
        reference_native, reference_stats = strict_hll_rollout(
            initial_condition(case, REFERENCE_CELLS), periodic
        )
        reference = conservative_restrict(reference_native, TARGET_CELLS)
        native, native_stats = strict_hll_rollout(reference[0], periodic)
        native_row = {
            "boundary": "periodic" if periodic else "transmissive",
            "case": case,
            "display": display,
            "method": "native_hll_512",
            "family": "native_hll",
            "stencil_cells": None,
            **diagnostics(reference, native, state_std, periodic),
            **{f"run_{key}": value for key, value in native_stats.items()},
        }
        rows.append(native_row)
        predictions: dict[str, np.ndarray] = {}
        for arm in base.ARMS:
            prediction, run_stats = learned_rollout(
                models[arm.name], reference[0], periodic
            )
            predictions[arm.name] = prediction
            rows.append(
                {
                    "boundary": "periodic" if periodic else "transmissive",
                    "case": case,
                    "display": display,
                    "method": arm.name,
                    "family": arm.model_name,
                    "stencil_cells": arm.stencil_cells,
                    **diagnostics(reference, prediction, state_std, periodic),
                    **proposal_diagnostics(
                        models[arm.name], reference, periodic
                    ),
                    **{f"run_{key}": value for key, value in run_stats.items()},
                }
            )
        trajectories[case] = {
            "reference_native": reference_native,
            "reference": reference,
            "native": native,
            "predictions": predictions,
            "reference_stats": reference_stats,
        }
    return rows, trajectories


def plot_profiles(
    boundary: str,
    cases: dict[str, str],
    trajectories: dict[str, Any],
    family: str,
    output: Path,
) -> None:
    family_arms = [arm for arm in base.ARMS if arm.model_name == family]
    figure, axes = plt.subplots(2, len(cases), figsize=(14.8, 6.3), squeeze=False)
    fine_x = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    coarse_x = (np.arange(TARGET_CELLS) + 0.5) / TARGET_CELLS
    for column, (case, display) in enumerate(cases.items()):
        item = trajectories[case]
        reference = primitive_np(item["reference_native"][-1])
        native = primitive_np(item["native"][-1])
        for row, (component, label) in enumerate(((0, "h"), (1, "u"))):
            axis = axes[row, column]
            axis.plot(
                fine_x,
                reference[:, component],
                color="#B8B8B8",
                linewidth=1.2,
                label="HLL-2048 reference" if column == 0 and row == 0 else None,
                zorder=1,
            )
            axis.plot(
                coarse_x,
                native[:, component],
                color="#222222",
                linestyle="--",
                linewidth=1.15,
                label="HLL-512" if column == 0 and row == 0 else None,
                zorder=2,
            )
            for arm in family_arms:
                values = primitive_np(item["predictions"][arm.name][-1])
                axis.plot(
                    coarse_x,
                    values[:, component],
                    color=COLORS[arm.stencil_cells],
                    linewidth=1.15,
                    label=(
                        f"HCFL {arm.stencil_cells}-cell"
                        if column == 0 and row == 0
                        else None
                    ),
                    zorder=3,
                )
            axis.grid(alpha=0.18)
            axis.set_xlim(0.0, 1.0)
            if row == 0:
                axis.set_title(display)
            if column == 0:
                axis.set_ylabel(label)
            if row == 1:
                axis.set_xlabel("x")
    figure.suptitle(
        f"{FAMILIES[family]} — {boundary}, t={base.DT_SNAPSHOT * (EVAL_SNAPSHOTS - 1):.4f}",
        fontsize=13,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    figure.tight_layout(rect=(0.0, 0.07, 1.0, 0.94))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_summary(rows: list[dict[str, Any]], output: Path) -> None:
    boundaries = ("periodic", "transmissive")
    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    for axis, boundary in zip(axes, boundaries):
        selected = [row for row in rows if row["boundary"] == boundary]
        aggregate_rows = aggregate(selected)
        labels = ["HLL-512"]
        values = [aggregate_rows["native_hll_512"]["mean_rollout_nrmse"]]
        colors = ["#333333"]
        method_colors = {
            ("hll_roe", 4): "#0072B2",
            ("hll_roe", 6): "#56B4E9",
            ("central_nonnegative", 4): "#D55E00",
            ("central_nonnegative", 6): "#E69F00",
        }
        for family in FAMILIES:
            for cells in (4, 6):
                arm = next(
                    candidate
                    for candidate in base.ARMS
                    if candidate.model_name == family
                    and candidate.stencil_cells == cells
                )
                labels.append(f"{family.split('_')[0]}-{cells}")
                values.append(aggregate_rows[arm.name]["mean_rollout_nrmse"])
                colors.append(method_colors[(family, cells)])
        axis.bar(np.arange(len(values)), values, color=colors, alpha=0.88)
        axis.set_xticks(np.arange(len(values)), labels, rotation=25, ha="right")
        axis.set_ylabel("mean rollout NRMSE")
        axis.set_title(boundary)
        axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_stability_summary(rows: list[dict[str, Any]], output: Path) -> None:
    methods = ["native_hll_512", *(arm.name for arm in base.ARMS)]
    labels = ["HLL-512", "hll-4", "hll-6", "central-4", "central-6"]
    colors = ["#333333", "#0072B2", "#56B4E9", "#D55E00", "#E69F00"]
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 7.5))
    for column, boundary in enumerate(("periodic", "transmissive")):
        selected = [row for row in rows if row["boundary"] == boundary]
        grouped = {
            method: [row for row in selected if row["method"] == method]
            for method in methods
        }
        tv_values = [
            float(
                np.mean(
                    [row["mean_positive_final_tv_excess"] for row in grouped[method]]
                )
            )
            for method in methods
        ]
        extrema_values = [
            float(
                np.sum(
                    [
                        row["total_excess_significant_extrema"]
                        for row in grouped[method]
                    ]
                )
            )
            for method in methods
        ]
        for row_index, (values, ylabel) in enumerate(
            (
                (tv_values, "mean positive final TV excess"),
                (extrema_values, "total excess significant extrema"),
            )
        ):
            axis = axes[row_index, column]
            axis.bar(np.arange(len(methods)), values, color=colors, alpha=0.88)
            axis.set_xticks(
                np.arange(len(methods)), labels, rotation=25, ha="right"
            )
            axis.set_ylabel(ylabel)
            axis.grid(axis="y", alpha=0.2)
            if row_index == 0:
                axis.set_title(boundary)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_convergence(results: Path, seed: int, output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.1), sharey=True)
    for axis, family in zip(axes, FAMILIES):
        for arm in (item for item in base.ARMS if item.model_name == family):
            curve_path = results / f"training_curve_{arm.name}_seed{seed}.csv"
            with curve_path.open(newline="", encoding="utf-8") as handle:
                records = list(csv.DictReader(handle))
            updates = np.array([int(row["update"]) for row in records])
            metrics = np.array(
                [float(row["validation_rollout_nrmse"]) for row in records]
            )
            axis.plot(
                updates,
                metrics,
                color=COLORS[arm.stencil_cells],
                linewidth=1.35,
                label=f"{arm.stencil_cells}-cell",
            )
            best = int(np.argmin(metrics))
            axis.scatter(
                [updates[best]],
                [metrics[best]],
                color=COLORS[arm.stencil_cells],
                s=24,
                zorder=3,
            )
        axis.set_title(FAMILIES[family])
        axis.set_xlabel("optimizer updates")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    axes[0].set_ylabel("validation rollout NRMSE")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def boundary_operator_self_test(model: base.Solver) -> dict[str, Any]:
    cells = 24
    x = (torch.arange(cells, dtype=torch.float32) + 0.5) / cells
    h = 1.0 + 0.08 * torch.sin(2.0 * math.pi * x)
    u = 0.2 * torch.cos(2.0 * math.pi * x)
    state = torch.stack([h, h * u], dim=-1)[None]
    flux = nonperiodic_model_flux(model, state)
    left, right = boundary_physical_flux(state)
    boundary_error = max(
        float((flux[:, :1] - left).detach().abs().max()),
        float((flux[:, -1:] - right).detach().abs().max()),
    )
    modified = state.clone()
    modified[:, -1, 0] *= 1.7
    modified[:, -1, 1] *= -0.8
    modified_flux = nonperiodic_model_flux(model, modified)
    circular_leakage = float(
        (modified_flux[:, 1] - flux[:, 1]).detach().abs().max()
    )
    left_state = state[:, :-1].double()
    right_state = state[:, 1:].double()
    residual = (
        (base.entropy_variables(right_state) - base.entropy_variables(left_state))
        * flux[:, 1:-1].double()
    ).sum(-1) - (
        base.entropy_potential(right_state) - base.entropy_potential(left_state)
    )
    result = {
        "boundary_flux_max_error": boundary_error,
        "circular_wrap_leakage_at_first_interior_interface": circular_leakage,
        "maximum_interior_entropy_residual": float(residual.detach().max()),
    }
    if boundary_error > 1.0e-7 or circular_leakage > 1.0e-7:
        raise RuntimeError(f"Boundary operator self-test failed: {result}")
    if result["maximum_interior_entropy_residual"] > 1.0e-8:
        raise RuntimeError(f"Interior hard projection self-test failed: {result}")
    result["status"] = "pass"
    return result


def self_test() -> None:
    mean = np.array([1.0, 0.0], dtype=np.float32)
    std = np.ones(2, dtype=np.float32)
    for arm in base.ARMS:
        model = base.Solver(
            arm.model_name,
            mean,
            std,
            width=16,
            stencil_cells=arm.stencil_cells,
        )
        boundary_operator_self_test(model)
        initial = initial_condition("left_exit", 64)
        trajectory, _ = learned_rollout(model, initial, periodic=False)
        if not np.isfinite(trajectory).all() or trajectory[..., 0].min() < base.H_FLOOR:
            raise RuntimeError(f"Nonperiodic smoke rollout failed for {arm.name}")
    print("deployment self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--results-dir", type=Path, default=base.HERE / "results")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    train_data, _ = base.load_or_make_data(output, args.seed)
    mean, std, state_std_tensor = base.prepare_statistics(train_data)
    state_std = state_std_tensor.numpy()
    models = {
        arm.name: base.load_model(
            arm,
            base.checkpoint_path(output, arm, args.seed),
            mean,
            std,
            args.width,
        )
        for arm in base.ARMS
    }
    boundary_tests = {
        name: boundary_operator_self_test(model) for name, model in models.items()
    }
    periodic_rows, periodic_trajectories = evaluate_cases(
        PERIODIC_CASES, True, models, state_std
    )
    nonperiodic_rows, nonperiodic_trajectories = evaluate_cases(
        NONPERIODIC_CASES, False, models, state_std
    )
    rows = periodic_rows + nonperiodic_rows
    write_csv(output / f"deployment_metrics_seed{args.seed}.csv", rows)
    for family in FAMILIES:
        plot_profiles(
            "periodic",
            PERIODIC_CASES,
            periodic_trajectories,
            family,
            output / f"periodic512_{family}_s4_vs_s6_seed{args.seed}.png",
        )
        plot_profiles(
            "transmissive nonperiodic",
            NONPERIODIC_CASES,
            nonperiodic_trajectories,
            family,
            output / f"nonperiodic512_{family}_s4_vs_s6_seed{args.seed}.png",
        )
    plot_summary(rows, output / f"deployment_summary_seed{args.seed}.png")
    plot_stability_summary(
        rows, output / f"deployment_stability_seed{args.seed}.png"
    )
    plot_convergence(
        output,
        args.seed,
        output / f"validation_convergence_s4_vs_s6_seed{args.seed}.png",
    )
    summary = {
        "scope": "single-seed controlled 4-vs-6-cell screening",
        "seed": args.seed,
        "training_boundary_condition": "periodic",
        "deployment_boundaries": ["periodic", "transmissive"],
        "learned_boundary_fluxes": 0,
        "final_time": base.DT_SNAPSHOT * (EVAL_SNAPSHOTS - 1),
        "reference": (
            "HLL + SSP-RK2 on 2048 cells, conservatively restricted to 512 "
            "cells for metrics; raw 2048-cell curve retained in figures"
        ),
        "stencils": {
            "4": list(base.stencil_shifts(4)),
            "6": list(base.stencil_shifts(6)),
        },
        "training_wave_coverage": wave_coverage(train_data),
        "boundary_operator_self_tests": boundary_tests,
        "periodic_aggregate": aggregate(periodic_rows),
        "nonperiodic_aggregate": aggregate(nonperiodic_rows),
    }
    (output / f"deployment_summary_seed{args.seed}.json").write_text(
        json.dumps(json_ready(summary), indent=2), encoding="utf-8"
    )
    print(json.dumps(json_ready(summary), indent=2), flush=True)


if __name__ == "__main__":
    main()
