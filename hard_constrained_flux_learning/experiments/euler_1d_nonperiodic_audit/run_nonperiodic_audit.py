"""Zero-shot nonperiodic audit for the retained 1D Euler HCFL methods.

The neural network is evaluated only on interior interfaces.  The two domain
boundary fluxes use constant-extrapolation (transmissive) ghost states, for
which the numerical boundary flux is the physical Euler flux.  No boundary
condition or boundary flux is learned.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
CONVERGENCE_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_convergence_audit"
FULL_FLUX_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_full_flux_baseline"
ROE_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_roe_upwind_ablation"
for module_path in (
    CONVERGENCE_EXPERIMENT,
    FULL_FLUX_EXPERIMENT,
    ROE_EXPERIMENT,
):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import plot_full_flux512 as full_plot  # noqa: E402
import plot_hllc512_vs_hcfl512 as periodic_eval  # noqa: E402


base = periodic_eval.base
shared = periodic_eval.shared
comparison = periodic_eval.comparison
REFERENCE_CELLS = 2048
TARGET_CELLS = 512
SNAPSHOTS = shared.CANONICAL_NSNAP
VARIABLES = ("density", "velocity", "pressure")

CASES: dict[str, dict[str, Any]] = {
    "sod": {
        "display": "Sod",
        "left": (1.0, 0.0, 1.0),
        "right": (0.125, 0.0, 0.1),
        "cut": 0.5,
        "group": "centered",
    },
    "lax": {
        "display": "Lax",
        "left": (0.445, 0.698, 3.528),
        "right": (0.5, 0.0, 0.571),
        "cut": 0.5,
        "group": "centered",
    },
    "collision": {
        "display": "Collision",
        "left": (1.0, 2.0, 1.0),
        "right": (1.0, -2.0, 1.0),
        "cut": 0.5,
        "group": "centered",
    },
    "strong_pressure": {
        "display": "Strong pressure",
        "left": (1.0, 0.0, 5.0),
        "right": (1.0, 0.0, 0.05),
        "cut": 0.5,
        "group": "centered",
    },
    "near_vacuum_expansion": {
        "display": "Near-vacuum",
        "left": (1.0, -2.0, 0.4),
        "right": (1.0, 2.0, 0.4),
        "cut": 0.5,
        "group": "centered",
    },
    "contact_right_exit": {
        "display": "Contact exits right",
        "left": (1.0, 1.0, 1.0),
        "right": (0.25, 1.0, 1.0),
        "cut": 0.985,
        "group": "boundary_interaction",
    },
    "contact_left_exit": {
        "display": "Contact exits left",
        "left": (0.25, -1.0, 1.0),
        "right": (1.0, -1.0, 1.0),
        "cut": 0.015,
        "group": "boundary_interaction",
    },
    "pressure_right_exit": {
        "display": "Pressure wave exits right",
        "left": (1.0, 0.0, 5.0),
        "right": (1.0, 0.0, 0.05),
        "cut": 0.96,
        "group": "boundary_interaction",
    },
}

MODEL_SPECS = {
    "hllc_roe_correction": (
        "dissipation",
        CONVERGENCE_EXPERIMENT
        / "results"
        / "dissipation_broad_converged_best_seed{seed}.pt",
    ),
    "nonnegative_control": (
        "central_roe_upwind",
        ROE_EXPERIMENT
        / "results"
        / "central_roe_upwind_broad_converged_best_seed{seed}.pt",
    ),
    "nonnegative_feasibility": (
        "central_roe_upwind",
        ROE_EXPERIMENT
        / "results"
        / "central_roe_upwind_feas_broad_converged_best_seed{seed}.pt",
    ),
}
METHODS = (
    "native_hllc_512",
    "hllc_roe_correction",
    "nonnegative_control",
    "nonnegative_feasibility",
)
LABELS = {
    "reference": "nonperiodic HLLC-2048 reference",
    "native_hllc_512": "nonperiodic HLLC-512",
    "hllc_roe_correction": "HLLC + Roe correction",
    "nonnegative_control": r"nonnegative Roe, $\lambda_{feas}=0$",
    "nonnegative_feasibility": r"nonnegative Roe, $\lambda_{feas}=10^{-3}$",
}
COLORS = {
    "reference": "#CBD5E1",
    "native_hllc_512": "#0072B2",
    "hllc_roe_correction": "#009E73",
    "nonnegative_control": "#E69F00",
    "nonnegative_feasibility": "#7B2CBF",
}
LINESTYLES = {
    "native_hllc_512": "-",
    "hllc_roe_correction": "-",
    "nonnegative_control": "--",
    "nonnegative_feasibility": "-.",
}


def initial_condition(name: str, cells: int) -> torch.Tensor:
    spec = CASES[name]
    cut = min(max(int(round(float(spec["cut"]) * cells)), 1), cells - 1)
    left = torch.tensor(spec["left"], dtype=torch.float64)
    right = torch.tensor(spec["right"], dtype=torch.float64)
    primitive = torch.empty(1, cells, 3, dtype=torch.float64)
    primitive[:, :cut] = left
    primitive[:, cut:] = right
    return torch.from_numpy(
        base.prim_to_cons(
            primitive[..., 0].numpy(),
            primitive[..., 1].numpy(),
            primitive[..., 2].numpy(),
        )
    )


def update_from_flux(
    state: torch.Tensor,
    flux: torch.Tensor,
    lam: float,
) -> torch.Tensor:
    return state - lam * (flux[:, 1:] - flux[:, :-1])


def pair_rusanov(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    flux_left = base.t_flux(left)
    flux_right = base.t_flux(right)
    rho_left = left[..., 0].clamp_min(1.0e-10)
    rho_right = right[..., 0].clamp_min(1.0e-10)
    velocity_left = left[..., 1] / rho_left
    velocity_right = right[..., 1] / rho_right
    pressure_left = base.t_pressure(left).clamp_min(1.0e-10)
    pressure_right = base.t_pressure(right).clamp_min(1.0e-10)
    sound_left = torch.sqrt(base.GAMMA * pressure_left / rho_left)
    sound_right = torch.sqrt(base.GAMMA * pressure_right / rho_right)
    speed = torch.maximum(
        torch.abs(velocity_left) + sound_left,
        torch.abs(velocity_right) + sound_right,
    )
    return 0.5 * (flux_left + flux_right) - 0.5 * speed[..., None] * (
        right - left
    )


def assemble_flux(
    state: torch.Tensor,
    interior: torch.Tensor,
) -> torch.Tensor:
    left_boundary = base.t_flux(state[:, :1])
    right_boundary = base.t_flux(state[:, -1:])
    return torch.cat([left_boundary, interior, right_boundary], dim=1)


def native_hllc_flux(state: torch.Tensor) -> torch.Tensor:
    return assemble_flux(
        state,
        base.hllc_pair(state[:, :-1], state[:, 1:]),
    )


def low_order_flux(state: torch.Tensor) -> torch.Tensor:
    return assemble_flux(
        state,
        pair_rusanov(state[:, :-1], state[:, 1:]),
    )


def entropy_flux(state: torch.Tensor) -> torch.Tensor:
    density = state[..., 0].clamp_min(1.0e-12)
    velocity = state[..., 1] / density
    return velocity * base.entropy(state)


def entropy_target(state: torch.Tensor, lam: float) -> torch.Tensor:
    boundary_change = entropy_flux(state[:, -1]) - entropy_flux(state[:, 0])
    return base.entropy(state.double()).sum(dim=-1) - lam * boundary_change.double()


def admissible(state: torch.Tensor) -> torch.Tensor:
    return (state[..., 0] >= shared.RHO_FLOOR) & (
        base.t_pressure(state) >= shared.PRESSURE_FLOOR
    )


def assert_strict_admissible(state: torch.Tensor, context: str) -> None:
    if (
        not bool(torch.isfinite(state).all())
        or not bool((state[..., 0] > 0.0).all())
        or not bool((base.t_pressure(state) > 0.0).all())
    ):
        raise RuntimeError(f"Nonperiodic HLLC lost admissibility: {context}")


def project_entropy_pair(
    proposal: torch.Tensor,
    left: torch.Tensor,
    right: torch.Tensor,
) -> torch.Tensor:
    """Roundoff-safe Tadmor projection on interior interfaces only."""
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
        correction = torch.zeros_like(residual)
        correction[mask] = (
            residual[mask] - target[mask]
        ) / norm_squared[mask]
        output = (flux64 - correction[..., None] * normal).to(proposal.dtype)
    return output


def model_flux(
    model: base.Solver,
    state: torch.Tensor,
) -> tuple[torch.Tensor, dict[str, float]]:
    """Return boundary-classical/interior-neural fluxes without circular wrap."""
    left_ghost = state[:, :1].expand(-1, 2, -1)
    right_ghost = state[:, -1:].expand(-1, 2, -1)
    padded = torch.cat([left_ghost, state, right_ghost], dim=1)
    raw_padded = model.flux_net(padded)
    cells = state.shape[1]
    raw_interior = raw_padded[:, 2 : cells + 1]
    projected_interior = project_entropy_pair(
        raw_interior,
        state[:, :-1],
        state[:, 1:],
    )
    delta = projected_interior - raw_interior
    tolerance = 1.0e-7 * (1.0 + raw_interior.abs().amax(dim=-1))
    stats = {
        "projection_active": float(
            (
                torch.linalg.vector_norm(delta, dim=-1) > tolerance
            ).sum()
        ),
        "projection_total": float(delta.shape[0] * delta.shape[1]),
        "projection_squared": float((delta.double() ** 2).sum()),
        "raw_squared": float((raw_interior.double() ** 2).sum()),
        "projection_maximum": float(delta.abs().max()),
    }
    return assemble_flux(state, projected_interior), stats


@torch.no_grad()
def boundary_operator_self_test(model: base.Solver) -> dict[str, float | str]:
    """Verify consistency, no circular leakage, and flux-form conservation."""
    cells = 32
    primitive = np.empty((1, cells, 3), dtype=np.float32)
    primitive[..., 0] = 1.1
    primitive[..., 1] = 0.2
    primitive[..., 2] = 0.9
    constant = torch.from_numpy(
        base.prim_to_cons(
            primitive[..., 0], primitive[..., 1], primitive[..., 2]
        ).astype(np.float32)
    )
    flux, _ = model_flux(model, constant)
    physical = base.t_flux(constant)
    consistency_error = float(
        (flux - physical[:, :1].expand_as(flux)).abs().max()
    )
    constant_update_error = float(
        (update_from_flux(constant, flux, 0.0256) - constant).abs().max()
    )

    probe = initial_condition("sod", cells).float()
    changed_far_right = probe.clone()
    changed_far_right[:, -1] = torch.from_numpy(
        base.prim_to_cons(
            np.array([0.7]), np.array([-0.4]), np.array([1.3])
        ).astype(np.float32)
    )[0]
    probe_flux, _ = model_flux(model, probe)
    changed_flux, _ = model_flux(model, changed_far_right)
    left_interface_wrap_leak = float(
        (probe_flux[:, 1] - changed_flux[:, 1]).abs().max()
    )
    boundary_flux_error = max(
        float((probe_flux[:, 0] - base.t_flux(probe[:, 0])).abs().max()),
        float((probe_flux[:, -1] - base.t_flux(probe[:, -1])).abs().max()),
    )

    lam = 0.01
    next_state = update_from_flux(probe, probe_flux, lam)
    conservation_residual = (
        (next_state - probe).sum(dim=1)
        + lam * (probe_flux[:, -1] - probe_flux[:, 0])
    )
    conservation_error = float(conservation_residual.abs().max())
    tolerances = {
        "equal_state_consistency_max_error": 2.0e-6,
        "constant_update_max_error": 2.0e-6,
        "left_interface_wrap_leak_max_error": 1.0e-7,
        "boundary_flux_max_error": 1.0e-7,
        "flux_form_conservation_max_error": 2.0e-6,
    }
    measured = {
        "equal_state_consistency_max_error": consistency_error,
        "constant_update_max_error": constant_update_error,
        "left_interface_wrap_leak_max_error": left_interface_wrap_leak,
        "boundary_flux_max_error": boundary_flux_error,
        "flux_form_conservation_max_error": conservation_error,
    }
    if any(measured[key] > limit for key, limit in tolerances.items()):
        raise RuntimeError(
            f"Nonperiodic boundary-operator self-test failed: {measured}"
        )
    return {"status": "pass", **measured}


def local_admissibility_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
    max_outer: int = 12,
    n_bisect: int = 30,
) -> tuple[torch.Tensor, torch.Tensor]:
    delta_flux = high_flux - low_flux
    low_state = update_from_flux(state, low_flux, lam)
    if not bool(admissible(low_state).all()):
        raise RuntimeError("Nonperiodic low-order state is not admissible")

    batch, cells, _ = state.shape
    alpha = torch.ones(batch, cells + 1, dtype=state.dtype)
    for _ in range(max_outer):
        flux = low_flux + alpha[..., None] * delta_flux
        current = update_from_flux(state, flux, lam)
        good = admissible(current)
        if bool(good.all()):
            return flux, alpha

        bad = ~good
        direction = current - low_state
        count = int(bad.sum())
        lo = torch.zeros(count, dtype=torch.float64)
        hi = torch.ones(count, dtype=torch.float64)
        origin = low_state[bad].double()
        vector = direction[bad].double()
        for _ in range(n_bisect):
            midpoint = 0.5 * (lo + hi)
            okay = admissible(origin + midpoint[:, None] * vector)
            lo = torch.where(okay, midpoint, lo)
            hi = torch.where(okay, hi, midpoint)

        theta = torch.ones(batch, cells, dtype=state.dtype)
        theta[bad] = torch.clamp(lo.to(state.dtype) - 1.0e-5, min=0.0)
        interface_factor = torch.ones_like(alpha)
        interface_factor[:, 1:cells] = torch.minimum(
            theta[:, :-1], theta[:, 1:]
        )
        alpha = alpha * interface_factor

    flux = low_flux + alpha[..., None] * delta_flux
    current = update_from_flux(state, flux, lam)
    failed = (~admissible(current)).any(dim=-1)
    if bool(failed.any()):
        alpha[failed] = 0.0
        flux = low_flux + alpha[..., None] * delta_flux
    if not bool(admissible(update_from_flux(state, flux, lam)).all()):
        raise RuntimeError("Nonperiodic admissibility limiter failed")
    return flux, alpha


def global_entropy_limiter(
    state: torch.Tensor,
    high_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
    target: torch.Tensor,
    n_bisect: int = 40,
) -> tuple[torch.Tensor, torch.Tensor]:
    low_state = update_from_flux(state, low_flux, lam)
    high_state = update_from_flux(state, high_flux, lam)
    high_entropy = base.entropy(high_state.double()).sum(dim=-1)
    beta = torch.ones(state.shape[0], dtype=state.dtype)
    need = high_entropy > target

    if bool(need.any()):
        count = int(need.sum())
        lo = torch.zeros(count, dtype=torch.float64)
        hi = torch.ones(count, dtype=torch.float64)
        origin = low_state[need].double()
        direction = (high_state - low_state)[need].double()
        target_need = target[need]
        for _ in range(n_bisect):
            midpoint = 0.5 * (lo + hi)
            candidate = origin + midpoint[:, None, None] * direction
            okay = base.entropy(candidate).sum(dim=-1) <= target_need
            lo = torch.where(okay, midpoint, lo)
            hi = torch.where(okay, hi, midpoint)
        beta[need] = torch.clamp(lo.to(state.dtype) - 1.0e-6, min=0.0)

    delta_flux = high_flux - low_flux
    flux = low_flux + beta[:, None, None] * delta_flux
    final_state = update_from_flux(state, flux, lam)
    still_bad = (
        base.entropy(final_state.double()).sum(dim=-1)
        > target + shared.FD_ENTROPY_TOLERANCE
    )
    if bool(still_bad.any()):
        beta[still_bad] = 0.0
        flux = low_flux + beta[:, None, None] * delta_flux
    return flux, beta


@torch.no_grad()
def strict_hllc_rollout(
    name: str,
    cells: int,
    cfl: float = 0.25,
) -> tuple[torch.Tensor, dict[str, Any]]:
    state = initial_condition(name, cells)
    dx = 1.0 / cells
    snapshots = [state.clone()]
    substeps = 0
    for snapshot in range(1, SNAPSHOTS):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            density = state[..., 0]
            pressure = base.t_pressure(state)
            velocity = state[..., 1] / density
            sound = torch.sqrt(base.GAMMA * pressure / density)
            max_speed = float((torch.abs(velocity) + sound).max())
            dt = min(remaining, cfl * dx / max(max_speed, 1.0e-12))
            lam = dt / dx
            stage_one = update_from_flux(state, native_hllc_flux(state), lam)
            assert_strict_admissible(
                stage_one, f"{name}, snapshot {snapshot}, stage one"
            )
            stage_two = update_from_flux(
                stage_one, native_hllc_flux(stage_one), lam
            )
            state = 0.5 * state + 0.5 * stage_two
            assert_strict_admissible(
                state, f"{name}, snapshot {snapshot}, completion"
            )
            remaining -= dt
            substeps += 1
        snapshots.append(state.clone())
    return torch.stack(snapshots, dim=1).float(), {
        "completed": True,
        "ssprk2_substeps": substeps,
        "flux_evaluations": 2 * substeps,
    }


@torch.no_grad()
def learned_rollout(
    model: base.Solver,
    name: str,
    cells: int,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, Any]]:
    state = initial_condition(name, cells).float()
    dx = 1.0 / cells
    nominal_dt = base.DT_SNAPSHOT * base.NCOARSE / cells
    snapshots = [state.clone()]
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["max_entropy_balance_violation"] = -float("inf")
    totals["max_interior_entropy_residual"] = -float("inf")
    totals["projection_maximum"] = 0.0

    for _ in range(1, SNAPSHOTS):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            density = state[..., 0].clamp_min(1.0e-10)
            pressure = base.t_pressure(state).clamp_min(1.0e-10)
            velocity = state[..., 1] / density
            sound = torch.sqrt(base.GAMMA * pressure / density)
            max_speed = float((torch.abs(velocity) + sound).max())
            dt = min(
                remaining,
                nominal_dt,
                cfl * dx / max(max_speed, 1.0e-12),
            )

            for _ in range(24):
                lam = dt / dx
                low_flux = low_order_flux(state)
                low_state = update_from_flux(state, low_flux, lam)
                target = entropy_target(state, lam)
                entropy_ok = bool(
                    (
                        base.entropy(low_state.double()).sum(dim=-1)
                        <= target + 1.0e-10
                    ).all()
                )
                if bool(admissible(low_state).all()) and entropy_ok:
                    break
                dt *= 0.5
                totals["low_order_dt_halvings"] += 1.0
            else:
                raise RuntimeError(
                    "Could not establish nonperiodic low-order safety premise"
                )

            high_flux, projection = model_flux(model, state)
            local_flux, alpha = local_admissibility_limiter(
                state, high_flux, low_flux, lam
            )
            final_flux, beta = global_entropy_limiter(
                state, local_flux, low_flux, lam, target
            )
            next_state = update_from_flux(state, final_flux, lam)
            if (
                not bool(torch.isfinite(next_state).all())
                or not bool(admissible(next_state).all())
            ):
                raise RuntimeError("Nonperiodic learned rollout became inadmissible")

            left = state[:, :-1].double()
            right = state[:, 1:].double()
            interior_flux = final_flux[:, 1:-1].double()
            normal = base.entropy_variables(right) - base.entropy_variables(left)
            bound = base.entropy_potential(right) - base.entropy_potential(left)
            residual = (normal * interior_flux).sum(dim=-1) - bound
            balance_violation = (
                base.entropy(next_state.double()).sum(dim=-1) - target
            )

            for key, value in projection.items():
                if key == "projection_maximum":
                    totals[key] = max(totals[key], value)
                else:
                    totals[key] += value
            interior_alpha = alpha[:, 1:-1]
            totals["local_active"] += float(
                (interior_alpha < 1.0 - 1.0e-7).sum()
            )
            totals["local_total"] += float(interior_alpha.numel())
            totals["fd_active"] += float((beta < 1.0 - 1.0e-7).sum())
            totals["fd_total"] += float(beta.numel())
            totals["fd_beta_sum"] += float(beta.sum())
            totals["fd_beta_min"] = min(totals["fd_beta_min"], float(beta.min()))
            totals["max_entropy_balance_violation"] = max(
                totals["max_entropy_balance_violation"],
                float(balance_violation.max()),
            )
            totals["max_interior_entropy_residual"] = max(
                totals["max_interior_entropy_residual"], float(residual.max())
            )
            totals["substeps"] += 1.0
            state = next_state
            remaining -= dt
        snapshots.append(state.clone())

    trajectory = torch.stack(snapshots, dim=1)
    return trajectory, {
        "completed": True,
        "internal_substeps": int(totals["substeps"]),
        "hard_projection_intervention_rate": totals["projection_active"]
        / max(totals["projection_total"], 1.0),
        "hard_projection_relative_flux_rms": float(
            np.sqrt(
                totals["projection_squared"]
                / max(totals["raw_squared"], 1.0e-30)
            )
        ),
        "hard_projection_max_absolute_change": totals["projection_maximum"],
        "local_limiter_intervention_rate": totals["local_active"]
        / max(totals["local_total"], 1.0),
        "fd_entropy_intervention_rate": totals["fd_active"]
        / max(totals["fd_total"], 1.0),
        "mean_fd_beta": totals["fd_beta_sum"]
        / max(totals["fd_total"], 1.0),
        "min_fd_beta": totals["fd_beta_min"],
        "max_boundary_aware_entropy_balance_violation": totals[
            "max_entropy_balance_violation"
        ],
        "max_interior_tadmor_residual": totals[
            "max_interior_entropy_residual"
        ],
        "low_order_dt_halvings": int(totals["low_order_dt_halvings"]),
    }


def restrict_reference(trajectory: torch.Tensor, cells: int) -> torch.Tensor:
    factor = trajectory.shape[-2] // cells
    shape = trajectory.shape
    return trajectory.reshape(
        shape[0], shape[1], cells, factor, shape[-1]
    ).mean(dim=3)


def integrity_metrics(trajectory: torch.Tensor) -> dict[str, float]:
    return {
        "minimum_density": float(trajectory[..., 0].min()),
        "minimum_pressure": float(base.t_pressure(trajectory).min()),
    }


def total_variation(values: np.ndarray) -> float:
    """Non-circular total variation; the endpoints are not adjacent."""
    return float(np.abs(np.diff(values)).sum())


def significant_extrema(values: np.ndarray, scale: float) -> int:
    left = values[1:-1] - values[:-2]
    right = values[2:] - values[1:-1]
    threshold = 1.0e-3 * max(scale, 1.0e-8)
    return int(
        (
            (left * right < 0.0)
            & (np.minimum(np.abs(left), np.abs(right)) > threshold)
        ).sum()
    )


def oscillation_metrics(
    reference: torch.Tensor,
    candidate: torch.Tensor,
) -> dict[str, float | int]:
    truth = base.primitive(reference)[0, -1].numpy()
    prediction = base.primitive(candidate)[0, -1].numpy()
    normalized_tv_excesses: list[float] = []
    range_violations: list[float] = []
    excess_extrema = 0
    for index in range(3):
        exact = truth[:, index]
        estimated = prediction[:, index]
        scale = max(
            float(np.ptp(exact)),
            float(np.max(np.abs(exact))),
            1.0,
        )
        normalized_tv_excesses.append(
            max(total_variation(estimated) - total_variation(exact), 0.0)
            / scale
        )
        range_violations.append(
            (
                max(float(exact.min() - estimated.min()), 0.0)
                + max(float(estimated.max() - exact.max()), 0.0)
            )
            / scale
        )
        excess_extrema += max(
            significant_extrema(estimated, scale)
            - significant_extrema(exact, scale),
            0,
        )
    return {
        "mean_normalized_final_tv_excess": float(
            np.mean(normalized_tv_excesses)
        ),
        "mean_normalized_global_range_violation": float(
            np.mean(range_violations)
        ),
        "total_excess_significant_extrema": excess_extrema,
    }


def aggregate(
    cases: dict[str, dict[str, Any]],
    method: str,
    group: str | None = None,
) -> dict[str, Any]:
    names = [
        name
        for name in cases
        if group is None or CASES[name]["group"] == group
    ]
    rows = [cases[name][method] for name in names]
    result = {
        "case_count": len(rows),
        "mean_rollout_nrmse": float(
            np.mean([row["rollout_nrmse"] for row in rows])
        ),
        "mean_final_snapshot_nrmse": float(
            np.mean([row["final_snapshot_nrmse"] for row in rows])
        ),
        "mean_normalized_final_tv_excess": float(
            np.mean([
                row["mean_normalized_final_tv_excess"] for row in rows
            ])
        ),
        "mean_normalized_global_range_violation": float(
            np.mean([
                row["mean_normalized_global_range_violation"] for row in rows
            ])
        ),
        "total_excess_significant_extrema": int(
            sum(row["total_excess_significant_extrema"] for row in rows)
        ),
        "minimum_density": float(min(row["minimum_density"] for row in rows)),
        "minimum_pressure": float(
            min(row["minimum_pressure"] for row in rows)
        ),
    }
    for key in (
        "hard_projection_intervention_rate",
        "local_limiter_intervention_rate",
        "fd_entropy_intervention_rate",
    ):
        if all(key in row for row in rows):
            result[f"mean_{key}"] = float(
                np.mean([row[key] for row in rows])
            )
    for key in (
        "max_boundary_aware_entropy_balance_violation",
        "max_interior_tadmor_residual",
    ):
        if all(key in row for row in rows):
            result[key] = float(max(row[key] for row in rows))
    return result


def plot_profiles(
    trajectories: dict[str, dict[str, torch.Tensor]],
    output: Path,
) -> None:
    names = list(CASES)
    x_reference = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    x_target = (np.arange(TARGET_CELLS) + 0.5) / TARGET_CELLS
    figure, axes = plt.subplots(
        3, len(names), figsize=(23, 8.8), sharex=True
    )
    figure.subplots_adjust(
        left=0.045, right=0.995, bottom=0.09, top=0.79,
        wspace=0.23, hspace=0.14,
    )
    for column, name in enumerate(names):
        primitive = {
            method: base.primitive(value)[0, -1].numpy()
            for method, value in trajectories[name].items()
        }
        axes[0, column].set_title(
            CASES[name]["display"], fontsize=9.2, fontweight="semibold"
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                primitive["reference"][:, row],
                color=COLORS["reference"],
                linewidth=2.0,
                alpha=0.9,
                zorder=1,
            )
            for zorder, method in enumerate(METHODS, start=2):
                axis.plot(
                    x_target,
                    primitive[method][:, row],
                    color=COLORS[method],
                    linestyle=LINESTYLES[method],
                    linewidth=(
                        1.5 if method == "nonnegative_feasibility" else 1.15
                    ),
                    alpha=0.94,
                    zorder=zorder,
                )
            displayed = np.concatenate([
                primitive[method][:, row]
                for method in ("reference", *METHODS)
            ])
            data_min = float(displayed.min())
            data_max = float(displayed.max())
            reference_scale = max(
                float(np.max(np.abs(primitive["reference"][:, row]))),
                1.0,
            )
            span = max(data_max - data_min, 0.05 * reference_scale)
            center = 0.5 * (data_min + data_max)
            axis.set_ylim(center - 0.55 * span, center + 0.55 * span)
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(VARIABLES[row], fontsize=9.5)
            if row == 2:
                axis.set_xlabel("x", fontsize=8.5)
    handles = [
        Line2D(
            [0], [0], color=COLORS["reference"], linewidth=2.0,
            label=LABELS["reference"],
        )
    ]
    handles.extend(
        Line2D(
            [0], [0], color=COLORS[method],
            linestyle=LINESTYLES[method], linewidth=1.4,
            label=LABELS[method],
        )
        for method in METHODS
    )
    final_time = (SNAPSHOTS - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        "Zero-shot nonperiodic Euler: transmissive boundaries, "
        f"t = {final_time:.4f}",
        y=0.975,
        fontsize=16,
        fontweight="bold",
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.91),
        ncol=5,
        frameon=False,
        fontsize=8.4,
    )
    figure.text(
        0.995,
        0.015,
        "NN flux is used only on interior interfaces. Boundary fluxes are "
        "physical Euler fluxes from constant-extrapolation ghost states.",
        ha="right",
        fontsize=8.2,
        color="#64748B",
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    for name in cases:
        for method in METHODS:
            row = cases[name][method]
            rows.append(
                {
                    "case": name,
                    "group": CASES[name]["group"],
                    "method": method,
                    "rollout_nrmse": row["rollout_nrmse"],
                    "final_snapshot_nrmse": row["final_snapshot_nrmse"],
                    "mean_normalized_final_tv_excess": row[
                        "mean_normalized_final_tv_excess"
                    ],
                    "mean_normalized_global_range_violation": row[
                        "mean_normalized_global_range_violation"
                    ],
                    "total_excess_significant_extrema": row[
                        "total_excess_significant_extrema"
                    ],
                    "minimum_density": row["minimum_density"],
                    "minimum_pressure": row["minimum_pressure"],
                    "hard_projection_intervention_rate": row.get(
                        "hard_projection_intervention_rate"
                    ),
                    "max_boundary_aware_entropy_balance_violation": row.get(
                        "max_boundary_aware_entropy_balance_violation"
                    ),
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    _, _, state_std = periodic_eval.training_statistics(args.seed)
    models = {
        method: full_plot.load_solver(
            model_name,
            checkpoint.with_name(checkpoint.name.format(seed=args.seed)),
            args.width,
        )
        for method, (model_name, checkpoint) in MODEL_SPECS.items()
    }
    operator_self_tests = {
        method: boundary_operator_self_test(model)
        for method, model in models.items()
    }

    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    cases: dict[str, dict[str, Any]] = {}
    for name in CASES:
        print(f"Evaluating {CASES[name]['display']}...", flush=True)
        reference, reference_stats = strict_hllc_rollout(name, REFERENCE_CELLS)
        scoring_reference = restrict_reference(reference, TARGET_CELLS)
        hllc_512, hllc_stats = strict_hllc_rollout(name, TARGET_CELLS)
        trajectories[name] = {
            "reference": reference,
            "native_hllc_512": hllc_512,
        }
        case: dict[str, Any] = {
            "reference_hllc_2048": reference_stats,
            "native_hllc_512": hllc_stats,
        }
        for method, model in models.items():
            trajectory, stats = learned_rollout(model, name, TARGET_CELLS)
            trajectories[name][method] = trajectory
            case[method] = stats

        for method in METHODS:
            trajectory = trajectories[name][method]
            case[method].update(
                comparison.diagnostics(
                    scoring_reference, trajectory, state_std
                )
            )
            case[method].update(integrity_metrics(trajectory))
            case[method].update(
                oscillation_metrics(scoring_reference, trajectory)
            )
        cases[name] = case

    aggregates = {
        group: {
            method: aggregate(
                cases, method, None if group == "all" else group
            )
            for method in METHODS
        }
        for group in ("all", "centered", "boundary_interaction")
    }
    summary = {
        "seed": args.seed,
        "training_boundary_condition": "periodic",
        "deployment_boundary_condition": (
            "transmissive constant extrapolation; physical flux at both "
            "domain boundaries"
        ),
        "learned_boundary_interfaces": 0,
        "learned_interior_interfaces": TARGET_CELLS - 1,
        "reference": (
            "strict nonperiodic native HLLC + SSP-RK2 on 2048 cells; "
            "conservatively restricted to 512 cells only for metrics"
        ),
        "snapshots": SNAPSHOTS,
        "snapshot_dt": base.DT_SNAPSHOT,
        "final_time": (SNAPSHOTS - 1) * base.DT_SNAPSHOT,
        "boundary_operator_self_tests": operator_self_tests,
        "cases": cases,
        "aggregates": aggregates,
    }
    stem = f"nonperiodic_transmissive512_seed{args.seed}"
    plot_profiles(trajectories, output / f"{stem}.png")
    (output / f"{stem}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(output / f"{stem}.csv", cases)
    print(json.dumps(aggregates, indent=2))


if __name__ == "__main__":
    main()
