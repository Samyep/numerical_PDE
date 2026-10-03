"""Matched 64-cell benchmark for Euler, SWE, HCFL, FVM, MUSCL, and FNO.

Every candidate starts from exactly the same 64 finite-volume cell averages.
The scoring reference is obtained by piecewise-constant prolongation of those
same averages to a fine grid, native fine-grid evolution, and conservative
restriction back to 64 cells.  Thus no method receives hidden subcell initial
condition information.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
EULER_DIR = EXPERIMENTS / "euler_1d_stencil_ablation"
SWE_DIR = EXPERIMENTS / "swe_1d_consistent_hcfl"
SWE_DEPLOY_DIR = SWE_DIR
for module_dir in (HERE, EULER_DIR, SWE_DIR):
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import run_stencil_ablation as euler  # noqa: E402
import run_swe_consistent as swe  # noqa: E402
import evaluate_deployment as swe_deploy  # noqa: E402
import train_hcfl_replicates as replicates  # noqa: E402
import train_operator_baselines as operator  # noqa: E402
import train_learned_flux_baselines as learned_flux  # noqa: E402


CELLS = 64
EULER_RANDOM_SPLITS = ("ordinary", "broad", "extreme")
EULER_CANONICAL = (
    "sod",
    "lax",
    "collision",
    "strong_pressure",
    "near_vacuum_expansion",
)
SWE_RANDOM_SPLITS = ("ordinary", "broad", "froude")
SWE_CANONICAL = tuple(swe_deploy.PERIODIC_CASES)
METHOD_LABELS = {
    "native_fvm64": "native first-order FVM",
    "muscl_fvm64": "TVD MUSCL FVM",
    "hcfl64": "HCFL",
    "hcfl64_no_low": "HCFL without low-order blending",
    "learned_flux64": "learned neural FV flux",
    "fno64": "FNO baseline",
}
METHOD_COLORS = {
    "native_fvm64": "#0072B2",
    "muscl_fvm64": "#009E73",
    "hcfl64": "#D55E00",
    "hcfl64_no_low": "#E69F00",
    "learned_flux64": "#56B4E9",
    "fno64": "#CC79A7",
}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def restrict(values: torch.Tensor, target_cells: int = CELLS) -> torch.Tensor:
    source_cells = values.shape[-2]
    if source_cells % target_cells:
        raise ValueError("Source cells must be divisible by target cells")
    factor = source_cells // target_cells
    return values.reshape(
        *values.shape[:-2], target_cells, factor, values.shape[-1]
    ).mean(dim=-2)


def minmod_three(
    first: torch.Tensor,
    second: torch.Tensor,
    third: torch.Tensor,
) -> torch.Tensor:
    same = (
        (torch.sign(first) == torch.sign(second))
        & (torch.sign(second) == torch.sign(third))
    )
    magnitude = torch.minimum(
        torch.minimum(first.abs(), second.abs()), third.abs()
    )
    return torch.where(same, torch.sign(first) * magnitude, torch.zeros_like(first))


def limited_slope(values: torch.Tensor, theta: float = 1.5) -> torch.Tensor:
    backward = values - torch.roll(values, 1, dims=-2)
    forward = torch.roll(values, -1, dims=-2) - values
    centered = 0.5 * (torch.roll(values, -1, dims=-2) - torch.roll(values, 1, dims=-2))
    return minmod_three(theta * backward, centered, theta * forward)


def euler_primitive_to_conservative(values: torch.Tensor) -> torch.Tensor:
    rho = values[..., 0]
    velocity = values[..., 1]
    pressure = values[..., 2]
    return torch.stack(
        [
            rho,
            rho * velocity,
            pressure / (euler.base.GAMMA - 1.0)
            + 0.5 * rho * velocity.square(),
        ],
        dim=-1,
    )


def euler_muscl_flux(state: torch.Tensor) -> torch.Tensor:
    values = euler.base.primitive(state)
    slope = limited_slope(values)
    left_values = values + 0.5 * slope
    right_values = torch.roll(values, -1, dims=-2) - 0.5 * torch.roll(
        slope, -1, dims=-2
    )
    # A componentwise TVD limiter does not, by itself, guarantee positive
    # reconstructed density and pressure.  Revert only the affected interface
    # to its two adjacent cell averages, which is the standard robust fallback
    # rather than letting HLLC receive an unphysical state.
    reconstruction_ok = (
        torch.isfinite(left_values).all(dim=-1)
        & torch.isfinite(right_values).all(dim=-1)
        & (left_values[..., 0] > 0.0)
        & (right_values[..., 0] > 0.0)
        & (left_values[..., 2] > 0.0)
        & (right_values[..., 2] > 0.0)
    )
    left_values = torch.where(
        reconstruction_ok[..., None], left_values, values
    )
    right_values = torch.where(
        reconstruction_ok[..., None],
        right_values,
        torch.roll(values, -1, dims=-2),
    )
    left = euler_primitive_to_conservative(left_values)
    right = euler_primitive_to_conservative(right_values)
    return euler.base.hllc_pair(left, right)


def swe_primitive_to_conservative(values: torch.Tensor) -> torch.Tensor:
    return torch.stack(
        [values[..., 0], values[..., 0] * values[..., 1]], dim=-1
    )


def swe_hll_pair(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    flux_left = swe.t_flux(left)
    flux_right = swe.t_flux(right)
    h_left = left[..., 0]
    h_right = right[..., 0]
    u_left = left[..., 1] / h_left
    u_right = right[..., 1] / h_right
    c_left = torch.sqrt(swe.G * h_left)
    c_right = torch.sqrt(swe.G * h_right)
    speed_left = torch.minimum(u_left - c_left, u_right - c_right)
    speed_right = torch.maximum(u_left + c_left, u_right + c_right)
    middle = (
        speed_right[..., None] * flux_left
        - speed_left[..., None] * flux_right
        + (speed_left * speed_right)[..., None] * (right - left)
    ) / (speed_right - speed_left).clamp_min(1.0e-14)[..., None]
    return torch.where(
        (speed_left >= 0.0)[..., None],
        flux_left,
        torch.where((speed_right <= 0.0)[..., None], flux_right, middle),
    )


def swe_first_order_flux(state: torch.Tensor) -> torch.Tensor:
    return swe_hll_pair(state, torch.roll(state, -1, dims=-2))


def swe_muscl_flux(state: torch.Tensor) -> torch.Tensor:
    values = swe.primitive(state)
    slope = limited_slope(values)
    left_values = values + 0.5 * slope
    right_values = torch.roll(values, -1, dims=-2) - 0.5 * torch.roll(
        slope, -1, dims=-2
    )
    reconstruction_ok = (
        torch.isfinite(left_values).all(dim=-1)
        & torch.isfinite(right_values).all(dim=-1)
        & (left_values[..., 0] > 0.0)
        & (right_values[..., 0] > 0.0)
    )
    left_values = torch.where(
        reconstruction_ok[..., None], left_values, values
    )
    right_values = torch.where(
        reconstruction_ok[..., None],
        right_values,
        torch.roll(values, -1, dims=-2),
    )
    left = swe_primitive_to_conservative(left_values)
    right = swe_primitive_to_conservative(right_values)
    return swe_hll_pair(left, right)


def is_admissible(system: str, state: torch.Tensor) -> torch.Tensor:
    finite = torch.isfinite(state).all(dim=-1)
    if system == "euler":
        density = state[..., 0]
        pressure = (euler.base.GAMMA - 1.0) * (
            state[..., 2]
            - 0.5 * state[..., 1].square() / density.clamp_min(1.0e-300)
        )
        return finite & (density > 0.0) & (pressure > 0.0)
    if system == "swe":
        return finite & (state[..., 0] > 0.0)
    raise ValueError(system)


def maximum_speed(system: str, state: torch.Tensor) -> float:
    if system == "euler":
        rho = state[..., 0]
        velocity = state[..., 1] / rho
        pressure = (euler.base.GAMMA - 1.0) * (
            state[..., 2] - 0.5 * state[..., 1].square() / rho
        )
        speed = velocity.abs() + torch.sqrt(euler.base.GAMMA * pressure / rho)
    else:
        h = state[..., 0]
        speed = (state[..., 1] / h).abs() + torch.sqrt(swe.G * h)
    return float(speed.max())


def spatial_rhs(
    state: torch.Tensor,
    dx: float,
    flux_function: Callable[[torch.Tensor], torch.Tensor],
) -> torch.Tensor:
    flux = flux_function(state)
    return -(flux - torch.roll(flux, 1, dims=-2)) / dx


@torch.no_grad()
def ssprk2_rollout(
    system: str,
    initial: torch.Tensor,
    snapshots: int,
    dt_snapshot: float,
    flux_function: Callable[[torch.Tensor], torch.Tensor],
    cfl: float,
) -> tuple[torch.Tensor, dict[str, Any]]:
    state = initial.double().clone()
    cells = state.shape[-2]
    dx = 1.0 / cells
    saved = [state.float().clone()]
    substeps = 0
    halvings = 0
    started = time.perf_counter()
    for _ in range(1, snapshots):
        remaining = dt_snapshot
        while remaining > 1.0e-14:
            dt = min(
                remaining,
                cfl * dx / max(maximum_speed(system, state), 1.0e-12),
            )
            accepted = False
            for _ in range(24):
                first = state + dt * spatial_rhs(state, dx, flux_function)
                if bool(is_admissible(system, first).all()):
                    second = first + dt * spatial_rhs(first, dx, flux_function)
                    candidate = 0.5 * state + 0.5 * second
                    if bool(is_admissible(system, candidate).all()):
                        accepted = True
                        break
                dt *= 0.5
                halvings += 1
            if not accepted:
                raise RuntimeError(
                    f"{system} SSP-RK2 could not find an admissible step"
                )
            state = candidate
            remaining -= dt
            substeps += 1
        saved.append(state.float().clone())
    return torch.stack(saved, dim=1), {
        "substeps": substeps,
        "admissibility_dt_halvings": halvings,
        "seconds": time.perf_counter() - started,
    }


def reference_rollout(
    system: str,
    initial64: torch.Tensor,
    snapshots: int,
    reference_cells: int,
) -> tuple[torch.Tensor, dict[str, Any]]:
    if reference_cells % CELLS:
        raise ValueError("reference_cells must be divisible by 64")
    fine = initial64.double().repeat_interleave(
        reference_cells // CELLS, dim=-2
    )
    if system == "euler":
        trajectory, stats = ssprk2_rollout(
            system,
            fine,
            snapshots,
            euler.base.DT_SNAPSHOT,
            euler.base.t_hllc,
            0.20,
        )
    else:
        trajectory, stats = ssprk2_rollout(
            system,
            fine,
            snapshots,
            swe.DT_SNAPSHOT,
            swe_first_order_flux,
            0.20,
        )
    return restrict(trajectory), stats


def euler_initial_suites(seed: int, count: int) -> dict[str, torch.Tensor]:
    suites = {
        "ordinary": torch.from_numpy(
            euler.base.generate_ic(count, CELLS, 91000 + seed, ood=False)
        ).float(),
        "broad": torch.from_numpy(
            euler.base.generate_ic(count, CELLS, 92000 + seed, ood=True)
        ).float(),
        "extreme": torch.from_numpy(
            euler.base.generate_extreme_ic(count, CELLS, 93000 + seed)
        ).float(),
    }
    canonical = euler.shared.canonical_initial_conditions()
    for name in EULER_CANONICAL:
        value = torch.from_numpy(canonical[name]).float()
        suites[name] = restrict(value)
    return suites


def swe_initial_suites(seed: int, count: int) -> dict[str, torch.Tensor]:
    suites = {
        "ordinary": torch.from_numpy(
            swe.legacy.generate_ic(count, CELLS, 94000 + seed, ood=False)
        ).float(),
        "broad": torch.from_numpy(
            swe.legacy.generate_ic(count, CELLS, 95000 + seed, ood=True)
        ).float(),
        "froude": torch.from_numpy(
            swe.generate_froude_coverage_ic(count, CELLS, 96000 + seed)
        ).float(),
    }
    for name in SWE_CANONICAL:
        suites[name] = torch.from_numpy(
            swe_deploy.initial_condition(name, CELLS)
        )[None].float()
    return suites


def suite_snapshots(system: str, split: str) -> int:
    if system == "euler":
        return euler.shared.CANONICAL_NSNAP if split in EULER_CANONICAL else euler.base.NSNAP
    return swe_deploy.EVAL_SNAPSHOTS if split in SWE_CANONICAL else swe.NSNAP


def load_or_make_references(
    output: Path,
    seed: int,
    count: int,
    reference_cells: int,
) -> tuple[dict[str, dict[str, torch.Tensor]], dict[str, Any]]:
    cache = output / (
        f"reference_cache_seed{seed}_n{count}_cells{reference_cells}.pt"
    )
    if cache.exists():
        payload = torch.load(cache, map_location="cpu", weights_only=True)
        return payload["suites"], payload["stats"]
    initial = {
        "euler": euler_initial_suites(seed, count),
        "swe": swe_initial_suites(seed, count),
    }
    suites: dict[str, dict[str, torch.Tensor]] = {"euler": {}, "swe": {}}
    stats: dict[str, Any] = {"euler": {}, "swe": {}}
    for system in ("euler", "swe"):
        for split, initial64 in initial[system].items():
            print(
                json.dumps(
                    {
                        "stage": "reference",
                        "system": system,
                        "split": split,
                        "trajectories": int(initial64.shape[0]),
                        "reference_cells": reference_cells,
                    }
                ),
                flush=True,
            )
            reference, run_stats = reference_rollout(
                system,
                initial64,
                suite_snapshots(system, split),
                reference_cells,
            )
            suites[system][split] = torch.cat(
                [initial64[:, None], reference[:, 1:]], dim=1
            )
            stats[system][split] = run_stats
    torch.save({"suites": suites, "stats": stats}, cache)
    return suites, stats


def training_statistics(
    system: str, seed: int, output: Path
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if system == "euler":
        train, _ = replicates.euler_cache(output, seed)
    else:
        cache_dir = SWE_DIR / "results" if seed == 0 else output
        train, _ = swe.load_or_make_data(cache_dir, seed)
    state_mean = train.mean(dim=(0, 1, 2))
    state_std = train.std(dim=(0, 1, 2))
    return train, state_mean, state_std


def load_hcfl(system: str, seed: int, output: Path) -> torch.nn.Module:
    train, _, _ = training_statistics(system, seed, output)
    if system == "euler":
        values = euler.base.primitive(train)
        mean = values.mean(dim=(0, 1, 2)).numpy()
        std = values.std(dim=(0, 1, 2)).numpy()
        if seed == 0:
            checkpoint = (
                EULER_DIR
                / "results"
                / "central_nonnegative_feas_s4_converged_best_seed0.pt"
            )
        else:
            checkpoint = output / f"euler_hcfl64_s4_converged_best_seed{seed}.pt"
        return euler.load_model(
            euler.ARM_BY_NAME["central_nonnegative_feas_s4"],
            checkpoint,
            mean,
            std,
            72,
        )
    mean, std, _ = swe.prepare_statistics(train)
    if seed == 0:
        checkpoint = (
            SWE_DIR
            / "results"
            / "central_nonnegative_feas_s4_converged_best_seed0.pt"
        )
    else:
        checkpoint = output / f"swe_hcfl64_s4_converged_best_seed{seed}.pt"
    return swe.load_model(
        next(arm for arm in swe.ARMS if arm.name == "central_nonnegative_feas_s4"),
        checkpoint,
        mean,
        std,
        72,
    )


def load_fno(system: str, output: Path, seed: int = 0) -> operator.ResidualFNO1d:
    report_path = output / f"{system}_fno64_report_seed{seed}.json"
    checkpoint = output / f"{system}_fno64_converged_best_seed{seed}.pt"
    if not report_path.exists() or not checkpoint.exists():
        raise FileNotFoundError(
            f"Missing {system} FNO artifacts; run train_operator_baselines.py"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    _, state_mean, state_std = training_statistics(system, seed, output)
    model = operator.ResidualFNO1d(
        state_mean.numel(),
        state_mean,
        state_std,
        width=int(report["width"]),
        modes=int(report["modes"]),
        layers=int(report["layers"]),
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


def load_learned_flux(
    system: str, output: Path, seed: int = 0
) -> learned_flux.UnconstrainedLearnedFluxFV:
    report_path = output / f"{system}_learned_flux64_report_seed{seed}.json"
    checkpoint = output / (
        f"{system}_learned_flux64_converged_best_seed{seed}.pt"
    )
    if not report_path.exists() or not checkpoint.exists():
        raise FileNotFoundError(
            f"Missing {system} learned-flux artifacts; run "
            "train_learned_flux_baselines.py"
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    train, _, _ = training_statistics(system, seed, output)
    primitives = (
        euler.base.primitive(train)
        if system == "euler"
        else swe.primitive(train)
    )
    model = learned_flux.UnconstrainedLearnedFluxFV(
        system,
        primitives.mean(dim=(0, 1, 2)),
        primitives.std(dim=(0, 1, 2)),
        width=int(report["width"]),
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


@torch.no_grad()
def hcfl_rollout(
    system: str,
    model: torch.nn.Module,
    initial: torch.Tensor,
    snapshots: int,
) -> tuple[torch.Tensor, dict[str, Any]]:
    state = initial.clone()
    saved = [state.clone()]
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["max_entropy_residual"] = -float("inf")
    totals["maximum_entropy_residual"] = -float("inf")
    totals["max_total_entropy_change"] = -float("inf")
    totals["maximum_total_entropy_change"] = -float("inf")
    for _ in range(1, snapshots):
        if system == "euler":
            state, step = euler.shared.advance_safe_snapshot(model, state)
            euler.shared.merge_step_stats(totals, step)
        else:
            state, step = swe.advance_safe_snapshot(model, state)
            swe.merge_stats(totals, step)
        saved.append(state.clone())
    entropy_total = totals.get("entropy_total", 0.0)
    entropy_violations = totals.get("entropy_violations", 0.0)
    maximum_residual = max(
        totals.get("max_entropy_residual", -float("inf")),
        totals.get("maximum_entropy_residual", -float("inf")),
    )
    maximum_change = max(
        totals.get("max_total_entropy_change", -float("inf")),
        totals.get("maximum_total_entropy_change", -float("inf")),
    )
    return torch.stack(saved, dim=1), {
        "interface_entropy_violation_rate": entropy_violations
        / max(entropy_total, 1.0),
        "maximum_interface_entropy_residual": maximum_residual,
        "maximum_internal_total_entropy_change": maximum_change,
        "local_limiter_intervention_rate": totals.get("local_active", 0.0)
        / max(totals.get("local_total", 0.0), 1.0),
        "fd_entropy_intervention_rate": totals.get("fd_active", 0.0)
        / max(totals.get("fd_total", 0.0), 1.0),
        "minimum_fd_beta": totals["fd_beta_min"],
        "internal_substeps": totals.get("substeps", 0.0),
    }


@torch.no_grad()
def hcfl_no_low_rollout(
    system: str,
    model: torch.nn.Module,
    initial: torch.Tensor,
    snapshots: int,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, Any]]:
    """Deploy the hard-projected proposal without an F_low convex anchor.

    The accepted update is the proposal flux itself.  Positivity and the
    fully-discrete entropy inequality are enforced by rejecting the proposed
    time step and halving ``dt``; no low-order flux is computed or blended.
    This is the clean ablation requested for deciding whether F_low is useful
    enough to justify its extra design complexity.
    """
    state = initial.clone()
    saved = [state.clone()]
    dx = 1.0 / state.shape[-2]
    dt_snapshot = (
        euler.base.DT_SNAPSHOT if system == "euler" else swe.DT_SNAPSHOT
    )
    entropy_tolerance = 1.0e-8
    totals: defaultdict[str, float] = defaultdict(float)
    totals["maximum_entropy_residual"] = -float("inf")
    totals["maximum_total_entropy_change"] = -float("inf")

    for _ in range(1, snapshots):
        remaining = dt_snapshot
        while remaining > 1.0e-14:
            dt = min(
                remaining,
                cfl * dx / max(maximum_speed(system, state), 1.0e-12),
            )
            entropy_before = entropy_values(system, state[:, None])[:, 0]
            accepted = False
            for _ in range(30):
                lam = dt / dx
                if system == "euler":
                    flux = euler.shared.strict_entropy_projection(
                        model.flux(state), state
                    )
                else:
                    flux = model.flux(state)
                candidate = state - lam * (
                    flux - torch.roll(flux, 1, dims=-2)
                )
                entropy_change = (
                    entropy_values(system, candidate[:, None])[:, 0]
                    - entropy_before
                )
                if (
                    bool(is_admissible(system, candidate).all())
                    and bool((entropy_change <= entropy_tolerance).all())
                ):
                    accepted = True
                    break
                dt *= 0.5
                totals["direct_dt_halvings"] += 1.0
            if not accepted:
                raise RuntimeError(
                    f"{system} no-F_low proposal could not find an "
                    "admissible entropy-nonincreasing time step"
                )

            if system == "euler":
                residual = euler.shared.entropy_residual64(flux, state)
            else:
                residual = swe.entropy_residual(
                    flux.double(), state.double()
                )
            totals["entropy_violations"] += float(
                (residual > 1.0e-5).sum()
            )
            totals["entropy_total"] += float(residual.numel())
            totals["maximum_entropy_residual"] = max(
                totals["maximum_entropy_residual"], float(residual.max())
            )
            totals["maximum_total_entropy_change"] = max(
                totals["maximum_total_entropy_change"],
                float(entropy_change.max()),
            )
            totals["internal_substeps"] += 1.0
            state = candidate
            remaining -= dt
        saved.append(state.clone())

    return torch.stack(saved, dim=1), {
        "interface_entropy_violation_rate": totals["entropy_violations"]
        / max(totals["entropy_total"], 1.0),
        "maximum_interface_entropy_residual": totals[
            "maximum_entropy_residual"
        ],
        "maximum_internal_total_entropy_change": totals[
            "maximum_total_entropy_change"
        ],
        "direct_dt_halvings": totals["direct_dt_halvings"],
        "internal_substeps": totals["internal_substeps"],
        "uses_low_order_flux": False,
    }


@torch.no_grad()
def fno_rollout(
    model: operator.ResidualFNO1d,
    initial: torch.Tensor,
    snapshots: int,
) -> tuple[torch.Tensor, dict[str, Any]]:
    state = initial.clone()
    saved = [state.clone()]
    for _ in range(1, snapshots):
        state = model(state)
        saved.append(state.clone())
    return torch.stack(saved, dim=1), {}


def entropy_values(system: str, trajectory: torch.Tensor) -> torch.Tensor:
    if system == "euler":
        return euler.base.entropy(trajectory.double()).sum(dim=-1)
    return swe.entropy(trajectory.double()).sum(dim=-1)


def physical_minima(system: str, trajectory: torch.Tensor) -> tuple[float, float]:
    primary = float(trajectory[..., 0].min())
    if system == "euler":
        pressure = (euler.base.GAMMA - 1.0) * (
            trajectory[..., 2]
            - 0.5
            * trajectory[..., 1].square()
            / trajectory[..., 0].clamp_min(1.0e-30)
        )
        return primary, float(pressure.min())
    return primary, float("nan")


def trajectory_metrics(
    system: str,
    reference: torch.Tensor,
    candidate: torch.Tensor,
    state_std: torch.Tensor,
) -> dict[str, Any]:
    finite_trajectory = torch.isfinite(candidate).all(dim=(-1, -2, -3))
    alive = is_admissible(system, candidate).all(dim=(-1, -2))
    completion_rate = float(alive.float().mean())
    if bool(finite_trajectory.all()):
        normalized_all = (candidate - reference) / state_std
        rollout_nrmse_all = float(
            torch.sqrt(normalized_all[:, 1:].double().square().mean())
        )
    else:
        rollout_nrmse_all = float("inf")

    # Conservation is meaningful for finite nonphysical states too: a
    # negative-pressure learned-FV trajectory may have failed positivity while
    # still conserving exactly.  Do not conflate these two properties.
    finite_selected = candidate[finite_trajectory]
    if finite_selected.shape[0]:
        finite_initial_sum = finite_selected[:, 0].double().sum(dim=-2)
        finite_sums = finite_selected.double().sum(dim=-2)
        finite_denominator = (
            finite_selected[:, 0]
            .double()
            .abs()
            .sum(dim=(-2, -1))
            .clamp_min(1.0e-14)
        )
        finite_conservation = (
            (finite_sums - finite_initial_sum[:, None]).abs().sum(dim=-1)
            / finite_denominator[:, None]
        )
        maximum_conservation_drift = float(finite_conservation.max())
        minimum_primary, minimum_pressure = physical_minima(
            system, finite_selected
        )
    else:
        maximum_conservation_drift = float("inf")
        minimum_primary = float("nan")
        minimum_pressure = float("nan")

    if bool(alive.any()):
        normalized = (candidate[alive] - reference[alive]) / state_std
        rollout_nrmse = float(
            torch.sqrt(normalized[:, 1:].double().square().mean())
        )
        final_nrmse = float(
            torch.sqrt(normalized[:, -1].double().square().mean())
        )
        selected = candidate[alive]
        entropy = entropy_values(system, selected)
        changes = entropy[:, 1:] - entropy[:, :-1]
        entropy_tolerance = 1.0e-8
        saved_violation_rate = float(
            (changes > entropy_tolerance).double().mean()
        )
        maximum_saved_entropy_increase = float(changes.max())
    else:
        rollout_nrmse = float("inf")
        final_nrmse = float("inf")
        saved_violation_rate = float("nan")
        maximum_saved_entropy_increase = float("nan")
    return {
        "trajectory_count": int(candidate.shape[0]),
        "completed_trajectories": int(alive.sum()),
        "completion_rate": completion_rate,
        "rollout_nrmse_all_finite_outputs": rollout_nrmse_all,
        "rollout_nrmse_completed_only": rollout_nrmse,
        "final_snapshot_nrmse_completed_only": final_nrmse,
        "maximum_relative_conservation_drift": maximum_conservation_drift,
        "saved_snapshot_entropy_violation_rate": saved_violation_rate,
        "maximum_saved_snapshot_entropy_increase": maximum_saved_entropy_increase,
        "minimum_density_or_depth": minimum_primary,
        "minimum_pressure": minimum_pressure,
    }


def evaluate(
    args: argparse.Namespace,
    output: Path,
    references: dict[str, dict[str, torch.Tensor]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    models = {
        system: {seed: load_hcfl(system, seed, output) for seed in args.hcfl_seeds}
        for system in ("euler", "swe")
    }
    fnos = {system: load_fno(system, output, args.fno_seed) for system in ("euler", "swe")}
    learned_flux_models = {
        system: load_learned_flux(system, output, args.learned_flux_seed)
        for system in ("euler", "swe")
    }
    state_stds = {
        system: training_statistics(system, 0, output)[2]
        for system in ("euler", "swe")
    }
    rows: list[dict[str, Any]] = []
    trajectories: dict[str, Any] = {"euler": {}, "swe": {}}

    for system in ("euler", "swe"):
        random_names = EULER_RANDOM_SPLITS if system == "euler" else SWE_RANDOM_SPLITS
        canonical_names = EULER_CANONICAL if system == "euler" else SWE_CANONICAL
        for split in (*random_names, *canonical_names):
            reference = references[system][split]
            initial = reference[:, 0]
            snapshots = reference.shape[1]
            print(json.dumps({"stage": "evaluate", "system": system, "split": split}), flush=True)
            if system == "euler":
                native, native_stats = ssprk2_rollout(
                    system,
                    initial,
                    snapshots,
                    euler.base.DT_SNAPSHOT,
                    euler.base.t_hllc,
                    0.20,
                )
                muscl, muscl_stats = ssprk2_rollout(
                    system,
                    initial,
                    snapshots,
                    euler.base.DT_SNAPSHOT,
                    euler_muscl_flux,
                    0.20,
                )
            else:
                native, native_stats = ssprk2_rollout(
                    system,
                    initial,
                    snapshots,
                    swe.DT_SNAPSHOT,
                    swe_first_order_flux,
                    0.20,
                )
                muscl, muscl_stats = ssprk2_rollout(
                    system,
                    initial,
                    snapshots,
                    swe.DT_SNAPSHOT,
                    swe_muscl_flux,
                    0.20,
                )
            method_trajectories: dict[str, Any] = {
                "reference": reference,
                "native_fvm64": native,
                "muscl_fvm64": muscl,
            }
            for method, candidate, run_stats in (
                ("native_fvm64", native, native_stats),
                ("muscl_fvm64", muscl, muscl_stats),
            ):
                row = {
                    "system": system,
                    "split_group": "random" if split in random_names else "canonical",
                    "split": split,
                    "method": method,
                    "seed": -1,
                    "parameter_count": 0,
                    **trajectory_metrics(system, reference, candidate, state_stds[system]),
                    **run_stats,
                }
                rows.append(row)

            for seed, model in models[system].items():
                candidate, run_stats = hcfl_rollout(system, model, initial, snapshots)
                method_trajectories[f"hcfl64_seed{seed}"] = candidate
                rows.append(
                    {
                        "system": system,
                        "split_group": "random" if split in random_names else "canonical",
                        "split": split,
                        "method": "hcfl64",
                        "seed": seed,
                        "parameter_count": sum(
                            parameter.numel() for parameter in model.parameters()
                        ),
                        **trajectory_metrics(system, reference, candidate, state_stds[system]),
                        **run_stats,
                    }
                )

                no_low_candidate, no_low_stats = hcfl_no_low_rollout(
                    system, model, initial, snapshots
                )
                method_trajectories[
                    f"hcfl64_no_low_seed{seed}"
                ] = no_low_candidate
                rows.append(
                    {
                        "system": system,
                        "split_group": (
                            "random" if split in random_names else "canonical"
                        ),
                        "split": split,
                        "method": "hcfl64_no_low",
                        "seed": seed,
                        "parameter_count": sum(
                            parameter.numel()
                            for parameter in model.parameters()
                        ),
                        **trajectory_metrics(
                            system,
                            reference,
                            no_low_candidate,
                            state_stds[system],
                        ),
                        **no_low_stats,
                    }
                )

            learned_flux_candidate = learned_flux.raw_rollout(
                system,
                learned_flux_models[system],
                initial,
                snapshots,
            )
            method_trajectories["learned_flux64"] = learned_flux_candidate
            rows.append(
                {
                    "system": system,
                    "split_group": (
                        "random" if split in random_names else "canonical"
                    ),
                    "split": split,
                    "method": "learned_flux64",
                    "seed": args.learned_flux_seed,
                    "parameter_count": sum(
                        parameter.numel()
                        for parameter in learned_flux_models[system].parameters()
                    ),
                    **trajectory_metrics(
                        system,
                        reference,
                        learned_flux_candidate,
                        state_stds[system],
                    ),
                }
            )

            fno_candidate, fno_stats = fno_rollout(fnos[system], initial, snapshots)
            method_trajectories["fno64"] = fno_candidate
            rows.append(
                {
                    "system": system,
                    "split_group": "random" if split in random_names else "canonical",
                    "split": split,
                    "method": "fno64",
                    "seed": args.fno_seed,
                    "parameter_count": sum(
                        parameter.numel()
                        for parameter in fnos[system].parameters()
                    ),
                    **trajectory_metrics(system, reference, fno_candidate, state_stds[system]),
                    **fno_stats,
                }
            )
            trajectories[system][split] = method_trajectories
    return rows, trajectories


def aggregate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    aggregates: list[dict[str, Any]] = []
    keys = sorted(
        {
            (row["system"], row["split_group"], row["method"], row["seed"])
            for row in rows
        }
    )
    for system, group, method, seed in keys:
        selected = [
            row
            for row in rows
            if (
                row["system"],
                row["split_group"],
                row["method"],
                row["seed"],
            )
            == (system, group, method, seed)
        ]
        complete = [row for row in selected if row["completion_rate"] == 1.0]
        aggregates.append(
            {
                "system": system,
                "split_group": group,
                "method": method,
                "seed": seed,
                "parameter_count": selected[0]["parameter_count"],
                "case_count": len(selected),
                "fully_completed_case_count": len(complete),
                "mean_completion_rate": float(
                    np.mean([row["completion_rate"] for row in selected])
                ),
                "mean_rollout_nrmse_over_completed_cases": (
                    float(
                        np.mean(
                            [
                                row["rollout_nrmse_completed_only"]
                                for row in complete
                            ]
                        )
                    )
                    if complete
                    else float("inf")
                ),
                "mean_rollout_nrmse_all_finite_outputs": float(
                    np.mean(
                        [
                            row["rollout_nrmse_all_finite_outputs"]
                            for row in selected
                        ]
                    )
                ),
                "maximum_relative_conservation_drift": float(
                    max(row["maximum_relative_conservation_drift"] for row in selected)
                ),
                "maximum_saved_snapshot_entropy_increase": float(
                    np.nanmax(
                        [row["maximum_saved_snapshot_entropy_increase"] for row in selected]
                    )
                ),
                "minimum_density_or_depth": float(
                    min(row["minimum_density_or_depth"] for row in selected)
                ),
                "minimum_pressure": (
                    float(min(row["minimum_pressure"] for row in selected))
                    if system == "euler"
                    else float("nan")
                ),
                "maximum_interface_entropy_residual": (
                    float(
                        max(
                            row.get("maximum_interface_entropy_residual", -float("inf"))
                            for row in selected
                        )
                    )
                    if method in ("hcfl64", "hcfl64_no_low")
                    else float("nan")
                ),
                "maximum_internal_total_entropy_change": (
                    float(
                        max(
                            row.get("maximum_internal_total_entropy_change", -float("inf"))
                            for row in selected
                        )
                    )
                    if method in ("hcfl64", "hcfl64_no_low")
                    else float("nan")
                ),
                "mean_internal_substeps": float(
                    np.mean(
                        [row.get("internal_substeps", float("nan")) for row in selected]
                    )
                ),
                "mean_local_limiter_intervention_rate": float(
                    np.nanmean(
                        [
                            row.get(
                                "local_limiter_intervention_rate",
                                float("nan"),
                            )
                            for row in selected
                        ]
                    )
                )
                if any("local_limiter_intervention_rate" in row for row in selected)
                else float("nan"),
                "mean_fd_entropy_intervention_rate": float(
                    np.nanmean(
                        [
                            row.get(
                                "fd_entropy_intervention_rate",
                                float("nan"),
                            )
                            for row in selected
                        ]
                    )
                )
                if any("fd_entropy_intervention_rate" in row for row in selected)
                else float("nan"),
                "mean_direct_dt_halvings": float(
                    np.mean(
                        [row.get("direct_dt_halvings", 0.0) for row in selected]
                    )
                ),
            }
        )
    return aggregates


def plot_accuracy(aggregates: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.0), constrained_layout=True)
    methods = (
        "native_fvm64",
        "muscl_fvm64",
        "hcfl64",
        "learned_flux64",
        "fno64",
    )
    for row_index, system in enumerate(("euler", "swe")):
        for column_index, group in enumerate(("random", "canonical")):
            axis = axes[row_index, column_index]
            positions = np.arange(len(methods))
            means: list[float] = []
            errors: list[float] = []
            completion: list[float] = []
            for method in methods:
                selected = [
                    row
                    for row in aggregates
                    if row["system"] == system
                    and row["split_group"] == group
                    and row["method"] == method
                ]
                values = np.array(
                    [row["mean_rollout_nrmse_all_finite_outputs"] for row in selected],
                    dtype=np.float64,
                )
                means.append(float(values.mean()))
                errors.append(float(values.std(ddof=1)) if len(values) > 1 else 0.0)
                completion.append(float(np.mean([row["mean_completion_rate"] for row in selected])))
            bars = axis.bar(
                positions,
                means,
                yerr=errors,
                capsize=3,
                color=[METHOD_COLORS[method] for method in methods],
                alpha=0.84,
            )
            for bar, rate in zip(bars, completion):
                if rate < 1.0:
                    axis.annotate(
                        f"{100 * rate:.0f}% complete",
                        xy=(
                            bar.get_x() + bar.get_width() / 2,
                            bar.get_height(),
                        ),
                        xytext=(0, 4),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        rotation=90,
                        clip_on=False,
                    )
            axis.set_xticks(positions, [METHOD_LABELS[m] for m in methods], rotation=18, ha="right")
            axis.set_ylabel("rollout NRMSE")
            axis.set_title(f"{system.upper()} — {group} tests")
            axis.grid(axis="y", alpha=0.2)
            axis.margins(y=0.16)
    fig.suptitle("Matched 64-cell accuracy (HCFL bars show mean ± sample std over seeds)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_low_order_ablation(
    aggregates: list[dict[str, Any]], output: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    configurations = (
        ("hcfl64", "with $F_{low}$ limiting"),
        ("hcfl64_no_low", "no $F_{low}$; reject/halve $dt$"),
    )
    group_positions = np.arange(2)
    width = 0.34
    for axis, system in zip(axes, ("euler", "swe")):
        for offset, (method, label) in enumerate(configurations):
            means: list[float] = []
            errors: list[float] = []
            for group in ("random", "canonical"):
                selected = [
                    row
                    for row in aggregates
                    if row["system"] == system
                    and row["split_group"] == group
                    and row["method"] == method
                ]
                values = np.asarray(
                    [
                        row["mean_rollout_nrmse_all_finite_outputs"]
                        for row in selected
                    ],
                    dtype=np.float64,
                )
                means.append(float(values.mean()))
                errors.append(
                    float(values.std(ddof=1)) if len(values) > 1 else 0.0
                )
            axis.bar(
                group_positions + (offset - 0.5) * width,
                means,
                width,
                yerr=errors,
                capsize=3,
                color=METHOD_COLORS[method],
                alpha=0.84,
                hatch=None if method == "hcfl64" else "//",
                label=label,
            )
        axis.set_xticks(group_positions, ("random", "canonical"))
        axis.set_ylabel("rollout NRMSE")
        axis.set_title(system.upper())
        axis.grid(axis="y", alpha=0.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.subplots_adjust(top=0.72, bottom=0.16, left=0.08, right=0.98, wspace=0.20)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.87),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        "Low-order-anchor ablation at 64 cells (mean ± seed std)", y=0.98
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_profiles(trajectories: dict[str, Any], output: Path) -> None:
    selections = [
        ("euler", "lax", ("density", "velocity", "pressure")),
        ("swe", "counterflow_collision", ("depth", "velocity")),
    ]
    fig, axes = plt.subplots(5, 1, figsize=(10.2, 12.0))
    axis_index = 0
    x = (np.arange(CELLS) + 0.5) / CELLS
    for system, split, names in selections:
        data = trajectories[system][split]
        candidates = {
            "reference": data["reference"],
            "native_fvm64": data["native_fvm64"],
            "muscl_fvm64": data["muscl_fvm64"],
            "hcfl64": data["hcfl64_seed0"],
            "learned_flux64": data["learned_flux64"],
            "fno64": data["fno64"],
        }
        primitives: dict[str, np.ndarray] = {}
        for method, trajectory in candidates.items():
            final = trajectory[0, -1]
            if system == "euler":
                primitives[method] = euler.base.primitive(final).detach().numpy()
            else:
                primitives[method] = swe.primitive(final).detach().numpy()
        for variable, name in enumerate(names):
            axis = axes[axis_index]
            axis_index += 1
            axis.plot(x, primitives["reference"][:, variable], color="black", linewidth=2.0, label="fine reference")
            axis.plot(x, primitives["native_fvm64"][:, variable], color=METHOD_COLORS["native_fvm64"], linestyle="--", linewidth=1.6, label=METHOD_LABELS["native_fvm64"])
            axis.plot(x, primitives["muscl_fvm64"][:, variable], color=METHOD_COLORS["muscl_fvm64"], linestyle="-.", linewidth=1.6, label=METHOD_LABELS["muscl_fvm64"])
            axis.plot(x, primitives["hcfl64"][:, variable], color=METHOD_COLORS["hcfl64"], linestyle=(0, (6, 2)), linewidth=1.8, label=METHOD_LABELS["hcfl64"])
            if bool(is_admissible(system, candidates["learned_flux64"]).all()):
                axis.plot(x, primitives["learned_flux64"][:, variable], color=METHOD_COLORS["learned_flux64"], linestyle=(0, (4, 2, 1, 2)), linewidth=1.5, label=METHOD_LABELS["learned_flux64"])
            if bool(is_admissible(system, candidates["fno64"]).all()):
                axis.plot(x, primitives["fno64"][:, variable], color=METHOD_COLORS["fno64"], linestyle=(0, (2, 2)), linewidth=1.5, label=METHOD_LABELS["fno64"])
            axis.set_ylabel(name)
            axis.set_title(f"{system.upper()} {split.replace('_', ' ')}")
            axis.grid(alpha=0.18)
    axes[-1].set_xlabel("x")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.subplots_adjust(
        top=0.94, bottom=0.055, left=0.075, right=0.99, hspace=0.38
    )
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=6,
        frameon=False,
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def self_test() -> None:
    euler_state = torch.from_numpy(
        euler.base.generate_ic(2, CELLS, 7, ood=False)
    ).double()
    swe_state = torch.from_numpy(
        swe.legacy.generate_ic(2, CELLS, 8, ood=False)
    ).double()
    for system, state, first, second in (
        ("euler", euler_state, euler.base.t_hllc, euler_muscl_flux),
        ("swe", swe_state, swe_first_order_flux, swe_muscl_flux),
    ):
        for flux in (first, second):
            value = flux(state)
            if value.shape != state.shape or not bool(torch.isfinite(value).all()):
                raise RuntimeError(f"Invalid {system} flux in self-test")
        trajectory, _ = ssprk2_rollout(
            system,
            state,
            3,
            euler.base.DT_SNAPSHOT if system == "euler" else swe.DT_SNAPSHOT,
            second,
            0.20,
        )
        if not bool(is_admissible(system, trajectory).all()):
            raise RuntimeError(f"{system} MUSCL self-test lost admissibility")
    print("self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--random-trajectories", type=int, default=24)
    parser.add_argument("--reference-cells", type=int, default=2048)
    parser.add_argument("--hcfl-seeds", nargs="+", type=int, default=[0, 1, 2])
    parser.add_argument("--fno-seed", type=int, default=0)
    parser.add_argument("--learned-flux-seed", type=int, default=0)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    references, reference_stats = load_or_make_references(
        output,
        args.seed,
        args.random_trajectories,
        args.reference_cells,
    )
    rows, trajectories = evaluate(args, output, references)
    aggregates = aggregate_rows(rows)
    write_csv(output / "benchmark64_metrics.csv", rows)
    write_csv(output / "benchmark64_aggregate.csv", aggregates)
    payload = {
        "protocol": {
            "candidate_cells": CELLS,
            "reference_cells": args.reference_cells,
            "reference_initialization": (
                "piecewise-constant prolongation of the exact same 64-cell "
                "initial averages; no hidden subcell initial information"
            ),
            "euler_reference": "HLLC + SSP-RK2",
            "swe_reference": "HLL + SSP-RK2",
            "hcfl_seeds": args.hcfl_seeds,
            "fno_scope": (
                "compact matched periodic residual FNO; architecture baseline, "
                "not an exact reproduction of published FNO experiments"
            ),
            "learned_flux_scope": (
                "matched conservative four-cell unconstrained learned-flux "
                "baseline; family comparison, not an exact reproduction of "
                "one published neural-FV architecture"
            ),
        },
        "reference_generation": reference_stats,
        "aggregate": aggregates,
        "rows": rows,
    }
    (output / "benchmark64_summary.json").write_text(
        json.dumps(json_ready(payload), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    plot_accuracy(aggregates, output / "benchmark64_accuracy.png")
    plot_low_order_ablation(
        aggregates, output / "benchmark64_low_order_ablation.png"
    )
    plot_profiles(trajectories, output / "benchmark64_profiles.png")
    print(json.dumps({"stage": "complete", "aggregate": aggregates}, indent=2), flush=True)


if __name__ == "__main__":
    main()
