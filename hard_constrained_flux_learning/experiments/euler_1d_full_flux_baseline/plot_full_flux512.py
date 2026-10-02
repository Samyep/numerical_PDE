"""Test direct full-flux transfer from 64-cell training to 512 cells."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
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
AUDIT_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_convergence_audit"
if str(AUDIT_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(AUDIT_EXPERIMENT))

import plot_hllc512_vs_hcfl512 as hllc  # noqa: E402


comparison = hllc.comparison
base = hllc.base
shared = hllc.shared
REFERENCE_COLOR = comparison.BACKGROUND_REFERENCE_COLOR
HLLC_COLOR = "#0072B2"
DISSIPATION_COLOR = "#009E73"
FULL_COLOR = "#D55E00"
MAX_INTERNAL_SUBSTEPS_PER_SAFE_CALL = 64


def load_solver(
    model_name: str,
    checkpoint: Path,
    width: int,
) -> base.Solver:
    model = base.Solver(
        model_name, np.zeros(3), np.ones(3), width=width
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


@torch.no_grad()
def guarded_advance_safe_snapshot(
    model: base.Solver,
    state: torch.Tensor,
    cfl: float = 0.42,
    max_substeps: int = MAX_INTERNAL_SUBSTEPS_PER_SAFE_CALL,
) -> tuple[torch.Tensor, dict[str, float], bool]:
    """Shared safe update with an explicit computational-cost guard."""
    dx = 1.0 / base.NCOARSE
    remaining = base.DT_SNAPSHOT
    stats = {
        "local_active": 0.0,
        "local_total": 0.0,
        "fd_active": 0.0,
        "fd_total": 0.0,
        "fd_beta_sum": 0.0,
        "fd_beta_min": 1.0,
        "entropy_violations": 0.0,
        "entropy_total": 0.0,
        "max_entropy_residual": -float("inf"),
        "max_total_entropy_change": -float("inf"),
        "substeps": 0.0,
        "low_order_dt_halvings": 0.0,
        "maximum_characteristic_speed": 0.0,
        "minimum_dt": float("inf"),
    }

    while remaining > 1.0e-14:
        if stats["substeps"] >= max_substeps:
            stats["remaining_nominal_time"] = remaining
            return state, stats, False

        density = state[..., 0].clamp_min(1.0e-10)
        pressure = base.t_pressure(state).clamp_min(1.0e-10)
        velocity = state[..., 1] / density
        sound_speed = torch.sqrt(base.GAMMA * pressure / density)
        max_speed = float((torch.abs(velocity) + sound_speed).max())
        dt = min(remaining, cfl * dx / max(max_speed, 1.0e-12))
        stats["maximum_characteristic_speed"] = max(
            stats["maximum_characteristic_speed"], max_speed
        )

        for _ in range(24):
            lam = dt / dx
            low_flux = shared.strict_entropy_projection(
                base.t_rusanov(state), state
            )
            low_state = state - lam * (
                low_flux - torch.roll(low_flux, 1, dims=-2)
            )
            entropy_ok = bool(
                (
                    shared.total_entropy(low_state)
                    <= shared.total_entropy(state) + 1.0e-10
                ).all()
            )
            if bool(shared.admissible(low_state).all()) and entropy_ok:
                break
            dt *= 0.5
            stats["low_order_dt_halvings"] += 1.0
        else:
            raise RuntimeError(
                "Could not establish the guarded low-order safety premise"
            )

        stats["minimum_dt"] = min(stats["minimum_dt"], dt)
        high_flux = shared.strict_entropy_projection(model.flux(state), state)
        local_flux, alpha = shared.local_admissibility_limiter(
            state, high_flux, low_flux, lam
        )
        final_flux, beta = shared.global_entropy_limiter(
            state, local_flux, low_flux, lam
        )
        residual = shared.entropy_residual64(final_flux, state)
        next_state = state - lam * (
            final_flux - torch.roll(final_flux, 1, dims=-2)
        )
        entropy_change = (
            shared.total_entropy(next_state) - shared.total_entropy(state)
        )
        if (
            not bool(torch.isfinite(next_state).all())
            or not bool(shared.admissible(next_state).all())
        ):
            raise RuntimeError("Guarded safe update became inadmissible")

        stats["local_active"] += float((alpha < 1.0 - 1.0e-7).sum())
        stats["local_total"] += float(alpha.numel())
        stats["fd_active"] += float((beta < 1.0 - 1.0e-7).sum())
        stats["fd_total"] += float(beta.numel())
        stats["fd_beta_sum"] += float(beta.sum())
        stats["fd_beta_min"] = min(
            stats["fd_beta_min"], float(beta.min())
        )
        stats["entropy_violations"] += float(
            (residual > shared.ENTROPY_RESIDUAL_TOLERANCE).sum()
        )
        stats["entropy_total"] += float(residual.numel())
        stats["max_entropy_residual"] = max(
            stats["max_entropy_residual"], float(residual.max())
        )
        stats["max_total_entropy_change"] = max(
            stats["max_total_entropy_change"], float(entropy_change.max())
        )
        stats["substeps"] += 1.0
        state = next_state
        remaining -= dt

    stats["remaining_nominal_time"] = 0.0
    return state, stats, True


@torch.no_grad()
def guarded_matched_safe_rollout(
    model: base.Solver,
    name: str,
    cells: int,
) -> tuple[torch.Tensor, dict[str, Any]]:
    """Matched rollout that records a deterministic CFL-cost failure."""
    updates_per_snapshot = cells // base.NCOARSE
    state = torch.from_numpy(comparison.initial_condition(name, cells)).float()
    snapshots = [state.clone()]
    totals: defaultdict[str, float] = hllc.empty_safety_totals()
    update_calls = 0
    projection_active = 0
    projection_total = 0
    projection_squared = 0.0
    raw_squared = 0.0
    projection_maximum = 0.0
    maximum_speed = 0.0
    minimum_dt = float("inf")
    low_order_dt_halvings = 0
    started = time.perf_counter()

    for snapshot in range(1, shared.CANONICAL_NSNAP):
        for _ in range(updates_per_snapshot):
            raw_flux = model.flux_net(state)
            projected_flux = model.flux(state)
            projection_delta = projected_flux - raw_flux
            tolerance = 1.0e-7 * (
                1.0 + raw_flux.abs().amax(dim=-1)
            )
            projection_active += int(
                (
                    torch.linalg.vector_norm(projection_delta, dim=-1)
                    > tolerance
                ).sum()
            )
            projection_total += projection_delta.shape[0] * (
                projection_delta.shape[1]
            )
            projection_squared += float(
                (projection_delta.double() ** 2).sum()
            )
            raw_squared += float((raw_flux.double() ** 2).sum())
            projection_maximum = max(
                projection_maximum, float(projection_delta.abs().max())
            )

            state, step, completed = guarded_advance_safe_snapshot(
                model, state
            )
            shared.merge_step_stats(totals, step)
            update_calls += 1
            maximum_speed = max(
                maximum_speed, step["maximum_characteristic_speed"]
            )
            minimum_dt = min(minimum_dt, step["minimum_dt"])
            low_order_dt_halvings += int(step["low_order_dt_halvings"])
            if not completed:
                pressure = hllc.physical_pressure(state)
                return torch.stack(snapshots, dim=1), {
                    "completed": False,
                    "failure_reason": "safe_internal_substep_limit_exceeded",
                    "internal_substep_limit_per_safe_call": (
                        MAX_INTERNAL_SUBSTEPS_PER_SAFE_CALL
                    ),
                    "failed_saved_snapshot": snapshot,
                    "completed_saved_snapshots": len(snapshots),
                    "safe_update_calls_before_failure": update_calls,
                    "internal_substeps_before_failure": int(
                        totals["substeps"]
                    ),
                    "maximum_characteristic_speed_before_failure": (
                        maximum_speed
                    ),
                    "minimum_dt_before_failure": minimum_dt,
                    "low_order_dt_halvings_before_failure": (
                        low_order_dt_halvings
                    ),
                    "minimum_density_before_failure": float(
                        state[..., 0].min()
                    ),
                    "minimum_pressure_before_failure": float(pressure.min()),
                    "wall_seconds_before_failure": (
                        time.perf_counter() - started
                    ),
                }
        snapshots.append(state.clone())

    trajectory = torch.stack(snapshots, dim=1)
    return trajectory, {
        "completed": True,
        "effective_physical_dt_per_call": (
            base.DT_SNAPSHOT / updates_per_snapshot
        ),
        "safe_update_calls": update_calls,
        "internal_substeps": int(totals["substeps"]),
        "maximum_characteristic_speed": maximum_speed,
        "minimum_internal_dt": minimum_dt,
        "low_order_dt_halvings": low_order_dt_halvings,
        "wall_seconds": time.perf_counter() - started,
        "hard_projection_intervention_rate": (
            projection_active / max(projection_total, 1)
        ),
        "hard_projection_relative_flux_rms": np.sqrt(
            projection_squared / max(raw_squared, 1.0e-30)
        ),
        "hard_projection_max_absolute_change": projection_maximum,
        "local_limiter_intervention_rate": (
            totals["local_active"] / max(totals["local_total"], 1.0)
        ),
        "fd_entropy_intervention_rate": (
            totals["fd_active"] / max(totals["fd_total"], 1.0)
        ),
        "mean_fd_beta": (
            totals["fd_beta_sum"] / max(totals["fd_total"], 1.0)
        ),
        "min_fd_beta": totals["fd_beta_min"],
        "entropy_violation_rate": (
            totals["entropy_violations"]
            / max(totals["entropy_total"], 1.0)
        ),
        "max_entropy_residual": totals["max_entropy_residual"],
        "max_total_entropy_change": totals["max_total_entropy_change"],
    }


def plot_profiles(
    trajectories: dict[str, dict[str, torch.Tensor]],
    cases: dict[str, dict[str, Any]],
    output: Path,
) -> None:
    x_reference = (
        np.arange(hllc.REFERENCE_CELLS) + 0.5
    ) / hllc.REFERENCE_CELLS
    x_target = (
        np.arange(hllc.TARGET_CELLS) + 0.5
    ) / hllc.TARGET_CELLS
    row_labels = [r"Density $\rho$", r"Velocity $u$", r"Pressure $p$"]
    figure, axes = plt.subplots(
        3,
        len(comparison.CASES),
        figsize=(18, 8.8),
        sharex=True,
    )
    figure.subplots_adjust(
        left=0.065,
        right=0.99,
        bottom=0.09,
        top=0.78,
        wspace=0.22,
        hspace=0.14,
    )

    for column, name in enumerate(comparison.CASES):
        primitive = {
            method: base.primitive(trajectory).numpy()[0, -1]
            for method, trajectory in trajectories[name].items()
        }
        dissipation_error = cases[name]["roe_dissipation_hcfl_512"][
            "rollout_nrmse"
        ]
        if cases[name]["direct_full_flux_512"]["completed"]:
            full_error = cases[name]["direct_full_flux_512"][
                "rollout_nrmse"
            ]
            subtitle = f"full / Roe NRMSE = {full_error / dissipation_error:.2f}x"
        else:
            subtitle = "full flux: safe-step $\\Delta t$ collapse"
        axes[0, column].set_title(
            f"{comparison.DISPLAY[name]}\n{subtitle}",
            fontsize=9.8,
            fontweight="semibold",
            color=comparison.TEXT_COLOR,
        )

        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                primitive["reference"][:, row],
                color=REFERENCE_COLOR,
                linewidth=1.7,
                alpha=0.72,
                zorder=1,
            )
            axis.plot(
                x_target,
                primitive["hllc_512"][:, row],
                color=HLLC_COLOR,
                linewidth=1.3,
                zorder=2,
            )
            axis.plot(
                x_target,
                primitive["dissipation_512"][:, row],
                color=DISSIPATION_COLOR,
                linewidth=1.45,
                zorder=3,
            )
            if "full_512" in primitive:
                axis.plot(
                    x_target,
                    primitive["full_512"][:, row],
                    color=FULL_COLOR,
                    linewidth=1.35,
                    linestyle="--",
                    zorder=4,
                )
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    handles = [
        Line2D(
            [0], [0], color=REFERENCE_COLOR, linewidth=1.7, alpha=0.72,
            label="native HLLC-2048 reference",
        ),
        Line2D(
            [0], [0], color=HLLC_COLOR, linewidth=1.3,
            label="native HLLC-512",
        ),
        Line2D(
            [0], [0], color=DISSIPATION_COLOR, linewidth=1.45,
            label="Roe-dissipation HCFL-512",
        ),
        Line2D(
            [0], [0], color=FULL_COLOR, linewidth=1.35, linestyle="--",
            label="direct full-flux HCFL-512",
        ),
    ]
    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"Direct full flux versus Roe-dissipation HCFL at t = {final_time:.4f}",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.895),
        ncol=4,
        frameon=False,
        fontsize=9.1,
    )
    figure.text(
        0.99,
        0.018,
        "Both neural models were trained on 64-cell trajectories and deployed "
        "on 512 cells. Curves are native-grid values; NRMSE uses conservative "
        "HLLC-2048 to 512 restriction only for the metric.",
        ha="right",
        fontsize=8.3,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    for name, case in cases.items():
        for method in (
            "native_hllc_512",
            "roe_dissipation_hcfl_512",
            "direct_full_flux_512",
        ):
            metrics = case[method]
            rows.append(
                {
                    "case": name,
                    "method": method,
                    "completed": metrics["completed"],
                    "failure_reason": metrics.get("failure_reason"),
                    "rollout_nrmse": metrics.get("rollout_nrmse"),
                    "final_snapshot_nrmse": metrics.get(
                        "final_snapshot_nrmse"
                    ),
                    "minimum_density": metrics.get("minimum_density"),
                    "minimum_pressure": metrics.get("minimum_pressure"),
                    "hard_projection_intervention_rate": metrics.get(
                        "hard_projection_intervention_rate"
                    ),
                    "local_limiter_intervention_rate": metrics.get(
                        "local_limiter_intervention_rate"
                    ),
                    "fd_entropy_intervention_rate": metrics.get(
                        "fd_entropy_intervention_rate"
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
    parser.add_argument(
        "--results-dir", type=Path, default=HERE / "results"
    )
    args = parser.parse_args()
    results = args.results_dir.resolve()
    results.mkdir(parents=True, exist_ok=True)

    _mean, _std, state_std = hllc.training_statistics(args.seed)
    full = load_solver(
        "full",
        results / f"full_flux_broad_converged_best_seed{args.seed}.pt",
        args.width,
    )
    dissipation = load_solver(
        "dissipation",
        AUDIT_EXPERIMENT
        / "results"
        / f"dissipation_broad_converged_best_seed{args.seed}.pt",
        args.width,
    )

    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    cases: dict[str, dict[str, Any]] = {}
    for name in comparison.CASES:
        print(f"Evaluating {comparison.DISPLAY[name]}...", flush=True)
        reference, reference_stats = hllc.strict_native_hllc_rollout(
            name, hllc.REFERENCE_CELLS
        )
        restricted = comparison.precision.conservative_restrict(
            reference.numpy(), hllc.TARGET_CELLS
        )
        reference_512 = torch.from_numpy(restricted.astype(np.float32))
        hllc_512, hllc_stats = hllc.strict_native_hllc_rollout(
            name, hllc.TARGET_CELLS
        )
        dissipation_512, dissipation_stats = hllc.matched_safe_rollout(
            dissipation, name, hllc.TARGET_CELLS
        )
        full_512, full_stats = guarded_matched_safe_rollout(
            full, name, hllc.TARGET_CELLS
        )

        trajectories[name] = {
            "reference": reference,
            "hllc_512": hllc_512,
            "dissipation_512": dissipation_512,
        }
        if full_stats["completed"]:
            trajectories[name]["full_512"] = full_512
            full_metrics = hllc.enrich_metrics(
                reference_512, full_512, state_std, full_stats
            )
        else:
            full_metrics = full_stats
        cases[name] = {
            "reference_hllc_2048": reference_stats,
            "native_hllc_512": hllc.enrich_metrics(
                reference_512, hllc_512, state_std, hllc_stats
            ),
            "roe_dissipation_hcfl_512": hllc.enrich_metrics(
                reference_512,
                dissipation_512,
                state_std,
                dissipation_stats,
            ),
            "direct_full_flux_512": full_metrics,
        }

    methods = (
        "native_hllc_512",
        "roe_dissipation_hcfl_512",
        "direct_full_flux_512",
    )
    mean_metrics: dict[str, dict[str, Any]] = {
        method: {
            "mean_rollout_nrmse": float(np.mean([
                cases[name][method]["rollout_nrmse"] for name in cases
            ])),
            "mean_final_snapshot_nrmse": float(np.mean([
                cases[name][method]["final_snapshot_nrmse"]
                for name in cases
            ])),
        }
        for method in methods
        if method != "direct_full_flux_512"
    }
    completed_full_cases = [
        name
        for name in cases
        if cases[name]["direct_full_flux_512"]["completed"]
    ]
    failed_full_cases = [
        name for name in cases if name not in completed_full_cases
    ]
    full_values = [
        cases[name]["direct_full_flux_512"]["rollout_nrmse"]
        for name in completed_full_cases
    ]
    mean_metrics["direct_full_flux_512"] = {
        "completed_cases": len(completed_full_cases),
        "failed_cases": failed_full_cases,
        "mean_rollout_nrmse_over_completed_cases_only": (
            float(np.mean(full_values)) if full_values else None
        ),
        "all_five_case_mean_is_valid": not failed_full_cases,
    }

    summary = {
        "seed": args.seed,
        "full_flux_checkpoint": (
            f"full_flux_broad_converged_best_seed{args.seed}.pt"
        ),
        "roe_dissipation_checkpoint": (
            f"dissipation_broad_converged_best_seed{args.seed}.pt"
        ),
        "checkpoint_training_cells": base.NCOARSE,
        "deployment_cells": hllc.TARGET_CELLS,
        "safe_step_cost_failure_rule": (
            "fail a nominal safe update if it needs more than "
            f"{MAX_INTERNAL_SUBSTEPS_PER_SAFE_CALL} internal substeps; "
            "the matched expected count is one"
        ),
        "reference": (
            "strict native HLLC + SSP-RK2 on 2048 cells, conservatively "
            "restricted to 512 cells only for metrics"
        ),
        "mean_metrics": mean_metrics,
        "cases": cases,
    }
    stem = f"full_flux512_vs_dissipation512_with_hllc2048_seed{args.seed}"
    plot_profiles(trajectories, cases, results / f"{stem}.png")
    (results / f"{stem}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(results / f"{stem}.csv", cases)
    print(json.dumps(mean_metrics, indent=2))


if __name__ == "__main__":
    main()
