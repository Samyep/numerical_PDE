"""Compare only Central + nonnegative Roe feasibility-loss variants at 512 cells."""

from __future__ import annotations

import argparse
import csv
import json
import sys
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
FULL_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_full_flux_baseline"
for module_path in (HERE, AUDIT_EXPERIMENT, FULL_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import plot_full_flux512 as full_plot  # noqa: E402
import plot_hllc512_vs_hcfl512 as hllc  # noqa: E402
import run_roe_upwind_ablation as training  # noqa: E402


comparison = hllc.comparison
base = hllc.base
METHODS = {
    "control": (
        "central_roe_upwind_broad_converged_best_seed{seed}.pt",
        0.0,
    ),
    "feas_1e4": (
        "central_roe_upwind_feas_light_broad_converged_best_seed{seed}.pt",
        1.0e-4,
    ),
    "feas_1e3": (
        "central_roe_upwind_feas_broad_converged_best_seed{seed}.pt",
        1.0e-3,
    ),
}
LABELS = {
    "reference": "native HLLC-2048 reference",
    "control": r"nonnegative Roe, $\lambda_{feas}=0$",
    "feas_1e4": r"nonnegative Roe, $\lambda_{feas}=10^{-4}$",
    "feas_1e3": r"nonnegative Roe, $\lambda_{feas}=10^{-3}$",
}
COLORS = {
    "reference": comparison.BACKGROUND_REFERENCE_COLOR,
    "control": "#0072B2",
    "feas_1e4": "#D55E00",
    "feas_1e3": "#7B2CBF",
}
LINESTYLES = {"control": "-", "feas_1e4": "--", "feas_1e3": ":"}
VARIABLES = ("density", "velocity", "pressure")


def circular_total_variation(values: np.ndarray) -> float:
    return float(np.abs(np.roll(values, -1) - values).sum())


def significant_extrema(values: np.ndarray, scale: float) -> int:
    left = values - np.roll(values, 1)
    right = np.roll(values, -1) - values
    threshold = 1.0e-3 * max(scale, 1.0e-8)
    return int(
        (
            (left * right < 0.0)
            & (np.minimum(np.abs(left), np.abs(right)) > threshold)
        ).sum()
    )


def oscillation_diagnostics(
    reference: torch.Tensor,
    candidate: torch.Tensor,
) -> dict[str, Any]:
    reference_primitive = base.primitive(reference).numpy()[0, -1]
    candidate_primitive = base.primitive(candidate).numpy()[0, -1]
    variables: dict[str, dict[str, float | int]] = {}
    for index, name in enumerate(VARIABLES):
        truth = reference_primitive[:, index]
        prediction = candidate_primitive[:, index]
        scale = max(float(np.ptp(truth)), 1.0e-8)
        reference_tv = circular_total_variation(truth)
        candidate_tv = circular_total_variation(prediction)
        range_violation = (
            max(float(truth.min() - prediction.min()), 0.0)
            + max(float(prediction.max() - truth.max()), 0.0)
        ) / scale
        variables[name] = {
            "reference_total_variation": reference_tv,
            "candidate_total_variation": candidate_tv,
            "total_variation_ratio": candidate_tv / max(reference_tv, 1.0e-10),
            "normalized_global_range_violation": range_violation,
            "reference_significant_extrema": significant_extrema(truth, scale),
            "candidate_significant_extrema": significant_extrema(
                prediction, scale
            ),
            "excess_significant_extrema": max(
                significant_extrema(prediction, scale)
                - significant_extrema(truth, scale),
                0,
            ),
        }
    return {"variables": variables}


def aggregate(
    cases: dict[str, dict[str, Any]], method: str
) -> dict[str, Any]:
    completed = [
        name for name in cases if cases[name][method]["completed"]
    ]
    failed = [name for name in cases if name not in completed]
    profiles = [
        values
        for name in completed
        for values in cases[name][method]["oscillation"]["variables"].values()
    ]
    metrics = [cases[name][method] for name in completed]
    return {
        "completed_cases": len(completed),
        "failed_cases": failed,
        "all_five_case_mean_is_valid": not failed,
        "mean_rollout_nrmse": float(
            np.mean([row["rollout_nrmse"] for row in metrics])
        ),
        "mean_final_snapshot_nrmse": float(
            np.mean([row["final_snapshot_nrmse"] for row in metrics])
        ),
        "mean_hard_projection_intervention_rate": float(
            np.mean([
                row["hard_projection_intervention_rate"] for row in metrics
            ])
        ),
        "mean_local_limiter_intervention_rate": float(
            np.mean([
                row["local_limiter_intervention_rate"] for row in metrics
            ])
        ),
        "mean_fd_entropy_intervention_rate": float(
            np.mean([
                row["fd_entropy_intervention_rate"] for row in metrics
            ])
        ),
        "minimum_density": float(
            min(row["minimum_density"] for row in metrics)
        ),
        "minimum_pressure": float(
            min(row["minimum_pressure"] for row in metrics)
        ),
        "mean_final_total_variation_ratio": float(
            np.mean([row["total_variation_ratio"] for row in profiles])
        ),
        "maximum_final_total_variation_ratio": float(
            max(row["total_variation_ratio"] for row in profiles)
        ),
        "mean_normalized_global_range_violation": float(
            np.mean([
                row["normalized_global_range_violation"] for row in profiles
            ])
        ),
        "maximum_normalized_global_range_violation": float(
            max(row["normalized_global_range_violation"] for row in profiles)
        ),
        "total_excess_significant_extrema": int(
            sum(row["excess_significant_extrema"] for row in profiles)
        ),
    }


def plot_profiles(
    trajectories: dict[str, dict[str, torch.Tensor]],
    cases: dict[str, dict[str, Any]],
    output: Path,
) -> None:
    x_reference = (
        np.arange(hllc.REFERENCE_CELLS) + 0.5
    ) / hllc.REFERENCE_CELLS
    x_target = (np.arange(hllc.TARGET_CELLS) + 0.5) / hllc.TARGET_CELLS
    figure, axes = plt.subplots(
        3, len(comparison.CASES), figsize=(19, 9.2), sharex=True
    )
    figure.subplots_adjust(
        left=0.06, right=0.992, bottom=0.09, top=0.77,
        wspace=0.22, hspace=0.14,
    )

    for column, case_name in enumerate(comparison.CASES):
        primitive = {
            method: base.primitive(trajectory).numpy()[0, -1]
            for method, trajectory in trajectories[case_name].items()
        }
        control = cases[case_name]["control"]["rollout_nrmse"]
        ratios = [
            cases[case_name][method]["rollout_nrmse"] / control
            for method in ("feas_1e4", "feas_1e3")
        ]
        axes[0, column].set_title(
            f"{comparison.DISPLAY[case_name]}\n"
            f"1e-4 / 1e-3 vs control = {ratios[0]:.2f} / {ratios[1]:.2f}",
            fontsize=9.3,
            fontweight="semibold",
            color=comparison.TEXT_COLOR,
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                primitive["reference"][:, row],
                color=COLORS["reference"],
                linewidth=1.7,
                alpha=0.72,
                zorder=1,
            )
            for zorder, method in enumerate(METHODS, start=2):
                axis.plot(
                    x_target,
                    primitive[method][:, row],
                    color=COLORS[method],
                    linestyle=LINESTYLES[method],
                    linewidth=1.35,
                    alpha=0.94,
                    zorder=zorder,
                )
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(VARIABLES[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    handles = [
        Line2D(
            [0], [0], color=COLORS["reference"], linewidth=1.7,
            alpha=0.72, label=LABELS["reference"],
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
    final_time = (hllc.shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        "Central + nonnegative Roe: raw-proposal feasibility learning "
        f"at t = {final_time:.4f}",
        y=0.978,
        fontsize=16.5,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    figure.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.525, 0.905),
        ncol=4, frameon=False, fontsize=9,
    )
    figure.text(
        0.99,
        0.018,
        "All three models have the same architecture and hard-projected PDE "
        "update; only the raw-proposal training penalty differs.",
        ha="right",
        fontsize=8.5,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    for case_name, case in cases.items():
        for method in METHODS:
            metrics = case[method]
            rows.append({
                "case": case_name,
                "method": method,
                "feasibility_weight": METHODS[method][1],
                "rollout_nrmse": metrics["rollout_nrmse"],
                "final_snapshot_nrmse": metrics["final_snapshot_nrmse"],
                "hard_projection_intervention_rate": metrics[
                    "hard_projection_intervention_rate"
                ],
                "local_limiter_intervention_rate": metrics[
                    "local_limiter_intervention_rate"
                ],
                "fd_entropy_intervention_rate": metrics[
                    "fd_entropy_intervention_rate"
                ],
                "minimum_density": metrics["minimum_density"],
                "minimum_pressure": metrics["minimum_pressure"],
            })
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
    results = args.results_dir.resolve()
    results.mkdir(parents=True, exist_ok=True)

    _, _, state_std = hllc.training_statistics(args.seed)
    models = {
        method: full_plot.load_solver(
            "central_roe_upwind",
            results / pattern.format(seed=args.seed),
            args.width,
        )
        for method, (pattern, _weight) in METHODS.items()
    }
    validation_data = training.audit.make_validation_data(args.seed)

    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    cases: dict[str, dict[str, Any]] = {}
    for case_name in comparison.CASES:
        print(f"Evaluating {comparison.DISPLAY[case_name]}...", flush=True)
        reference, reference_stats = hllc.strict_native_hllc_rollout(
            case_name, hllc.REFERENCE_CELLS
        )
        restricted = comparison.precision.conservative_restrict(
            reference.numpy(), hllc.TARGET_CELLS
        )
        reference_512 = torch.from_numpy(restricted.astype(np.float32))
        trajectories[case_name] = {"reference": reference}
        case: dict[str, Any] = {"reference_hllc_2048": reference_stats}
        for method, model in models.items():
            trajectory, stats = full_plot.guarded_matched_safe_rollout(
                model, case_name, hllc.TARGET_CELLS
            )
            if not stats["completed"]:
                raise RuntimeError(
                    f"{method} failed {case_name}: {stats['failure_reason']}"
                )
            trajectories[case_name][method] = trajectory
            metrics = hllc.enrich_metrics(
                reference_512, trajectory, state_std, stats
            )
            metrics["oscillation"] = oscillation_diagnostics(
                reference_512, trajectory
            )
            case[method] = metrics
        cases[case_name] = case

    mean_metrics = {method: aggregate(cases, method) for method in METHODS}
    validation_reports: dict[str, Any] = {}
    report_names = {
        "control": "report_central_roe_upwind_broad_seed{seed}.json",
        "feas_1e4": (
            "report_central_roe_upwind_feas_light_broad_seed{seed}.json"
        ),
        "feas_1e3": "report_central_roe_upwind_feas_broad_seed{seed}.json",
    }
    for method, pattern in report_names.items():
        report = json.loads(
            (results / pattern.format(seed=args.seed)).read_text(
                encoding="utf-8"
            )
        )
        validation_reports[method] = {
            "best_validation_rollout_nrmse": report["convergence"][
                "best_validation_rollout_nrmse"
            ],
            "proposal_diagnostics": report[
                "proposal_diagnostics_on_validation_ground_truth_states"
            ],
            "raw_feasibility": training.proposal_feasibility_diagnostics(
                models[method], validation_data
            ),
        }

    output = {
        "seed": args.seed,
        "model": "Central + nonnegative Roe",
        "training_cells": base.NCOARSE,
        "deployment_cells": hllc.TARGET_CELLS,
        "reference": (
            "strict native HLLC-2048, conservatively restricted to 512 "
            "only for metrics"
        ),
        "validation": validation_reports,
        "mean_metrics": mean_metrics,
        "cases": cases,
    }
    stem = f"upwind_feasibility512_seed{args.seed}"
    (results / f"{stem}.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )
    write_csv(results / f"{stem}.csv", cases)
    plot_profiles(trajectories, cases, results / f"{stem}.png")
    print(json.dumps(mean_metrics, indent=2))


if __name__ == "__main__":
    main()
