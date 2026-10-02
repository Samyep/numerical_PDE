"""Test the 64-cell-trained consistent central model on 512 cells."""

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
for module_path in (AUDIT_EXPERIMENT, FULL_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import plot_full_flux512 as full_plot  # noqa: E402
import plot_hllc512_vs_hcfl512 as hllc  # noqa: E402


comparison = hllc.comparison
base = hllc.base
shared = hllc.shared
REFERENCE_COLOR = comparison.BACKGROUND_REFERENCE_COLOR
HLLC_COLOR = "#0072B2"
DISSIPATION_COLOR = "#009E73"
CENTRAL_COLOR = "#D55E00"


def completed_mean(
    cases: dict[str, dict[str, Any]], method: str
) -> dict[str, Any]:
    completed = [
        name for name in cases if cases[name][method]["completed"]
    ]
    failed = [name for name in cases if name not in completed]
    return {
        "completed_cases": len(completed),
        "failed_cases": failed,
        "mean_rollout_nrmse": (
            float(np.mean([
                cases[name][method]["rollout_nrmse"] for name in completed
            ]))
            if completed
            else None
        ),
        "mean_final_snapshot_nrmse": (
            float(np.mean([
                cases[name][method]["final_snapshot_nrmse"]
                for name in completed
            ]))
            if completed
            else None
        ),
        "all_five_case_mean_is_valid": not failed,
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
        roe_error = cases[name]["roe_dissipation_hcfl_512"][
            "rollout_nrmse"
        ]
        central_metrics = cases[name]["central_consistent_hcfl_512"]
        if central_metrics["completed"]:
            ratio = central_metrics["rollout_nrmse"] / roe_error
            subtitle = f"central / Roe NRMSE = {ratio:.2f}x"
        else:
            subtitle = "central model: safe-step cost limit"
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
            if "central_512" in primitive:
                axis.plot(
                    x_target,
                    primitive["central_512"][:, row],
                    color=CENTRAL_COLOR,
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
            [0], [0], color=CENTRAL_COLOR, linewidth=1.35,
            linestyle="--", label="consistent-central HCFL-512",
        ),
    ]
    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"Consistent central flux transfer to 512 cells at t = {final_time:.4f}",
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
        "unchanged on 512 cells. Curves are native-grid values; HLLC-2048 is "
        "conservatively restricted only when computing NRMSE.",
        ha="right",
        fontsize=8.3,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    methods = (
        "native_hllc_512",
        "roe_dissipation_hcfl_512",
        "central_consistent_hcfl_512",
    )
    for name, case in cases.items():
        for method in methods:
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
    central = full_plot.load_solver(
        "central_consistent",
        results
        / f"central_consistent_broad_converged_best_seed{args.seed}.pt",
        args.width,
    )
    dissipation = full_plot.load_solver(
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
        central_512, central_stats = full_plot.guarded_matched_safe_rollout(
            central, name, hllc.TARGET_CELLS
        )

        trajectories[name] = {
            "reference": reference,
            "hllc_512": hllc_512,
            "dissipation_512": dissipation_512,
        }
        if central_stats["completed"]:
            trajectories[name]["central_512"] = central_512
            central_metrics = hllc.enrich_metrics(
                reference_512, central_512, state_std, central_stats
            )
        else:
            central_metrics = central_stats
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
            "central_consistent_hcfl_512": central_metrics,
        }

    methods = (
        "native_hllc_512",
        "roe_dissipation_hcfl_512",
        "central_consistent_hcfl_512",
    )
    mean_metrics = {
        method: completed_mean(cases, method) for method in methods
    }
    summary = {
        "seed": args.seed,
        "central_consistent_checkpoint": (
            f"central_consistent_broad_converged_best_seed{args.seed}.pt"
        ),
        "roe_dissipation_checkpoint": (
            f"dissipation_broad_converged_best_seed{args.seed}.pt"
        ),
        "checkpoint_training_cells": base.NCOARSE,
        "deployment_cells": hllc.TARGET_CELLS,
        "safe_step_cost_failure_rule": (
            "fail a nominal safe update if it needs more than "
            f"{full_plot.MAX_INTERNAL_SUBSTEPS_PER_SAFE_CALL} internal "
            "substeps; the matched expected count is one"
        ),
        "reference": (
            "strict native HLLC + SSP-RK2 on 2048 cells, conservatively "
            "restricted to 512 cells only for metrics"
        ),
        "mean_metrics": mean_metrics,
        "cases": cases,
    }
    stem = (
        f"central_consistent512_vs_hllc_and_roe_with_hllc2048_seed"
        f"{args.seed}"
    )
    plot_profiles(trajectories, cases, results / f"{stem}.png")
    (results / f"{stem}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(results / f"{stem}.csv", cases)
    print(json.dumps(mean_metrics, indent=2))


if __name__ == "__main__":
    main()
