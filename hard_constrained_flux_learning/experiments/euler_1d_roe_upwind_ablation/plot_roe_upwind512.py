"""Evaluate all Roe-coordinate/upwind ablation arms on 512 cells."""

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
NEW_MODELS = {
    "roe_complete_hcfl_512": (
        "roe_complete",
        "roe_complete_broad_converged_best_seed{seed}.pt",
    ),
    "central_roe_signed_hcfl_512": (
        "central_roe_signed",
        "central_roe_signed_broad_converged_best_seed{seed}.pt",
    ),
    "central_roe_upwind_hcfl_512": (
        "central_roe_upwind",
        "central_roe_upwind_broad_converged_best_seed{seed}.pt",
    ),
}
METHODS = (
    "native_hllc_512",
    "hllc_roe_correction_hcfl_512",
    *NEW_MODELS,
)
COLORS = {
    "reference": comparison.BACKGROUND_REFERENCE_COLOR,
    "native_hllc_512": "#0072B2",
    "hllc_roe_correction_hcfl_512": "#009E73",
    "roe_complete_hcfl_512": "#CC79A7",
    "central_roe_signed_hcfl_512": "#7B2CBF",
    "central_roe_upwind_hcfl_512": "#D55E00",
}
LABELS = {
    "reference": "native HLLC-2048 reference",
    "native_hllc_512": "native HLLC-512",
    "hllc_roe_correction_hcfl_512": "HLLC + signed Roe correction",
    "roe_complete_hcfl_512": "direct complete flux in Roe coordinates",
    "central_roe_signed_hcfl_512": "central + signed Roe multipliers",
    "central_roe_upwind_hcfl_512": "central + nonnegative auto-upwind",
}
LINESTYLES = {
    "native_hllc_512": "-",
    "hllc_roe_correction_hcfl_512": "-",
    "roe_complete_hcfl_512": ":",
    "central_roe_signed_hcfl_512": "--",
    "central_roe_upwind_hcfl_512": "-.",
}


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
        figsize=(19, 9.2),
        sharex=True,
    )
    figure.subplots_adjust(
        left=0.06,
        right=0.992,
        bottom=0.09,
        top=0.76,
        wspace=0.22,
        hspace=0.14,
    )

    for column, name in enumerate(comparison.CASES):
        primitive = {
            method: base.primitive(trajectory).numpy()[0, -1]
            for method, trajectory in trajectories[name].items()
        }
        baseline = cases[name]["hllc_roe_correction_hcfl_512"][
            "rollout_nrmse"
        ]
        signed = cases[name]["central_roe_signed_hcfl_512"]
        upwind = cases[name]["central_roe_upwind_hcfl_512"]
        if signed["completed"] and upwind["completed"]:
            subtitle = (
                "signed / upwind / old Roe = "
                f"{signed['rollout_nrmse'] / baseline:.2f} / "
                f"{upwind['rollout_nrmse'] / baseline:.2f} / 1.00"
            )
        else:
            subtitle = "one or more new arms hit the safe-step cost limit"
        axes[0, column].set_title(
            f"{comparison.DISPLAY[name]}\n{subtitle}",
            fontsize=9.4,
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
                alpha=0.7,
                zorder=1,
            )
            for zorder, method in enumerate(METHODS, start=2):
                if method not in primitive:
                    continue
                axis.plot(
                    x_target,
                    primitive[method][:, row],
                    color=COLORS[method],
                    linewidth=1.2 if method != "central_roe_upwind_hcfl_512" else 1.5,
                    linestyle=LINESTYLES[method],
                    alpha=0.92,
                    zorder=zorder,
                )
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    handles = [
        Line2D(
            [0], [0], color=COLORS["reference"], linewidth=1.7,
            alpha=0.7, label=LABELS["reference"],
        )
    ]
    handles.extend(
        Line2D(
            [0], [0], color=COLORS[method], linewidth=1.35,
            linestyle=LINESTYLES[method], label=LABELS[method],
        )
        for method in METHODS
    )
    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"Roe-coordinate and automatic-upwind transfer at t = {final_time:.4f}",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.525, 0.895),
        ncol=3,
        frameon=False,
        fontsize=8.8,
    )
    figure.text(
        0.99,
        0.018,
        "All neural checkpoints were trained on 64-cell trajectories and "
        "deployed unchanged on 512 cells. HLLC-2048 is restricted only for "
        "metrics; its plotted curve remains native-grid.",
        ha="right",
        fontsize=8.3,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    for name, case in cases.items():
        for method in METHODS:
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
                    "minimum_internal_dt": metrics.get(
                        "minimum_internal_dt"
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
    models = {
        method: full_plot.load_solver(
            model_name,
            results / pattern.format(seed=args.seed),
            args.width,
        )
        for method, (model_name, pattern) in NEW_MODELS.items()
    }
    existing_roe = full_plot.load_solver(
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
        old_roe_512, old_roe_stats = hllc.matched_safe_rollout(
            existing_roe, name, hllc.TARGET_CELLS
        )

        trajectories[name] = {
            "reference": reference,
            "native_hllc_512": hllc_512,
            "hllc_roe_correction_hcfl_512": old_roe_512,
        }
        case: dict[str, Any] = {
            "reference_hllc_2048": reference_stats,
            "native_hllc_512": hllc.enrich_metrics(
                reference_512, hllc_512, state_std, hllc_stats
            ),
            "hllc_roe_correction_hcfl_512": hllc.enrich_metrics(
                reference_512,
                old_roe_512,
                state_std,
                old_roe_stats,
            ),
        }
        for method, model in models.items():
            trajectory, stats = full_plot.guarded_matched_safe_rollout(
                model, name, hllc.TARGET_CELLS
            )
            if stats["completed"]:
                trajectories[name][method] = trajectory
                case[method] = hllc.enrich_metrics(
                    reference_512, trajectory, state_std, stats
                )
            else:
                case[method] = stats
        cases[name] = case

    mean_metrics = {
        method: completed_mean(cases, method) for method in METHODS
    }
    summary = {
        "seed": args.seed,
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
    stem = f"roe_upwind512_with_hllc2048_seed{args.seed}"
    plot_profiles(trajectories, cases, results / f"{stem}.png")
    (results / f"{stem}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(results / f"{stem}.csv", cases)
    print(json.dumps(mean_metrics, indent=2))


if __name__ == "__main__":
    main()
