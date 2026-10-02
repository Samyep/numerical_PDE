"""Plot convergence and 64-cell accuracy of the consistent central-flux arm."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_RESULTS = HERE.parent / "euler_1d_convergence_audit" / "results"
FULL_RESULTS = HERE.parent / "euler_1d_full_flux_baseline" / "results"
ARMS = (
    "direct_broad",
    "dissipation_broad",
    "full_flux_broad",
    "central_consistent_broad",
)
COLORS = {
    "direct_broad": "#0072B2",
    "dissipation_broad": "#009E73",
    "full_flux_broad": "#CC79A7",
    "central_consistent_broad": "#D55E00",
}
LABELS = {
    "direct_broad": "HLLC + direct vector correction",
    "dissipation_broad": "HLLC + Roe-dissipation correction",
    "full_flux_broad": "direct complete flux",
    "central_consistent_broad": "central physical flux + consistent correction",
}
DISPLAY = {
    "ordinary_id": "Ordinary ID",
    "broad_random_in_support": "Broad in-support",
    "moderate_ood_high_frequency": "Moderate OOD",
    "sod": "Sod",
    "lax": "Lax",
    "collision": "Collision",
    "strong_pressure": "Strong pressure",
    "near_vacuum_expansion": "Near-vacuum",
}
GRID = "#D8DEE9"
TEXT = "#1F2933"
MUTED = "#66788A"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def style_axis(axis: plt.Axes) -> None:
    axis.grid(True, color=GRID, linewidth=0.7, alpha=0.65, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#9AA5B1")
    axis.spines["bottom"].set_color("#9AA5B1")
    axis.tick_params(colors="#465362", labelsize=9)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--results-dir", type=Path, default=HERE / "results"
    )
    args = parser.parse_args()
    results = args.results_dir.resolve()

    audit_curve = read_csv(
        AUDIT_RESULTS / f"training_curve_seed{args.seed}.csv"
    )
    curves = {
        "direct_broad": [
            row for row in audit_curve if row["arm"] == "direct_broad"
        ],
        "dissipation_broad": [
            row
            for row in audit_curve
            if row["arm"] == "dissipation_broad"
        ],
        "full_flux_broad": read_csv(
            FULL_RESULTS / f"training_curve_seed{args.seed}.csv"
        ),
        "central_consistent_broad": read_csv(
            results / f"training_curve_seed{args.seed}.csv"
        ),
    }

    audit_metrics = read_csv(AUDIT_RESULTS / f"metrics_seed{args.seed}.csv")
    metrics = {
        "direct_broad": {
            row["split"]: row
            for row in audit_metrics
            if row["arm"] == "direct_broad"
            and row["stage"] == "converged_best"
        },
        "dissipation_broad": {
            row["split"]: row
            for row in audit_metrics
            if row["arm"] == "dissipation_broad"
            and row["stage"] == "converged_best"
        },
        "full_flux_broad": {
            row["split"]: row
            for row in read_csv(FULL_RESULTS / f"metrics_seed{args.seed}.csv")
        },
        "central_consistent_broad": {
            row["split"]: row
            for row in read_csv(results / f"metrics_seed{args.seed}.csv")
        },
    }
    report = json.loads(
        (results / f"report_seed{args.seed}.json").read_text(
            encoding="utf-8"
        )
    )

    figure, (curve_axis, metric_axis) = plt.subplots(
        2,
        1,
        figsize=(14.8, 9.2),
        gridspec_kw={"height_ratios": (1.0, 1.15)},
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.16,
        top=0.88,
        hspace=0.58,
    )

    for arm in ARMS:
        updates = np.array([int(row["update"]) for row in curves[arm]])
        values = np.array([
            float(row["validation_rollout_nrmse"]) for row in curves[arm]
        ])
        curve_axis.plot(
            updates,
            values,
            color=COLORS[arm],
            linewidth=1.65,
            label=LABELS[arm],
            zorder=3,
        )
        best_index = int(np.argmin(values))
        curve_axis.scatter(
            [updates[best_index]],
            [values[best_index]],
            color=COLORS[arm],
            edgecolor="white",
            marker="*",
            s=125,
            linewidth=0.7,
            zorder=4,
        )
    curve_axis.set_yscale("log")
    curve_axis.set_xlabel("optimizer updates", fontsize=10)
    curve_axis.set_ylabel("independent validation rollout NRMSE", fontsize=10)
    curve_axis.set_title(
        "A. Validation convergence (stars are selected checkpoints)",
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=TEXT,
    )
    curve_axis.legend(frameon=False, fontsize=8.7, ncol=2)
    style_axis(curve_axis)

    splits = list(DISPLAY)
    positions = np.arange(len(splits))
    width = 0.19
    offsets = (-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width)
    for offset, arm in zip(offsets, ARMS, strict=True):
        values = [
            float(metrics[arm][split]["rollout_nrmse"])
            for split in splits
        ]
        metric_axis.bar(
            positions + offset,
            values,
            width=width,
            color=COLORS[arm],
            label=LABELS[arm],
            zorder=3,
        )
    metric_axis.axvline(2.5, color="#8793A0", linewidth=1.0)
    metric_axis.text(
        1.0,
        1.03,
        "random-distribution tests",
        transform=metric_axis.get_xaxis_transform(),
        ha="center",
        fontsize=9,
        color=MUTED,
    )
    metric_axis.text(
        5.0,
        1.03,
        "canonical Riemann tests (not used for selection)",
        transform=metric_axis.get_xaxis_transform(),
        ha="center",
        fontsize=9,
        color=MUTED,
    )
    metric_axis.set_xticks(
        positions,
        [DISPLAY[split] for split in splits],
        rotation=24,
        ha="right",
    )
    metric_axis.set_ylabel("64-cell rollout NRMSE (lower is better)", fontsize=10)
    metric_axis.set_title(
        "B. Exact consistency helps some shocks, but is not uniformly best",
        loc="left",
        y=1.12,
        fontsize=12,
        fontweight="semibold",
        color=TEXT,
    )
    style_axis(metric_axis)

    proposal = report[
        "proposal_diagnostics_on_validation_ground_truth_states"
    ]
    correction = report["learned_correction_diagnostics"]
    consistency = report[
        "equal_interface_with_varied_outer_stencil_consistency"
    ]
    figure.suptitle(
        "Exactly consistent central flux + learned correction",
        y=0.97,
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    figure.text(
        0.985,
        0.025,
        "All neural models: 6,627 parameters, identical data and validation "
        "stopping. Equal-interface error: "
        f"{consistency['maximum_projected_absolute_error']:.1e}; correction "
        f"RMS / central-flux RMS: "
        f"{100.0 * correction['correction_relative_flux_rms']:.2f}%; "
        "hard projection active on "
        f"{100.0 * proposal['hard_projection_intervention_rate']:.2f}% of "
        "validation interfaces.",
        ha="right",
        fontsize=8.35,
        color=MUTED,
    )
    output = results / f"central_consistent_64cell_comparison_seed{args.seed}.png"
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
