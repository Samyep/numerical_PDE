"""Plot the validation-converged unconstrained DNN flux baseline."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


HERE = Path(__file__).resolve().parent
AUDIT_RESULTS = (
    HERE.parent / "euler_1d_convergence_audit" / "results"
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
RED = "#B23A48"
GRAY = "#8A96A3"
GRID = "#D8DEE9"
TEXT = "#1F2933"
MUTED = "#66788A"

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
        "--results-dir",
        type=Path,
        default=HERE / "results",
    )
    args = parser.parse_args()
    results = args.results_dir.resolve()

    raw_curve = read_csv(results / f"training_curve_seed{args.seed}.csv")
    constrained_curve = [
        row
        for row in read_csv(
            AUDIT_RESULTS / f"training_curve_seed{args.seed}.csv"
        )
        if row["arm"] == "direct_broad"
    ]
    comparison = read_csv(
        results / f"comparison_vs_constrained_direct_seed{args.seed}.csv"
    )
    raw_metrics = [
        row
        for row in read_csv(results / f"metrics_seed{args.seed}.csv")
        if row["stage"] == "converged_best"
    ]

    figure = plt.figure(figsize=(14.5, 9.2), constrained_layout=False)
    grid = figure.add_gridspec(
        2,
        2,
        height_ratios=(1.18, 1.0),
        width_ratios=(1.08, 1.0),
        left=0.07,
        right=0.985,
        bottom=0.10,
        top=0.86,
        hspace=0.42,
        wspace=0.25,
    )
    validation_axis = figure.add_subplot(grid[0, :])
    change_axis = figure.add_subplot(grid[1, 0])
    completion_axis = figure.add_subplot(grid[1, 1])

    constrained_updates = np.array(
        [int(row["update"]) for row in constrained_curve]
    )
    constrained_values = np.array(
        [float(row["validation_rollout_nrmse"]) for row in constrained_curve]
    )
    validation_axis.plot(
        constrained_updates,
        constrained_values,
        color=BLUE,
        linewidth=2.0,
        label="Constrained direct HCFL",
        zorder=3,
    )

    raw_updates = np.array([int(row["update"]) for row in raw_curve])
    raw_partial = np.array(
        [float(row["validation_partial_rollout_nrmse"]) for row in raw_curve]
    )
    raw_complete = np.array(
        [float(row["validation_completion_rate"]) == 1.0 for row in raw_curve]
    )
    validation_axis.plot(
        raw_updates,
        raw_partial,
        color=GRAY,
        linewidth=1.3,
        linestyle="--",
        label="Unconstrained partial NRMSE",
        zorder=2,
    )
    validation_axis.scatter(
        raw_updates[raw_complete],
        raw_partial[raw_complete],
        color=ORANGE,
        s=22,
        label="Unconstrained: all validation trajectories complete",
        zorder=4,
    )
    validation_axis.scatter(
        raw_updates[~raw_complete],
        raw_partial[~raw_complete],
        color=RED,
        marker="x",
        s=30,
        linewidths=1.4,
        label="Unconstrained: at least one validation failure",
        zorder=5,
    )
    best_row = min(
        (row for row in raw_curve if float(row["validation_completion_rate"]) == 1.0),
        key=lambda row: float(row["validation_rollout_nrmse"]),
    )
    validation_axis.scatter(
        [int(best_row["update"])],
        [float(best_row["validation_rollout_nrmse"])],
        color=ORANGE,
        edgecolor="white",
        marker="*",
        s=180,
        linewidth=0.8,
        zorder=6,
    )
    validation_axis.annotate(
        "selected: 6,500",
        (
            int(best_row["update"]),
            float(best_row["validation_rollout_nrmse"]),
        ),
        xytext=(12, 12),
        textcoords="offset points",
        fontsize=9,
        color=ORANGE,
    )
    validation_axis.set_title(
        "A. Validation convergence and raw-rollout feasibility",
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=TEXT,
    )
    validation_axis.set_xlabel("optimizer updates", fontsize=10)
    validation_axis.set_ylabel("rollout NRMSE", fontsize=10)
    validation_axis.legend(
        frameon=False,
        fontsize=9,
        ncol=2,
        loc="upper right",
    )
    style_axis(validation_axis)

    successful = [
        row for row in comparison if row["candidate_rollout_nrmse"]
    ]
    names = [DISPLAY[row["split"]] for row in successful]
    changes = np.array(
        [float(row["relative_change_percent"]) for row in successful]
    )
    colors = np.where(changes <= 0.0, BLUE, ORANGE)
    positions = np.arange(len(successful))
    change_axis.bar(
        positions,
        changes,
        color=colors,
        width=0.68,
        zorder=3,
    )
    change_axis.axhline(0.0, color="#67727E", linewidth=0.9, zorder=2)
    for position, value in zip(positions, changes, strict=True):
        change_axis.text(
            position,
            value + (0.08 if value >= 0.0 else -0.08),
            f"{value:+.2f}%",
            ha="center",
            va="bottom" if value >= 0.0 else "top",
            fontsize=8.5,
            color=TEXT,
        )
    change_axis.set_xticks(positions, names, rotation=25, ha="right")
    change_axis.set_ylim(-1.55, 1.45)
    change_axis.set_ylabel("NRMSE change vs constrained direct (%)", fontsize=9.5)
    change_axis.set_title(
        "B. Accuracy changes are all within 1.3%",
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=TEXT,
    )
    style_axis(change_axis)

    canonical_names = [
        "sod",
        "lax",
        "collision",
        "strong_pressure",
        "near_vacuum_expansion",
    ]
    indexed_metrics = {row["split"]: row for row in raw_metrics}
    y = np.arange(len(canonical_names))
    constrained_end = np.full(len(canonical_names), 63.0)
    raw_end = []
    failed = []
    for name in canonical_names:
        row = indexed_metrics[name]
        failure = row["first_failure_snapshot"]
        if failure:
            raw_end.append(float(int(failure) - 1))
            failed.append(True)
        else:
            raw_end.append(63.0)
            failed.append(False)
    raw_end_values = np.array(raw_end)
    completion_axis.barh(
        y - 0.16,
        constrained_end,
        height=0.28,
        color=BLUE,
        label="Constrained direct HCFL",
        zorder=2,
    )
    completion_axis.barh(
        y + 0.16,
        raw_end_values,
        height=0.28,
        color=[RED if value else ORANGE for value in failed],
        label="Unconstrained DNN flux",
        zorder=3,
    )
    for index, (name, is_failed) in enumerate(
        zip(canonical_names, failed, strict=True)
    ):
        if not is_failed:
            continue
        failure = int(indexed_metrics[name]["first_failure_snapshot"])
        completion_axis.scatter(
            [failure],
            [index + 0.16],
            color=RED,
            marker="x",
            s=48,
            linewidths=1.6,
            zorder=5,
        )
        completion_axis.text(
            failure + 1.5,
            index + 0.16,
            f"negative p at {failure}",
            va="center",
            fontsize=8.5,
            color=RED,
        )
    completion_axis.set_yticks(
        y,
        [DISPLAY[name] for name in canonical_names],
    )
    completion_axis.invert_yaxis()
    completion_axis.set_xlim(0.0, 66.0)
    completion_axis.set_xticks([0, 16, 32, 48, 63])
    completion_axis.set_xlabel("last admissible saved snapshot (of 63)", fontsize=9.5)
    completion_axis.set_title(
        "C. Hard constraints prevent severe-rollout failure",
        loc="left",
        fontsize=12,
        fontweight="semibold",
        color=TEXT,
    )
    completion_axis.legend(
        handles=[
            Patch(color=BLUE, label="Constrained direct HCFL"),
            Patch(color=ORANGE, label="Unconstrained: completed"),
            Patch(color=RED, label="Unconstrained: failed"),
        ],
        frameon=False,
        fontsize=9,
        loc="lower right",
    )
    style_axis(completion_axis)

    figure.suptitle(
        "Unconstrained DNN flux: negligible accuracy gain, incomplete severe rollouts",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
        y=0.965,
    )
    figure.text(
        0.985,
        0.025,
        "Seed 0; same DirectFlux proposal, data, optimizer, and convergence rule. "
        "Only Tadmor/admissibility/fully-discrete output safeguards are removed.",
        ha="right",
        fontsize=8.5,
        color=MUTED,
    )
    output = results / f"unconstrained_flux_baseline_seed{args.seed}.png"
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
