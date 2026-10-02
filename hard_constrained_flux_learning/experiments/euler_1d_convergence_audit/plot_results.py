"""Plot the seed-0 validation-convergence audit.

The figure deliberately separates the quantity used for checkpoint selection
(validation rollout NRMSE) from held-out test changes.  Negative percentages
in the heatmap mean that convergence improved the test metric.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm


ARM_ORDER = [
    "direct_broad",
    "invariant_broad",
    "characteristic_broad",
    "dissipation_broad",
    "conv_broad",
    "direct_wave",
]

ARM_LABELS = {
    "direct_broad": "Direct / broad",
    "invariant_broad": "Invariant / broad",
    "characteristic_broad": "Characteristic / broad",
    "dissipation_broad": "Dissipation / broad",
    "conv_broad": "CNN / broad",
    "direct_wave": "Direct / wave",
}

SPLIT_ORDER = [
    "ordinary_id",
    "broad_random_in_support",
    "moderate_ood_high_frequency",
    "sod",
    "lax",
    "collision",
    "strong_pressure",
    "near_vacuum_expansion",
]

SPLIT_LABELS = {
    "ordinary_id": "Ordinary ID",
    "broad_random_in_support": "Broad random",
    "moderate_ood_high_frequency": "Moderate OOD",
    "sod": "Sod",
    "lax": "Lax",
    "collision": "Collision",
    "strong_pressure": "Strong pressure",
    "near_vacuum_expansion": "Near vacuum",
}

# Okabe-Ito-derived, color-vision-friendly ordering.
COLORS = ["#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "results",
    )
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    curve_rows = read_csv(results_dir / f"training_curve_seed{args.seed}.csv")
    convergence_rows = read_csv(results_dir / f"convergence_seed{args.seed}.csv")
    metric_rows = read_csv(results_dir / f"metrics_seed{args.seed}.csv")

    convergence = {row["arm"]: row for row in convergence_rows}
    metrics: dict[tuple[str, str, str], float] = {}
    for row in metric_rows:
        metrics[(row["arm"], row["stage"], row["split"])] = float(
            row["rollout_nrmse"]
        )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(13.5, 10.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0])

    ax_curve = fig.add_subplot(grid[0, 0])
    for arm, color in zip(ARM_ORDER, COLORS, strict=True):
        rows = sorted(
            (row for row in curve_rows if row["arm"] == arm),
            key=lambda row: int(row["update"]),
        )
        updates = np.asarray([int(row["update"]) for row in rows])
        values = np.asarray([float(row["validation_rollout_nrmse"]) for row in rows])
        ax_curve.plot(
            updates,
            values,
            color=color,
            linewidth=1.8,
            label=ARM_LABELS[arm],
        )
        best_update = int(convergence[arm]["best_update"])
        best_value = float(convergence[arm]["best_validation_rollout_nrmse"])
        ax_curve.scatter(
            [best_update],
            [best_value],
            marker="*",
            s=90,
            color=color,
            edgecolor="black",
            linewidth=0.45,
            zorder=4,
        )

    ax_curve.axvline(1100, color="0.25", linestyle="--", linewidth=1.2)
    ax_curve.text(
        1100,
        ax_curve.get_ylim()[1],
        "  old fixed budget",
        ha="left",
        va="top",
        fontsize=9,
        color="0.25",
    )
    ax_curve.set_title("A. Validation rollout error until convergence (stars: selected checkpoints)")
    ax_curve.set_xlabel("Adam updates")
    ax_curve.set_ylabel("Validation rollout NRMSE")
    ax_curve.legend(ncol=3, frameon=False, fontsize=9)
    ax_curve.grid(True, color="0.88", linewidth=0.7)

    changes = np.empty((len(ARM_ORDER), len(SPLIT_ORDER)), dtype=float)
    for row_index, arm in enumerate(ARM_ORDER):
        for column_index, split in enumerate(SPLIT_ORDER):
            fixed = metrics[(arm, "fixed_1100", split)]
            converged = metrics[(arm, "converged_best", split)]
            changes[row_index, column_index] = 100.0 * (converged / fixed - 1.0)

    ax_heatmap = fig.add_subplot(grid[1, 0])
    norm = TwoSlopeNorm(vmin=-50.0, vcenter=0.0, vmax=90.0)
    image = ax_heatmap.imshow(changes, aspect="auto", cmap="RdBu_r", norm=norm)
    ax_heatmap.set_title(
        "B. Held-out test change from 1,100 updates to validation-selected convergence"
    )
    ax_heatmap.set_xticks(range(len(SPLIT_ORDER)))
    ax_heatmap.set_xticklabels(
        [SPLIT_LABELS[split] for split in SPLIT_ORDER], rotation=25, ha="right"
    )
    ax_heatmap.set_yticks(range(len(ARM_ORDER)))
    ax_heatmap.set_yticklabels([ARM_LABELS[arm] for arm in ARM_ORDER])
    ax_heatmap.set_xlabel("Final test split (not used for checkpoint selection)")
    ax_heatmap.set_ylabel("Training arm")
    ax_heatmap.grid(False)

    for row_index in range(changes.shape[0]):
        for column_index in range(changes.shape[1]):
            value = changes[row_index, column_index]
            text_color = "white" if value < -32.0 or value > 48.0 else "black"
            ax_heatmap.text(
                column_index,
                row_index,
                f"{value:+.1f}%",
                ha="center",
                va="center",
                fontsize=8.2,
                color=text_color,
            )

    colorbar = fig.colorbar(image, ax=ax_heatmap, fraction=0.025, pad=0.02)
    colorbar.set_label("NRMSE change (%) — negative is better")

    output = results_dir / f"convergence_audit_seed{args.seed}.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
