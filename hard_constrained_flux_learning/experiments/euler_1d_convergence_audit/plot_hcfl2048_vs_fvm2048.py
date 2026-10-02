"""Compare HCFL and FVM directly on the same 2048-cell grid.

The learned flux is the validation-selected checkpoint trained on 64-cell
states.  It is deployed without retraining on 2048 cells, using 32 updates per
saved interval so that its update ratio dt/dx matches the training grid.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

import plot_best_vs_fvm as comparison


base = comparison.base
shared = comparison.shared
TARGET_CELLS = comparison.precision.HIGH_REFERENCE_CELLS
SELECTED_ARM = "dissipation_broad"
VARIABLE_NAMES = ("density", "velocity", "pressure")


def final_difference_metrics(
    reference: torch.Tensor,
    prediction: torch.Tensor,
    state_std: torch.Tensor,
) -> dict[str, object]:
    diagnostics = comparison.diagnostics(reference, prediction, state_std)
    primitive_reference = base.primitive(reference)[0, -1]
    primitive_prediction = base.primitive(prediction)[0, -1]
    primitive_delta = primitive_prediction - primitive_reference
    primitive_absolute_error = torch.abs(primitive_delta)
    final_normalized_error = (
        prediction[0, -1] - reference[0, -1]
    ) / state_std
    per_cell_nrmse = torch.sqrt(torch.mean(final_normalized_error**2, dim=-1))
    worst_cell = int(torch.argmax(per_cell_nrmse))

    return {
        **diagnostics,
        "final_primitive_mae": {
            name: float(primitive_absolute_error[:, index].mean())
            for index, name in enumerate(VARIABLE_NAMES)
        },
        "final_primitive_rmse": {
            name: float(torch.sqrt(torch.mean(primitive_delta[:, index] ** 2)))
            for index, name in enumerate(VARIABLE_NAMES)
        },
        "worst_final_cell": {
            "index": worst_cell,
            "x": (worst_cell + 0.5) / TARGET_CELLS,
            "normalized_rmse": float(per_cell_nrmse[worst_cell]),
        },
    }


def plot_final_profiles(
    references: dict[str, torch.Tensor],
    predictions: dict[str, torch.Tensor],
    metrics: dict[str, dict[str, object]],
    output: Path,
) -> None:
    x = (np.arange(TARGET_CELLS) + 0.5) / TARGET_CELLS
    row_labels = [r"Density $\rho$", r"Velocity $u$", r"Pressure $p$"]
    figure, axes = plt.subplots(
        3,
        len(comparison.CASES),
        figsize=(18, 8.8),
        sharex=True,
        constrained_layout=False,
    )
    figure.subplots_adjust(
        left=0.065,
        right=0.99,
        bottom=0.09,
        top=0.80,
        wspace=0.22,
        hspace=0.14,
    )

    for column, name in enumerate(comparison.CASES):
        reference = base.primitive(references[name]).numpy()[0, -1]
        prediction = base.primitive(predictions[name]).numpy()[0, -1]
        axes[0, column].set_title(
            f"{comparison.DISPLAY[name]}\n"
            f"final NRMSE = {metrics[name]['final_snapshot_nrmse']:.4f}",
            fontsize=10.5,
            fontweight="semibold",
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x,
                reference[:, row],
                color=comparison.FVM_COLOR,
                linewidth=1.8,
                label="FVM-2048",
                zorder=2,
            )
            axis.plot(
                x,
                prediction[:, row],
                color=comparison.HCFL_COLOR,
                linewidth=1.45,
                linestyle="--",
                label="HCFL-2048 (64-grid-trained checkpoint)",
                zorder=3,
            )
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"HCFL-2048 versus FVM-2048 at t = {final_time:.4f}",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.89),
        ncol=2,
        frameon=False,
        fontsize=10,
    )
    figure.text(
        0.99,
        0.018,
        "Same 2048 cell centers; no averaging or point subsampling. HCFL uses "
        "the 64-grid-trained checkpoint with 32 CFL-matched updates per saved "
        "interval.",
        ha="right",
        fontsize=8.5,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=comparison.HERE / "results",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    state_std = comparison.training_state_std(args.seed)
    model = comparison.load_model(
        results_dir,
        SELECTED_ARM,
        args.seed,
        args.width,
    )

    references: dict[str, torch.Tensor] = {}
    predictions: dict[str, torch.Tensor] = {}
    metrics: dict[str, dict[str, object]] = {}
    for name in comparison.CASES:
        references[name] = torch.from_numpy(
            comparison.strict_native_rollout(
                comparison.initial_condition(name, TARGET_CELLS),
                TARGET_CELLS,
            ).astype(np.float32)
        )
        predictions[name] = comparison.hcfl_rollout_on_grid(
            model,
            name,
            TARGET_CELLS,
        )
        metrics[name] = final_difference_metrics(
            references[name],
            predictions[name],
            state_std,
        )

    mean_final_nrmse = float(
        np.mean([row["final_snapshot_nrmse"] for row in metrics.values()])
    )
    worst_case = max(
        metrics,
        key=lambda name: float(metrics[name]["final_snapshot_nrmse"]),
    )
    summary = {
        "seed": args.seed,
        "selected_model": f"{SELECTED_ARM}_converged_best",
        "training_cells": base.NCOARSE,
        "deployment_cells": TARGET_CELLS,
        "fvm_reference_cells": TARGET_CELLS,
        "cfl_matched_updates_per_saved_interval": (
            TARGET_CELLS // base.NCOARSE
        ),
        "final_time": (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT,
        "comparison": (
            "same-grid cell-centered values; no averaging or subsampling"
        ),
        "mean_final_snapshot_nrmse": mean_final_nrmse,
        "worst_case_by_final_snapshot_nrmse": worst_case,
        "cases": metrics,
    }

    plot_path = results_dir / f"hcfl2048_vs_fvm2048_seed{args.seed}.png"
    summary_path = results_dir / f"hcfl2048_vs_fvm2048_seed{args.seed}.json"
    plot_final_profiles(references, predictions, metrics, plot_path)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(plot_path)
    print(summary_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
