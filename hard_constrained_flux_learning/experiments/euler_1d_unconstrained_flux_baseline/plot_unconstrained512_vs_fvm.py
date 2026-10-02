"""Compare the raw unconstrained flux model with native FVM on 512 cells.

The learned checkpoint was trained on 64-cell trajectories.  At deployment it
is advanced eight times per saved-output interval so that each learned update
has the same dt/dx ratio as in training.  No admissibility repair, fallback, or
hard output constraint is applied.  Native FVM-2048 is retained only as a
light visual background reference.
"""

from __future__ import annotations

import argparse
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
AUDIT_EXPERIMENT = HERE.parent / "euler_1d_convergence_audit"
if str(AUDIT_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(AUDIT_EXPERIMENT))

import plot_best_vs_fvm as comparison  # noqa: E402
import run_unconstrained_baseline as experiment  # noqa: E402


base = comparison.base
shared = comparison.shared
TARGET_CELLS = 512
FAILURE_COLOR = "#B42318"


def load_model(
    results_dir: Path,
    seed: int,
    width: int,
) -> experiment.UnconstrainedDirectSolver:
    checkpoint = (
        results_dir
        / f"{experiment.ARM}_converged_best_seed{seed}.pt"
    )
    if not checkpoint.is_file():
        raise FileNotFoundError(
            f"Missing converged checkpoint: {checkpoint}. "
            "Run run_unconstrained_baseline.py first."
        )
    model = experiment.UnconstrainedDirectSolver(
        np.zeros(3, dtype=np.float32),
        np.ones(3, dtype=np.float32),
        width=width,
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


def finite_minimum(values: torch.Tensor) -> float | None:
    finite = values[torch.isfinite(values)]
    if finite.numel() == 0:
        return None
    return float(finite.min())


@torch.no_grad()
def raw_rollout_on_grid(
    model: experiment.UnconstrainedDirectSolver,
    name: str,
    cells: int,
) -> tuple[torch.Tensor | None, dict[str, Any]]:
    """Deploy the unmodified learned flux and stop at first inadmissibility."""
    if cells % base.NCOARSE:
        raise ValueError(
            f"Target grid ({cells}) must be a multiple of the training grid "
            f"({base.NCOARSE})"
        )

    updates_per_saved_interval = cells // base.NCOARSE
    physical_dt_per_update = base.DT_SNAPSHOT / updates_per_saved_interval
    state = torch.from_numpy(comparison.initial_condition(name, cells)).float()
    snapshots = [state.clone()]
    completed_updates = 0
    minimum_density = finite_minimum(state[..., 0])
    minimum_pressure = finite_minimum(base.t_pressure(state))

    for saved_snapshot in range(1, shared.CANONICAL_NSNAP):
        for substep in range(1, updates_per_saved_interval + 1):
            candidate = model.one_step(state)
            density = candidate[..., 0]
            pressure = base.t_pressure(candidate)
            finite = bool(torch.isfinite(candidate).all())
            density_ok = finite and bool(
                (density >= experiment.RHO_FLOOR).all()
            )
            pressure_ok = finite and bool(
                (pressure >= experiment.PRESSURE_FLOOR).all()
            )

            if not (finite and density_ok and pressure_ok):
                attempted_update = completed_updates + 1
                if not finite:
                    reason = "nonfinite_state"
                elif not density_ok:
                    reason = "density_below_floor"
                else:
                    reason = "pressure_below_floor"
                return None, {
                    "completed": False,
                    "failure_reason": reason,
                    "failure_saved_interval": saved_snapshot,
                    "failure_substep_within_interval": substep,
                    "failure_attempted_update": attempted_update,
                    "failure_time": attempted_update * physical_dt_per_update,
                    "last_admissible_time": (
                        completed_updates * physical_dt_per_update
                    ),
                    "candidate_min_density": finite_minimum(density),
                    "candidate_min_pressure": finite_minimum(pressure),
                    "minimum_accepted_density": minimum_density,
                    "minimum_accepted_pressure": minimum_pressure,
                }

            state = candidate
            completed_updates += 1
            candidate_min_density = finite_minimum(density)
            candidate_min_pressure = finite_minimum(pressure)
            if candidate_min_density is not None:
                minimum_density = min(
                    float(minimum_density), candidate_min_density
                )
            if candidate_min_pressure is not None:
                minimum_pressure = min(
                    float(minimum_pressure), candidate_min_pressure
                )

        snapshots.append(state.clone())

    trajectory = torch.stack(snapshots, dim=1)
    return trajectory, {
        "completed": True,
        "completed_updates": completed_updates,
        "final_time": completed_updates * physical_dt_per_update,
        "minimum_accepted_density": minimum_density,
        "minimum_accepted_pressure": minimum_pressure,
    }


def plot_profiles(
    fvm_2048: dict[str, torch.Tensor],
    fvm_512: dict[str, torch.Tensor],
    predictions: dict[str, torch.Tensor | None],
    records: dict[str, dict[str, Any]],
    output: Path,
) -> None:
    x_2048 = (
        np.arange(comparison.precision.HIGH_REFERENCE_CELLS) + 0.5
    ) / comparison.precision.HIGH_REFERENCE_CELLS
    x_512 = (np.arange(TARGET_CELLS) + 0.5) / TARGET_CELLS
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
        top=0.79,
        wspace=0.22,
        hspace=0.14,
    )

    for column, name in enumerate(comparison.CASES):
        reference_2048 = base.primitive(fvm_2048[name]).numpy()[0, -1]
        native_512 = base.primitive(fvm_512[name]).numpy()[0, -1]
        prediction = predictions[name]
        record = records[name]
        if prediction is None:
            subtitle = f"failed at t = {record['failure_time']:.5f}"
            learned_512 = None
        else:
            subtitle = (
                "final NRMSE = "
                f"{record['final_snapshot_nrmse']:.4f}"
            )
            learned_512 = base.primitive(prediction).numpy()[0, -1]
        axes[0, column].set_title(
            f"{comparison.DISPLAY[name]}\n{subtitle}",
            fontsize=9.8,
            fontweight="semibold",
            color=(
                comparison.TEXT_COLOR
                if prediction is not None
                else FAILURE_COLOR
            ),
        )

        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_2048,
                reference_2048[:, row],
                color=comparison.BACKGROUND_REFERENCE_COLOR,
                linewidth=1.7,
                alpha=0.72,
                zorder=1,
            )
            axis.plot(
                x_512,
                native_512[:, row],
                color=comparison.FVM_512_COLOR,
                linewidth=1.45,
                zorder=3,
            )
            if learned_512 is not None:
                axis.plot(
                    x_512,
                    learned_512[:, row],
                    color=comparison.HCFL_COLOR,
                    linewidth=1.45,
                    linestyle="--",
                    zorder=4,
                )
            else:
                axis.text(
                    0.5,
                    0.5,
                    "no admissible final state",
                    transform=axis.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    color=FAILURE_COLOR,
                    fontweight="semibold",
                    bbox={
                        "boxstyle": "round,pad=0.25",
                        "facecolor": "white",
                        "edgecolor": FAILURE_COLOR,
                        "alpha": 0.88,
                        "linewidth": 0.8,
                    },
                    zorder=6,
                )
            comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=comparison.BACKGROUND_REFERENCE_COLOR,
            linewidth=1.7,
            alpha=0.72,
            label="FVM-2048 (background reference)",
        ),
        Line2D(
            [0],
            [0],
            color=comparison.FVM_512_COLOR,
            linewidth=1.45,
            label="FVM-512",
        ),
        Line2D(
            [0],
            [0],
            color=comparison.HCFL_COLOR,
            linewidth=1.45,
            linestyle="--",
            label="Unconstrained DNN-512 (64-grid-trained)",
        ),
    ]
    if any(not record["completed"] for record in records.values()):
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=FAILURE_COLOR,
                marker="x",
                linestyle="None",
                markersize=7,
                label="raw rollout failed before final time",
            )
        )

    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"Unconstrained DNN-512 versus FVM-512 at t = {final_time:.4f}",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.895),
        ncol=len(legend_handles),
        frameon=False,
        fontsize=9.5,
    )
    figure.text(
        0.99,
        0.018,
        "FVM-2048 is the light background. FVM-512 and the raw DNN flux use "
        "their native 512-cell states; no averaging, sampling, constraint, "
        "repair, or fallback. The 64-grid-trained DNN takes 8 CFL-matched "
        "updates per saved interval.",
        ha="right",
        fontsize=8.3,
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
        default=HERE / "results",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    model = load_model(results_dir, args.seed, args.width)
    state_std = comparison.training_state_std(args.seed)

    fvm_2048: dict[str, torch.Tensor] = {}
    fvm_512: dict[str, torch.Tensor] = {}
    predictions: dict[str, torch.Tensor | None] = {}
    records: dict[str, dict[str, Any]] = {}

    for name in comparison.CASES:
        print(f"Generating {comparison.DISPLAY[name]} references and raw rollout...")
        fvm_2048[name] = torch.from_numpy(
            comparison.strict_native_rollout(
                comparison.initial_condition(
                    name, comparison.precision.HIGH_REFERENCE_CELLS
                ),
                comparison.precision.HIGH_REFERENCE_CELLS,
            ).astype(np.float32)
        )
        fvm_512[name] = torch.from_numpy(
            comparison.strict_native_rollout(
                comparison.initial_condition(name, TARGET_CELLS),
                TARGET_CELLS,
            ).astype(np.float32)
        )
        prediction, record = raw_rollout_on_grid(model, name, TARGET_CELLS)
        predictions[name] = prediction
        if prediction is not None:
            record.update(
                comparison.diagnostics(
                    fvm_512[name], prediction, state_std
                )
            )
        records[name] = record

    completed_records = [
        record for record in records.values() if record["completed"]
    ]
    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    summary = {
        "seed": args.seed,
        "checkpoint": f"{experiment.ARM}_converged_best_seed{args.seed}.pt",
        "training_cells": base.NCOARSE,
        "deployment_cells": TARGET_CELLS,
        "updates_per_saved_interval": TARGET_CELLS // base.NCOARSE,
        "saved_interval": base.DT_SNAPSHOT,
        "final_time": final_time,
        "comparison": (
            "native FVM-512 versus raw unconstrained DNN on the same 512 "
            "cell centers; native FVM-2048 is a visual background only"
        ),
        "data_reduction": "none; no averaging or stride sampling",
        "raw_rollout_policy": (
            "stop at first nonfinite state, density below 1e-5, or pressure "
            "below 1e-5; no constraint, repair, or fallback"
        ),
        "completed_cases": len(completed_records),
        "total_cases": len(records),
        "completed_case_mean_final_snapshot_nrmse": (
            float(
                np.mean(
                    [
                        record["final_snapshot_nrmse"]
                        for record in completed_records
                    ]
                )
            )
            if completed_records
            else None
        ),
        "cases": records,
    }

    stem = f"unconstrained512_vs_fvm512_with_fvm2048_seed{args.seed}"
    figure_path = results_dir / f"{stem}.png"
    summary_path = results_dir / f"{stem}.json"
    plot_profiles(fvm_2048, fvm_512, predictions, records, figure_path)
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(figure_path)
    print(summary_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
