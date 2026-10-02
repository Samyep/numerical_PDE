"""Compare the validation-selected best HCFL checkpoint with FVM references.

The checkpoint is selected by the convergence audit, not by these canonical
tests.  The comparison reference is a strict 2048-cell periodic Rusanov +
SSP-RK2 finite-volume rollout, conservatively restricted to the 64-cell grid.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
PRECISION_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_reference_precision_ablation"
if str(PRECISION_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(PRECISION_EXPERIMENT))

import run_precision_ablation as precision  # noqa: E402


shared = precision.shared
base = precision.base

CASES = {
    "sod": ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
    "lax": ((0.445, 0.698, 3.528), (0.5, 0.0, 0.571)),
    "collision": ((1.0, 2.0, 1.0), (1.0, -2.0, 1.0)),
    "strong_pressure": ((1.0, 0.0, 5.0), (1.0, 0.0, 0.05)),
    "near_vacuum_expansion": ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4)),
}

ARM_MODELS = {
    "direct_broad": "direct",
    "invariant_broad": "invariant",
    "characteristic_broad": "characteristic",
    "dissipation_broad": "dissipation",
    "conv_broad": "conv",
    "direct_wave": "direct",
}

DISPLAY = {
    "sod": "Sod",
    "lax": "Lax",
    "collision": "Collision",
    "strong_pressure": "Strong pressure",
    "near_vacuum_expansion": "Near-vacuum expansion",
}

FVM_COLOR = "#202A35"
HCFL_COLOR = "#D55E00"
GRID_COLOR = "#D8DEE9"
SPINE_COLOR = "#9AA5B1"
TEXT_COLOR = "#1F2933"
MUTED_COLOR = "#66788A"


def load_model(
    results_dir: Path,
    arm: str,
    seed: int,
    width: int,
) -> base.Solver:
    model = base.Solver(ARM_MODELS[arm], np.zeros(3), np.ones(3), width=width)
    checkpoint = results_dir / f"{arm}_converged_best_seed{seed}.pt"
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


def training_state_std(seed: int) -> torch.Tensor:
    """Recreate the normalization used by the convergence audit."""
    training_data = shared.make_baseline_training_data(seed)
    return training_data.std(dim=(0, 1, 2))


def fvm_reference(name: str) -> torch.Tensor:
    left, right = CASES[name]
    initial = shared._two_state_periodic(
        precision.HIGH_REFERENCE_CELLS,
        precision.HIGH_REFERENCE_CELLS // 2,
        left,
        right,
        0,
    )[None, ...]
    return torch.from_numpy(
        precision.strict_rollout_reference(
            initial,
            precision.HIGH_REFERENCE_CELLS,
            nsnap=shared.CANONICAL_NSNAP,
        )
    )


@torch.no_grad()
def hcfl_rollout(model: base.Solver, reference: torch.Tensor) -> torch.Tensor:
    state = reference[:, 0].clone()
    snapshots = [state.clone()]
    for _ in range(1, reference.shape[1]):
        state, _ = shared.advance_safe_snapshot(model, state)
        snapshots.append(state.clone())
    return torch.stack(snapshots, dim=1)


def diagnostics(
    reference: torch.Tensor,
    prediction: torch.Tensor,
    state_std: torch.Tensor,
) -> dict[str, object]:
    normalized = (prediction[:, 1:] - reference[:, 1:]) / state_std
    per_snapshot = torch.sqrt(torch.mean(normalized**2, dim=(0, 2, 3)))
    primitive_reference = base.primitive(reference)[0]
    primitive_prediction = base.primitive(prediction)[0]
    primitive_error = torch.abs(primitive_prediction - primitive_reference)
    return {
        "rollout_nrmse": float(torch.sqrt(torch.mean(normalized**2))),
        "final_snapshot_nrmse": float(per_snapshot[-1]),
        "peak_snapshot_nrmse": float(per_snapshot.max()),
        "peak_snapshot_index": int(per_snapshot.argmax()) + 1,
        "max_final_primitive_absolute_error": {
            "density": float(primitive_error[-1, :, 0].max()),
            "velocity": float(primitive_error[-1, :, 1].max()),
            "pressure": float(primitive_error[-1, :, 2].max()),
        },
    }


def style_axis(axis: plt.Axes) -> None:
    axis.grid(True, color=GRID_COLOR, linewidth=0.7, alpha=0.65)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(SPINE_COLOR)
    axis.spines["bottom"].set_color(SPINE_COLOR)
    axis.tick_params(labelsize=8, colors="#465362")


def plot_final_profiles(
    references: dict[str, torch.Tensor],
    predictions: dict[str, torch.Tensor],
    summary: dict[str, dict[str, object]],
    output: Path,
) -> None:
    x = (np.arange(base.NCOARSE) + 0.5) / base.NCOARSE
    row_labels = [r"Density $\rho$", r"Velocity $u$", r"Pressure $p$"]
    figure, axes = plt.subplots(
        3,
        len(CASES),
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

    for column, name in enumerate(CASES):
        reference = base.primitive(references[name]).numpy()[0, -1]
        prediction = base.primitive(predictions[name]).numpy()[0, -1]
        axes[0, column].set_title(
            f"{DISPLAY[name]}\nrollout NRMSE {summary[name]['rollout_nrmse']:.4f}",
            fontsize=10.5,
            fontweight="semibold",
        )
        for row in range(3):
            axis = axes[row, column]
            axis.step(
                x,
                reference[:, row],
                where="mid",
                color=FVM_COLOR,
                linewidth=2.1,
                label="FVM-2048 reference",
            )
            axis.plot(
                x,
                prediction[:, row],
                color=HCFL_COLOR,
                linewidth=1.55,
                linestyle="--",
                marker="o",
                markersize=2.0,
                markevery=4,
                label="best HCFL (dissipation)",
            )
            style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    figure.suptitle(
        "Validation-selected HCFL versus high-resolution FVM at t = 0.0252",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.895),
        ncol=2,
        frameon=False,
        fontsize=10,
    )
    figure.text(
        0.99,
        0.018,
        "Periodic domain; FVM uses 2048 cells and is conservatively restricted to 64 cell averages.",
        ha="right",
        fontsize=8.5,
        color=MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def plot_error_maps(
    references: dict[str, torch.Tensor],
    predictions: dict[str, torch.Tensor],
    state_std: torch.Tensor,
    output: Path,
) -> None:
    maps: dict[str, np.ndarray] = {}
    for name in CASES:
        normalized = (
            (predictions[name][:, 1:] - references[name][:, 1:]) / state_std
        ).numpy()[0]
        maps[name] = np.sqrt(np.mean(normalized**2, axis=-1))

    positive = np.concatenate([value.ravel() for value in maps.values()])
    positive = positive[positive > 0.0]
    lower = max(float(np.quantile(positive, 0.01)), 1.0e-5)
    upper = float(np.quantile(positive, 0.995))
    norm = LogNorm(vmin=lower, vmax=upper)

    figure, axes = plt.subplots(
        len(CASES),
        1,
        figsize=(12.5, 10.2),
        sharex=True,
        constrained_layout=False,
    )
    figure.subplots_adjust(
        left=0.14,
        right=0.91,
        bottom=0.09,
        top=0.88,
        hspace=0.22,
    )
    image = None
    for axis, name in zip(axes, CASES, strict=True):
        image = axis.imshow(
            np.maximum(maps[name], lower),
            origin="lower",
            aspect="auto",
            extent=(
                0.0,
                1.0,
                base.DT_SNAPSHOT,
                (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT,
            ),
            cmap="magma",
            norm=norm,
            interpolation="nearest",
        )
        axis.set_ylabel(f"{DISPLAY[name]}\ntime", fontsize=9)
        axis.tick_params(labelsize=8)
    axes[-1].set_xlabel("x", fontsize=10)
    assert image is not None
    color_axis = figure.add_axes((0.93, 0.12, 0.018, 0.71))
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("normalized per-cell RMSE", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)
    figure.suptitle(
        "Where the best HCFL rollout departs from FVM-2048",
        y=0.97,
        fontsize=17,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    figure.text(
        0.5,
        0.92,
        "Error combines density, momentum, and energy; color uses a logarithmic scale.",
        ha="center",
        fontsize=9.5,
        color=MUTED_COLOR,
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
    state_std = training_state_std(args.seed)
    references = {name: fvm_reference(name) for name in CASES}
    models = {
        arm: load_model(results_dir, arm, args.seed, args.width)
        for arm in ARM_MODELS
    }
    all_predictions = {
        arm: {
            name: hcfl_rollout(model, references[name])
            for name in CASES
        }
        for arm, model in models.items()
    }
    all_diagnostics = {
        arm: {
            name: diagnostics(references[name], predictions[name], state_std)
            for name in CASES
        }
        for arm, predictions in all_predictions.items()
    }
    canonical_means = {
        arm: float(
            np.mean(
                [row["rollout_nrmse"] for row in case_summary.values()]
            )
        )
        for arm, case_summary in all_diagnostics.items()
    }

    selected_arm = "dissipation_broad"
    predictions = all_predictions[selected_arm]
    case_summary = all_diagnostics[selected_arm]
    post_hoc_best_arm = min(canonical_means, key=canonical_means.get)

    summary = {
        "seed": args.seed,
        "selected_model": f"{selected_arm}_converged_best",
        "selection_metric": "independent validation rollout NRMSE",
        "fvm_reference": (
            "2048-cell periodic Rusanov + SSP-RK2, strict/no repair, "
            "conservatively restricted to 64 cells"
        ),
        "final_time": (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT,
        "cases": case_summary,
        "canonical_mean_rollout_nrmse": canonical_means[selected_arm],
        "all_converged_arm_canonical_means": canonical_means,
        "post_hoc_best_arm_on_fvm_2048": post_hoc_best_arm,
    }
    summary_path = results_dir / f"best_vs_fvm_summary_seed{args.seed}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    comparison_rows: list[dict[str, object]] = []
    for arm, case_rows in all_diagnostics.items():
        for name, row in case_rows.items():
            comparison_rows.append(
                {
                    "seed": args.seed,
                    "arm": arm,
                    "model": ARM_MODELS[arm],
                    "split": name,
                    "rollout_nrmse": row["rollout_nrmse"],
                    "final_snapshot_nrmse": row["final_snapshot_nrmse"],
                    "peak_snapshot_nrmse": row["peak_snapshot_nrmse"],
                    "peak_snapshot_index": row["peak_snapshot_index"],
                }
            )
        comparison_rows.append(
            {
                "seed": args.seed,
                "arm": arm,
                "model": ARM_MODELS[arm],
                "split": "canonical_mean",
                "rollout_nrmse": canonical_means[arm],
                "final_snapshot_nrmse": "",
                "peak_snapshot_nrmse": "",
                "peak_snapshot_index": "",
            }
        )
    comparison_path = results_dir / f"all_arms_vs_fvm_seed{args.seed}.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0]))
        writer.writeheader()
        writer.writerows(comparison_rows)

    profiles_path = results_dir / f"best_vs_fvm_profiles_seed{args.seed}.png"
    error_maps_path = results_dir / f"best_vs_fvm_error_maps_seed{args.seed}.png"
    plot_final_profiles(references, predictions, case_summary, profiles_path)
    plot_error_maps(references, predictions, state_std, error_maps_path)

    print(summary_path)
    print(comparison_path)
    print(profiles_path)
    print(error_maps_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
