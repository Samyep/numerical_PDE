"""Compare the validation-selected best HCFL checkpoint with FVM references.

The checkpoint is selected by the convergence audit, not by these canonical
tests.  Profiles show the strict 2048-cell periodic Rusanov + SSP-RK2 rollout
on its native grid, a separately evolved 64-cell FVM baseline, and HCFL-64.
Metrics compare both coarse solvers with conservative 64-cell averages of the
2048-cell reference.
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
FVM_512_COLOR = "#0072B2"
BACKGROUND_REFERENCE_COLOR = "#B8C0CA"
HCFL_COLOR = "#D55E00"
GRID_COLOR = "#D8DEE9"
SPINE_COLOR = "#9AA5B1"
TEXT_COLOR = "#1F2933"
MUTED_COLOR = "#66788A"
COMPARISON_CELLS = 512


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


def initial_condition(name: str, cells: int) -> np.ndarray:
    left, right = CASES[name]
    return shared._two_state_periodic(
        cells,
        cells // 2,
        left,
        right,
        0,
    )[None, ...]


def strict_native_rollout(initial: np.ndarray, cells: int) -> np.ndarray:
    """Return every snapshot on the solver's native finite-volume grid."""
    if initial.shape[-2:] != (cells, 3):
        raise ValueError(f"Expected initial shape (..., {cells}, 3), got {initial.shape}")

    state = np.asarray(initial, dtype=np.float64).copy()
    cell_width = 1.0 / cells
    snapshots = [state.copy()]
    for _ in range(1, shared.CANONICAL_NSNAP):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            rho = state[..., 0]
            pressure = base.np_pressure(state)
            velocity = state[..., 1] / rho
            sound_speed = np.sqrt(base.GAMMA * pressure / rho)
            max_speed = float(np.max(np.abs(velocity) + sound_speed))
            dt = min(remaining, 0.25 * cell_width / max(max_speed, 1.0e-12))
            state = base.np_ssprk2(state, dt, cell_width)
            if (
                not np.isfinite(state).all()
                or (state[..., 0] <= 0).any()
                or (base.np_pressure(state) <= 0).any()
            ):
                raise RuntimeError(
                    f"Strict {cells}-cell FVM rollout lost admissibility; "
                    "no state repair was applied."
                )
            remaining -= dt
        snapshots.append(state.copy())
    return np.stack(snapshots, axis=1)


def fvm_trajectories(
    name: str,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return native FVM-2048, its 64-cell averages, and native FVM-64."""
    native_high_values = strict_native_rollout(
        initial_condition(name, precision.HIGH_REFERENCE_CELLS),
        precision.HIGH_REFERENCE_CELLS,
    )
    projected_high = torch.from_numpy(
        precision.conservative_restrict(
            native_high_values, base.NCOARSE
        ).astype(np.float32)
    )
    native_high = torch.from_numpy(native_high_values.astype(np.float32))
    native_coarse = torch.from_numpy(
        strict_native_rollout(
            initial_condition(name, base.NCOARSE),
            base.NCOARSE,
        ).astype(np.float32)
    )
    initial_mismatch = float(
        torch.max(torch.abs(projected_high[:, 0] - native_coarse[:, 0]))
    )
    if initial_mismatch > 2.0e-6:
        raise RuntimeError(
            f"Projected and native 64-cell initial states differ by "
            f"{initial_mismatch:.3e}"
        )
    return native_high, projected_high, native_coarse


@torch.no_grad()
def hcfl_rollout(model: base.Solver, reference: torch.Tensor) -> torch.Tensor:
    state = reference[:, 0].clone()
    snapshots = [state.clone()]
    for _ in range(1, reference.shape[1]):
        state, _ = shared.advance_safe_snapshot(model, state)
        snapshots.append(state.clone())
    return torch.stack(snapshots, dim=1)


@torch.no_grad()
def hcfl_rollout_on_grid(
    model: base.Solver,
    name: str,
    cells: int,
) -> torch.Tensor:
    """Deploy the 64-cell-trained flux on a finer grid at the same CFL."""
    if cells % base.NCOARSE:
        raise ValueError(
            f"Target grid ({cells}) must be a multiple of the training grid "
            f"({base.NCOARSE})"
        )
    substeps_per_snapshot = cells // base.NCOARSE
    state = torch.from_numpy(initial_condition(name, cells)).float()
    snapshots = [state.clone()]
    for _ in range(1, shared.CANONICAL_NSNAP):
        for _ in range(substeps_per_snapshot):
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
    native_references: dict[str, torch.Tensor],
    predictions: dict[str, torch.Tensor],
    output: Path,
) -> None:
    sampling_stride = precision.HIGH_REFERENCE_CELLS // base.NCOARSE
    if sampling_stride * base.NCOARSE != precision.HIGH_REFERENCE_CELLS:
        raise ValueError("Fine and coarse grids must be nested")
    sampling_offset = sampling_stride // 2
    sample_indices = (
        np.arange(base.NCOARSE) * sampling_stride + sampling_offset
    )
    x_coarse = (np.arange(base.NCOARSE) + 0.5) / base.NCOARSE
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
        native_reference = base.primitive(native_references[name]).numpy()[0, -1]
        sampled_reference = native_reference[sample_indices]
        prediction = base.primitive(predictions[name]).numpy()[0, -1]
        axes[0, column].set_title(
            DISPLAY[name],
            fontsize=10.5,
            fontweight="semibold",
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_coarse,
                sampled_reference[:, row],
                color=FVM_COLOR,
                linewidth=1.7,
                marker="o",
                markersize=2.2,
                label="FVM-2048 (one central sample per 32 cells)",
                zorder=2,
            )
            axis.plot(
                x_coarse,
                prediction[:, row],
                color=HCFL_COLOR,
                linewidth=1.55,
                linestyle="--",
                label="HCFL-64 (dissipation)",
                zorder=3,
            )
            style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    figure.suptitle(
        "HCFL-64 versus stride-sampled FVM-2048 at t = 0.0252",
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
        "FVM-2048 is not averaged: one near-center fine cell is selected from "
        "each consecutive block of 32 and the resulting 64 values are connected.",
        ha="right",
        fontsize=8.5,
        color=MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def plot_fvm512_vs_hcfl512(
    native_references: dict[str, torch.Tensor],
    fvm_512: dict[str, torch.Tensor],
    hcfl_512: dict[str, torch.Tensor],
    output: Path,
) -> None:
    x_2048 = (
        np.arange(precision.HIGH_REFERENCE_CELLS) + 0.5
    ) / precision.HIGH_REFERENCE_CELLS
    x_512 = (np.arange(COMPARISON_CELLS) + 0.5) / COMPARISON_CELLS
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
        reference_2048 = base.primitive(native_references[name]).numpy()[0, -1]
        native_512 = base.primitive(fvm_512[name]).numpy()[0, -1]
        learned_512 = base.primitive(hcfl_512[name]).numpy()[0, -1]
        axes[0, column].set_title(
            DISPLAY[name],
            fontsize=10.5,
            fontweight="semibold",
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_2048,
                reference_2048[:, row],
                color=BACKGROUND_REFERENCE_COLOR,
                linewidth=1.7,
                alpha=0.72,
                label="FVM-2048 (background reference)",
                zorder=1,
            )
            axis.plot(
                x_512,
                native_512[:, row],
                color=FVM_512_COLOR,
                linewidth=1.45,
                label="FVM-512",
                zorder=3,
            )
            axis.plot(
                x_512,
                learned_512[:, row],
                color=HCFL_COLOR,
                linewidth=1.45,
                linestyle="--",
                label="HCFL-512 (64-grid-trained checkpoint)",
                zorder=4,
            )
            style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(row_labels[row], fontsize=10)
            if row == 2:
                axis.set_xlabel("x", fontsize=9)

    figure.suptitle(
        "FVM-512 versus HCFL-512 at t = 0.0252",
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
        ncol=3,
        frameon=False,
        fontsize=10,
    )
    figure.text(
        0.99,
        0.018,
        "FVM-2048 is the light background. HCFL uses the validation-selected "
        "64-grid-trained checkpoint on 512 cells with 8 CFL-matched substeps "
        "per snapshot.",
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
    native_references: dict[str, torch.Tensor] = {}
    references: dict[str, torch.Tensor] = {}
    coarse_fvm: dict[str, torch.Tensor] = {}
    for name in CASES:
        native, projected, coarse = fvm_trajectories(name)
        native_references[name] = native
        references[name] = projected
        coarse_fvm[name] = coarse
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
    coarse_diagnostics = {
        name: diagnostics(references[name], coarse_fvm[name], state_std)
        for name in CASES
    }
    coarse_canonical_mean = float(
        np.mean(
            [row["rollout_nrmse"] for row in coarse_diagnostics.values()]
        )
    )

    selected_arm = "dissipation_broad"
    predictions = all_predictions[selected_arm]
    case_summary = all_diagnostics[selected_arm]
    post_hoc_best_arm = min(canonical_means, key=canonical_means.get)
    fvm_512 = {
        name: torch.from_numpy(
            strict_native_rollout(
                initial_condition(name, COMPARISON_CELLS),
                COMPARISON_CELLS,
            ).astype(np.float32)
        )
        for name in CASES
    }
    hcfl_512 = {
        name: hcfl_rollout_on_grid(
            models[selected_arm],
            name,
            COMPARISON_CELLS,
        )
        for name in CASES
    }

    summary = {
        "seed": args.seed,
        "selected_model": f"{selected_arm}_converged_best",
        "selection_metric": "independent validation rollout NRMSE",
        "profile_reference": (
            "2048-cell periodic Rusanov + SSP-RK2 shown on its native grid"
        ),
        "profile_comparison": (
            "HCFL-64 and 64 point samples from native FVM-2048; sample index "
            "16 + 32*i in each fine-grid trajectory; no averaging"
        ),
        "fvm512_hcfl512_profile": (
            "native FVM-512 versus the 64-grid-trained dissipation checkpoint "
            "deployed on 512 cells with 8 CFL-matched substeps per snapshot; "
            "native FVM-2048 retained as a light background reference"
        ),
        "metric_reference": (
            "conservative 64-cell averages of the 2048-cell periodic "
            "Rusanov + SSP-RK2 trajectory"
        ),
        "coarse_fvm_baseline": (
            "64-cell periodic Rusanov + SSP-RK2 evolved directly on the "
            "coarse grid"
        ),
        "final_time": (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT,
        "cases": case_summary,
        "fvm_64_cases": coarse_diagnostics,
        "canonical_mean_rollout_nrmse": canonical_means[selected_arm],
        "fvm_64_canonical_mean_rollout_nrmse": coarse_canonical_mean,
        "hcfl_improvement_over_fvm_64_percent": (
            100.0
            * (coarse_canonical_mean - canonical_means[selected_arm])
            / coarse_canonical_mean
        ),
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
    for name, row in coarse_diagnostics.items():
        comparison_rows.append(
            {
                "seed": args.seed,
                "arm": "fvm_64",
                "model": "Rusanov + SSP-RK2",
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
            "arm": "fvm_64",
            "model": "Rusanov + SSP-RK2",
            "split": "canonical_mean",
            "rollout_nrmse": coarse_canonical_mean,
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
    fvm512_profiles_path = (
        results_dir
        / f"fvm512_vs_hcfl512_with_fvm2048_seed{args.seed}.png"
    )
    error_maps_path = results_dir / f"best_vs_fvm_error_maps_seed{args.seed}.png"
    plot_final_profiles(
        native_references,
        predictions,
        profiles_path,
    )
    plot_fvm512_vs_hcfl512(
        native_references,
        fvm_512,
        hcfl_512,
        fvm512_profiles_path,
    )
    plot_error_maps(references, predictions, state_std, error_maps_path)

    print(summary_path)
    print(comparison_path)
    print(profiles_path)
    print(fvm512_profiles_path)
    print(error_maps_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
