"""Compare learned HCFL with strong and exactly matched HLLC baselines.

The scientific control is the zero-correction model: it retains the same
HLLC proposal, hard safety stack, grid, and update protocol as learned HCFL,
but its neural correction is identically zero.  A classical native-grid
HLLC + SSP-RK2 solver is included as a stronger conventional baseline.
All quantitative errors use conservative 512-cell averages of a strict
native HLLC-2048 trajectory as the common reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import plot_best_vs_fvm as comparison  # noqa: E402


base = comparison.base
shared = comparison.shared
TARGET_CELLS = 512
REFERENCE_CELLS = 2048
HLLC_512_COLOR = "#0072B2"
MATCHED_HLLC_COLOR = "#009E73"
HCFL_COLOR = "#D55E00"


def physical_pressure(state: torch.Tensor) -> torch.Tensor:
    density = state[..., 0]
    momentum = state[..., 1]
    energy = state[..., 2]
    return (base.GAMMA - 1.0) * (
        energy - 0.5 * momentum * momentum / density
    )


def assert_admissible(state: torch.Tensor, context: str) -> None:
    pressure = physical_pressure(state)
    if (
        not bool(torch.isfinite(state).all())
        or not bool((state[..., 0] > 0.0).all())
        or not bool((pressure > 0.0).all())
    ):
        raise RuntimeError(
            f"Strict HLLC rollout lost admissibility during {context}; "
            "no repair was applied."
        )


@torch.no_grad()
def strict_native_hllc_rollout(
    name: str,
    cells: int,
    cfl: float = 0.25,
) -> tuple[torch.Tensor, dict[str, Any]]:
    """Evolve native-grid HLLC with adaptive SSP-RK2 and no repair."""
    state = torch.from_numpy(
        comparison.initial_condition(name, cells)
    ).double()
    cell_width = 1.0 / cells
    snapshots = [state.clone()]
    substeps = 0

    for snapshot in range(1, shared.CANONICAL_NSNAP):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            density = state[..., 0]
            pressure = physical_pressure(state)
            velocity = state[..., 1] / density
            sound_speed = torch.sqrt(
                base.GAMMA * pressure / density
            )
            max_speed = float((torch.abs(velocity) + sound_speed).max())
            dt = min(
                remaining,
                cfl * cell_width / max(max_speed, 1.0e-12),
            )
            lam = dt / cell_width

            flux = base.t_hllc(state)
            stage_one = state - lam * (
                flux - torch.roll(flux, 1, dims=-2)
            )
            assert_admissible(
                stage_one,
                f"{name}, snapshot {snapshot}, SSP-RK2 stage one",
            )

            stage_flux = base.t_hllc(stage_one)
            stage_two = stage_one - lam * (
                stage_flux - torch.roll(stage_flux, 1, dims=-2)
            )
            candidate = 0.5 * state + 0.5 * stage_two
            assert_admissible(
                candidate,
                f"{name}, snapshot {snapshot}, SSP-RK2 completion",
            )
            state = candidate
            remaining -= dt
            substeps += 1
        snapshots.append(state.clone())

    trajectory = torch.stack(snapshots, dim=1).float()
    return trajectory, {
        "completed": True,
        "cfl": cfl,
        "ssprk2_substeps": substeps,
        "flux_evaluations": 2 * substeps,
    }


def training_statistics(
    seed: int,
) -> tuple[np.ndarray, np.ndarray, torch.Tensor]:
    data = shared.make_baseline_training_data(seed)
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    return mean, std, state_std


def zero_correction_solver(
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
    seed: int,
) -> base.Solver:
    torch.manual_seed(12000 + seed)
    model = base.Solver("dissipation", mean, std, width=width)
    final_layer = model.flux_net.net[-1]
    if (
        float(final_layer.weight.detach().abs().max()) != 0.0
        or float(final_layer.bias.detach().abs().max()) != 0.0
    ):
        raise RuntimeError("The matched HLLC control is not zero initialized")
    model.eval()
    return model


def empty_safety_totals() -> defaultdict[str, float]:
    totals: defaultdict[str, float] = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["max_entropy_residual"] = -float("inf")
    totals["max_total_entropy_change"] = -float("inf")
    return totals


@torch.no_grad()
def matched_safe_rollout(
    model: base.Solver,
    name: str,
    cells: int,
) -> tuple[torch.Tensor, dict[str, Any]]:
    """Use the exact HCFL safety/update protocol on a target grid."""
    if cells % base.NCOARSE:
        raise ValueError(
            f"Target cells ({cells}) must be a multiple of "
            f"training cells ({base.NCOARSE})"
        )
    updates_per_saved_interval = cells // base.NCOARSE
    state = torch.from_numpy(
        comparison.initial_condition(name, cells)
    ).float()
    snapshots = [state.clone()]
    totals = empty_safety_totals()
    update_calls = 0
    hard_projection_active = 0
    hard_projection_total = 0
    hard_projection_squared = 0.0
    raw_flux_squared = 0.0
    hard_projection_maximum = 0.0

    for _ in range(1, shared.CANONICAL_NSNAP):
        for _ in range(updates_per_saved_interval):
            raw_flux = model.flux_net(state)
            projected_flux = model.flux(state)
            projection_delta = projected_flux - raw_flux
            interface_tolerance = 1.0e-7 * (
                1.0 + raw_flux.abs().amax(dim=-1)
            )
            hard_projection_active += int(
                (
                    torch.linalg.vector_norm(projection_delta, dim=-1)
                    > interface_tolerance
                ).sum()
            )
            hard_projection_total += projection_delta.shape[0] * (
                projection_delta.shape[1]
            )
            hard_projection_squared += float(
                (projection_delta.double() ** 2).sum()
            )
            raw_flux_squared += float((raw_flux.double() ** 2).sum())
            hard_projection_maximum = max(
                hard_projection_maximum,
                float(projection_delta.abs().max()),
            )
            state, step = shared.advance_safe_snapshot(model, state)
            shared.merge_step_stats(totals, step)
            update_calls += 1
        snapshots.append(state.clone())

    trajectory = torch.stack(snapshots, dim=1)
    return trajectory, {
        "completed": True,
        "effective_physical_dt_per_call": (
            base.DT_SNAPSHOT / updates_per_saved_interval
        ),
        "safe_update_calls": update_calls,
        "internal_substeps": int(totals["substeps"]),
        "hard_projection_intervention_rate": (
            hard_projection_active / max(hard_projection_total, 1)
        ),
        "hard_projection_relative_flux_rms": (
            np.sqrt(hard_projection_squared / max(raw_flux_squared, 1.0e-30))
        ),
        "hard_projection_max_absolute_change": hard_projection_maximum,
        "local_limiter_intervention_rate": (
            totals["local_active"] / max(totals["local_total"], 1.0)
        ),
        "fd_entropy_intervention_rate": (
            totals["fd_active"] / max(totals["fd_total"], 1.0)
        ),
        "mean_fd_beta": (
            totals["fd_beta_sum"] / max(totals["fd_total"], 1.0)
        ),
        "min_fd_beta": totals["fd_beta_min"],
        "entropy_violation_rate": (
            totals["entropy_violations"]
            / max(totals["entropy_total"], 1.0)
        ),
        "max_entropy_residual": totals["max_entropy_residual"],
        "max_total_entropy_change": totals[
            "max_total_entropy_change"
        ],
    }


def integrity_metrics(trajectory: torch.Tensor) -> dict[str, Any]:
    pressure = physical_pressure(trajectory)
    initial_integral = trajectory[:, :1].mean(dim=-2)
    integral = trajectory.mean(dim=-2)
    drift = torch.abs(integral - initial_integral)
    return {
        "minimum_density": float(trajectory[..., 0].min()),
        "minimum_pressure": float(pressure.min()),
        "max_conservation_drift": {
            "density": float(drift[..., 0].max()),
            "momentum": float(drift[..., 1].max()),
            "energy": float(drift[..., 2].max()),
        },
    }


def enrich_metrics(
    reference: torch.Tensor,
    prediction: torch.Tensor,
    state_std: torch.Tensor,
    run_stats: dict[str, Any],
) -> dict[str, Any]:
    output = dict(run_stats)
    output.update(comparison.diagnostics(reference, prediction, state_std))
    output.update(integrity_metrics(prediction))
    return output


def improvement_percent(candidate: float, comparator: float) -> float:
    return 100.0 * (comparator - candidate) / comparator


def plot_profiles(
    native_reference: dict[str, torch.Tensor],
    classical_hllc: dict[str, torch.Tensor],
    matched_hllc: dict[str, torch.Tensor],
    learned_hcfl: dict[str, torch.Tensor],
    cases: dict[str, dict[str, Any]],
    output: Path,
) -> None:
    x_reference = (
        np.arange(REFERENCE_CELLS) + 0.5
    ) / REFERENCE_CELLS
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
        top=0.78,
        wspace=0.22,
        hspace=0.14,
    )

    for column, name in enumerate(comparison.CASES):
        reference = base.primitive(native_reference[name]).numpy()[0, -1]
        classical = base.primitive(classical_hllc[name]).numpy()[0, -1]
        matched = base.primitive(matched_hllc[name]).numpy()[0, -1]
        learned = base.primitive(learned_hcfl[name]).numpy()[0, -1]
        reduction = cases[name][
            "learned_reduction_vs_matched_hllc_rollout_percent"
        ]
        axes[0, column].set_title(
            f"{comparison.DISPLAY[name]}\n"
            f"learned reduction: {reduction:.1f}%",
            fontsize=9.8,
            fontweight="semibold",
            color=comparison.TEXT_COLOR,
        )

        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                reference[:, row],
                color=comparison.BACKGROUND_REFERENCE_COLOR,
                linewidth=1.7,
                alpha=0.72,
                zorder=1,
            )
            axis.plot(
                x_512,
                classical[:, row],
                color=HLLC_512_COLOR,
                linewidth=1.35,
                zorder=2,
            )
            axis.plot(
                x_512,
                matched[:, row],
                color=MATCHED_HLLC_COLOR,
                linewidth=1.25,
                linestyle=":",
                zorder=3,
            )
            axis.plot(
                x_512,
                learned[:, row],
                color=HCFL_COLOR,
                linewidth=1.45,
                linestyle="--",
                zorder=4,
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
            label="HLLC + SSP-RK2, 2048 cells (reference)",
        ),
        Line2D(
            [0],
            [0],
            color=HLLC_512_COLOR,
            linewidth=1.35,
            label="HLLC + SSP-RK2, 512 cells",
        ),
        Line2D(
            [0],
            [0],
            color=MATCHED_HLLC_COLOR,
            linewidth=1.25,
            linestyle=":",
            label="matched hard-safe HLLC-512 (zero NN)",
        ),
        Line2D(
            [0],
            [0],
            color=HCFL_COLOR,
            linewidth=1.45,
            linestyle="--",
            label="learned HCFL-512 (trained on 64-cell states)",
        ),
    ]
    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    figure.suptitle(
        f"HLLC-512 versus learned HCFL-512 at t = {final_time:.4f}",
        y=0.975,
        fontsize=17,
        fontweight="bold",
        color=comparison.TEXT_COLOR,
    )
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.895),
        ncol=4,
        frameon=False,
        fontsize=9.1,
    )
    figure.text(
        0.99,
        0.018,
        "Profiles are native-grid values. Reported NRMSE uses conservative "
        "HLLC-2048 to 512 restriction. The matched zero-NN control retains "
        "the exact HCFL safety stack and update schedule.",
        ha="right",
        fontsize=8.3,
        color=comparison.MUTED_COLOR,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def write_csv(path: Path, cases: dict[str, dict[str, Any]]) -> None:
    rows: list[dict[str, Any]] = []
    for name, case in cases.items():
        for method in (
            "classical_hllc_512",
            "matched_safe_hllc_512_zero_nn",
            "learned_hcfl_512",
        ):
            metrics = case[method]
            rows.append(
                {
                    "case": name,
                    "method": method,
                    "rollout_nrmse": metrics["rollout_nrmse"],
                    "final_snapshot_nrmse": metrics[
                        "final_snapshot_nrmse"
                    ],
                    "minimum_density": metrics["minimum_density"],
                    "minimum_pressure": metrics["minimum_pressure"],
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
        "--results-dir",
        type=Path,
        default=HERE / "results",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    mean, std, state_std = training_statistics(args.seed)
    learned = comparison.load_model(
        results_dir,
        "dissipation_broad",
        args.seed,
        args.width,
    )
    matched_zero = zero_correction_solver(
        mean,
        std,
        args.width,
        args.seed,
    )

    native_reference: dict[str, torch.Tensor] = {}
    reference_512: dict[str, torch.Tensor] = {}
    classical_hllc: dict[str, torch.Tensor] = {}
    matched_hllc: dict[str, torch.Tensor] = {}
    learned_hcfl: dict[str, torch.Tensor] = {}
    cases: dict[str, dict[str, Any]] = {}

    for name in comparison.CASES:
        print(f"Evaluating {comparison.DISPLAY[name]}...", flush=True)
        native_reference[name], reference_stats = (
            strict_native_hllc_rollout(name, REFERENCE_CELLS)
        )
        restricted = comparison.precision.conservative_restrict(
            native_reference[name].numpy(), TARGET_CELLS
        )
        reference_512[name] = torch.from_numpy(
            restricted.astype(np.float32)
        )
        classical_hllc[name], classical_stats = (
            strict_native_hllc_rollout(name, TARGET_CELLS)
        )
        matched_hllc[name], matched_stats = matched_safe_rollout(
            matched_zero,
            name,
            TARGET_CELLS,
        )
        learned_hcfl[name], learned_stats = matched_safe_rollout(
            learned,
            name,
            TARGET_CELLS,
        )

        classical_metrics = enrich_metrics(
            reference_512[name],
            classical_hllc[name],
            state_std,
            classical_stats,
        )
        matched_metrics = enrich_metrics(
            reference_512[name],
            matched_hllc[name],
            state_std,
            matched_stats,
        )
        learned_metrics = enrich_metrics(
            reference_512[name],
            learned_hcfl[name],
            state_std,
            learned_stats,
        )
        cases[name] = {
            "reference_hllc_2048": reference_stats,
            "classical_hllc_512": classical_metrics,
            "matched_safe_hllc_512_zero_nn": matched_metrics,
            "learned_hcfl_512": learned_metrics,
            "learned_reduction_vs_classical_hllc_rollout_percent": (
                improvement_percent(
                    learned_metrics["rollout_nrmse"],
                    classical_metrics["rollout_nrmse"],
                )
            ),
            "learned_reduction_vs_matched_hllc_rollout_percent": (
                improvement_percent(
                    learned_metrics["rollout_nrmse"],
                    matched_metrics["rollout_nrmse"],
                )
            ),
        }

    method_names = (
        "classical_hllc_512",
        "matched_safe_hllc_512_zero_nn",
        "learned_hcfl_512",
    )
    mean_metrics = {
        method: {
            "mean_rollout_nrmse": float(
                np.mean(
                    [cases[name][method]["rollout_nrmse"] for name in cases]
                )
            ),
            "mean_final_snapshot_nrmse": float(
                np.mean(
                    [
                        cases[name][method]["final_snapshot_nrmse"]
                        for name in cases
                    ]
                )
            ),
        }
        for method in method_names
    }
    learned_mean = mean_metrics["learned_hcfl_512"][
        "mean_rollout_nrmse"
    ]
    mean_metrics["learned_hcfl_512"][
        "reduction_vs_classical_hllc_percent"
    ] = improvement_percent(
        learned_mean,
        mean_metrics["classical_hllc_512"]["mean_rollout_nrmse"],
    )
    mean_metrics["learned_hcfl_512"][
        "reduction_vs_matched_zero_nn_percent"
    ] = improvement_percent(
        learned_mean,
        mean_metrics["matched_safe_hllc_512_zero_nn"][
            "mean_rollout_nrmse"
        ],
    )

    final_time = (shared.CANONICAL_NSNAP - 1) * base.DT_SNAPSHOT
    summary = {
        "seed": args.seed,
        "checkpoint": (
            f"dissipation_broad_converged_best_seed{args.seed}.pt"
        ),
        "checkpoint_training_cells": base.NCOARSE,
        "deployment_cells": TARGET_CELLS,
        "reference": (
            "strict native HLLC + SSP-RK2 on 2048 cells, conservatively "
            "restricted to 512 cells for every metric"
        ),
        "final_time": final_time,
        "metric": (
            "RMSE of conserved variables standardized by the 64-cell "
            "training-set component standard deviations"
        ),
        "matched_control": (
            "zero neural correction with identical HLLC proposal, hard "
            "safety stack, 512-cell state, and update schedule"
        ),
        "mean_metrics": mean_metrics,
        "cases": cases,
    }

    stem = f"hllc512_vs_hcfl512_with_hllc2048_seed{args.seed}"
    figure_path = results_dir / f"{stem}.png"
    summary_path = results_dir / f"{stem}.json"
    csv_path = results_dir / f"{stem}.csv"
    plot_profiles(
        native_reference,
        classical_hllc,
        matched_hllc,
        learned_hcfl,
        cases,
        figure_path,
    )
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    write_csv(csv_path, cases)
    print(figure_path)
    print(summary_path)
    print(csv_path)
    print(json.dumps(summary["mean_metrics"], indent=2))


if __name__ == "__main__":
    main()
