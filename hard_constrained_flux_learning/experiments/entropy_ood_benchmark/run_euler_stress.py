"""Exploratory 64-cell Euler entropy/OOD stress evaluation.

The exploratory design record is stored beside this script.  The evaluator
loads supplied checkpoints and performs no checkpoint or hyperparameter
selection itself; RoeNet training/selection is handled by the separate
``train_roenet_adaptation.py`` script.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
PAPER_DIR = EXPERIMENTS / "paper_64_benchmark"
for module_dir in (PAPER_DIR,):
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import run_64_benchmark as paper  # noqa: E402
from roenet_adaptation import RoeNetEuler1d  # noqa: E402


CELLS = 64
SNAPSHOTS = 64
REFERENCE_CELLS = 2048
TRAINING_INTERVALS = 15
CASES: dict[str, dict[str, Any]] = {
    "sod": {
        "display": "Sod control",
        "left": (1.0, 0.0, 1.0),
        "right": (0.125, 0.0, 0.1),
        "role": "ordinary_control",
    },
    "transonic_rarefaction": {
        "display": "Transonic rarefaction",
        "left": (0.1, -2.0, 0.1),
        "right": (1.0, -1.0, 1.0),
        "role": "entropy_stress",
    },
    "near_vacuum_expansion": {
        "display": "Near-vacuum expansion",
        "left": (1.0, -2.0, 0.4),
        "right": (1.0, 2.0, 0.4),
        "role": "positivity_entropy_stress",
    },
    "collision": {
        "display": "Collision",
        "left": (1.0, 2.0, 1.0),
        "right": (1.0, -2.0, 1.0),
        "role": "compression_stress",
    },
    "strong_pressure": {
        "display": "Strong pressure ratio",
        "left": (1.0, 0.0, 5.0),
        "right": (1.0, 0.0, 0.05),
        "role": "strong_shock_stress",
    },
}

LABELS = {
    "hllc64": "HLLC-64",
    "muscl64": "MUSCL-HLLC-64",
    "roenet64": "RoeNet-64 adaptation",
    "fno64": "vanilla residual FNO-64",
    "hcfl64_seed0": "HCFL-64 seed 0",
    "hcfl64_seed1": "HCFL-64 seed 1",
    "hcfl64_seed2": "HCFL-64 seed 2",
}
COLORS = {
    "reference": "#B8C2CC",
    "hllc64": "#0072B2",
    "muscl64": "#009E73",
    "roenet64": "#56B4E9",
    "fno64": "#CC79A7",
    "hcfl64_seed0": "#D55E00",
    "hcfl64_seed1": "#E69F00",
    "hcfl64_seed2": "#7B2CBF",
}
LINESTYLES = {
    "reference": "-",
    "hllc64": "--",
    "muscl64": "-.",
    "roenet64": (0, (4, 2)),
    "fno64": (0, (1, 2)),
    "hcfl64_seed0": "-",
    "hcfl64_seed1": (0, (5, 2)),
    "hcfl64_seed2": (0, (3, 1, 1, 1)),
}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def initial_condition(spec: dict[str, Any], cells: int = CELLS) -> torch.Tensor:
    primitive = np.empty((1, cells, 3), dtype=np.float64)
    cut = cells // 2
    primitive[:, :cut] = np.asarray(spec["left"], dtype=np.float64)
    primitive[:, cut:] = np.asarray(spec["right"], dtype=np.float64)
    conservative = paper.euler.base.prim_to_cons(
        primitive[..., 0], primitive[..., 1], primitive[..., 2]
    )
    return torch.from_numpy(conservative).float()


def make_references(output: Path, reference_cells: int) -> dict[str, torch.Tensor]:
    cache = output / f"reference_named_cells{reference_cells}.pt"
    if cache.exists():
        return torch.load(cache, map_location="cpu", weights_only=True)
    references: dict[str, torch.Tensor] = {}
    for name, spec in CASES.items():
        print(json.dumps({"stage": "reference", "case": name}), flush=True)
        initial = initial_condition(spec)
        trajectory, _ = paper.reference_rollout(
            "euler", initial, SNAPSHOTS, reference_cells
        )
        references[name] = torch.cat([initial[:, None], trajectory[:, 1:]], dim=1)
    torch.save(references, cache)
    return references


def state_is_physical(trajectory: torch.Tensor) -> torch.Tensor:
    finite = torch.isfinite(trajectory).all(dim=(-1, -2))
    safe = torch.where(torch.isfinite(trajectory), trajectory, torch.zeros_like(trajectory))
    density = safe[..., 0]
    pressure = paper.euler.base.t_pressure(safe)
    return (
        finite
        & (density >= paper.euler.shared.RHO_FLOOR).all(dim=-1)
        & (pressure >= paper.euler.shared.PRESSURE_FLOOR).all(dim=-1)
    )


def primitive(trajectory: torch.Tensor) -> torch.Tensor:
    safe = torch.where(torch.isfinite(trajectory), trajectory, torch.zeros_like(trajectory))
    density = safe[..., 0].clamp_min(1.0e-12)
    velocity = safe[..., 1] / density
    pressure = paper.euler.base.t_pressure(safe)
    return torch.stack((density, velocity, pressure), dim=-1)


def total_variation(values: torch.Tensor) -> torch.Tensor:
    return torch.abs(values - torch.roll(values, 1, dims=-2)).sum(dim=-2)


def metrics(
    name: str,
    method: str,
    trajectory: torch.Tensor,
    reference: torch.Tensor,
    state_std: torch.Tensor,
    run_stats: dict[str, Any],
    seconds: float,
) -> dict[str, Any]:
    physical_by_snapshot = state_is_physical(trajectory)
    completed = bool(physical_by_snapshot.all())
    finite = bool(torch.isfinite(trajectory).all())
    first_failure = None
    if not completed:
        failed = torch.nonzero(~physical_by_snapshot[0], as_tuple=False)
        first_failure = int(failed[0]) if failed.numel() else 0

    safe = torch.where(torch.isfinite(trajectory), trajectory, torch.zeros_like(trajectory))
    state_error = (safe - reference) / state_std.reshape(1, 1, 1, -1)
    rollout_nrmse = float(torch.sqrt(state_error.double().square().mean())) if completed else None
    rollout_nmae = float(state_error.double().abs().mean()) if completed else None
    final_nrmse = (
        float(torch.sqrt(state_error[:, -1].double().square().mean()))
        if completed
        else None
    )
    final_nmae = float(state_error[:, -1].double().abs().mean()) if completed else None
    predicted_primitive = primitive(safe)
    reference_primitive = primitive(reference)
    primitive_absolute_error = (predicted_primitive - reference_primitive).abs().double()
    primitive_rollout_mae = (
        primitive_absolute_error.mean(dim=(0, 1, 2)) if completed else None
    )
    primitive_final_mae = (
        primitive_absolute_error[:, -1].mean(dim=(0, 1)) if completed else None
    )
    primitive_scale = reference_primitive.abs().mean(dim=(0, 1, 2)).clamp_min(1.0e-8)
    relative_l1 = (
        float(
            (
                (predicted_primitive - reference_primitive).abs()
                / primitive_scale.reshape(1, 1, 1, -1)
            ).mean()
        )
        if completed
        else None
    )

    initial_state = safe[:, 0].double()
    initial_sum = initial_state.sum(dim=-2)
    sums = safe.double().sum(dim=-2)
    # Net momentum can be exactly or nearly zero.  Normalizing its drift by
    # the signed conserved total would therefore turn roundoff into an
    # arbitrarily large number.  Use the extensive componentwise L1 scale.
    component_scale = initial_state.abs().sum(dim=-2)
    mass_scale = initial_state[..., 0].abs().sum(dim=-1, keepdim=True)
    conservation_scale = torch.maximum(component_scale, mass_scale).clamp_min(1.0e-12)
    conservation_drift = float(
        ((sums - initial_sum[:, None]).abs() / conservation_scale[:, None]).max()
    )

    if completed:
        entropy = paper.euler.base.entropy(safe.double()).sum(dim=-1)
        entropy_change = entropy[:, 1:] - entropy[:, :-1]
        entropy_violation_rate = float((entropy_change > 1.0e-8).double().mean())
        maximum_entropy_increase = float(entropy_change.max())
    else:
        # Thermodynamic entropy is undefined after density or pressure becomes
        # nonphysical.  Keep the failure/minimum-state diagnostics, but never
        # turn clamped invalid states into a seemingly quantitative entropy.
        entropy_violation_rate = None
        maximum_entropy_increase = None

    pressure = paper.euler.base.t_pressure(safe)
    minimum_density = float(safe[..., 0].min())
    minimum_pressure = float(pressure.min())
    predicted_tv = total_variation(predicted_primitive[:, -1])
    reference_tv = total_variation(reference_primitive[:, -1]).clamp_min(1.0e-12)
    normalized_tv_excess = (
        float(((predicted_tv - reference_tv) / reference_tv).mean()) if completed else None
    )

    row: dict[str, Any] = {
        "case": name,
        "role": CASES[name]["role"],
        "method": method,
        "label": LABELS[method],
        "cells": CELLS,
        "snapshots": SNAPSHOTS,
        "horizon_multiple": (SNAPSHOTS - 1) / TRAINING_INTERVALS,
        "completed": completed,
        "finite": finite,
        "first_failure_snapshot": first_failure,
        "rollout_nrmse_completed_only": rollout_nrmse,
        "rollout_nmae_completed_only": rollout_nmae,
        "final_nrmse_completed_only": final_nrmse,
        "final_nmae_completed_only": final_nmae,
        "rollout_primitive_mae_density_completed_only": (
            float(primitive_rollout_mae[0]) if primitive_rollout_mae is not None else None
        ),
        "rollout_primitive_mae_velocity_completed_only": (
            float(primitive_rollout_mae[1]) if primitive_rollout_mae is not None else None
        ),
        "rollout_primitive_mae_pressure_completed_only": (
            float(primitive_rollout_mae[2]) if primitive_rollout_mae is not None else None
        ),
        "final_primitive_mae_density_completed_only": (
            float(primitive_final_mae[0]) if primitive_final_mae is not None else None
        ),
        "final_primitive_mae_velocity_completed_only": (
            float(primitive_final_mae[1]) if primitive_final_mae is not None else None
        ),
        "final_primitive_mae_pressure_completed_only": (
            float(primitive_final_mae[2]) if primitive_final_mae is not None else None
        ),
        "primitive_relative_l1_completed_only": relative_l1,
        "maximum_relative_conservation_drift": conservation_drift,
        "saved_snapshot_entropy_violation_rate": entropy_violation_rate,
        "maximum_saved_snapshot_entropy_increase": maximum_entropy_increase,
        "minimum_density": minimum_density,
        "minimum_pressure": minimum_pressure,
        "mean_normalized_final_tv_excess": normalized_tv_excess,
        "rollout_seconds": seconds,
    }
    row.update(run_stats)
    return row


def run_method(
    method: str,
    initial: torch.Tensor,
    models: dict[str, torch.nn.Module],
) -> tuple[torch.Tensor, dict[str, Any]]:
    if method == "hllc64":
        return paper.ssprk2_rollout(
            "euler",
            initial,
            SNAPSHOTS,
            paper.euler.base.DT_SNAPSHOT,
            paper.euler.base.t_hllc,
            0.25,
        )
    if method == "muscl64":
        return paper.ssprk2_rollout(
            "euler",
            initial,
            SNAPSHOTS,
            paper.euler.base.DT_SNAPSHOT,
            paper.euler_muscl_flux,
            0.22,
        )
    if method == "fno64":
        return paper.fno_rollout(models[method], initial, SNAPSHOTS)
    if method == "roenet64":
        state = initial.clone()
        saved = [state.clone()]
        with torch.no_grad():
            for _ in range(1, SNAPSHOTS):
                state = models[method](state, paper.euler.base.DT_SNAPSHOT)
                saved.append(state.clone())
        return torch.stack(saved, dim=1), {}
    if method.startswith("hcfl64_seed"):
        return paper.hcfl_rollout("euler", models[method], initial, SNAPSHOTS)
    raise ValueError(method)


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for method in LABELS:
        selected = [row for row in rows if row["method"] == method]
        completed = [row for row in selected if row["completed"]]
        result[method] = {
            "label": LABELS[method],
            "completed_cases": len(completed),
            "total_cases": len(selected),
            "completion_rate": len(completed) / max(len(selected), 1),
            "mean_rollout_nrmse_completed_only": (
                float(np.mean([row["rollout_nrmse_completed_only"] for row in completed]))
                if completed
                else None
            ),
            "mean_rollout_nmae_completed_only": (
                float(np.mean([row["rollout_nmae_completed_only"] for row in completed]))
                if completed
                else None
            ),
            "worst_rollout_nrmse_completed_only": (
                float(np.max([row["rollout_nrmse_completed_only"] for row in completed]))
                if completed
                else None
            ),
            "worst_rollout_nmae_completed_only": (
                float(np.max([row["rollout_nmae_completed_only"] for row in completed]))
                if completed
                else None
            ),
            "maximum_entropy_increase_completed_only": (
                float(
                    np.max(
                        [row["maximum_saved_snapshot_entropy_increase"] for row in completed]
                    )
                )
                if completed
                else None
            ),
            "maximum_conservation_drift": float(
                np.max([row["maximum_relative_conservation_drift"] for row in selected])
            ),
            "minimum_density": float(np.min([row["minimum_density"] for row in selected])),
            "minimum_pressure": float(np.min([row["minimum_pressure"] for row in selected])),
        }
    return result


def plot_profiles(
    trajectories: dict[str, dict[str, torch.Tensor]],
    references: dict[str, torch.Tensor],
    output: Path,
) -> None:
    x = (np.arange(CELLS) + 0.5) / CELLS
    variables = ("density", "velocity", "pressure")
    figure, axes = plt.subplots(
        len(CASES), 3, figsize=(12.0, 2.35 * len(CASES)), constrained_layout=False
    )
    figure.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.06,
        top=0.88,
        wspace=0.16,
        hspace=0.18,
    )
    plotted_methods = (
        "hllc64",
        "muscl64",
        "roenet64",
        "fno64",
        "hcfl64_seed0",
    )
    for row_index, (name, spec) in enumerate(CASES.items()):
        reference = primitive(references[name])[0, -1].numpy()
        completed: dict[str, bool] = {}
        for method in plotted_methods:
            values = primitive(trajectories[name][method])
            completed[method] = bool(
                torch.isfinite(values).all()
                and (values[..., 0] > 0.0).all()
                and (values[..., 2] > 0.0).all()
            )
        for column, variable in enumerate(variables):
            axis = axes[row_index, column]
            axis.plot(
                x,
                reference[:, column],
                color=COLORS["reference"],
                linewidth=3.2,
                label="HLLC-2048 reference",
                zorder=1,
            )
            for method in plotted_methods:
                if not completed[method]:
                    continue
                values = primitive(trajectories[name][method])[0, -1, :, column].numpy()
                axis.plot(
                    x,
                    values,
                    color=COLORS[method],
                    linestyle=LINESTYLES[method],
                    linewidth=1.45,
                    label=LABELS[method],
                    zorder=3,
                )
            if row_index == 0:
                axis.set_title(variable)
            if column == 0:
                axis.set_ylabel(spec["display"])
                failures = [LABELS[method] for method in plotted_methods if not completed[method]]
                if failures:
                    axis.text(
                        0.02,
                        0.04,
                        "failed: " + ", ".join(failures),
                        transform=axis.transAxes,
                        fontsize=7.2,
                        color="#8B0000",
                        bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none"},
                    )
            if row_index == len(CASES) - 1:
                axis.set_xlabel("x")
            axis.grid(alpha=0.16)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.988),
        ncol=3,
        fontsize=8.5,
        frameon=False,
    )
    figure.suptitle(
        f"64-cell Euler stress suite, t={(SNAPSHOTS - 1) * paper.euler.base.DT_SNAPSHOT:.4f}",
        y=0.925,
    )
    figure.savefig(output / "frozen64_profiles.png", dpi=240)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-cells", type=int, default=REFERENCE_CELLS)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results" / "frozen64")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    paper_results = PAPER_DIR / "results"
    references = make_references(output, args.reference_cells)
    training, _, state_std = paper.training_statistics("euler", 0, paper_results)
    del training
    models: dict[str, torch.nn.Module] = {
        "fno64": paper.load_fno("euler", paper_results, 0),
    }
    roenet_dir = HERE / "results" / "roenet"
    roenet_report = json.loads(
        (roenet_dir / "roenet_adaptation64_report_seed0.json").read_text(
            encoding="utf-8"
        )
    )
    roenet = RoeNetEuler1d(
        paper.training_statistics("euler", 0, paper_results)[1],
        state_std,
        cells=CELLS,
        hidden_waves=int(roenet_report["hidden_waves"]),
        internal_steps=int(roenet_report["internal_steps"]),
    )
    roenet.load_state_dict(
        torch.load(
            roenet_dir / roenet_report["checkpoint"],
            map_location="cpu",
            weights_only=True,
        )
    )
    roenet.eval()
    models["roenet64"] = roenet
    for seed in (0, 1, 2):
        models[f"hcfl64_seed{seed}"] = paper.load_hcfl(
            "euler", seed, paper_results
        )
    methods = tuple(LABELS)
    rows: list[dict[str, Any]] = []
    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    for name, spec in CASES.items():
        initial = initial_condition(spec)
        trajectories[name] = {}
        for method in methods:
            print(json.dumps({"stage": "evaluate", "case": name, "method": method}), flush=True)
            started = time.perf_counter()
            try:
                trajectory, stats = run_method(method, initial, models)
            except Exception as error:  # preserve failures in the output table
                trajectory = torch.full(
                    (1, SNAPSHOTS, CELLS, 3), float("nan"), dtype=torch.float32
                )
                trajectory[:, 0] = initial
                stats = {"runtime_error": f"{type(error).__name__}: {error}"}
            elapsed = time.perf_counter() - started
            trajectories[name][method] = trajectory.cpu()
            rows.append(
                metrics(
                    name,
                    method,
                    trajectory.cpu(),
                    references[name],
                    state_std,
                    stats,
                    elapsed,
                )
            )
    summary = {
        "protocol": json.loads((HERE / "protocol.json").read_text(encoding="utf-8")),
        "reference_cells": args.reference_cells,
        "final_time": (SNAPSHOTS - 1) * paper.euler.base.DT_SNAPSHOT,
        "aggregate": aggregate(rows),
        "rows": rows,
    }
    write_csv(output / "frozen64_metrics.csv", rows)
    (output / "frozen64_summary.json").write_text(
        json.dumps(json_ready(summary), indent=2), encoding="utf-8"
    )
    torch.save(trajectories, output / "frozen64_trajectories.pt")
    plot_profiles(trajectories, references, output)
    print(json.dumps(json_ready(summary["aggregate"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
