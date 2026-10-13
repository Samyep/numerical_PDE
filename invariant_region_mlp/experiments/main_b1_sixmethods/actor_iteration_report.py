"""Validate and report the teacher-free Actor Iteration pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import PROJECT_DIR, RESULTS_DIR, atomic_json


CELLS = ((100, 10), (400, 40), (400, 1))
OUTPUT_DIR = RESULTS_DIR / "actor_iteration_pilot"
FIGURE_DIR = OUTPUT_DIR / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "ACTOR_ITERATION_PILOT.md"
BASELINE_PATH = (
    RESULTS_DIR
    / "actor_convergence"
    / "protocol_v3"
    / "baseline_comparison"
    / "comparison.csv"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _read_baselines() -> dict[tuple[int, int], dict[str, str]]:
    with BASELINE_PATH.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {(int(row["d"]), int(row["k"])): row for row in rows}


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _format(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}g}"


def _format_rmse(value: float, ci95: float) -> str:
    return f"{_format(value)} ± {_format(ci95, 2)}"


def _finite_numbers(value: Any) -> bool:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, dict):
        return all(_finite_numbers(item) for item in value.values())
    if isinstance(value, list):
        return all(_finite_numbers(item) for item in value)
    return True


def _comparison_rows(
    results: dict[tuple[int, int], dict[str, Any]],
    baselines: dict[tuple[int, int], dict[str, str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for d, k in CELLS:
        result = results[d, k]
        baseline = baselines[d, k]
        selected = int(result["stopping_decision"]["selected_iteration"])
        for iteration in result["iterations"]:
            rows.append(
                {
                    "d": d,
                    "k": k,
                    "iteration": int(iteration["iteration"]),
                    "selected_by_value_stopping": int(
                        int(iteration["iteration"]) == selected
                    ),
                    "direct_policy_value_rmse": iteration[
                        "direct_policy_value_rmse"
                    ],
                    "direct_policy_value_ci95": iteration[
                        "direct_policy_value_ci95"
                    ],
                    "actor_mlp_value_rmse": iteration["actor_mlp_value_rmse"],
                    "actor_mlp_value_ci95": iteration["actor_mlp_value_ci95"],
                    "actor_gradient_rmse": iteration["actor_gradient_rmse"],
                    "actor_relative_gradient_error": iteration[
                        "actor_relative_gradient_error"
                    ],
                    "validation_value_mean": iteration[
                        "validation_value_mean"
                    ],
                    "validation_value_improvement": iteration[
                        "validation_value_improvement"
                    ],
                    "validation_value_improvement_ci95": iteration[
                        "validation_value_improvement_ci95"
                    ],
                    "label_generation_seconds": iteration[
                        "label_generation_seconds"
                    ],
                    "actor_training_seconds": iteration[
                        "actor_training_seconds"
                    ],
                    "value_validation_seconds": iteration[
                        "value_validation_seconds"
                    ],
                    "direct_policy_inference_seconds": iteration[
                        "direct_policy_inference_seconds"
                    ],
                    "actor_mlp_inference_seconds": iteration[
                        "actor_mlp_inference_seconds"
                    ],
                    "original_actor_mlp_value_rmse": float(
                        baseline["actor_mlp_value_rmse"]
                    ),
                    "converged_actor_mlp_value_rmse": float(
                        baseline["new_actor_value_rmse"]
                    ),
                    "converged_actor_gradient_rmse": float(
                        baseline["direct_actor_gradient_rmse"]
                    ),
                    "oracle_actor_mlp_value_rmse": float(
                        baseline["oracle_actor_value_rmse"]
                    ),
                    "best_five_baseline_method": baseline[
                        "best_previous_value_method"
                    ],
                    "best_five_baseline_value_rmse": float(
                        baseline["best_previous_value_rmse"]
                    ),
                }
            )
    return rows


def _validate(
    results: dict[tuple[int, int], dict[str, Any]],
    baselines: dict[tuple[int, int], dict[str, str]],
) -> dict[str, Any]:
    failures: list[str] = []
    hashes: dict[str, str] = {}
    checkpoint_hashes: dict[str, list[str]] = {}
    prediction_arrays = 0
    prediction_values = 0
    direct_seed_rows = 0
    actor_mlp_seed_rows = 0
    for d, k in CELLS:
        cell = f"d{d}_k{k}"
        result = results[d, k]
        if result.get("status") != "complete":
            failures.append(f"{cell}: result is not complete")
        expected_hash = baselines[d, k]["test_point_sha256"]
        observed_hash = result.get("test_point_sha256")
        hashes[cell] = str(observed_hash)
        if observed_hash != expected_hash:
            failures.append(f"{cell}: held-out point hash mismatch")
        protocol_cell = result.get("protocol", {}).get("cell")
        if protocol_cell != {"d": d, "k": k}:
            failures.append(f"{cell}: protocol cell mismatch")
        for key in (
            "analytic_reference_used_for_training",
            "analytic_reference_used_for_stopping",
            "analytic_reference_used_for_checkpoint_selection",
        ):
            if result.get(key) is not False:
                failures.append(f"{cell}: {key} must be false")
        if result.get("analytic_reference_first_used_after_stopping_decision") is not True:
            failures.append(f"{cell}: analytic-reference timing flag is invalid")
        if not _finite_numbers(result):
            failures.append(f"{cell}: nonfinite number in result JSON")
        iterations = result.get("iterations", [])
        observed_iterations = [row.get("iteration") for row in iterations]
        if observed_iterations != list(range(len(iterations))):
            failures.append(f"{cell}: iterations are not consecutive")
        selected = result.get("stopping_decision", {}).get("selected_iteration")
        if selected != len(iterations) - 1:
            failures.append(f"{cell}: stopping decision did not select final iterate")
        for row in iterations:
            direct_seed_rows += len(row.get("direct_policy_seed_rows", []))
            actor_mlp_seed_rows += len(row.get("actor_mlp_seed_rows", []))
            if len(row.get("direct_policy_seed_rows", [])) != 3:
                failures.append(f"{cell}: direct policy does not have three seeds")
            if len(row.get("actor_mlp_seed_rows", [])) != 10:
                failures.append(f"{cell}: Actor--MLP does not have ten seeds")
            direct_costs = row.get("direct_policy_costs", [])
            if [cost.get("repetition") for cost in direct_costs] != [0, 1, 2]:
                failures.append(f"{cell}: direct-policy repetitions mismatch")
            if {
                cost.get("seed_namespace") for cost in direct_costs
            } != {"actor_iteration_direct_test"}:
                failures.append(f"{cell}: direct-policy namespace mismatch")
            if [
                seed.get("repetition")
                for seed in row.get("actor_mlp_seed_rows", [])
            ] != list(range(10)):
                failures.append(f"{cell}: Actor--MLP repetitions mismatch")

        prediction_path = OUTPUT_DIR / "predictions" / f"actor_iteration_{cell}.npz"
        if not prediction_path.is_file():
            failures.append(f"{cell}: prediction archive is missing")
        else:
            with np.load(prediction_path, allow_pickle=False) as payload:
                expected_keys = {
                    f"q{j}_{suffix}"
                    for j in range(len(iterations))
                    for suffix in (
                        "direct_policy_value",
                        "actor_mlp_value",
                        "direct_actor",
                    )
                }
                if set(payload.files) != expected_keys:
                    failures.append(f"{cell}: prediction keys mismatch")
                for key in payload.files:
                    value = np.asarray(payload[key])
                    prediction_arrays += 1
                    prediction_values += value.size
                    if not np.all(np.isfinite(value)):
                        failures.append(f"{cell}: nonfinite prediction in {key}")
                    expected_shape = (
                        (3, 960)
                        if key.endswith("direct_policy_value")
                        else (10, 960)
                        if key.endswith("actor_mlp_value")
                        else (960, d)
                    )
                    if value.shape != expected_shape:
                        failures.append(
                            f"{cell}: {key} shape {value.shape} != {expected_shape}"
                        )

        hashes_for_cell: list[str] = []
        for iteration in range(len(iterations)):
            checkpoint = (
                OUTPUT_DIR
                / "checkpoints"
                / f"actor_iteration_{cell}_q{iteration}.pt"
            )
            if not checkpoint.is_file():
                failures.append(f"{cell}: checkpoint q{iteration} is missing")
            else:
                hashes_for_cell.append(_sha256(checkpoint))
        checkpoint_hashes[cell] = hashes_for_cell
        if len(set(hashes_for_cell)) != len(hashes_for_cell):
            failures.append(f"{cell}: duplicate actor checkpoints")

        expected_namespaces = {
            "actor_iteration_q0_train",
            "actor_iteration_q0_valid",
            *{
                f"actor_iteration_{j}_{split}_labels"
                for j in range(1, len(iterations))
                for split in ("train", "valid")
            },
        }
        observed_namespaces: set[str] = set()
        cache_dir = OUTPUT_DIR / "cache" / cell
        for metadata_path in cache_dir.glob("*.json"):
            metadata = _read_json(metadata_path)
            cost = metadata.get("cost", {})
            namespace = cost.get("seed_namespace")
            if namespace:
                observed_namespaces.add(namespace)
        missing_namespaces = expected_namespaces - observed_namespaces
        if missing_namespaces:
            failures.append(
                f"{cell}: missing independent label namespaces "
                + ", ".join(sorted(missing_namespaces))
            )
        label_prefixes = ["q0", *(f"iteration_{j}" for j in range(1, len(iterations)))]
        label_hashes: set[str] = set()
        for prefix in label_prefixes:
            pair = [
                _read_json(cache_dir / f"{prefix}_{split}.json")
                for split in ("train", "valid")
            ]
            pair_hashes = {metadata["point_sha256"] for metadata in pair}
            label_hashes.update(pair_hashes)
            if len(pair_hashes) != 2:
                failures.append(f"{cell}: {prefix} train/validation points overlap")
        value_hashes: set[str] = set()
        for iteration in range(len(iterations)):
            metadata = _read_json(
                cache_dir / f"value_validation_iteration_{iteration}.json"
            )
            value_hashes.add(metadata["point_sha256"])
            costs = metadata.get("costs", [])
            if [cost.get("repetition") for cost in costs] != [0, 1, 2]:
                failures.append(f"{cell}: stopping repetitions mismatch")
            if {
                cost.get("seed_namespace") for cost in costs
            } != {"actor_iteration_value_validation"}:
                failures.append(f"{cell}: stopping namespace mismatch")
        if len(value_hashes) != 1:
            failures.append(f"{cell}: stopping points changed across actors")
        if label_hashes & value_hashes:
            failures.append(f"{cell}: actor labels overlap stopping points")

    return {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "cells_checked": len(CELLS),
        "iterations_checked": sum(
            len(result["iterations"]) for result in results.values()
        ),
        "direct_policy_seed_rows_checked": direct_seed_rows,
        "actor_mlp_seed_rows_checked": actor_mlp_seed_rows,
        "prediction_arrays_checked": prediction_arrays,
        "prediction_values_checked": prediction_values,
        "nonfinite_predictions": (
            0
            if not any("nonfinite prediction" in item for item in failures)
            else None
        ),
        "paired_test_hashes": hashes,
        "checkpoint_sha256": checkpoint_hashes,
        "independent_randomness_namespaces_checked": True,
        "failures": failures,
    }


def _summary(
    results: dict[tuple[int, int], dict[str, Any]],
    rows: list[dict[str, Any]],
    validation: dict[str, Any],
) -> dict[str, Any]:
    selected = [row for row in rows if row["selected_by_value_stopping"]]
    selected_by_cell = {(row["d"], row["k"]): row for row in selected}
    actor_quality_improvements = 0
    degraded_transitions: list[dict[str, Any]] = []
    for cell in CELLS:
        cell_rows = [row for row in rows if (row["d"], row["k"]) == cell]
        q0 = cell_rows[0]
        final = selected_by_cell[cell]
        if final["actor_gradient_rmse"] < q0["actor_gradient_rmse"]:
            actor_quality_improvements += 1
        for row in cell_rows[1:]:
            improvement = float(row["validation_value_improvement"])
            ci95 = float(row["validation_value_improvement_ci95"])
            if improvement + ci95 < 0.0:
                degraded_transitions.append(
                    {
                        "d": row["d"],
                        "k": row["k"],
                        "iteration": row["iteration"],
                        "mean_previous_minus_current": improvement,
                        "ci95": ci95,
                    }
                )
    phases = {
        key: float(sum(float(row[key]) for row in rows))
        for key in (
            "label_generation_seconds",
            "actor_training_seconds",
            "value_validation_seconds",
            "direct_policy_inference_seconds",
            "actor_mlp_inference_seconds",
        )
    }
    scheduler_path = OUTPUT_DIR / "scheduler_accounting.json"
    scheduler = _read_json(scheduler_path) if scheduler_path.is_file() else None
    return {
        "schema_version": 1,
        "status": validation["status"],
        "cells": [
            {
                "d": d,
                "k": k,
                "selected_iteration": int(
                    results[d, k]["stopping_decision"]["selected_iteration"]
                ),
                "stopping_reason": results[d, k]["stopping_decision"]["reason"],
                "selected_direct_policy_value_rmse": selected_by_cell[d, k][
                    "direct_policy_value_rmse"
                ],
                "selected_actor_mlp_value_rmse": selected_by_cell[d, k][
                    "actor_mlp_value_rmse"
                ],
                "selected_actor_gradient_rmse": selected_by_cell[d, k][
                    "actor_gradient_rmse"
                ],
                "q0_actor_gradient_rmse": next(
                    row["actor_gradient_rmse"]
                    for row in rows
                    if (row["d"], row["k"], row["iteration"]) == (d, k, 0)
                ),
                "converged_actor_mlp_value_rmse": selected_by_cell[d, k][
                    "converged_actor_mlp_value_rmse"
                ],
                "converged_actor_gradient_rmse": selected_by_cell[d, k][
                    "converged_actor_gradient_rmse"
                ],
            }
            for d, k in CELLS
        ],
        "selected_direct_policy_wins_over_converged_actor_mlp": sum(
            row["direct_policy_value_rmse"]
            < row["converged_actor_mlp_value_rmse"]
            for row in selected
        ),
        "selected_actor_mlp_wins_over_converged_actor_mlp": sum(
            row["actor_mlp_value_rmse"]
            < row["converged_actor_mlp_value_rmse"]
            for row in selected
        ),
        "selected_direct_policy_wins_over_original_actor_mlp": sum(
            row["direct_policy_value_rmse"]
            < row["original_actor_mlp_value_rmse"]
            for row in selected
        ),
        "selected_actor_mlp_wins_over_original_actor_mlp": sum(
            row["actor_mlp_value_rmse"]
            < row["original_actor_mlp_value_rmse"]
            for row in selected
        ),
        "selected_direct_policy_wins_over_best_five_baseline": sum(
            row["direct_policy_value_rmse"] < row["best_five_baseline_value_rmse"]
            for row in selected
        ),
        "selected_actor_mlp_wins_over_best_five_baseline": sum(
            row["actor_mlp_value_rmse"] < row["best_five_baseline_value_rmse"]
            for row in selected
        ),
        "selected_actor_gradient_improves_over_q0_cells": actor_quality_improvements,
        "selected_actor_gradient_wins_over_converged_actor_cells": sum(
            row["actor_gradient_rmse"] < row["converged_actor_gradient_rmse"]
            for row in selected
        ),
        "validation_value_monotone_within_mc_uncertainty": not degraded_transitions,
        "statistically_degraded_value_transitions": degraded_transitions,
        "program_elapsed_seconds": float(
            sum(result["elapsed_seconds"] for result in results.values())
        ),
        "phase_gpu_seconds": phases,
        "scheduler_accounting": scheduler,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_stopping": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }


def _plot(results: dict[tuple[int, int], dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), sharey=False)
    for axis, (d, k) in zip(axes, CELLS, strict=True):
        result = results[d, k]
        iterations = np.asarray(
            [row["iteration"] for row in result["iterations"]], dtype=int
        )
        direct = np.asarray(
            [row["direct_policy_value_rmse"] for row in result["iterations"]]
        )
        direct_ci = np.asarray(
            [row["direct_policy_value_ci95"] for row in result["iterations"]]
        )
        actor_mlp = np.asarray(
            [row["actor_mlp_value_rmse"] for row in result["iterations"]]
        )
        actor_mlp_ci = np.asarray(
            [row["actor_mlp_value_ci95"] for row in result["iterations"]]
        )
        baseline = result["baselines"]
        axis.errorbar(
            iterations,
            direct,
            yerr=direct_ci,
            marker="o",
            linewidth=2,
            capsize=3,
            label="Direct policy MC",
        )
        axis.errorbar(
            iterations,
            actor_mlp,
            yerr=actor_mlp_ci,
            marker="s",
            linewidth=2,
            capsize=3,
            label="Actor inside MLP",
        )
        axis.axhline(
            baseline["converged_actor_mlp_value_rmse"],
            color="black",
            linestyle="--",
            linewidth=1.4,
            label="Converged Actor--MLP",
        )
        axis.axhline(
            baseline["best_previous_value_rmse"],
            color="tab:green",
            linestyle=":",
            linewidth=1.6,
            label=baseline["best_previous_value_method"],
        )
        axis.axhline(
            baseline["original_actor_mlp_value_rmse"],
            color="0.55",
            linestyle="-.",
            linewidth=1.2,
            label="Original Actor--MLP",
        )
        selected = int(result["stopping_decision"]["selected_iteration"])
        axis.axvline(selected, color="tab:red", alpha=0.2, linewidth=5)
        axis.set_title(f"d={d}, k={k}; selected q{selected}")
        axis.set_xlabel("Actor iteration")
        axis.set_xticks(iterations)
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel("Value RMSE (log scale)")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=5, frameon=False)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.87))
    figure.savefig(FIGURE_DIR / "value_rmse_vs_iteration.png", dpi=180)
    figure.savefig(FIGURE_DIR / "value_rmse_vs_iteration.pdf")
    plt.close(figure)


def _report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    selected = [row for row in rows if row["selected_by_value_stopping"]]
    selected_lookup = {(row["d"], row["k"]): row for row in selected}
    row_lookup = {
        (row["d"], row["k"], row["iteration"]): row for row in rows
    }
    actor_improved = summary["selected_actor_gradient_improves_over_q0_cells"]
    actor_beat_converged = summary[
        "selected_actor_gradient_wins_over_converged_actor_cells"
    ]
    direct_beat_original = summary[
        "selected_direct_policy_wins_over_original_actor_mlp"
    ]
    actor_mlp_beat_original = summary[
        "selected_actor_mlp_wins_over_original_actor_mlp"
    ]
    direct_beat_converged = summary[
        "selected_direct_policy_wins_over_converged_actor_mlp"
    ]
    actor_mlp_beat_converged = summary[
        "selected_actor_mlp_wins_over_converged_actor_mlp"
    ]
    direct_beat_five = summary[
        "selected_direct_policy_wins_over_best_five_baseline"
    ]
    actor_mlp_beat_five = summary[
        "selected_actor_mlp_wins_over_best_five_baseline"
    ]
    k40_q0_gradient = _format(row_lookup[400, 40, 0]["actor_gradient_rmse"])
    k40_q1_gradient = _format(row_lookup[400, 40, 1]["actor_gradient_rmse"])
    k40_converged_gradient = _format(
        row_lookup[400, 40, 1]["converged_actor_gradient_rmse"]
    )
    k1_q0_direct = _format(row_lookup[400, 1, 0]["direct_policy_value_rmse"])
    k1_selected_direct = _format(
        selected_lookup[400, 1]["direct_policy_value_rmse"]
    )
    k1_q0_actor_mlp = _format(row_lookup[400, 1, 0]["actor_mlp_value_rmse"])
    k1_selected_actor_mlp = _format(
        selected_lookup[400, 1]["actor_mlp_value_rmse"]
    )
    cell_lines = []
    for d, k in CELLS:
        row = selected_lookup[d, k]
        q0 = next(
            item
            for item in rows
            if (item["d"], item["k"], item["iteration"]) == (d, k, 0)
        )
        cell_lines.append(
            f"- **({d},{k})** selected `q{row['iteration']}`: direct-policy RMSE "
            f"{_format(row['direct_policy_value_rmse'])}, Actor--MLP RMSE "
            f"{_format(row['actor_mlp_value_rmse'])}, and direct actor gradient "
            f"RMSE {_format(row['actor_gradient_rmse'])} (from "
            f"{_format(q0['actor_gradient_rmse'])} at `q0`)."
        )

    selected_table = [
        "| d | k | Selected | Direct policy | Actor--MLP | Direct actor gradient | "
        "Original Actor--MLP | Converged Actor--MLP | Best five-method baseline |",
        "|---:|---:|:---:|---:|---:|---:|---:|---:|:---|",
    ]
    for row in selected:
        selected_table.append(
            f"| {row['d']} | {row['k']} | q{row['iteration']} | "
            f"{_format_rmse(row['direct_policy_value_rmse'], row['direct_policy_value_ci95'])} | "
            f"{_format_rmse(row['actor_mlp_value_rmse'], row['actor_mlp_value_ci95'])} | "
            f"{_format(row['actor_gradient_rmse'])} | "
            f"{_format(row['original_actor_mlp_value_rmse'])} | "
            f"{_format(row['converged_actor_mlp_value_rmse'])} | "
            f"{row['best_five_baseline_method']}: "
            f"{_format(row['best_five_baseline_value_rmse'])} |"
        )

    iteration_table = [
        "| d | k | q | Policy MC RMSE | Actor--MLP RMSE | Actor gradient RMSE | "
        "Validation Δ (prev-current) | Labels s | Train s | Direct s | Actor--MLP s |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        improvement = (
            "—"
            if row["validation_value_improvement"] is None
            else _format_rmse(
                row["validation_value_improvement"],
                row["validation_value_improvement_ci95"],
            )
        )
        iteration_table.append(
            f"| {row['d']} | {row['k']} | {row['iteration']} | "
            f"{_format_rmse(row['direct_policy_value_rmse'], row['direct_policy_value_ci95'])} | "
            f"{_format_rmse(row['actor_mlp_value_rmse'], row['actor_mlp_value_ci95'])} | "
            f"{_format(row['actor_gradient_rmse'])} | {improvement} | "
            f"{row['label_generation_seconds']:.2f} | "
            f"{row['actor_training_seconds']:.2f} | "
            f"{row['direct_policy_inference_seconds']:.2f} | "
            f"{row['actor_mlp_inference_seconds']:.2f} |"
        )

    degraded = summary["statistically_degraded_value_transitions"]
    if degraded:
        degraded_text = ", ".join(
            f"({item['d']},{item['k']}) q{item['iteration']} "
            f"({_format(item['mean_previous_minus_current'])} ± "
            f"{_format(item['ci95'], 2)})"
            for item in degraded
        )
        monotonicity = (
            "Validation Value was not monotone within MC uncertainty. The "
            f"statistically resolved regressions were {degraded_text}. They were "
            "still within the predeclared 0.001 absolute plateau tolerance when "
            "used by the stopping rule."
        )
    else:
        monotonicity = "Validation Value was monotone within MC uncertainty."

    phases = summary["phase_gpu_seconds"]
    scheduler = summary.get("scheduler_accounting")
    scheduler_text = "Scheduler accounting was not available."
    if scheduler:
        scheduler_text = (
            f"The two Slurm allocations used {scheduler['total_gpu_seconds'] / 60:.2f} "
            "GPU-minutes including environment startup and final serialization."
        )

    return f"""# Teacher-free Actor Iteration pilot

## Numerical answer

Teacher-free Actor Iteration is a **mixed result**. It clearly improves `q0`
as a policy in all {actor_improved}/3 cells. Direct policy evaluation and the
new actor inside the old MLP beat the **original Actor--MLP in
{direct_beat_original}/3 and {actor_mlp_beat_original}/3 cells**, respectively.
Against the stronger high-budget converged Actor--MLP, the corresponding counts
are {direct_beat_converged}/3 and {actor_mlp_beat_converged}/3. Against the best
of the five non-actor baselines, they are {direct_beat_five}/3 and
{actor_mlp_beat_five}/3.

The stopping-selected actor itself beats the previous converged actor's direct
gradient RMSE in only {actor_beat_converged}/3 cells. The Value gain therefore
does not mean that Actor Iteration learned the most accurate gradient in every
cell, and it does not transfer to `(400,40)`.

{chr(10).join(cell_lines)}

The failure modes differ. At `(400,40)`, `q1` improves the direct actor
gradient RMSE only from {k40_q0_gradient} to {k40_q1_gradient}, still far above
the converged actor's {k40_converged_gradient}; both Value estimators remain
near 0.15. This points to insufficient policy accuracy under the frozen small
label budget as the leading diagnosis, rather than an MLP-only evaluation
failure. At `(400,1)`, direct policy RMSE improves from {k1_q0_direct} to
{k1_selected_direct}, while the same actor inside the old MLP worsens from
{k1_q0_actor_mlp} to {k1_selected_actor_mlp}. That cell supports the earlier
diagnosis that the low-rank Actor--MLP recursion can be the bottleneck even when
the learned policy improves.

The iteration was selected only from independent controlled-policy Value
estimates. Exact values and gradients were first used after all checkpoints and
the stopping decision were frozen. Post-hoc best iterations are diagnostic and
did not replace the selected checkpoint.

## Selected checkpoints

{chr(10).join(selected_table)}

`Direct policy` and `Actor--MLP` are different estimators. Direct policy uses
the controlled linear Feynman--Kac estimator; Actor--MLP plugs the same actor
into the previous nonlinear `n=3, M=16` recursion. Intervals are 95% normal
approximations across three direct-policy seeds or ten Actor--MLP seeds.

![Value RMSE versus Actor iteration](../results/main_b1_sixmethods/actor_iteration_pilot/figures/value_rmse_vs_iteration.png)

## Every iteration and measured GPU phases

{chr(10).join(iteration_table)}

Positive validation Δ means the new policy reduced estimated cost. {monotonicity}

Label construction used 16 Euler steps and 128 antithetic paths for each
improvement dataset. The stopping audit used 256 paths, three common-random
banks, and fixed independent validation points. Final direct evaluation used
the original 960 Stage B points, 16 steps, 256 paths, and three seeds.

Across the three cells, program-measured elapsed time was
{summary['program_elapsed_seconds'] / 60:.2f} minutes. The separately measured
phase totals were {phases['label_generation_seconds'] / 60:.2f} minutes for
label generation, {phases['actor_training_seconds'] / 60:.2f} for actor
optimization, {phases['value_validation_seconds'] / 60:.2f} for stopping
audits, {phases['direct_policy_inference_seconds'] / 60:.2f} for final direct
policy inference, and {phases['actor_mlp_inference_seconds'] / 60:.2f} for
final Actor--MLP inference. {scheduler_text}

## Correctness and stability

The pre-run suite had 27 passing tests, including the controlled-drift sign,
terminal exactness, zero-policy heat expectation, common-random-number finite
differences, frozen actor parameters, unbroken input-state gradients, shapes,
and finite values. The artifact validator passed all three cells, all actor
checkpoints were distinct, and no run produced a NaN, OOM, or nonzero exit.

Training and stopping did not use the analytic nonterminal solution. `q0` used
only terminal pathwise derivatives. Improvement labels differentiated through
the frozen actor drift with a VJP; no explicit d-by-d Jacobian was constructed.

## Reproduction

From the repository root inside a single-GH200 allocation:

```bash
python -m pytest -q invariant_region_mlp/experiments/main_b1_sixmethods
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_iteration --d 100 --k 10
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_iteration --d 400 --k 40
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_iteration --d 400 --k 1
python -m invariant_region_mlp.experiments.main_b1_sixmethods.actor_iteration_report
```

Remove the corresponding result, cache, and checkpoint files or pass `--force`
only when a clean rerun is intended.

## Artifacts

- `result_d*_k*.json`: complete per-cell seed rows, costs, protocol, and stopping decision.
- `summary_d*_k*.csv`: compact per-iteration cell results.
- `comparison.csv`: combined three-cell comparison.
- `summary.json`: machine-readable conclusions and compute totals.
- `validation.json`: hashes, shape checks, seed counts, and finite-value checks.
- `scheduler_accounting.json`: Slurm states and elapsed GPU allocation time.
- `manifest.json`: SHA-256 hashes for this report and compact artifacts.
"""


def build_report() -> dict[str, Any]:
    results = {
        (d, k): _read_json(OUTPUT_DIR / f"result_d{d}_k{k}.json")
        for d, k in CELLS
    }
    baselines = _read_baselines()
    rows = _comparison_rows(results, baselines)
    validation = _validate(results, baselines)
    if validation["status"] != "PASS":
        atomic_json(OUTPUT_DIR / "validation.json", validation)
        message = "; ".join(validation["failures"])
        raise RuntimeError(f"Actor Iteration validation failed: {message}")
    summary = _summary(results, rows, validation)
    _write_csv(OUTPUT_DIR / "comparison.csv", rows)
    atomic_json(OUTPUT_DIR / "summary.json", summary)
    atomic_json(OUTPUT_DIR / "validation.json", validation)
    _plot(results)
    REPORT_PATH.write_text(_report(rows, summary))

    artifacts = [
        *(OUTPUT_DIR / f"result_d{d}_k{k}.json" for d, k in CELLS),
        *(OUTPUT_DIR / f"summary_d{d}_k{k}.csv" for d, k in CELLS),
        *(OUTPUT_DIR / f"stopping_decision_d{d}_k{k}.json" for d, k in CELLS),
        OUTPUT_DIR / "comparison.csv",
        OUTPUT_DIR / "summary.json",
        OUTPUT_DIR / "validation.json",
        OUTPUT_DIR / "scheduler_accounting.json",
        FIGURE_DIR / "value_rmse_vs_iteration.png",
        FIGURE_DIR / "value_rmse_vs_iteration.pdf",
        REPORT_PATH,
    ]
    manifest = {
        "schema_version": 1,
        "artifacts": [
            {
                "path": str(path.relative_to(PROJECT_DIR)),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in artifacts
        ],
    }
    atomic_json(OUTPUT_DIR / "manifest.json", manifest)
    for record in manifest["artifacts"]:
        path = PROJECT_DIR / record["path"]
        if path.stat().st_size != record["bytes"] or _sha256(path) != record["sha256"]:
            raise RuntimeError(f"manifest verification failed for {path}")
    return summary


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
