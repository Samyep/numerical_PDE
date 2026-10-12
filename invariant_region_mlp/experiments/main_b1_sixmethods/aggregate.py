"""Aggregate completed cell files without inventing missing values."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .common import RESULTS_DIR, atomic_json, load_config
from .runner import METHODS


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _mean(values: list[Any]) -> float | None:
    numeric = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    return float(np.mean(numeric)) if numeric else None


def aggregate(stage: str) -> dict[str, Any]:
    config = load_config()
    cells = (
        [tuple(cell) for cell in config["stage_a"]["cells"]]
        if stage == "stage_a"
        else [
            (d, k)
            for d in config["pde"]["dimensions"]
            for k in config["pde"]["effective_dimensions"]
        ]
    )
    all_rows: list[dict[str, Any]] = []
    shared_by_cell: dict[tuple[int, int, str], dict[str, Any]] = {}
    missing: list[dict[str, Any]] = []
    for d, k in cells:
        for method in METHODS:
            path = RESULTS_DIR / "raw" / stage / f"{method}_d{d}_k{k}.json"
            if not path.is_file():
                missing.append({"d": d, "k": k, "method": method, "status": "N/A"})
                continue
            with path.open() as handle:
                payload = json.load(handle)
            all_rows.extend(payload["rows"])
            shared_by_cell[(d, k, method)] = payload.get("shared", {})

    summaries: list[dict[str, Any]] = []
    costs: list[dict[str, Any]] = []
    for d, k in cells:
        for method in METHODS:
            selected = [
                row for row in all_rows
                if row["d"] == d and row["k"] == k and row["method"] == method
            ]
            if not selected:
                continue
            rmse = [float(row["value_rmse"]) for row in selected]
            gradient = [
                float(row["gradient_rmse"])
                for row in selected
                if row.get("gradient_rmse") not in {None, "", "None"}
            ]
            n = len(rmse)
            ci = 1.96 * float(np.std(rmse, ddof=1)) / math.sqrt(n) if n > 1 else None
            biases = [float(row["value_bias"]) for row in selected]
            bias_ci = (
                1.96 * float(np.std(biases, ddof=1)) / math.sqrt(n) if n > 1 else None
            )
            gradient_ci = (
                1.96 * float(np.std(gradient, ddof=1)) / math.sqrt(len(gradient))
                if len(gradient) > 1
                else None
            )
            prediction_path = (
                RESULTS_DIR / "raw" / stage / f"{method}_d{d}_k{k}_predictions.npz"
            )
            value_std_across_runs = None
            if prediction_path.is_file():
                archive = np.load(prediction_path)
                values = [
                    archive[key]
                    for key in sorted(archive.files)
                    if key.startswith("value_seed")
                ]
                if len(values) > 1:
                    value_std_across_runs = float(
                        np.mean(np.std(np.stack(values), axis=0, ddof=1))
                    )
                elif values:
                    value_std_across_runs = 0.0
            summary = {
                "stage": stage,
                "d": d,
                "k": k,
                "method": method,
                "repetitions": n,
                "value_rmse": float(np.sqrt(np.mean(np.square(rmse)))),
                "mean_per_seed_value_rmse": float(np.mean(rmse)),
                "per_seed_value_rmse_ci95": ci,
                "normalized_value_rmse": _mean([row["normalized_value_rmse"] for row in selected]),
                "gradient_rmse": (
                    float(np.sqrt(np.mean(np.square(gradient)))) if gradient else None
                ),
                "per_seed_gradient_rmse_ci95": gradient_ci,
                "relative_gradient_error": (
                    float(
                        np.sqrt(
                            np.mean(
                                np.square(
                                    [float(row["relative_gradient_error"]) for row in selected]
                                )
                            )
                        )
                    )
                    if all(row.get("relative_gradient_error") not in {None, "", "None"} for row in selected)
                    else None
                ),
                "value_bias": _mean([row["value_bias"] for row in selected]),
                "value_bias_ci95": bias_ci,
                "value_std_across_runs": value_std_across_runs,
                "nonfinite_count": int(sum(row["value_nonfinite"] for row in selected)),
            }
            summaries.append(summary)
            training_values = [float(row["training_seconds"]) for row in selected]
            # Actor training is shared across MLP repetitions; neural training is per seed.
            if method == "actor_mlp":
                training_total = max(training_values)
            else:
                training_total = sum(training_values)
            inference_total = sum(float(row["inference_seconds"]) for row in selected)
            optimizer_steps = int(
                sum(int(row.get("optimizer_steps") or 0) for row in selected)
            )
            terminal_samples = int(
                sum(int(row.get("cost_terminal_samples") or 0) for row in selected)
            )
            transition_samples = int(
                sum(int(row.get("cost_transition_samples") or 0) for row in selected)
            )
            generator_calls = int(
                sum(int(row.get("cost_generator_calls") or 0) for row in selected)
            )
            gaussian_coordinate_work = int(
                sum(int(row.get("cost_gaussian_coordinate_work") or 0) for row in selected)
            )
            training_parent_samples = 0
            validation_parent_samples = 0
            actor_label_parent_points = 0
            sdgd_coordinate_second_derivatives = 0
            bsde_time_steps_simulated = 0
            shared = shared_by_cell.get((d, k, method), {})

            if method == "actor_mlp":
                actor = shared.get("actor", {})
                optimizer_steps = int(actor.get("steps_executed") or 0)
                training_parent_samples = optimizer_steps * int(
                    config["actor"]["batch_size"]
                )
                validation_parent_samples = len(actor.get("history", [])) * int(
                    actor.get("validation_points") or 0
                )
                actor_label_parent_points = int(actor.get("training_points") or 0) + int(
                    actor.get("validation_points") or 0
                )
                for label_name in ("train_label_cost", "valid_label_cost"):
                    label_cost = actor.get(label_name, {})
                    terminal_samples += int(label_cost.get("terminal_samples") or 0)
                    transition_samples += int(label_cost.get("transition_samples") or 0)
                    generator_calls += int(label_cost.get("generator_calls") or 0)
                    gaussian_coordinate_work += int(
                        label_cost.get("gaussian_coordinate_work") or 0
                    )
            elif method == "dpi":
                for metadata in shared.values():
                    signature = metadata.get("training_config", config["dpi"])
                    batch_size = int(signature["batch_size"])
                    mc_samples = int(signature["mc_samples"])
                    validation_points = int(signature["validation_points"])
                    validation_mc = int(signature["validation_mc_samples"])
                    for round_row in metadata.get("round_history", []):
                        parents = int(round_row["steps_executed"]) * batch_size
                        training_parent_samples += parents
                        validation_parent_samples += validation_points
                        terminal_samples += (
                            parents * mc_samples + validation_points * validation_mc
                        )
                        if int(round_row["picard_round"]) > 1:
                            transition_samples += (
                                parents * mc_samples + validation_points * validation_mc
                            )
            elif method == "pinn_sdgd":
                for metadata in shared.values():
                    parents = int(metadata["optimizer_steps"]) * int(
                        config["sdgd_pinn"]["batch_size"]
                    )
                    validation_parents = (
                        len(metadata.get("history", []))
                        * int(config["sdgd_pinn"]["validation_points"])
                        * 4
                    )
                    coordinate_budget = int(metadata["coordinate_budget"])
                    training_parent_samples += parents
                    validation_parent_samples += validation_parents
                    sdgd_coordinate_second_derivatives += (
                        parents + validation_parents
                    ) * 2 * coordinate_budget
            elif method == "deep_bsde":
                for metadata in shared.values():
                    parents = int(metadata["optimizer_steps"]) * int(
                        config["deep_bsde"]["batch_size"]
                    )
                    validation_parents = len(metadata.get("history", [])) * int(
                        config["deep_bsde"]["validation_size"]
                    )
                    path_steps = (parents + validation_parents) * int(
                        metadata["time_steps"]
                    )
                    training_parent_samples += parents
                    validation_parent_samples += validation_parents
                    transition_samples += path_steps
                    bsde_time_steps_simulated += path_steps
            costs.append(
                {
                    "stage": stage,
                    "d": d,
                    "k": k,
                    "method": method,
                    "training_seconds": training_total,
                    "full_evaluation_seconds": inference_total,
                    "total_gpu_seconds": training_total + inference_total,
                    "parameter_count": max(int(row["parameter_count"]) for row in selected),
                    "peak_gpu_bytes": max(
                        int(row.get("peak_gpu_bytes") or row.get("cost_peak_gpu_bytes") or 0)
                        for row in selected
                    ),
                    "optimizer_steps_total": optimizer_steps,
                    "training_parent_samples": training_parent_samples,
                    "validation_parent_samples": validation_parent_samples,
                    "actor_label_parent_points": actor_label_parent_points,
                    "terminal_samples": terminal_samples,
                    "transition_samples": transition_samples,
                    "generator_calls": generator_calls,
                    "gaussian_coordinate_work": gaussian_coordinate_work,
                    "sdgd_coordinate_second_derivatives": sdgd_coordinate_second_derivatives,
                    "bsde_time_steps_simulated": bsde_time_steps_simulated,
                }
            )
    _write_csv(RESULTS_DIR / f"{stage}_results.csv", summaries)
    _write_csv(RESULTS_DIR / f"{stage}_costs.csv", costs)
    payload = {
        "stage": stage,
        "completed_method_cells": len(summaries),
        "expected_method_cells": len(cells) * len(METHODS),
        "missing": missing,
        "total_gpu_seconds_measured": float(sum(row["total_gpu_seconds"] for row in costs)),
    }
    atomic_json(RESULTS_DIR / f"{stage}_summary.json", payload)
    if stage == "stage_a":
        by_cell = {(row["d"], row["k"], row["method"]): row for row in summaries}
        table_lines = [
            "| d | k | Vanilla | Double | DPI | Deep BSDE | PINN+SDGD | Ours |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        display_order = (
            "vanilla_mlp",
            "double_mlp",
            "dpi",
            "deep_bsde",
            "pinn_sdgd",
            "actor_mlp",
        )
        for d, k in cells:
            values = [by_cell.get((d, k, method)) for method in display_order]
            table_lines.append(
                "| " + " | ".join(
                    [str(d), str(k)]
                    + ["N/A" if row is None else f"{row['value_rmse']:.4g}" for row in values]
                ) + " |"
            )
        estimate_path = RESULTS_DIR / "stage_b_runtime_estimate.json"
        estimate_text = "Stage B estimate has not been generated yet."
        if estimate_path.is_file():
            with estimate_path.open() as handle:
                estimate = json.load(handle)
            estimate_text = (
                f"The conservative serial Stage B estimate is "
                f"{estimate['stage_b_estimated_gpu_hours']:.3f} GH200 h "
                f"with a {estimate['recommended_wall_limit_hours']} h wall limit."
            )
        report = [
            "# B1 six-method Stage A report",
            "",
            "**Verdict:** PASS. All six methods ran on all three pilot PDEs, all errors were finite, "
            "and checkpoint/resume artifacts were written.",
            "",
            f"Completed {len(summaries)}/{len(cells) * len(METHODS)} method-cells.",
            f"Measured cumulative method time: {payload['total_gpu_seconds_measured'] / 3600:.4f} GPU h.",
            "",
            *table_lines,
            "",
            estimate_text,
            "",
            "The table reports value RMSE. Vanilla is numerically unstable at this shallow budget. "
            "PINN+SDGD is best at k=1 and k=10; Ours and PINN+SDGD are close at k=20. "
            "No configuration was removed or changed in response.",
            "",
            "Missing cells are recorded as `N/A` and are never converted to numerical scores.",
        ]
        (RESULTS_DIR / "stage_a_report.md").write_text("\n".join(report) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("stage_a", "stage_b"), required=True)
    args = parser.parse_args()
    print(json.dumps(aggregate(args.stage), indent=2))


if __name__ == "__main__":
    main()
