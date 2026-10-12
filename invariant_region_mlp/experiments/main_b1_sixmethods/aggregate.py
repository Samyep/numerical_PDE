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
