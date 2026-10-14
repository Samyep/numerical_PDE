"""Audit and summarize frozen-Actor evaluation sensitivity runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .common import RESULTS_DIR, atomic_json, model_metrics, write_csv


ACTORS = ("q0", "q1", "q1_refresh", "q2", "q3")
NEW_BUDGETS = ("paths1024_steps64", "paths256_steps128")
BASE_BUDGET = "paths256_steps64"


def _audit_prediction(row: dict[str, Any]) -> float:
    path = Path(row["prediction_path"])
    with np.load(path, allow_pickle=False) as payload:
        independent = model_metrics(payload["value"], payload["truth_u"])
    differences = [
        abs(float(row[key]) - float(independent[key]))
        for key in ("relative_l2", "rmse", "bias", "rmse_over_truth_std")
    ]
    maximum = max(differences)
    if maximum > 5e-12:
        raise RuntimeError(f"prediction audit failed for {path}: {maximum}")
    return maximum


def _sample_summary(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(array)),
        "sample_sd": float(np.std(array, ddof=1)) if len(array) > 1 else 0.0,
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
    }


def aggregate() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    for seed in range(3):
        final_path = (
            RESULTS_DIR
            / "runs"
            / "hjb_log_d400"
            / "ours"
            / f"seed{seed}_final.json"
        )
        sensitivity_path = (
            RESULTS_DIR / "evaluation_sensitivity" / f"hjb_log_seed{seed}.json"
        )
        if not final_path.is_file() or not sensitivity_path.is_file():
            raise FileNotFoundError(
                f"missing final or sensitivity record for seed {seed}"
            )
        final = json.loads(final_path.read_text())
        sensitivity = json.loads(sensitivity_path.read_text())
        if sensitivity["source_run_sha256"] != hashlib.sha256(
            final_path.read_bytes()
        ).hexdigest():
            raise RuntimeError(f"source run hash changed for seed {seed}")
        for item in final["evaluations"]:
            row = {
                "seed": seed,
                "actor": item["actor"],
                "budget": BASE_BUDGET,
                "paths": 256,
                "time_steps": 64,
                "path_training_seconds": item["path_training_seconds"],
                "inference_seconds_960": item["inference_seconds_960"],
                "inference_seconds_2048": item["inference_seconds_2048"],
                "mean_mc_standard_error": item["mean_mc_standard_error"],
                "rms_mc_standard_error": "",
                "mean_estimator_variance": "",
                "mc_variance_fraction_of_observed_mse": "",
                "noise_corrected_rmse_diagnostic": "",
                "prediction_path": item["prediction_path"],
                **{
                    key: item[key]
                    for key in (
                        "relative_l2",
                        "rmse",
                        "bias",
                        "rmse_over_truth_std",
                    )
                },
            }
            maximum = _audit_prediction(row)
            audit.append(
                {"prediction_path": row["prediction_path"], "maximum_difference": maximum}
            )
            rows.append(row)
        for item in sensitivity["evaluations"]:
            maximum = _audit_prediction(item)
            audit.append(
                {
                    "prediction_path": item["prediction_path"],
                    "maximum_difference": maximum,
                }
            )
            rows.append({"seed": seed, **item})

    expected = 3 * len(ACTORS) * (1 + len(NEW_BUDGETS))
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} rows, found {len(rows)}")
    output_csv = RESULTS_DIR / "evaluation_sensitivity_summary.csv"
    fields = list(rows[0])
    fields += sorted({key for row in rows for key in row} - set(fields))
    write_csv(output_csv, rows, fields)

    grouped: dict[str, Any] = {}
    for budget in (BASE_BUDGET, *NEW_BUDGETS):
        grouped[budget] = {}
        for actor in ACTORS:
            members = [
                row
                for row in rows
                if row["budget"] == budget and row["actor"] == actor
            ]
            grouped[budget][actor] = {
                metric: _sample_summary([float(row[metric]) for row in members])
                for metric in (
                    "rmse",
                    "relative_l2",
                    "bias",
                    "inference_seconds_960",
                    "inference_seconds_2048",
                )
            }
            if budget != BASE_BUDGET:
                grouped[budget][actor].update(
                    {
                        metric: _sample_summary(
                            [float(row[metric]) for row in members]
                        )
                        for metric in (
                            "rms_mc_standard_error",
                            "mc_variance_fraction_of_observed_mse",
                            "noise_corrected_rmse_diagnostic",
                        )
                    }
                )

    selection: dict[str, Any] = {}
    for budget in (BASE_BUDGET, *NEW_BUDGETS):
        chosen: list[dict[str, Any]] = []
        for seed in range(3):
            candidates = [
                row
                for row in rows
                if row["budget"] == budget and int(row["seed"]) == seed
            ]
            chosen.append(min(candidates, key=lambda row: float(row["rmse"])))
        selection[budget] = {
            "actors_by_seed": [row["actor"] for row in chosen],
            "oracle_selected_rmse": _sample_summary(
                [float(row["rmse"]) for row in chosen]
            ),
            "explicit_test_oracle_selection": True,
        }

    result = {
        "status": "passed",
        "rows": len(rows),
        "prediction_artifacts_audited": len(audit),
        "maximum_metric_difference": max(
            item["maximum_difference"] for item in audit
        ),
        "budgets": grouped,
        "test_oracle_selection": selection,
        "summary_csv": str(output_csv),
        "audit": audit,
    }
    atomic_json(RESULTS_DIR / "evaluation_sensitivity_aggregate.json", result)
    return result


def main() -> None:
    print(json.dumps(aggregate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
