"""Independently recompute all Value metrics from committed prediction arrays."""

from __future__ import annotations

import json

import numpy as np

from .common import PREDICTION_DIR, RESULTS_DIR, atomic_json, write_csv


def independent_metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    """Deliberately independent from the metric helper used by the runner."""
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    truth = np.asarray(truth, dtype=np.float64).reshape(-1)
    error = prediction - truth
    rmse = float(np.sqrt(np.dot(error, error) / len(error)))
    truth_norm = float(np.sqrt(np.dot(truth, truth)))
    truth_std = float(np.sqrt(np.mean((truth - truth.mean()) ** 2)))
    return {
        "rmse": rmse,
        "relative_l2": float(np.sqrt(np.dot(error, error)) / truth_norm),
        "policy_regret": float(error.sum() / len(error)),
        "bias": float(error.sum() / len(error)),
        "truth_mean": float(truth.sum() / len(truth)),
        "truth_std": truth_std,
        "rmse_over_truth_std": rmse / truth_std,
    }


def main() -> None:
    rows = []
    comparisons = []
    for problem in ("lqg", "quartic", "bounded"):
        prediction_path = PREDICTION_DIR / f"{problem}.npz"
        raw_path = RESULTS_DIR / f"raw_{problem}.json"
        with np.load(prediction_path) as arrays, raw_path.open() as handle:
            raw = json.load(handle)
            recorded = {row["iteration"]: row for row in raw["rows"]}
            truth = arrays["reference_value"].astype(np.float64)
            for iteration in [row["iteration"] for row in raw["rows"]]:
                metrics = independent_metrics(arrays[f"{iteration}_value"], truth)
                row = {"pde": problem, "iteration": iteration, **metrics}
                rows.append(row)
                for key in ("rmse", "relative_l2", "policy_regret", "rmse_over_truth_std"):
                    comparisons.append(abs(metrics[key] - float(recorded[iteration][key])))
    write_csv(RESULTS_DIR / "recomputed_metrics.csv", rows)
    audit = {
        "method": "Metrics recomputed only from saved policy predictions and reference values in NPZ files.",
        "rows": len(rows),
        "maximum_absolute_difference_from_reported_scalar": float(max(comparisons)),
        "passed": bool(max(comparisons) < 1e-12),
    }
    atomic_json(RESULTS_DIR / "independent_recomputation.json", audit)
    if not audit["passed"]:
        raise AssertionError(audit)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

