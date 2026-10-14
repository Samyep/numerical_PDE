"""Independently recompute all Value metrics from committed prediction arrays."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .common import PREDICTION_DIR, RESULTS_DIR, atomic_json, model_metrics, write_csv


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
            for iteration in ("a0", "a1", "a2", "optimal_reference_policy"):
                metrics = model_metrics(arrays[f"{iteration}_value"], truth)
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

