"""Independently recompute reported metrics from the saved prediction arrays."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .common import RESULTS_DIR, atomic_json, model_metrics


def _jsons(prefix: str) -> list[dict]:
    result = []
    for path in sorted(RESULTS_DIR.glob(f"{prefix}_*.json")):
        with path.open() as handle:
            result.append(json.load(handle))
    return result


def main() -> None:
    project_dir = RESULTS_DIR.parent.parent
    training = _jsons("training")
    ground_truth = _jsons("ground_truth")
    max_mean_error = 0.0
    max_delta_error = 0.0
    max_ground_truth_metric_error = 0.0
    arrays_checked = 0
    all_finite = True
    all_variances_nonnegative = True
    all_training_reference_free = True
    all_ground_truth_posthoc = True

    for payload in training:
        with np.load(project_dir / payload["prediction_path"]) as arrays:
            values = [arrays[f"a{index}_working_value"] for index in range(4)]
            for index, (row, value) in enumerate(zip(payload["iteration_rows"], values, strict=True)):
                max_mean_error = max(
                    max_mean_error,
                    abs(float(np.mean(value)) - row["reference_free_working_value_mean"]),
                )
                if index:
                    max_delta_error = max(
                        max_delta_error,
                        abs(float(np.mean(value - values[index - 1])) - row["paired_delta_vs_previous"]),
                    )
                variance = arrays[f"a{index}_mc_variance_of_mean"]
                all_finite &= bool(np.all(np.isfinite(value)) and np.all(np.isfinite(variance)))
                all_variances_nonnegative &= bool(np.all(variance >= 0.0))
                arrays_checked += 2
        all_training_reference_free &= not payload["ground_truth_used_for_checkpoint_or_iteration_selection"]
        all_training_reference_free &= payload["terminal_signal"]["ambient_gradient_matrix_rank"] == payload["spec"]["k"]
        for phase in payload["phase_metadata"]:
            all_training_reference_free &= not phase["ground_truth_used"]
            all_training_reference_free &= not phase["label"]["oracle_subspace_projection_used"]
            all_training_reference_free &= not phase["label"]["analytic_nonterminal_reference_used"]

    for payload in ground_truth:
        with np.load(project_dir / payload["prediction_path"]) as arrays:
            truth = arrays["reference_original_value"]
            for row in payload["rows"]:
                prediction = arrays[f"{row['iteration']}_original_value"]
                recomputed = model_metrics(prediction, truth)
                for key in ("rmse", "relative_l2", "bias", "truth_mean", "truth_std"):
                    max_ground_truth_metric_error = max(
                        max_ground_truth_metric_error, abs(recomputed[key] - row[key])
                    )
                variance = arrays[f"{row['iteration']}_mc_variance_of_mean"]
                all_finite &= bool(
                    np.all(np.isfinite(prediction)) and np.all(np.isfinite(variance))
                )
                all_variances_nonnegative &= bool(np.all(variance >= 0.0))
                arrays_checked += 2
        all_ground_truth_posthoc &= payload["checkpoint_and_reference_free_evaluation_frozen_before_this_step"]
        all_ground_truth_posthoc &= not payload["ground_truth_used_for_selection"]

    passed = all(
        (
            len(training) == 27,
            len(ground_truth) == 9,
            all_finite,
            all_variances_nonnegative,
            all_training_reference_free,
            all_ground_truth_posthoc,
            max_mean_error < 1e-12,
            max_delta_error < 1e-12,
            max_ground_truth_metric_error < 1e-12,
        )
    )
    result = {
        "passed": passed,
        "training_configurations": len(training),
        "ground_truth_runs": len(ground_truth),
        "arrays_checked": arrays_checked,
        "all_arrays_finite": all_finite,
        "all_mc_variances_nonnegative": all_variances_nonnegative,
        "all_training_and_selection_reference_free": all_training_reference_free,
        "all_ground_truth_evaluation_posthoc": all_ground_truth_posthoc,
        "max_recomputed_value_mean_abs_difference": max_mean_error,
        "max_recomputed_paired_delta_abs_difference": max_delta_error,
        "max_recomputed_ground_truth_metric_abs_difference": max_ground_truth_metric_error,
        "prediction_arrays_committed": False,
        "prediction_arrays_regenerated_by": "run_extended.py",
    }
    atomic_json(RESULTS_DIR / "result_integrity.json", result)
    print(json.dumps(result, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
