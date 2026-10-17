"""Build the frozen comparison tables, audits, and publication figures.

The collector deliberately keeps native single-point and amortized function
tasks in separate tables.  It can be rerun while jobs are arriving; missing
methods are omitted rather than represented by invented numbers.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


EXPERIMENT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXPERIMENT_DIR.parents[1] / "results/original_pde_benchmark_comparison"
RAW_DIR = RESULTS_DIR / "raw"
FIGURE_DIR = EXPERIMENT_DIR / "figures"
BENCHMARKS = ("hjb_log_d100", "dpi_ou_gmm_d100", "allen_cahn_d100")
BENCHMARK_LABELS = {
    "hjb_log_d100": "HJB-Log (100D)",
    "dpi_ou_gmm_d100": "OU/GMM HJB (100D)",
    "allen_cahn_d100": "Allen–Cahn (100D)",
}
METHOD_LABELS = {
    "terminal_aware_tangent_actor_mcv": "Tangent-Actor + MCV",
    "deep_bsde_official_pytorch_compatibility_port": "Deep BSDE",
    "dpi": "DPI",
    "pinn_hte": "PINN-HTE",
    "d_dbsde": "D-DBSDE",
    "sdgd_pinn_official": "SDGD-PINN",
    "full_history_mlp_raw": "Full-History MLP",
}


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list, tuple))
                    else value
                    for key, value in row.items()
                }
            )


def _actor_training_records() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(RAW_DIR.glob("actor_*.json")):
        record = _read(path)
        if "phases" in record and record.get("method") == "terminal_aware_tangent_actor_mcv":
            records.append(record)
    return records


def _selection(benchmark: str) -> dict[str, Any] | None:
    path = RESULTS_DIR / f"final_actor_selection_{benchmark}.json"
    return _read(path) if path.exists() else None


def _screen_map(
    benchmark: str,
) -> tuple[dict[tuple[str, int], dict[str, Any]], int | None]:
    path = RESULTS_DIR / f"actor_postselection_screen_{benchmark}.json"
    if not path.exists():
        return {}, None
    data = _read(path)
    return {
        (str(row["trial"]), int(row["iteration"])): row
        for row in data["rows"]
    }, int(data["states"])


def actor_trial_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in _actor_training_records():
        benchmark = str(record["benchmark"])
        trial = str(record["trial"])
        config = record["configuration"]
        selection = _selection(benchmark)
        final_key = None
        if selection:
            selected = selection["selected"]
            final_key = (str(selected["trial"]), int(selected["iteration"]))
        stage_b_path = RESULTS_DIR / f"stage_b_selection_{benchmark}.json"
        stage_b = _read(stage_b_path) if stage_b_path.exists() else {}
        screen, screen_states = _screen_map(benchmark)
        for phase in record["phases"]:
            iteration = int(phase["iteration"])
            # The descriptive screen explicitly evaluates seed-0 checkpoints.
            # Never copy those test metrics onto independently trained seeds.
            screen_row = (
                screen.get((trial, iteration), {})
                if int(record["seed"]) == 0
                else {}
            )
            value_metrics = screen_row.get("value_metrics", {})
            validation = phase["reference_free_policy_validation"]
            training = phase["training"]
            labels = phase["labels"]
            train_label = labels.get("train", {})
            rows.append(
                {
                    "benchmark": benchmark,
                    "trial": trial,
                    "seed": int(record["seed"]),
                    "iteration": iteration,
                    "width": config["width"],
                    "depth": config["depth"],
                    "activation": config["activation"],
                    "training_points": config["training_points"],
                    "paths_per_label": config["paths_per_label"],
                    "label_time_steps": config["label_time_steps"],
                    "optimizer_steps_per_iteration": config["optimizer_steps"],
                    "label_estimator": config["label_estimator"],
                    "parameters": training["parameters"],
                    "cumulative_training_seconds": phase["cumulative_training_seconds"],
                    "phase_seconds": phase["phase_seconds"],
                    "label_seconds": phase["label_seconds_reported"],
                    "fit_seconds": phase["fit_seconds_reported"],
                    "heldout_tangent_label_rmse": training["selected_validation_rmse"],
                    "label_replication_rmse": labels[
                        "validation_label_replication_rmse"
                    ],
                    "reference_free_policy_value": validation["new_mean"],
                    "reference_free_source_minus_new": validation[
                        "source_minus_new_mean"
                    ],
                    "reference_free_policy_state_se": validation[
                        "source_minus_new_state_se"
                    ],
                    "training_label_weight_ess_fraction": train_label.get(
                        "weight_ess_fraction"
                    ),
                    "training_label_peak_gpu_bytes": train_label.get(
                        "peak_gpu_bytes"
                    ),
                    "stage_b_continued": trial
                    in stage_b.get("selected_for_continuation", []),
                    "reference_free_final_selected": final_key == (trial, iteration),
                    "descriptive_postselection_value_rmse": value_metrics.get("rmse"),
                    "descriptive_postselection_relative_l2": value_metrics.get(
                        "relative_l2"
                    ),
                    "descriptive_postselection_bias": value_metrics.get("bias"),
                    "descriptive_postselection_evaluation_seconds": screen_row.get(
                        "evaluation_seconds"
                    ),
                    "descriptive_postselection_states": (
                        screen_states if screen_row else None
                    ),
                    "descriptive_postselection_training_seed": (
                        0 if screen_row else None
                    ),
                    "descriptive_test_used_for_selection": False,
                    "ground_truth_used_for_training_or_selection": False,
                }
            )
    rows.sort(key=lambda row: (
        row["benchmark"], row["seed"], row["cumulative_training_seconds"], row["trial"]
    ))
    return rows


def _actor_function_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(RAW_DIR.glob("actor_evaluation_*.json")):
        data = _read(path)
        corrected = data["corrected_metrics_of_replicate_mean"]
        plain = data["plain_metrics_of_replicate_mean"]
        sampling_mse = _finite(data.get("sampling_variance_component_mse"))
        mse = float(corrected["rmse"]) ** 2
        selection = _selection(str(data["benchmark"]))
        selected = selection["selected"] if selection else {}
        is_selected = (
            str(data["trial"]) == str(selected.get("trial"))
            and int(data["iteration"]) == int(selected.get("iteration", -1))
        )
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": f"{data['trial']}/a{data['iteration']}",
                "seed": data["seed"],
                "protocol": "amortized_function_common_bank",
                "evaluation_role": (
                    "frozen_headline"
                    if is_selected and int(data["time_steps"]) == 64
                    else "frozen_discretization_diagnostic"
                    if is_selected
                    else "diagnostic"
                ),
                "states": data["states"],
                "time_steps": data["time_steps"],
                "paths_per_replicate": data["paths_per_replicate"],
                "replicates": data["replicates"],
                "value_rmse": corrected["rmse"],
                "value_relative_l2": corrected["relative_l2"],
                "value_bias": corrected["bias"],
                "value_mae": corrected["mae"],
                "plain_value_rmse": plain["rmse"],
                "gradient_vector_rmse": data["actor_diagnostics"].get("vector_rmse"),
                "gradient_component_rmse": data["actor_diagnostics"].get(
                    "component_rmse"
                ),
                "direct_value_actor_rmse": data["actor_diagnostics"].get("rmse"),
                "sampling_variance_component_mse": sampling_mse,
                "non_mc_rmse_proxy": math.sqrt(max(mse - (sampling_mse or 0.0), 0.0)),
                "variance_reduction": data["mean_variance_reduction"],
                "weight_ess_fraction": data["mean_weight_ess_fraction"],
                "training_seconds": data["training_seconds"],
                "inference_seconds": data.get(
                    "function_evaluation_seconds", data["evaluation_seconds"]
                ),
                "combined_function_and_native_evaluation_seconds": data[
                    "evaluation_seconds"
                ],
                "diagnostic_and_io_overhead_seconds": data.get(
                    "diagnostic_and_io_overhead_seconds"
                ),
                "total_seconds": data.get(
                    "function_total_seconds", data["total_seconds"]
                ),
                "allocation_total_seconds_including_diagnostics": data.get(
                    "total_seconds_including_diagnostics", data["total_seconds"]
                ),
                "parameter_count": data["parameter_count"],
                "peak_training_gpu_bytes": data.get("peak_training_gpu_bytes"),
                "peak_inference_gpu_bytes": data.get("peak_inference_gpu_bytes"),
                "peak_gpu_bytes": data.get("peak_gpu_bytes"),
                "nonfinite_count": int(data.get("training_nonfinite_count", 0))
                + int(data.get("evaluation_nonfinite_count", 0)),
                "frozen_before_test": data["frozen_before_test_reference"],
                "ground_truth_used_for_selection": False,
            }
        )
    return rows


def _official_function_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(RAW_DIR.glob("official_dpi_*_s*_official_evaluation.json")):
        data = _read(path)
        training_path = Path(str(path).removesuffix("_evaluation.json") + ".json")
        training_record = _read(training_path) if training_path.exists() else {}
        training_peak = training_record.get("peak_gpu_bytes")
        inference_peak = data.get("peak_gpu_bytes_inference")
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": "official",
                "seed": data["seed"],
                "protocol": "amortized_function_common_bank",
                "evaluation_role": "official_fixed",
                "states": data["states"],
                "value_rmse": data["value_rmse"],
                "value_relative_l2": data["value_relative_l2"],
                "value_bias": data["value_bias"],
                "gradient_vector_rmse": data["gradient_vector_rmse"],
                "gradient_component_rmse": data["gradient_component_rmse"],
                "training_seconds": data["training_seconds"],
                "inference_seconds": data["inference_seconds"],
                "diagnostic_gradient_seconds": data.get(
                    "gradient_evaluation_seconds_reported_separately"
                ),
                "total_seconds": data["total_seconds"],
                "parameter_count": data["parameter_count"],
                "peak_training_gpu_bytes": training_peak,
                "peak_inference_gpu_bytes": inference_peak,
                "peak_gpu_bytes": max(
                    [value for value in (training_peak, inference_peak) if value is not None],
                    default=None,
                ),
                "nonfinite_count": data.get("nonfinite_count", 0),
                "ground_truth_used_for_selection": False,
            }
        )
    for path in sorted(RAW_DIR.glob("official_sdgd_s*_official_evaluation.json")):
        data = _read(path)
        training_path = Path(str(path).removesuffix("_evaluation.json") + ".json")
        training_record = _read(training_path) if training_path.exists() else {}
        training_peak = training_record.get("peak_gpu_bytes")
        inference_peak = data.get("peak_gpu_bytes_inference")
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": "official",
                "seed": data["seed"],
                "protocol": "amortized_function_common_bank",
                "evaluation_role": "official_fixed",
                "states": 960,
                "value_rmse": data["common_960_value_rmse"],
                "value_relative_l2": data["common_960_value_relative_l2"],
                "value_bias": data["common_960_value_bias"],
                "gradient_vector_rmse": data["common_960_gradient_vector_rmse"],
                "gradient_component_rmse": data[
                    "common_960_gradient_component_rmse"
                ],
                "training_seconds": data["training_seconds"],
                "inference_seconds": data["inference_seconds"],
                "diagnostic_gradient_seconds": data.get(
                    "gradient_evaluation_seconds_reported_separately"
                ),
                "total_seconds": data["total_seconds"],
                "parameter_count": data["parameter_count"],
                "peak_training_gpu_bytes": training_peak,
                "peak_inference_gpu_bytes": inference_peak,
                "peak_gpu_bytes": max(
                    [value for value in (training_peak, inference_peak) if value is not None],
                    default=None,
                ),
                "nonfinite_count": data.get("nonfinite_count", 0),
                "ground_truth_used_for_selection": False,
                "official_test_labels_reported_each_epoch": True,
            }
        )
    for path in sorted(RAW_DIR.glob("full_history_mlp_*.json")):
        data = _read(path)
        metrics = data["function"]["metrics"]
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": f"n{data['n']}/M{data['M']}",
                "seed": data["seed"],
                "protocol": "amortized_function_common_bank",
                "evaluation_role": "frozen_control_grid",
                "states": data["function"]["states"],
                "value_rmse": metrics["rmse"],
                "value_relative_l2": metrics["relative_l2"],
                "value_bias": metrics["bias"],
                "value_mae": metrics["mae"],
                "training_seconds": 0.0,
                "inference_seconds": data.get(
                    "function_inference_seconds", data["inference_seconds"]
                ),
                "combined_function_and_native_evaluation_seconds": data[
                    "inference_seconds"
                ],
                "total_seconds": data.get(
                    "function_inference_seconds", data["total_seconds"]
                ),
                "parameter_count": 0,
                "peak_training_gpu_bytes": 0,
                "peak_inference_gpu_bytes": max(
                    data["function"].get("peak_gpu_bytes") or 0,
                    data["native"].get("peak_gpu_bytes") or 0,
                ),
                "peak_gpu_bytes": max(
                    data["function"].get("peak_gpu_bytes") or 0,
                    data["native"].get("peak_gpu_bytes") or 0,
                ),
                "nonfinite_count": data["function"]["stats"].get(
                    "nonfinite_count", 0
                )
                + data["native"]["stats"].get("nonfinite_count", 0),
                "ground_truth_used_for_selection": False,
            }
        )
    return rows


def function_rows() -> list[dict[str, Any]]:
    rows = _actor_function_rows() + _official_function_rows()
    rows.sort(key=lambda row: (
        row["benchmark"], row["method"], row["variant"], row["seed"], row.get("time_steps", 0)
    ))
    return rows


def native_rows(function: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(RAW_DIR.glob("actor_evaluation_*.json")):
        data = _read(path)
        selection = _selection(str(data["benchmark"]))
        selected = selection["selected"] if selection else {}
        is_selected = (
            str(data["trial"]) == str(selected.get("trial"))
            and int(data["iteration"]) == int(selected.get("iteration", -1))
        )
        native = data["native"]
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": f"{data['trial']}/a{data['iteration']}/n{data['time_steps']}",
                "seed": data["seed"],
                "protocol": "native_published_point"
                if native["is_published_native_point"]
                else "native_diagnostic_only",
                "published_native_point": native["is_published_native_point"],
                "prediction": native["corrected_mean"],
                "reference": native["reference"],
                "absolute_error": native["corrected_absolute_error"],
                "mc_standard_error": native["corrected_replicate_standard_error"],
                "training_seconds": data["training_seconds"],
                "inference_seconds": data.get(
                    "native_evaluation_seconds", data["evaluation_seconds"]
                ),
                "diagnostic_and_io_overhead_seconds": data.get(
                    "diagnostic_and_io_overhead_seconds"
                ),
                "total_seconds": data.get(
                    "native_total_seconds", data["total_seconds"]
                ),
                "allocation_total_seconds_including_diagnostics": data.get(
                    "total_seconds_including_diagnostics", data["total_seconds"]
                ),
                "parameter_count": data["parameter_count"],
                "peak_training_gpu_bytes": data.get("peak_training_gpu_bytes"),
                "peak_inference_gpu_bytes": data.get("peak_inference_gpu_bytes"),
                "peak_gpu_bytes": data.get("peak_gpu_bytes"),
                "nonfinite_count": int(data.get("training_nonfinite_count", 0))
                + int(data.get("evaluation_nonfinite_count", 0)),
                "evaluation_role": "frozen_selected" if is_selected else "diagnostic",
            }
        )
    for path in sorted(RAW_DIR.glob("deep_bsde_*.json")):
        data = _read(path)
        rows.append(
            {
                "benchmark": data["benchmark"],
                "method": data["method"],
                "method_label": METHOD_LABELS[data["method"]],
                "variant": "official_fixed_iterations",
                "seed": data["seed"],
                "protocol": "native_published_point",
                "published_native_point": True,
                "prediction": data["reported_y0"],
                "reference": data["native_reference"],
                "absolute_error": data["native_absolute_error"],
                "relative_error": data["native_relative_error"],
                "training_seconds": data["training_seconds"],
                "inference_seconds": 0.0,
                "total_seconds": data["total_process_seconds"],
                "parameter_count": data["parameter_count"],
                "nonfinite_count": data["nonfinite_count"],
                "peak_training_gpu_bytes": data.get("peak_gpu_bytes"),
                "peak_inference_gpu_bytes": 0,
                "peak_gpu_bytes": data.get("peak_gpu_bytes"),
                "evaluation_role": "official_fixed",
            }
        )
    for pattern in (
        "official_dpi_*_s*_official_evaluation.json",
        "official_sdgd_s*_official_evaluation.json",
        "full_history_mlp_*.json",
    ):
        for path in sorted(RAW_DIR.glob(pattern)):
            data = _read(path)
            training_record: dict[str, Any] = {}
            if method := data.get("method"):
                if method in {"dpi", "pinn_hte", "d_dbsde", "sdgd_pinn_official"}:
                    training_path = Path(str(path).removesuffix("_evaluation.json") + ".json")
                    if training_path.exists():
                        training_record = _read(training_path)
            native = data["native"]
            prediction = native.get("prediction", native.get("corrected_mean"))
            absolute_error = native.get(
                "absolute_error", native.get("corrected_absolute_error")
            )
            method = data["method"]
            if method == "full_history_mlp_raw":
                inference_peak = max(
                    data["function"].get("peak_gpu_bytes") or 0,
                    data["native"].get("peak_gpu_bytes") or 0,
                )
                training_peak = 0
                nonfinite_count = data["function"]["stats"].get(
                    "nonfinite_count", 0
                ) + data["native"]["stats"].get("nonfinite_count", 0)
            else:
                training_peak = training_record.get("peak_gpu_bytes")
                inference_peak = data.get("peak_gpu_bytes_inference")
                nonfinite_count = data.get("nonfinite_count", 0)
            rows.append(
                {
                    "benchmark": data["benchmark"],
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "variant": (
                        f"n{data['n']}/M{data['M']}"
                        if method == "full_history_mlp_raw"
                        else "official"
                    ),
                    "seed": data["seed"],
                    "protocol": "native_published_point"
                    if data["benchmark"] in {"hjb_log_d100", "allen_cahn_d100"}
                    else "native_diagnostic_only",
                    "published_native_point": data["benchmark"]
                    in {"hjb_log_d100", "allen_cahn_d100"},
                    "prediction": prediction,
                    "reference": native["reference"],
                    "absolute_error": absolute_error,
                    "training_seconds": data.get("training_seconds", 0.0),
                    "inference_seconds": (
                        data.get("native_inference_seconds", data["inference_seconds"])
                        if method == "full_history_mlp_raw"
                        else data.get("inference_seconds", 0.0)
                    ),
                    "diagnostic_gradient_seconds": data.get(
                        "gradient_evaluation_seconds_reported_separately"
                    ),
                    "total_seconds": (
                        data.get("native_inference_seconds", data["total_seconds"])
                        if method == "full_history_mlp_raw"
                        else data["total_seconds"]
                    ),
                    "parameter_count": data["parameter_count"],
                    "peak_training_gpu_bytes": training_peak,
                    "peak_inference_gpu_bytes": inference_peak,
                    "peak_gpu_bytes": max(
                        [
                            value
                            for value in (training_peak, inference_peak)
                            if value is not None
                        ],
                        default=None,
                    ),
                    "nonfinite_count": nonfinite_count,
                    "evaluation_role": "official_fixed"
                    if method != "full_history_mlp_raw"
                    else "frozen_control_grid",
                }
            )
    rows.sort(key=lambda row: (
        row["benchmark"], row["method"], row["variant"], row["seed"]
    ))
    return rows


def training_frontier_rows(
    trials: list[dict[str, Any]], function: list[dict[str, Any]], native: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in trials:
        if row.get("descriptive_postselection_value_rmse") is None:
            continue
        rows.append(
            {
                "benchmark": row["benchmark"],
                "method": "terminal_aware_tangent_actor_mcv",
                "method_label": METHOD_LABELS["terminal_aware_tangent_actor_mcv"],
                "variant": f"{row['trial']}/a{row['iteration']}",
                "seed": row["seed"],
                "protocol": "amortized_function_common_bank",
                "error_metric": "value_rmse",
                "error": row["descriptive_postselection_value_rmse"],
                "training_seconds": row["cumulative_training_seconds"],
                "selection_role": "reference_free_selected"
                if row["reference_free_final_selected"]
                else "postselection_descriptive_only",
            }
        )
    for row in function:
        if row["method"] == "terminal_aware_tangent_actor_mcv":
            if (
                row["evaluation_role"] == "frozen_headline"
                and row.get("time_steps") == 64
            ):
                rows.append(
                    {
                        "benchmark": row["benchmark"],
                        "method": row["method"],
                        "method_label": row["method_label"],
                        "variant": row["variant"],
                        "seed": row["seed"],
                        "protocol": row["protocol"],
                        "error_metric": "value_rmse",
                        "error": row["value_rmse"],
                        "training_seconds": row["training_seconds"],
                        "selection_role": "frozen_headline_actor",
                    }
                )
            continue
        rows.append(
            {
                "benchmark": row["benchmark"],
                "method": row["method"],
                "method_label": row["method_label"],
                "variant": row["variant"],
                "seed": row["seed"],
                "protocol": row["protocol"],
                "error_metric": "value_rmse",
                "error": row["value_rmse"],
                "training_seconds": row["training_seconds"],
                "selection_role": row["evaluation_role"],
            }
        )
    for row in native:
        if row["method"] != "deep_bsde_official_pytorch_compatibility_port":
            continue
        rows.append(
            {
                "benchmark": row["benchmark"],
                "method": row["method"],
                "method_label": row["method_label"],
                "variant": row["variant"],
                "seed": row["seed"],
                "protocol": row["protocol"],
                "error_metric": "native_absolute_error",
                "error": row["absolute_error"],
                "training_seconds": row["training_seconds"],
                "selection_role": row["evaluation_role"],
            }
        )
    rows.sort(key=lambda row: (row["benchmark"], row["training_seconds"], row["method"]))
    return rows


def total_compute_rows(
    trials: list[dict[str, Any]],
    function: list[dict[str, Any]],
    native: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    workloads = (1, 16, 128, 960)
    for row in function:
        states = max(int(row["states"]), 1)
        training = float(row["training_seconds"])
        inference = float(row["inference_seconds"])
        for queries in workloads:
            scaled_inference = inference * queries / states
            is_actor = row["method"] == "terminal_aware_tangent_actor_mcv"
            rows.append(
                {
                    "benchmark": row["benchmark"],
                    "method": row["method"],
                    "method_label": row["method_label"],
                    "variant": row["variant"],
                    "seed": row["seed"],
                    "protocol": row["protocol"],
                    "time_steps": row.get("time_steps"),
                    "queries": queries,
                    "error_metric": "value_rmse",
                    "error": row["value_rmse"],
                    "training_seconds": training,
                    "inference_seconds_scaled": scaled_inference,
                    "total_seconds": training + scaled_inference,
                    "selection_role": row["evaluation_role"],
                    "envelope_eligible": (
                        row["evaluation_role"] == "frozen_headline"
                        if is_actor
                        else True
                    ),
                    "row_kind": "individual_seed",
                    "cost_note": (
                        "measured protocol-specific batched inference scaled "
                        "linearly by query count"
                    ),
                }
            )
    for row in native:
        if row["method"] != "deep_bsde_official_pytorch_compatibility_port":
            continue
        for queries in workloads:
            rows.append(
                {
                    "benchmark": row["benchmark"],
                    "method": row["method"],
                    "method_label": row["method_label"],
                    "variant": row["variant"],
                    "seed": row["seed"],
                    "protocol": "native_single_point_retrained_per_query",
                    "queries": queries,
                    "error_metric": "native_absolute_error",
                    "error": row["absolute_error"],
                    "training_seconds": float(row["training_seconds"]) * queries,
                    "inference_seconds_scaled": 0.0,
                    "total_seconds": float(row["total_seconds"]) * queries,
                    "envelope_eligible": False,
                    "row_kind": "individual_seed",
                    "cost_note": "native Deep BSDE requires an independent training run per state",
                }
            )
    for row in trials:
        error = _finite(row.get("descriptive_postselection_value_rmse"))
        inference = _finite(row.get("descriptive_postselection_evaluation_seconds"))
        states = row.get("descriptive_postselection_states")
        if error is None or inference is None or not states:
            continue
        for queries in workloads:
            scaled_inference = inference * queries / int(states)
            rows.append(
                {
                    "benchmark": row["benchmark"],
                    "method": "terminal_aware_tangent_actor_mcv",
                    "method_label": METHOD_LABELS[
                        "terminal_aware_tangent_actor_mcv"
                    ],
                    "variant": f"{row['trial']}/a{row['iteration']}",
                    "seed": row["seed"],
                    "protocol": "postselection_descriptive_actor_envelope",
                    "queries": queries,
                    "error_metric": "value_rmse",
                    "error": error,
                    "training_seconds": row["cumulative_training_seconds"],
                    "inference_seconds_scaled": scaled_inference,
                    "total_seconds": float(row["cumulative_training_seconds"])
                    + scaled_inference,
                    "cost_note": (
                        "descriptive final-test screen after checkpoint freeze; "
                        "not used to select the reported Actor"
                    ),
                    "selection_role": (
                        "reference_free_selected"
                        if row["reference_free_final_selected"]
                        else "postselection_descriptive_only"
                    ),
                    "envelope_eligible": True,
                    "row_kind": "individual_seed",
                }
            )
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        if row["protocol"] not in {
            "amortized_function_common_bank",
            "postselection_descriptive_actor_envelope",
            "native_single_point_retrained_per_query",
        }:
            continue
        key = (
            row["benchmark"],
            row["method"],
            row["method_label"],
            row["variant"],
            row["protocol"],
            row["queries"],
            row.get("time_steps"),
            row.get("selection_role"),
            bool(row["envelope_eligible"]),
        )
        grouped.setdefault(key, []).append(row)
    configuration_means: list[dict[str, Any]] = []
    for key, members in grouped.items():
        (
            benchmark,
            method,
            method_label,
            variant,
            protocol,
            queries,
            time_steps,
            selection_role,
            envelope_eligible,
        ) = key
        configuration_means.append(
            {
                "benchmark": benchmark,
                "method": method,
                "method_label": method_label,
                "variant": variant,
                "seed": "mean",
                "seed_count": len(members),
                "protocol": protocol,
                "time_steps": time_steps,
                "queries": queries,
                "error_metric": members[0]["error_metric"],
                "error": float(np.mean([float(row["error"]) for row in members])),
                "training_seconds": float(
                    np.mean([float(row["training_seconds"]) for row in members])
                ),
                "inference_seconds_scaled": float(
                    np.mean(
                        [float(row["inference_seconds_scaled"]) for row in members]
                    )
                ),
                "total_seconds": float(
                    np.mean([float(row["total_seconds"]) for row in members])
                ),
                "cost_note": (
                    "configuration mean across independent seeds; budget envelope "
                    "never selects the best seed"
                ),
                "selection_role": selection_role,
                "envelope_eligible": envelope_eligible,
                "row_kind": "configuration_seed_mean",
            }
        )
    rows.extend(configuration_means)
    for benchmark in BENCHMARKS:
        for queries in workloads:
            for protocol in (
                "amortized_function_common_bank",
                "postselection_descriptive_actor_envelope",
            ):
                subset = [
                    row for row in rows
                    if row["benchmark"] == benchmark
                    and row["queries"] == queries
                    and row["protocol"] == protocol
                    and row.get("row_kind") == "configuration_seed_mean"
                    and row.get("envelope_eligible") is True
                ]
                subset.sort(key=lambda row: row["total_seconds"])
                best = float("inf")
                for row in subset:
                    best = min(best, float(row["error"]))
                    row["budget_envelope_error"] = best
    rows.sort(key=lambda row: (
        row["benchmark"], row["queries"], row["protocol"], row["total_seconds"]
    ))
    return rows


def accounting_rows() -> list[dict[str, Any]]:
    path = RESULTS_DIR / "slurm_accounting.json"
    if not path.exists():
        return []
    return _read(path)["records"]


def baseline_audit(function: list[dict[str, Any]], native: list[dict[str, Any]]) -> dict[str, Any]:
    source_path = RESULTS_DIR / "pde_definitions_and_references.json"
    sources = _read(source_path)["source_commits"] if source_path.exists() else {}
    methods: dict[str, Any] = {}
    for method in (
        "deep_bsde_official_pytorch_compatibility_port",
        "dpi",
        "pinn_hte",
        "d_dbsde",
        "sdgd_pinn_official",
    ):
        frows = [row for row in function if row["method"] == method]
        nrows = [row for row in native if row["method"] == method]
        methods[method] = {
            "completed_function_evaluations": len(frows),
            "completed_native_evaluations": len(nrows),
            "benchmarks": sorted({row["benchmark"] for row in frows + nrows}),
            "best_function_rmse": min(
                (float(row["value_rmse"]) for row in frows), default=None
            ),
            "best_native_absolute_error": min(
                (float(row["absolute_error"]) for row in nrows), default=None
            ),
        }
    calibration_trials: list[dict[str, Any]] = []
    for pattern in (
        "official_dpi_*_calibration_evaluation.json",
        "official_dpi_*_full_data_one_iteration_evaluation.json",
        "official_sdgd_*_calibration_evaluation.json",
    ):
        for path in sorted(RAW_DIR.glob(pattern)):
            data = _read(path)
            calibration_trials.append(
                {
                    "artifact": path.name,
                    "method": data["method"],
                    "benchmark": data["benchmark"],
                    "training_seconds": data.get("training_seconds"),
                    "value_rmse": data.get(
                        "value_rmse", data.get("common_960_value_rmse")
                    ),
                    "headline": False,
                }
            )
    accounting = accounting_rows()
    unsuccessful_allocations = [
        {
            "job_id": row.get("JobIDRaw"),
            "job_name": row.get("JobName"),
            "state": row.get("State"),
            "elapsed_seconds": row.get("elapsed_seconds"),
            "exit_code": row.get("ExitCode"),
        }
        for row in accounting
        if not str(row.get("State", "")).startswith(("COMPLETED", "RUNNING", "PENDING"))
    ]
    result = {
        "source_commits": sources,
        "methods": methods,
        "calibration_and_audit_trials": calibration_trials,
        "failed_or_cancelled_allocations": unsuccessful_allocations,
        "fidelity_controls": {
            "baseline_terminal_aware_feature_added": False,
            "baseline_mcv_added": False,
            "deep_bsde_framework_port_algorithmic_changes": [],
            "dpi_non_numerical_change": "Rich progress display disabled under Slurm",
            "dpi_checkpoint_compatibility": "restore pre-PyTorch-2.6 torch.load weights_only=False default for locally generated checkpoints",
            "d_dbsde_launcher_correction": "standalone solver executed once instead of inherited redundant outer N=40 loop",
            "sdgd_posttraining_dependency": "isolated openpyxl==3.1.5 for xlsx serialization",
        },
        "calibration_is_headline_result": False,
    }
    (EXPERIMENT_DIR / "baseline_reproduction_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result


def _plot_accuracy(training: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.1), constrained_layout=True)
    colors = plt.get_cmap("tab10")
    methods = sorted({row["method"] for row in training})
    color = {method: colors(index % 10) for index, method in enumerate(methods)}
    for axis, benchmark in zip(axes, BENCHMARKS, strict=True):
        subset = [row for row in training if row["benchmark"] == benchmark]
        for method in methods:
            local = [row for row in subset if row["method"] == method]
            if not local:
                continue
            for index, row in enumerate(local):
                if row["protocol"] == "native_published_point":
                    marker, size, face = "X", 48, color[method]
                elif row["selection_role"] == "frozen_headline_actor":
                    marker, size, face = "D", 42, color[method]
                elif row["selection_role"] == "reference_free_selected":
                    marker, size, face = "*", 85, color[method]
                elif row["selection_role"] == "postselection_descriptive_only":
                    marker, size, face = "o", 32, "none"
                else:
                    marker, size, face = "s", 34, color[method]
                axis.scatter(
                    [max(float(row["training_seconds"]), 1e-2)],
                    [float(row["error"])],
                    label=METHOD_LABELS.get(method, method) if index == 0 else None,
                    marker=marker,
                    s=size,
                    facecolors=face,
                    edgecolors=color[method],
                    alpha=0.85,
                )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.set_xlabel("cumulative training GPU seconds")
        axis.grid(True, which="both", alpha=0.22)
    axes[0].set_ylabel("error (function RMSE; X marks native absolute error)")
    handles, labels = axes[0].get_legend_handles_labels()
    for axis in axes[1:]:
        local_handles, local_labels = axis.get_legend_handles_labels()
        for handle, label in zip(local_handles, local_labels, strict=True):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    handles.extend(
        [
            Line2D([0], [0], marker="*", linestyle="none", color="0.25", markersize=9),
            Line2D([0], [0], marker="D", linestyle="none", color="0.25", markersize=6),
            Line2D(
                [0], [0], marker="o", linestyle="none", markerfacecolor="none",
                markeredgecolor="0.25", color="0.25", markersize=6,
            ),
            Line2D([0], [0], marker="X", linestyle="none", color="0.25", markersize=6),
        ]
    )
    labels.extend(
        [
            "reference-free selected Actor",
            "frozen headline Actor",
            "descriptive postselection Actor",
            "native single-point task",
        ]
    )
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.savefig(FIGURE_DIR / "accuracy_vs_training_gpu_time.png", dpi=220)
    fig.savefig(FIGURE_DIR / "accuracy_vs_training_gpu_time.pdf")
    plt.close(fig)


def _plot_budget(total: list[dict[str, Any]]) -> None:
    primary_rows = [
        row for row in total
        if row["queries"] == 960
        and row["protocol"] == "amortized_function_common_bank"
        and row.get("row_kind") == "configuration_seed_mean"
        and row.get("envelope_eligible") is True
        and row.get("budget_envelope_error") is not None
    ]
    descriptive_rows = [
        row for row in total
        if row["queries"] == 960
        and row["protocol"] == "postselection_descriptive_actor_envelope"
        and row.get("row_kind") == "configuration_seed_mean"
        and row.get("envelope_eligible") is True
        and row.get("budget_envelope_error") is not None
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.0), constrained_layout=True)
    for axis, benchmark in zip(axes, BENCHMARKS, strict=True):
        subset = sorted(
            (row for row in primary_rows if row["benchmark"] == benchmark),
            key=lambda row: float(row["total_seconds"]),
        )
        if subset:
            axis.step(
                [float(row["total_seconds"]) for row in subset],
                [float(row["budget_envelope_error"]) for row in subset],
                where="post",
                color="#1f4e79",
                linewidth=2.0,
                label="true budget envelope",
            )
            axis.scatter(
                [float(row["total_seconds"]) for row in subset],
                [float(row["error"]) for row in subset],
                c="#d6604d",
                s=25,
                alpha=0.75,
                label="measured configurations",
            )
        descriptive = sorted(
            (row for row in descriptive_rows if row["benchmark"] == benchmark),
            key=lambda row: float(row["total_seconds"]),
        )
        if descriptive:
            axis.step(
                [float(row["total_seconds"]) for row in descriptive],
                [float(row["budget_envelope_error"]) for row in descriptive],
                where="post",
                color="#f28e2b",
                linewidth=1.7,
                linestyle="--",
                label="descriptive Actor test envelope",
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.set_xlabel("total GPU seconds for K=960")
        axis.grid(True, which="both", alpha=0.22)
    axes[0].set_ylabel("function Value RMSE")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(FIGURE_DIR / "total_compute_budget_envelope.png", dpi=220)
    fig.savefig(FIGURE_DIR / "total_compute_budget_envelope.pdf")
    plt.close(fig)


def _plot_actor_scaling(trials: list[dict[str, Any]]) -> None:
    rows = [row for row in trials if row.get("descriptive_postselection_value_rmse")]
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.0), constrained_layout=True)
    trial_names = sorted({str(row["trial"]) for row in rows})
    palette = plt.get_cmap("tab20")
    colors = {
        trial: palette(index % 20) for index, trial in enumerate(trial_names)
    }
    for axis, benchmark in zip(axes, BENCHMARKS, strict=True):
        subset = [row for row in rows if row["benchmark"] == benchmark]
        for trial in trial_names:
            local = sorted(
                (row for row in subset if row["trial"] == trial),
                key=lambda row: int(row["iteration"]),
            )
            if not local:
                continue
            axis.plot(
                [float(row["cumulative_training_seconds"]) for row in local],
                [float(row["descriptive_postselection_value_rmse"]) for row in local],
                marker="o",
                markersize=4.2,
                linewidth=1.2 if len(local) > 1 else 0,
                linestyle="--" if trial == "focused_plain" else "-",
                color=colors[trial],
                alpha=0.78,
                label=trial,
            )
            for row in local:
                if not row["reference_free_final_selected"]:
                    continue
                axis.scatter(
                    row["cumulative_training_seconds"],
                    row["descriptive_postselection_value_rmse"],
                    s=100,
                    marker="*",
                    facecolor=colors[trial],
                    edgecolor="black",
                    linewidth=0.7,
                    zorder=4,
                )
                axis.annotate(
                    f"{row['trial']} a{row['iteration']}",
                    (row["cumulative_training_seconds"], row["descriptive_postselection_value_rmse"]),
                    xytext=(4, 5), textcoords="offset points", fontsize=7,
                )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(BENCHMARK_LABELS[benchmark])
        axis.set_xlabel("Actor training seconds")
        axis.grid(True, which="both", alpha=0.22)
    axes[0].set_ylabel("descriptive postselection Value RMSE")
    handles, labels = axes[0].get_legend_handles_labels()
    for axis in axes[1:]:
        local_handles, local_labels = axis.get_legend_handles_labels()
        for handle, label in zip(local_handles, local_labels, strict=True):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    handles.append(
        Line2D([0], [0], marker="*", linestyle="none", color="black", markersize=9)
    )
    labels.append("reference-free selected")
    fig.legend(handles, labels, loc="outside lower center", ncol=6, frameon=False)
    fig.savefig(FIGURE_DIR / "tangent_actor_hyperparameter_scaling.png", dpi=220)
    fig.savefig(FIGURE_DIR / "tangent_actor_hyperparameter_scaling.pdf")
    plt.close(fig)


def _plot_decomposition(function: list[dict[str, Any]]) -> None:
    rows = [
        row for row in function
        if row["method"] == "terminal_aware_tangent_actor_mcv"
        and row["evaluation_role"] == "frozen_headline"
        and row.get("time_steps") == 64
        and row["seed"] == 0
    ]
    if not rows:
        return
    rows.sort(key=lambda row: BENCHMARKS.index(row["benchmark"]))
    labels = [BENCHMARK_LABELS[row["benchmark"]] for row in rows]
    sampling = [math.sqrt(max(float(row.get("sampling_variance_component_mse") or 0), 0)) for row in rows]
    non_mc = [float(row["non_mc_rmse_proxy"]) for row in rows]
    observed = [float(row["value_rmse"]) for row in rows]
    actor_error = [
        float(
            row.get("gradient_vector_rmse")
            if row.get("gradient_vector_rmse") is not None
            else row.get("direct_value_actor_rmse")
        )
        for row in rows
    ]
    training = [float(row["training_seconds"]) for row in rows]
    inference = [float(row["inference_seconds"]) for row in rows]
    selected_variants = {
        row["benchmark"]: row["variant"]
        for row in rows
    }
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.1), constrained_layout=True)
    width = 0.24
    axes[0, 0].bar(x - width, observed, width, label="observed RMSE", color="#4c78a8")
    axes[0, 0].bar(x, non_mc, width, label="non-MC proxy", color="#72b7b2")
    axes[0, 0].bar(x + width, sampling, width, label="MC sampling component", color="#f58518")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xticks(x, labels, rotation=15, ha="right")
    axes[0, 0].set_ylabel("Value error")
    axes[0, 0].set_title("MC and non-MC error")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(True, axis="y", which="both", alpha=0.22)

    for benchmark in BENCHMARKS:
        variant = selected_variants.get(benchmark)
        local = sorted(
            (
                row for row in function
                if row["method"] == "terminal_aware_tangent_actor_mcv"
                and row["benchmark"] == benchmark
                and row["seed"] == 0
                and row["variant"] == variant
                and row.get("time_steps") in {32, 64, 128}
            ),
            key=lambda row: int(row["time_steps"]),
        )
        if local:
            axes[0, 1].plot(
                [int(row["time_steps"]) for row in local],
                [float(row["value_rmse"]) for row in local],
                marker="o",
                linewidth=1.8,
                label=BENCHMARK_LABELS[benchmark],
            )
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xticks((32, 64, 128))
    axes[0, 1].set_xlabel("Euler time steps")
    axes[0, 1].set_ylabel("Value RMSE")
    axes[0, 1].set_title("Euler sensitivity")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(True, which="both", alpha=0.22)

    axes[1, 0].bar(x, actor_error, color="#b279a2")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xticks(x, labels, rotation=15, ha="right")
    axes[1, 0].set_ylabel("direct Actor RMSE")
    axes[1, 0].set_title("Actor approximation (gradient vector; Allen value)")
    axes[1, 0].grid(True, axis="y", which="both", alpha=0.22)

    cost_width = 0.34
    axes[1, 1].bar(
        x - cost_width / 2,
        training,
        cost_width,
        label="training",
        color="#54a24b",
    )
    axes[1, 1].bar(
        x + cost_width / 2,
        inference,
        cost_width,
        label="frozen evaluation",
        color="#e45756",
    )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xticks(x, labels, rotation=15, ha="right")
    axes[1, 1].set_ylabel("GPU seconds")
    axes[1, 1].set_title("Training and inference cost")
    axes[1, 1].legend(frameon=False, fontsize=8)
    axes[1, 1].grid(True, axis="y", which="both", alpha=0.22)
    fig.savefig(FIGURE_DIR / "error_and_compute_decomposition.png", dpi=220)
    fig.savefig(FIGURE_DIR / "error_and_compute_decomposition.pdf")
    plt.close(fig)


def aggregate() -> dict[str, int]:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    trials = actor_trial_rows()
    function = function_rows()
    native = native_rows(function)
    training = training_frontier_rows(trials, function, native)
    total = total_compute_rows(trials, function, native)
    accounting = accounting_rows()
    _write_csv(EXPERIMENT_DIR / "hyperparameter_trials.csv", trials)
    _write_csv(EXPERIMENT_DIR / "function_approximation_comparison.csv", function)
    _write_csv(EXPERIMENT_DIR / "native_task_comparison.csv", native)
    _write_csv(EXPERIMENT_DIR / "training_compute_frontier.csv", training)
    _write_csv(EXPERIMENT_DIR / "total_compute_frontier.csv", total)
    _write_csv(EXPERIMENT_DIR / "runtime_resource_accounting.csv", accounting)
    baseline_audit(function, native)
    for name in ("pde_definitions_and_references.json", "reference_accuracy_audit.json"):
        source = RESULTS_DIR / name
        if source.exists():
            (EXPERIMENT_DIR / name).write_text(source.read_text())
    _plot_accuracy(training)
    _plot_budget(total)
    _plot_actor_scaling(trials)
    _plot_decomposition(function)
    return {
        "hyperparameter_rows": len(trials),
        "function_rows": len(function),
        "native_rows": len(native),
        "training_frontier_rows": len(training),
        "total_compute_rows": len(total),
        "accounting_rows": len(accounting),
    }


if __name__ == "__main__":
    print(json.dumps(aggregate(), sort_keys=True))
