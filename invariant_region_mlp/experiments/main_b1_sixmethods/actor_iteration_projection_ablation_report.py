"""Validate and report the paired rank-k Actor projection ablation."""

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
import torch

from .actor_convergence import load_converged_actor
from .common import PROJECT_DIR, RESULTS_DIR, array_sha256, atomic_json, fixed_test_points, load_config
from .equation import B1Equation
from .run_actor_iteration_projection_ablation import CELLS, _load_plain_actor


ROOT = RESULTS_DIR / "actor_iteration_projection_ablation"
ACCOUNTING_PATH = ROOT / "scheduler_accounting.json"
FIGURE_DIR = ROOT / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_ACTOR_ITERATION_PROJECTION_ABLATION.md"
STRUCTURES = ("projected", "unprojected")
LABELS = {"projected": "Rank-k projected", "unprojected": "No projection"}
SUFFIX_SHAPES = {
    "raw": lambda d: (3, 960),
    "corrected": lambda d: (3, 960),
    "policy_gap": lambda d: (3, 960),
    "raw_estimator_variance": lambda d: (3, 960),
    "corrected_estimator_variance": lambda d: (3, 960),
    "policy_gap_estimator_variance": lambda d: (3, 960),
    "direct_actor": lambda d: (960, d),
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _finite(value: Any) -> bool:
    if value is None or isinstance(value, (str, bool)):
        return True
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    if isinstance(value, list):
        return all(_finite(item) for item in value)
    if isinstance(value, dict):
        return all(_finite(item) for item in value.values())
    return True


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    return float(np.mean(array)), 1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array))


def _close(
    observed: float,
    expected: float,
    label: str,
    failures: list[str],
    *,
    tolerance: float = 2.0e-12,
) -> float:
    difference = abs(float(observed) - float(expected))
    if difference > tolerance:
        failures.append(f"{label}: difference {difference} exceeds {tolerance}")
    return difference


def _rows(result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["actor"]: row for row in result["rows"]}


def _validate_protocol(result: dict[str, Any], d: int, k: int, failures: list[str]) -> None:
    protocol = result.get("protocol", {})
    if protocol.get("cell") != {"d": d, "k": k}:
        failures.append(f"{d},{k}: protocol cell mismatch")
    expected = {
        "training_seed_repetitions": 1,
        "policy_improvements": 1,
        "q2_q3_or_q4_permitted": False,
        "analytic_reference_first_used_after_both_final_checkpoints_frozen": True,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            failures.append(f"{d},{k}: protocol mismatch for {key}")
    architecture = protocol.get("architecture", {})
    if architecture.get("hidden_width") != 384 or architecture.get("hidden_layers") != 4:
        failures.append(f"{d},{k}: Actor architecture is not 4x384")
    if architecture.get("trainable_parameter_count_is_matched") is not True:
        failures.append(f"{d},{k}: trainable parameter counts are not marked matched")
    structures = protocol.get("structures", {})
    if structures.get("unprojected", {}).get("label_basis_constructed_or_read") is not False:
        failures.append(f"{d},{k}: unprojected protocol reads a label basis")
    q0 = protocol.get("q0", {})
    for key, value in {
        "training_points": 4096,
        "validation_points": 1024,
        "terminal_paths": 256,
        "shared_raw_labels": True,
        "shared_model_initialization": True,
        "shared_optimizer_minibatches": True,
    }.items():
        if q0.get(key) != value:
            failures.append(f"{d},{k}: q0 mismatch for {key}")
    phases = protocol.get("phases", [])
    if [phase.get("name") for phase in phases] != ["improvement", "label_refresh"]:
        failures.append(f"{d},{k}: phase names mismatch")
    for phase in phases:
        for key, value in {
            "training_points": 8192,
            "validation_points": 512,
            "paths_per_label": 512,
            "sde_steps_per_label": 32,
            "shared_points_across_structures": True,
            "common_brownian_numbers_across_structures": True,
            "shared_optimizer_minibatches_across_structures": True,
        }.items():
            if phase.get(key) != value:
                failures.append(f"{d},{k}/{phase.get('name')}: mismatch for {key}")
        budget = phase.get("budget", {})
        if budget.get("max_steps") != 6000 or budget.get("minimum_steps") != 6000:
            failures.append(f"{d},{k}/{phase.get('name')}: not exactly 6000 optimizer steps")
    final = protocol.get("final_evaluation", {})
    for key, value in {
        "test_points": 960,
        "paths": 256,
        "seeds": 3,
        "time_steps": 64,
        "seed_namespace": "actor_iteration_direct_test",
        "martingale_control_variate": True,
        "common_brownian_numbers_across_structures": True,
    }.items():
        if final.get(key) != value:
            failures.append(f"{d},{k}: final evaluation mismatch for {key}")
    initialization = result.get("paired_q0_initialization", {})
    if initialization.get("parameter_states_identical") is not True:
        failures.append(f"{d},{k}: initial trainable states are not identical")
    if initialization.get("projected_parameter_sha256") != initialization.get("unprojected_parameter_sha256"):
        failures.append(f"{d},{k}: initial parameter hashes differ")
    if initialization.get("projected_trainable_parameters") != initialization.get("unprojected_trainable_parameters"):
        failures.append(f"{d},{k}: trainable parameter counts differ")
    if result.get("q2_q3_or_q4_run") is not False:
        failures.append(f"{d},{k}: result indicates q2/q3/q4 was run")
    if result.get("failed_runs") != 0 or result.get("nan_failures") != 0 or result.get("memory_failures") != 0:
        failures.append(f"{d},{k}: result records a failed, NaN, or memory run")


def _validate_training(result: dict[str, Any], d: int, k: int, failures: list[str]) -> dict[str, Any]:
    construction = result["construction"]
    checkpoint_hashes: dict[str, str] = {}
    cache_files = 0
    cache_values = 0
    shared_q0 = result["shared_q0_labels"]
    for split, expected_count in (("train", 4096), ("valid", 1024)):
        record = shared_q0[split]
        cache = RESULTS_DIR.parent / record["cache"]
        if not cache.is_file():
            failures.append(f"{d},{k}/q0/{split}: shared raw-label cache missing")
            continue
        if _sha256(cache) != record.get("cache_sha256"):
            failures.append(f"{d},{k}/q0/{split}: shared cache hash mismatch")
        with np.load(cache, allow_pickle=False) as payload:
            labels = np.asarray(payload["labels"], dtype=np.float64)
            point_hash = str(payload["point_hash"].item())
        if labels.shape != (expected_count, d) or not np.all(np.isfinite(labels)):
            failures.append(f"{d},{k}/q0/{split}: shared raw-label shape or finiteness failure")
        if array_sha256(labels) != record.get("raw_label_sha256"):
            failures.append(f"{d},{k}/q0/{split}: raw-label array hash mismatch")
        if point_hash != record.get("point_sha256"):
            failures.append(f"{d},{k}/q0/{split}: point hash mismatch")

    projected_q0_path = RESULTS_DIR.parent / construction["projected"]["q0"]["checkpoint"]
    projected_q0_metadata: dict[str, Any] | None = None
    if not projected_q0_path.is_file() or _sha256(projected_q0_path) != construction["projected"]["q0"].get("checkpoint_sha256"):
        failures.append(f"{d},{k}/projected/q0: checkpoint missing or hash mismatch")
    else:
        projected_q0, projected_q0_metadata = load_converged_actor(projected_q0_path, torch.device("cpu"))
        if tuple(projected_q0.label_basis.shape) != (d, k):
            failures.append(f"{d},{k}/projected/q0: invalid label basis")
        for split, key in (("train", "train_labels"), ("valid", "validation_labels")):
            if projected_q0_metadata[key].get("point_sha256") != shared_q0[split].get("point_sha256"):
                failures.append(f"{d},{k}/projected/q0/{split}: shared point provenance mismatch")
    unprojected_q0_path = ROOT / "checkpoints" / f"d{d}_k{k}" / "unprojected_q0.pt"
    if not unprojected_q0_path.is_file():
        failures.append(f"{d},{k}/unprojected/q0: checkpoint missing")
    else:
        config = load_config()
        equation = B1Equation(d=d, k=k, coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]))
        unprojected_q0, unprojected_q0_metadata = _load_plain_actor(unprojected_q0_path, equation, torch.device("cpu"))
        if hasattr(unprojected_q0, "label_basis"):
            failures.append(f"{d},{k}/unprojected/q0: checkpoint Actor has label_basis")
        if unprojected_q0_metadata.get("initialization") != result.get("paired_q0_initialization"):
            failures.append(f"{d},{k}/unprojected/q0: initialization record mismatch")
        for split, key in (("train", "train_labels"), ("valid", "validation_labels")):
            if unprojected_q0_metadata[key].get("raw_label_sha256") != shared_q0[split].get("raw_label_sha256"):
                failures.append(f"{d},{k}/unprojected/q0/{split}: raw-label provenance mismatch")
        if projected_q0_metadata is not None and projected_q0_metadata["training"].get("seed") != unprojected_q0_metadata["training"].get("seed"):
            failures.append(f"{d},{k}/q0: optimizer seeds differ")
    for structure in STRUCTURES:
        q0 = construction[structure]["q0"]
        if structure == "unprojected":
            training = q0["training"]
            if training.get("steps_executed", 0) > 6000:
                failures.append(f"{d},{k}/unprojected/q0 exceeds its optimizer cap")
            if q0.get("target_projection") is not None:
                failures.append(f"{d},{k}/unprojected/q0 has a target projection")
        for phase_name in ("improvement", "label_refresh"):
            record = construction[structure][phase_name]
            if record.get("reused_validated_current_checkpoint") is True:
                if (d, k, structure) != (400, 40, "projected"):
                    failures.append(f"{d},{k}/{structure}/{phase_name}: unexpected external reuse")
                continue
            training = record.get("training", {})
            if training.get("steps_executed") != 6000:
                failures.append(f"{d},{k}/{structure}/{phase_name}: did not execute 6000 steps")
            best_step = int(training.get("best_step", -1))
            if best_step < 1 or best_step > 6000:
                failures.append(f"{d},{k}/{structure}/{phase_name}: invalid selected validation step")
            if record.get("source_policy") != "the arm's frozen q0":
                failures.append(f"{d},{k}/{structure}/{phase_name}: wrong source policy")
            if record.get("labels_projected") != (structure == "projected"):
                failures.append(f"{d},{k}/{structure}/{phase_name}: target projection flag mismatch")
            if structure == "unprojected" and (
                record.get("train_target_projection") is not None
                or record.get("validation_target_projection") is not None
            ):
                failures.append(f"{d},{k}/{structure}/{phase_name}: projection metadata present")
            for split, expected_count in (("train", 8192), ("valid", 512)):
                cache = ROOT / "cache" / f"d{d}_k{k}" / f"{structure}_{phase_name}_{split}.npz"
                metadata_path = cache.with_suffix(".json")
                if not cache.is_file() or not metadata_path.is_file():
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: cache missing")
                    continue
                metadata = _load(metadata_path)
                with np.load(cache, allow_pickle=False) as payload:
                    values = np.asarray(payload["values"], dtype=np.float64)
                    labels = np.asarray(payload["labels"], dtype=np.float64)
                    point_hash = str(payload["point_hash"].item())
                    actor_hash = str(payload["actor_hash"].item())
                if values.shape != (expected_count,) or labels.shape != (expected_count, d):
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: cache shape mismatch")
                if not np.all(np.isfinite(values)) or not np.all(np.isfinite(labels)):
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: nonfinite cache")
                record_key = "train_labels" if split == "train" else "validation_labels"
                if point_hash != record[record_key].get("point_sha256"):
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: point hash mismatch")
                if actor_hash != record.get("source_q0_actor_state_sha256"):
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: source Actor hash mismatch")
                if metadata.get("common_brownian_numbers_across_structures") is not True:
                    failures.append(f"{d},{k}/{structure}/{phase_name}/{split}: common random numbers not recorded")
                cache_files += 1
                cache_values += values.size + labels.size

        final_record = construction[structure]["label_refresh"]
        checkpoint = RESULTS_DIR.parent / _rows(result)[structure]["checkpoint"]
        if not checkpoint.is_file():
            failures.append(f"{d},{k}/{structure}: final checkpoint missing")
            continue
        observed_hash = _sha256(checkpoint)
        if observed_hash != _rows(result)[structure]["checkpoint_sha256"]:
            failures.append(f"{d},{k}/{structure}: checkpoint hash mismatch")
        if final_record.get("checkpoint_sha256") not in (None, observed_hash):
            failures.append(f"{d},{k}/{structure}: construction checkpoint hash mismatch")
        if structure == "projected":
            actor, _ = load_converged_actor(checkpoint, torch.device("cpu"))
            if tuple(actor.label_basis.shape) != (d, k):
                failures.append(f"{d},{k}/projected: invalid label basis shape")
        else:
            config = load_config()
            equation = B1Equation(d=d, k=k, coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]))
            actor, _ = _load_plain_actor(checkpoint, equation, torch.device("cpu"))
            if hasattr(actor, "label_basis"):
                failures.append(f"{d},{k}/unprojected: checkpoint Actor has label_basis")
        checkpoint_hashes[f"d{d}_k{k}_{structure}"] = observed_hash
    for phase_name in ("improvement", "label_refresh"):
        projected = construction["projected"][phase_name]
        unprojected = construction["unprojected"][phase_name]
        if not projected.get("reused_validated_current_checkpoint"):
            for split_key in ("train_labels", "validation_labels"):
                if projected[split_key].get("point_sha256") != unprojected[split_key].get("point_sha256"):
                    failures.append(f"{d},{k}/{phase_name}: point banks differ")
                if projected[split_key].get("seed_namespace") != unprojected[split_key].get("seed_namespace"):
                    failures.append(f"{d},{k}/{phase_name}: Brownian namespaces differ")
            if projected.get("optimizer_seed") != unprojected.get("optimizer_seed"):
                failures.append(f"{d},{k}/{phase_name}: optimizer seeds differ")
        elif (d, k) == (400, 40):
            if phase_name == "improvement":
                prior = _load(ROOT.parent / "actor_quality_upgrade" / "upgrade" / "result.json")["training"]
            else:
                prior = _load(ROOT.parent / "actor_q2_control" / "result.json")["branch_training"]["continued_q1"]
            prior_train_hash = prior.get("train_point_sha256") or prior["train_labels"].get("point_sha256")
            prior_valid_hash = prior.get("validation_point_sha256") or prior["validation_labels"].get("point_sha256")
            prior_train_namespace = prior["train_labels"].get("seed_namespace") or prior["train_labels"].get("cost", {}).get("seed_namespace")
            prior_valid_namespace = prior["validation_labels"].get("seed_namespace") or prior["validation_labels"].get("cost", {}).get("seed_namespace")
            prior_optimizer_seed = prior.get("optimizer_seed") or prior["training"].get("seed")
            for observed, expected, label in (
                (unprojected["train_labels"].get("point_sha256"), prior_train_hash, "train point bank"),
                (unprojected["validation_labels"].get("point_sha256"), prior_valid_hash, "validation point bank"),
                (unprojected["train_labels"].get("seed_namespace"), prior_train_namespace, "train Brownian namespace"),
                (unprojected["validation_labels"].get("seed_namespace"), prior_valid_namespace, "validation Brownian namespace"),
                (unprojected.get("optimizer_seed"), prior_optimizer_seed, "optimizer seed"),
            ):
                if observed != expected:
                    failures.append(f"{d},{k}/{phase_name}: external projected {label} mismatch")
            if prior["training"].get("steps_executed") != 6000:
                failures.append(f"{d},{k}/{phase_name}: reused projected arm did not execute 6000 steps")
    return {
        "cache_files": cache_files,
        "cache_values": cache_values,
        "shared_q0_cache_files": 2,
        "checkpoint_sha256": checkpoint_hashes,
    }


def _validate_predictions(result: dict[str, Any], d: int, k: int, failures: list[str]) -> dict[str, Any]:
    path = ROOT / "predictions" / f"paired_predictions_d{d}_k{k}.npz"
    if not path.is_file():
        failures.append(f"{d},{k}: prediction archive missing")
        return {"arrays": 0, "values": 0, "nonfinite": 0, "maximum_metric_difference": math.inf}
    config = load_config()
    equation = B1Equation(d=d, k=k, coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]))
    points = fixed_test_points(d, k, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    if result.get("test_point_sha256") != array_sha256(points):
        failures.append(f"{d},{k}: test-point hash mismatch")
    rows = _rows(result)
    maximum_difference = 0.0
    value_count = 0
    nonfinite = 0
    with np.load(path, allow_pickle=False) as archive:
        if len(archive.files) != 14:
            failures.append(f"{d},{k}: expected 14 prediction arrays, found {len(archive.files)}")
        for structure in STRUCTURES:
            arrays: dict[str, np.ndarray] = {}
            for suffix, shape_fn in SUFFIX_SHAPES.items():
                key = f"{structure}_{suffix}"
                if key not in archive:
                    failures.append(f"{d},{k}: missing {key}")
                    continue
                value = np.asarray(archive[key], dtype=np.float64)
                arrays[suffix] = value
                value_count += value.size
                nonfinite += int(value.size - np.count_nonzero(np.isfinite(value)))
                if value.shape != shape_fn(d):
                    failures.append(f"{d},{k}: {key} shape is {value.shape}")
            if len(arrays) != len(SUFFIX_SHAPES):
                continue
            row = rows[structure]
            direct_error = arrays["direct_actor"] - exact_z
            direct_rmse = math.sqrt(float(np.mean(np.sum(direct_error * direct_error, axis=1))))
            maximum_difference = max(maximum_difference, _close(direct_rmse, row["direct_actor_gradient_rmse"], f"{d},{k}/{structure} direct RMSE", failures))
            reconstructed = []
            for repetition in range(3):
                raw_error = arrays["raw"][repetition] - exact_u
                corrected_error = arrays["corrected"][repetition] - exact_u
                gap = arrays["policy_gap"][repetition]
                identity = corrected_error - gap
                reconstructed.append({
                    "raw_value_rmse": float(np.sqrt(np.mean(raw_error**2))),
                    "corrected_value_rmse": float(np.sqrt(np.mean(corrected_error**2))),
                    "raw_value_bias": float(np.mean(raw_error)),
                    "corrected_value_bias": float(np.mean(corrected_error)),
                    "mean_analytic_policy_gap": float(np.mean(gap)),
                    "analytic_policy_gap_rms": float(np.sqrt(np.mean(gap**2))),
                    "mean_identity_residual": float(np.mean(identity)),
                    "identity_residual_rms": float(np.sqrt(np.mean(identity**2))),
                })
            if len(row.get("seed_rows", [])) != 3:
                failures.append(f"{d},{k}/{structure}: expected three evaluation rows")
            else:
                for repetition, (observed, stored) in enumerate(zip(reconstructed, row["seed_rows"], strict=True)):
                    for metric, value in observed.items():
                        maximum_difference = max(maximum_difference, _close(value, stored[metric], f"{d},{k}/{structure}/seed{repetition}/{metric}", failures))
            for metric in ("raw_value_rmse", "corrected_value_rmse", "corrected_value_bias", "mean_analytic_policy_gap", "analytic_policy_gap_rms", "mean_identity_residual"):
                mean, _ = _mean_ci95([item[metric] for item in reconstructed])
                maximum_difference = max(maximum_difference, _close(mean, row[metric], f"{d},{k}/{structure}/{metric}", failures))
            raw_variance = float(np.mean(arrays["raw_estimator_variance"]))
            corrected_variance = float(np.mean(arrays["corrected_estimator_variance"]))
            factor = raw_variance / max(corrected_variance, np.finfo(np.float64).tiny)
            maximum_difference = max(maximum_difference, _close(factor, row["variance_reduction_factor"], f"{d},{k}/{structure}/variance factor", failures))
            for repetition, cost in enumerate(row.get("costs", [])):
                if cost.get("seed_namespace") != "actor_iteration_direct_test" or cost.get("time_steps") != 64 or cost.get("paths") != 256 or cost.get("repetition") != repetition:
                    failures.append(f"{d},{k}/{structure}/seed{repetition}: evaluation budget mismatch")
    if nonfinite:
        failures.append(f"{d},{k}: {nonfinite} nonfinite predictions")
    return {
        "arrays": 14,
        "values": value_count,
        "nonfinite": nonfinite,
        "archive_sha256": _sha256(path),
        "maximum_metric_difference": maximum_difference,
    }


def _validate_prior_reproduction(results: dict[tuple[int, int], dict[str, Any]], failures: list[str]) -> float:
    current = ROOT / "predictions" / "paired_predictions_d400_k40.npz"
    prior = RESULTS_DIR / "actor_q2_control" / "predictions" / "paired_predictions.npz"
    if not current.is_file() or not prior.is_file():
        failures.append("cannot verify the reused (400,40) projected Actor against its prior archive")
        return math.inf
    maximum = 0.0
    with np.load(current, allow_pickle=False) as left, np.load(prior, allow_pickle=False) as right:
        for suffix in SUFFIX_SHAPES:
            left_key = f"projected_{suffix}"
            right_key = f"continued_q1_{suffix}"
            maximum = max(maximum, float(np.max(np.abs(np.asarray(left[left_key]) - np.asarray(right[right_key])))))
    if maximum != 0.0:
        failures.append(f"reused (400,40) projected Actor does not exactly reproduce prior predictions: {maximum}")
    expected = 0.05481227122676515
    observed = _rows(results[(400, 40)])["projected"]["corrected_value_rmse"]
    _close(observed, expected, "400,40 projected reference RMSE", failures)
    return maximum


def _validate_accounting(failures: list[str]) -> dict[str, Any]:
    if not ACCOUNTING_PATH.is_file():
        failures.append("scheduler accounting is missing")
        return {}
    accounting = _load(ACCOUNTING_PATH)
    jobs = accounting.get("jobs", [])
    if len(jobs) != 1:
        failures.append("expected one scheduler job")
    else:
        job = jobs[0]
        if job.get("state") != "COMPLETED" or job.get("exit_code") != "0:0" or job.get("gpu_count") != 1:
            failures.append("scheduler job did not complete successfully on one GPU")
    return accounting


def _vector_rms(value: np.ndarray) -> float:
    return math.sqrt(float(np.mean(np.sum(value * value, axis=1))))


def _subspace_diagnostic(result: dict[str, Any], d: int, k: int) -> dict[str, float]:
    rows = _rows(result)
    checkpoint = RESULTS_DIR.parent / rows["projected"]["checkpoint"]
    actor, _ = load_converged_actor(checkpoint, torch.device("cpu"))
    basis = actor.label_basis.detach().cpu().numpy().astype(np.float64)
    config = load_config()
    points = fixed_test_points(d, k, int(config["stage_b"]["test_points"]))
    equation = B1Equation(d=d, k=k, coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]))
    exact_z = equation.exact_z(np.zeros(len(points), dtype=np.float64), points)
    with np.load(ROOT / "predictions" / f"paired_predictions_d{d}_k{k}.npz", allow_pickle=False) as archive:
        unprojected = np.asarray(archive["unprojected_direct_actor"], dtype=np.float64)

    def project(value: np.ndarray) -> np.ndarray:
        return (value @ basis) @ basis.T

    exact_scale = _vector_rms(exact_z)
    actor_scale = _vector_rms(unprojected)
    exact_residual = _vector_rms(exact_z - project(exact_z))
    actor_residual = _vector_rms(unprojected - project(unprojected))
    return {
        "exact_t0_gradient_rms": exact_scale,
        "exact_t0_gradient_projection_residual_rms": exact_residual,
        "exact_t0_gradient_relative_projection_residual": exact_residual / max(exact_scale, np.finfo(np.float64).tiny),
        "unprojected_actor_t0_output_rms": actor_scale,
        "unprojected_actor_off_subspace_rms": actor_residual,
        "unprojected_actor_relative_off_subspace_rms": actor_residual / max(actor_scale, np.finfo(np.float64).tiny),
    }


def _validate(results: dict[tuple[int, int], dict[str, Any]]) -> dict[str, Any]:
    failures: list[str] = []
    arrays = values = nonfinite = cache_files = cache_values = q0_cache_files = 0
    maximum_difference = 0.0
    archive_hashes: dict[str, str] = {}
    checkpoint_hashes: dict[str, str] = {}
    for d, k in CELLS:
        result = results[(d, k)]
        if result.get("status") != "complete" or not _finite(result):
            failures.append(f"{d},{k}: result is incomplete or nonfinite")
        if set(_rows(result)) != set(STRUCTURES):
            failures.append(f"{d},{k}: expected exactly the projected and unprojected rows")
        _validate_protocol(result, d, k, failures)
        training = _validate_training(result, d, k, failures)
        predictions = _validate_predictions(result, d, k, failures)
        arrays += predictions["arrays"]
        values += predictions["values"]
        nonfinite += predictions["nonfinite"]
        maximum_difference = max(maximum_difference, predictions["maximum_metric_difference"])
        archive_hashes[f"d{d}_k{k}"] = predictions.get("archive_sha256", "")
        cache_files += training["cache_files"]
        cache_values += training["cache_values"]
        q0_cache_files += training["shared_q0_cache_files"]
        checkpoint_hashes.update(training["checkpoint_sha256"])
    prior_difference = _validate_prior_reproduction(results, failures)
    accounting = _validate_accounting(failures)
    return {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "cells_checked": len(CELLS),
        "actor_evaluations_checked": 2 * len(CELLS),
        "evaluation_seed_records_checked": 6 * len(CELLS),
        "prediction_arrays_checked": arrays,
        "prediction_values_checked": values,
        "nonfinite_predictions": nonfinite,
        "maximum_recomputed_metric_difference": maximum_difference,
        "prior_projected_d400_k40_max_abs_difference": prior_difference,
        "prediction_archive_sha256": archive_hashes,
        "cache_files_checked": cache_files,
        "cache_values_checked": cache_values,
        "shared_q0_cache_files_checked": q0_cache_files,
        "checkpoint_sha256": checkpoint_hashes,
        "scheduler_jobs_checked": len(accounting.get("jobs", [])),
    }


def _summary(results: dict[tuple[int, int], dict[str, Any]], validation: dict[str, Any]) -> dict[str, Any]:
    cells = []
    for d, k in CELLS:
        result = results[(d, k)]
        rows = _rows(result)
        cell = {
            "d": d,
            "k": k,
            "parameter_count": result["paired_q0_initialization"]["projected_trainable_parameters"],
            "projected": {key: rows["projected"][key] for key in (
                "direct_actor_gradient_rmse", "raw_value_rmse", "corrected_value_rmse", "corrected_value_ci95",
                "corrected_value_bias", "mean_analytic_policy_gap", "analytic_policy_gap_rms",
                "variance_reduction_factor", "total_evaluation_seconds",
            )},
            "unprojected": {key: rows["unprojected"][key] for key in (
                "direct_actor_gradient_rmse", "raw_value_rmse", "corrected_value_rmse", "corrected_value_ci95",
                "corrected_value_bias", "mean_analytic_policy_gap", "analytic_policy_gap_rms",
                "variance_reduction_factor", "total_evaluation_seconds",
            )},
            "comparison": result["comparison"],
            "subspace_diagnostic": _subspace_diagnostic(result, d, k),
            "elapsed_seconds": result["elapsed_seconds"],
        }
        cells.append(cell)
    wins = sum(cell["comparison"]["unprojected_beats_projected_on_value_rmse"] for cell in cells)
    return {
        "schema_version": 1,
        "status": validation["status"],
        "cells": cells,
        "headline": {
            "unprojected_value_rmse_wins": wins,
            "projected_value_rmse_wins": len(cells) - wins,
            "all_initial_parameter_states_identical": True,
            "single_training_seed_per_structure": True,
        },
        "scheduler_accounting": _load(ACCOUNTING_PATH),
    }


def _write_comparison(summary: dict[str, Any]) -> None:
    rows = []
    for cell in summary["cells"]:
        for structure in STRUCTURES:
            rows.append({"d": cell["d"], "k": cell["k"], "structure": structure, **cell[structure]})
    with (ROOT / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _plot(summary: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    cells = summary["cells"]
    labels = [f"({cell['d']},{cell['k']})" for cell in cells]
    x = np.arange(len(cells))
    width = 0.36
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
    metrics = (
        ("corrected_value_rmse", "Martingale-CV value RMSE"),
        ("mean_analytic_policy_gap", "Mean analytic policy gap"),
        ("direct_actor_gradient_rmse", "Direct t=0 gradient RMSE"),
    )
    for axis, (metric, ylabel) in zip(axes, metrics, strict=True):
        for index, structure in enumerate(STRUCTURES):
            values = [cell[structure][metric] for cell in cells]
            kwargs: dict[str, Any] = {}
            if metric == "corrected_value_rmse":
                kwargs["yerr"] = [cell[structure]["corrected_value_ci95"] for cell in cells]
                kwargs["capsize"] = 3
            axis.bar(x + (index - 0.5) * width, values, width, label=LABELS[structure], **kwargs)
        axis.set_xticks(x, labels)
        axis.set_ylabel(ylabel)
        axis.grid(True, axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    figure.suptitle("B1 one-improvement Actor: rank-k projection ablation")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "actor_iteration_projection_ablation.png", dpi=180)
    figure.savefig(FIGURE_DIR / "actor_iteration_projection_ablation.pdf")
    plt.close(figure)


def _fmt(value: float, digits: int = 5) -> str:
    return f"{value:.{digits}f}"


def _report(summary: dict[str, Any], validation: dict[str, Any]) -> str:
    lines = [
        "| Cell | Actor output | Direct grad. RMSE | Raw value RMSE | CV value RMSE | CV bias | Mean policy gap | Gap RMS | CV variance factor |",
        "|:---:|:---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    change_lines = []
    subspace_lines = []
    for cell in summary["cells"]:
        for structure in STRUCTURES:
            row = cell[structure]
            lines.append(
                f"| ({cell['d']},{cell['k']}) | {LABELS[structure]} | {_fmt(row['direct_actor_gradient_rmse'])} | "
                f"{_fmt(row['raw_value_rmse'])} | {_fmt(row['corrected_value_rmse'])} ± {_fmt(row['corrected_value_ci95'])} | "
                f"{_fmt(row['corrected_value_bias'])} | {_fmt(row['mean_analytic_policy_gap'])} | "
                f"{_fmt(row['analytic_policy_gap_rms'])} | {row['variance_reduction_factor']:.2f}x |"
            )
        comparison = cell["comparison"]
        paired = comparison["paired_differences"]["corrected_value_rmse"]
        direction = "lower" if comparison["unprojected_beats_projected_on_value_rmse"] else "higher"
        change_lines.append(
            f"- `({cell['d']},{cell['k']})`: removing projection gives {100 * abs(comparison['unprojected_relative_value_rmse_change']):.1f}% "
            f"{direction} CV RMSE; paired `no projection - projected` = {_fmt(paired['mean'])} ± {_fmt(paired['ci95'])}."
        )
        diagnostic = cell["subspace_diagnostic"]
        subspace_lines.append(
            f"| ({cell['d']},{cell['k']}) | {diagnostic['exact_t0_gradient_projection_residual_rms']:.3e} | "
            f"{diagnostic['unprojected_actor_off_subspace_rms']:.5f} | "
            f"{100 * diagnostic['unprojected_actor_relative_off_subspace_rms']:.3f}% |"
        )
    wins = summary["headline"]["unprojected_value_rmse_wins"]
    if wins == len(CELLS):
        conclusion = "Removing the rank-k Actor projection improves all three tested cells in this paired run."
    elif wins == 0:
        conclusion = "The rank-k Actor projection improves all three tested cells in this paired run."
    else:
        conclusion = f"The projection ablation is configuration-dependent: no projection wins {wins} of {len(CELLS)} tested cells."
    job = summary["scheduler_accounting"]["jobs"][0]
    return f"""# B1 Actor Iteration rank-k projection ablation

## Numerical answer

{conclusion}

{chr(10).join(change_lines)}

{chr(10).join(lines)}

![projection ablation](../results/main_b1_sixmethods/actor_iteration_projection_ablation/figures/actor_iteration_projection_ablation.png)

The analytic trajectory policy gap is the main check on the value ranking: it
uses the same controlled paths but does not rely on the martingale-CV value
noise. Direct `t=0` gradient RMSE is also reported as a post hoc diagnostic.

## What the no-projection Actor learned

| Cell | Exact z* residual outside fitted basis | No-projection Actor off-basis RMS | Relative off-basis RMS |
|:---:|---:|---:|---:|
{chr(10).join(subspace_lines)}

At the held-out starting states, the fitted terminal-gradient basis contains
the analytic gradient to numerical precision. The unconstrained Actor also
places only a small fraction of its output outside that basis. The final
differences at `t=0` therefore cannot be explained by useful true-gradient
signal that the fixed basis excluded. The trajectory policy gap remains the
relevant diagnostic away from the starting states.

## Fair comparison

For each cell, both arms use the same raw q0 terminal labels, point banks,
4x384 trainable network, initial trainable parameter state, optimizer seed and
minibatch sequence. Policy-evaluation labels use common Brownian random
numbers. Each improvement phase uses 8,192 training and 512 validation points,
512 paths per label, 32 Euler steps, and exactly 6,000 optimizer steps. Final
checkpoints use the same validation-based selection rule within those 6,000
steps. Final evaluation uses the same 960 points, three Brownian banks, 256
paths, 64 Euler steps, and martingale control variate.

The projected arm applies the terminal-label-derived rank-k map to its output
and targets. The no-projection arm has a plain ambient d-dimensional output,
uses raw targets, and does not read or construct that basis. Each arm evaluates
its own frozen q0 policy, so the simulated state paths and numerical labels can
differ after the matched random inputs are applied. That difference is part of
the end-to-end structural ablation.

The second training phase refreshes labels from the same frozen q0 and
continues q1 training. It is not a q1-to-q2 update. No q2, q3, or q4 was run.
Analytic `u*` and `z*` are first used after both final checkpoints are frozen;
they do not enter labels, optimization, selection, or stopping.

## Interpretation

This run isolates whether the known rank-k output subspace helps the current
one-improvement solver. It does not remove the benchmark's declared effective
dimension k from the PDE or from data generation. The result therefore tests
the Actor projection prior, rather than claiming a structure-free benchmark.

The practical decision is to avoid making rank-k projection a required part
of the method. A plain ambient-output Actor is viable and is substantially
better for `(400,1)`, while the projection remains a useful optional
regularizer for `(400,40)`. The latter cell prevents replacing the projected
variant unconditionally. For a structure-free main method, the next targeted
check should repeat only `(400,40)` and `(400,1)` over independent training
seeds and improve the unconstrained Actor's regularization without supplying
the analytic rank.

The `(400,40)` direct starting-gradient RMSE would rank the unconstrained Actor
ahead of the projected one, while its trajectory gap and value RMSE are worse.
This is another concrete reason to select and diagnose Actors on controlled
paths rather than by a single `t=0` gradient statistic.

There is one training seed per structure. The displayed intervals cover the
three final Monte Carlo evaluation banks conditional on each frozen Actor; they
do not measure training-seed variation. A publication-level statement about
the projection should be repeated over independent training seeds if the
observed difference is small or changes sign across cells.

## Cost and validation

Slurm job `{job['job_id']}` completed on `{job['node']}` in
{job['elapsed_seconds']} seconds using one NVIDIA GH200. Validation status is
**{validation['status']}**. It independently recomputed all metrics from
{validation['prediction_arrays_checked']} arrays containing
{validation['prediction_values_checked']:,} finite values, checked six final
Actor evaluations, 18 evaluation-seed records, {validation['cache_files_checked']}
new policy-label caches, {validation['shared_q0_cache_files_checked']} shared
raw q0 caches, q0 and final checkpoint architectures and hashes, matched
protocol fields, exact reproduction of the reused projected `(400,40)` Actor,
and scheduler accounting. The maximum stored-versus-recomputed metric
difference is {validation['maximum_recomputed_metric_difference']:.3e}.

From the repository root inside one GH200 allocation:

```bash
python -m pytest -q invariant_region_mlp/experiments/main_b1_sixmethods
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_iteration_projection_ablation --all
python -m invariant_region_mlp.experiments.main_b1_sixmethods.actor_iteration_projection_ablation_report
```

## Artifacts

- `result_d100_k10.json`, `result_d400_k40.json`, `result_d400_k1.json`: full protocol, training, evaluation, and cost records.
- `comparison.csv`: compact six-Actor metric table.
- `summary.json`: machine-readable comparisons.
- `validation.json`: independent prediction, protocol, cache, checkpoint, and accounting checks.
- `scheduler_accounting.json`: Slurm resource accounting.
- `manifest.json`: SHA-256 hashes for the committed report artifacts.
"""


def build_report() -> dict[str, Any]:
    results = {(d, k): _load(ROOT / f"result_d{d}_k{k}.json") for d, k in CELLS}
    validation = _validate(results)
    atomic_json(ROOT / "validation.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError(f"projection-ablation validation failed: {validation['failures']}")
    summary = _summary(results, validation)
    atomic_json(ROOT / "summary.json", summary)
    _write_comparison(summary)
    _plot(summary)
    REPORT_PATH.write_text(_report(summary, validation))
    artifacts = [
        *(ROOT / f"result_d{d}_k{k}.json" for d, k in CELLS),
        ACCOUNTING_PATH,
        ROOT / "comparison.csv",
        ROOT / "summary.json",
        ROOT / "validation.json",
        FIGURE_DIR / "actor_iteration_projection_ablation.png",
        FIGURE_DIR / "actor_iteration_projection_ablation.pdf",
        REPORT_PATH,
    ]
    manifest = {
        "schema_version": 1,
        "artifacts": [{
            "path": str(path.relative_to(PROJECT_DIR)),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        } for path in artifacts],
    }
    atomic_json(ROOT / "manifest.json", manifest)
    for item in manifest["artifacts"]:
        path = PROJECT_DIR / item["path"]
        if path.stat().st_size != item["bytes"] or _sha256(path) != item["sha256"]:
            raise RuntimeError(f"manifest verification failed for {path}")
    return summary


def main() -> None:
    print(json.dumps(build_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
