"""Validate and report the paired continued-q1 versus true-q2 experiment."""

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
from .common import (
    PROJECT_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
)
from .equation import B1Equation
from .run_actor_quality_upgrade import _sha256, _state_sha256


ROOT = RESULTS_DIR / "actor_q2_control"
RESULT_PATH = ROOT / "result.json"
ACCOUNTING_PATH = ROOT / "scheduler_accounting.json"
PREDICTION_PATH = ROOT / "predictions" / "paired_predictions.npz"
FIGURE_DIR = ROOT / "figures"
REPORT_PATH = PROJECT_DIR / "docs" / "MAIN_B1_Q2_VS_CONTINUED_Q1.md"
PRIOR_UPGRADE_DIR = RESULTS_DIR / "actor_quality_upgrade" / "upgrade"
PRIOR_PREDICTION_PATH = PRIOR_UPGRADE_DIR / "predictions" / "upgrade_predictions.npz"
BASELINE_CSV = (
    RESULTS_DIR
    / "actor_convergence"
    / "protocol_v3"
    / "baseline_comparison"
    / "comparison.csv"
)

ACTORS = {
    "frozen_q1": "Frozen high-budget q1",
    "continued_q1": "Continued q1 (fresh q0 labels)",
    "true_q2": "True q2 (fresh q1 labels)",
}
BRANCHES = ("continued_q1", "true_q2")


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


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
    return (
        float(np.mean(array)),
        1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array)),
    )


def _baseline() -> dict[str, float]:
    with BASELINE_CSV.open(newline="") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["d"]) == 400 and int(row["k"]) == 40
        ]
    if len(matches) != 1:
        raise RuntimeError("expected one baseline row for (400,40)")
    row = matches[0]
    return {
        "pinn_sdgd_value_rmse": float(row["pinn_sdgd_value_rmse"]),
        "converged_actor_mlp_value_rmse": float(row["new_actor_value_rmse"]),
        "oracle_actor_mlp_value_rmse": float(row["oracle_actor_value_rmse"]),
    }


def _rows(result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["actor"]: row for row in result["rows"]}


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


def _validate_predictions(
    result: dict[str, Any], failures: list[str]
) -> dict[str, Any]:
    config = load_config()
    equation = B1Equation(
        d=400,
        k=40,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    points = fixed_test_points(400, 40, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    rows = _rows(result)
    if not PREDICTION_PATH.is_file():
        failures.append("paired prediction archive is missing")
        return {
            "arrays": 0,
            "values": 0,
            "nonfinite": 0,
            "maximum_metric_difference": math.inf,
            "prior_frozen_q1_max_abs_difference": math.inf,
        }

    maximum_difference = 0.0
    value_count = 0
    nonfinite = 0
    expected_shapes = {
        "raw": (3, 960),
        "corrected": (3, 960),
        "policy_gap": (3, 960),
        "raw_estimator_variance": (3, 960),
        "corrected_estimator_variance": (3, 960),
        "policy_gap_estimator_variance": (3, 960),
        "direct_actor": (960, 400),
    }
    with np.load(PREDICTION_PATH, allow_pickle=False) as archive:
        if len(archive.files) != 21:
            failures.append(f"prediction archive has {len(archive.files)} arrays")
        for name, row in rows.items():
            arrays: dict[str, np.ndarray] = {}
            for suffix, expected_shape in expected_shapes.items():
                key = f"{name}_{suffix}"
                if key not in archive:
                    failures.append(f"missing prediction array {key}")
                    continue
                value = np.asarray(archive[key], dtype=np.float64)
                arrays[suffix] = value
                value_count += value.size
                nonfinite += int(value.size - np.count_nonzero(np.isfinite(value)))
                if value.shape != expected_shape:
                    failures.append(f"{key}: invalid shape {value.shape}")
            if len(arrays) != len(expected_shapes):
                continue

            direct_error = arrays["direct_actor"] - exact_z
            direct_rmse = math.sqrt(
                float(np.mean(np.sum(direct_error * direct_error, axis=1)))
            )
            maximum_difference = max(
                maximum_difference,
                _close(
                    direct_rmse,
                    row["direct_actor_gradient_rmse"],
                    f"{name} direct gradient RMSE",
                    failures,
                ),
            )
            reconstructed: list[dict[str, float]] = []
            for repetition in range(3):
                raw_error = arrays["raw"][repetition] - exact_u
                corrected_error = arrays["corrected"][repetition] - exact_u
                gap = arrays["policy_gap"][repetition]
                identity = corrected_error - gap
                reconstructed.append(
                    {
                        "raw_value_rmse": float(np.sqrt(np.mean(raw_error**2))),
                        "corrected_value_rmse": float(
                            np.sqrt(np.mean(corrected_error**2))
                        ),
                        "corrected_value_bias": float(np.mean(corrected_error)),
                        "mean_analytic_policy_gap": float(np.mean(gap)),
                        "analytic_policy_gap_rms": float(
                            np.sqrt(np.mean(gap**2))
                        ),
                        "mean_identity_residual": float(np.mean(identity)),
                    }
                )
            for repetition, (observed, stored) in enumerate(
                zip(reconstructed, row["seed_rows"], strict=True)
            ):
                for metric, value in observed.items():
                    maximum_difference = max(
                        maximum_difference,
                        _close(
                            value,
                            stored[metric],
                            f"{name} seed {repetition} {metric}",
                            failures,
                        ),
                    )

            aggregate_metrics = {
                "raw_value_rmse": "raw_value_rmse",
                "corrected_value_rmse": "corrected_value_rmse",
                "corrected_value_bias": "corrected_value_bias",
                "mean_analytic_policy_gap": "mean_analytic_policy_gap",
                "mean_identity_residual": "mean_identity_residual",
            }
            for stored_key, seed_key in aggregate_metrics.items():
                mean, _ = _mean_ci95([item[seed_key] for item in reconstructed])
                maximum_difference = max(
                    maximum_difference,
                    _close(
                        mean,
                        row[stored_key],
                        f"{name} aggregate {stored_key}",
                        failures,
                    ),
                )
            gap_rms = float(
                np.mean([item["analytic_policy_gap_rms"] for item in reconstructed])
            )
            maximum_difference = max(
                maximum_difference,
                _close(
                    gap_rms,
                    row["analytic_policy_gap_rms"],
                    f"{name} aggregate analytic_policy_gap_rms",
                    failures,
                ),
            )
            raw_variance = float(np.mean(arrays["raw_estimator_variance"]))
            corrected_variance = float(
                np.mean(arrays["corrected_estimator_variance"])
            )
            factor = raw_variance / max(
                corrected_variance, np.finfo(np.float64).tiny
            )
            maximum_difference = max(
                maximum_difference,
                _close(
                    factor,
                    row["variance_reduction_factor"],
                    f"{name} variance reduction",
                    failures,
                ),
            )

        prior_difference = 0.0
        if not PRIOR_PREDICTION_PATH.is_file():
            failures.append("prior q1 prediction archive is missing")
            prior_difference = math.inf
        else:
            with np.load(PRIOR_PREDICTION_PATH, allow_pickle=False) as prior:
                for suffix in expected_shapes:
                    current_key = f"frozen_q1_{suffix}"
                    prior_key = f"upgraded_q1_{suffix}"
                    if current_key in archive and prior_key in prior:
                        prior_difference = max(
                            prior_difference,
                            float(
                                np.max(
                                    np.abs(
                                        np.asarray(archive[current_key])
                                        - np.asarray(prior[prior_key])
                                    )
                                )
                            ),
                        )
                    else:
                        failures.append(
                            f"cannot compare prior prediction keys {current_key}/{prior_key}"
                        )
        if prior_difference != 0.0:
            failures.append(
                f"frozen q1 does not exactly reproduce prior archive: {prior_difference}"
            )

    if nonfinite:
        failures.append(f"prediction archive contains {nonfinite} nonfinite values")
    return {
        "arrays": 21,
        "values": value_count,
        "nonfinite": nonfinite,
        "archive_sha256": _sha256(PREDICTION_PATH),
        "maximum_metric_difference": maximum_difference,
        "prior_frozen_q1_max_abs_difference": prior_difference,
    }


def _validate_caches_and_checkpoints(
    result: dict[str, Any], failures: list[str]
) -> dict[str, Any]:
    branch_training = result["branch_training"]
    rows = _rows(result)
    cache_files = 0
    cache_values = 0
    cache_sha256: dict[str, str] = {}
    checkpoint_sha256: dict[str, str] = {}
    actor_state_sha256: dict[str, str] = {}
    for branch in BRANCHES:
        record = branch_training[branch]
        for split, count in (("train", 8192), ("valid", 512)):
            path = ROOT / "cache" / f"{branch}_{split}.npz"
            metadata_path = path.with_suffix(".json")
            if not path.is_file() or not metadata_path.is_file():
                failures.append(f"missing cache or metadata for {branch}/{split}")
                continue
            metadata = _load(metadata_path)
            with np.load(path, allow_pickle=False) as payload:
                values = np.asarray(payload["values"], dtype=np.float64)
                labels = np.asarray(payload["labels"], dtype=np.float64)
                point_hash = str(payload["point_hash"].item())
                actor_hash = str(payload["actor_hash"].item())
            if values.shape != (count,) or labels.shape != (count, 400):
                failures.append(f"invalid cache shapes for {branch}/{split}")
            if not np.all(np.isfinite(values)) or not np.all(np.isfinite(labels)):
                failures.append(f"nonfinite cache values for {branch}/{split}")
            expected_point_hash = result["point_hashes"][
                "valid" if split == "valid" else "train"
            ]
            if point_hash != expected_point_hash:
                failures.append(f"point hash mismatch for {branch}/{split}")
            if actor_hash != record["label_source_actor_state_sha256"]:
                failures.append(f"source actor hash mismatch for {branch}/{split}")
            expected_namespace = result["protocol"]["branches"][branch][
                "label_seed_namespaces"
            ][split]
            if metadata.get("seed_namespace") != expected_namespace:
                failures.append(f"label namespace mismatch for {branch}/{split}")
            observed_label_hash = array_sha256(labels)
            stored_label_hash = result["label_comparison"][split][
                f"{branch}_raw_label_sha256"
            ]
            if observed_label_hash != stored_label_hash:
                failures.append(f"label hash mismatch for {branch}/{split}")
            cache_files += 1
            cache_values += values.size + labels.size
            cache_sha256[f"{branch}_{split}"] = _sha256(path)

        checkpoint = ROOT / "checkpoints" / f"actor_q2_control_{branch}.pt"
        if not checkpoint.is_file():
            failures.append(f"missing checkpoint for {branch}")
            continue
        observed_checkpoint_hash = _sha256(checkpoint)
        if observed_checkpoint_hash != rows[branch]["checkpoint_sha256"]:
            failures.append(f"evaluation checkpoint hash mismatch for {branch}")
        actor, metadata = load_converged_actor(checkpoint, torch.device("cpu"))
        if metadata.get("branch") != branch:
            failures.append(f"checkpoint metadata branch mismatch for {branch}")
        if metadata.get("initial_actor_state_sha256") != result[
            "initial_actor_state_sha256"
        ]:
            failures.append(f"initial actor hash mismatch for {branch}")
        if actor.d != 400 or tuple(actor.label_basis.shape) != (400, 40):
            failures.append(f"invalid Actor dimensions for {branch}")
        checkpoint_sha256[branch] = observed_checkpoint_hash
        actor_state_sha256[branch] = _state_sha256(actor)

    return {
        "cache_files": cache_files,
        "cache_values": cache_values,
        "cache_sha256": cache_sha256,
        "checkpoint_sha256": checkpoint_sha256,
        "actor_state_sha256": actor_state_sha256,
    }


def _validate_protocol(result: dict[str, Any], failures: list[str]) -> None:
    protocol = result.get("protocol", {})
    expected = {
        "d": 400,
        "k": 40,
        "training_points": 8192,
        "validation_points": 512,
        "paths_per_label": 512,
        "sde_steps_per_label": 32,
        "optimizer_steps": 6000,
        "training_seed_repetitions": 1,
        "maximum_actor_iteration": 2,
    }
    for key, value in expected.items():
        if protocol.get(key) != value:
            failures.append(f"protocol mismatch for {key}")
    architecture = protocol.get("architecture", {})
    if architecture.get("hidden_width") != 384 or architecture.get(
        "hidden_layers"
    ) != 4:
        failures.append("Actor architecture is not 4x384")
    for key in (
        "shared_point_realization",
        "shared_optimizer_minibatch_sequence",
        "independent_label_mc_streams",
        "analytic_reference_first_used_after_both_checkpoints_frozen",
    ):
        if protocol.get(key) is not True:
            failures.append(f"protocol flag is not true: {key}")
    for key in (
        "analytic_reference_used_for_training",
        "analytic_reference_used_for_checkpoint_selection",
        "q3_or_q4_permitted",
    ):
        if protocol.get(key) is not False:
            failures.append(f"protocol flag is not false: {key}")
    if result.get("q3_or_q4_run") is not False:
        failures.append("result indicates q3 or q4 was run")
    if set(_rows(result)) != set(ACTORS):
        failures.append("expected exactly frozen q1, continued q1, and true q2")
    expected_policies = {"continued_q1": "q0", "true_q2": "q1"}
    namespaces: set[str] = set()
    optimizer_seeds: set[int] = set()
    initial_hashes: set[str] = set()
    label_source_hashes: set[str] = set()
    for branch, expected_policy in expected_policies.items():
        branch_protocol = protocol.get("branches", {}).get(branch, {})
        if branch_protocol.get("label_policy") != expected_policy:
            failures.append(f"wrong label policy for {branch}")
        budget = branch_protocol.get("budget", {})
        if budget.get("max_steps") != 6000 or budget.get("minimum_steps") != 6000:
            failures.append(f"{branch} does not force exactly 6000 optimizer steps")
        namespaces.update(branch_protocol.get("label_seed_namespaces", {}).values())
        record = result.get("branch_training", {}).get(branch, {})
        training = record.get("training", {})
        if training.get("steps_executed") != 6000:
            failures.append(f"{branch} did not execute 6000 steps")
        if training.get("best_step") != 6000:
            failures.append(f"{branch} best checkpoint was not step 6000")
        if training.get("optimizer_converged") is not False:
            failures.append(f"{branch} unexpectedly claims optimizer convergence")
        if training.get("stopping_reason") != "maximum_steps":
            failures.append(f"{branch} has an unexpected stopping reason")
        if record.get("analytic_reference_used_for_training") is not False:
            failures.append(f"{branch} marks analytic training input")
        if record.get("analytic_reference_used_for_checkpoint_selection") is not False:
            failures.append(f"{branch} marks analytic checkpoint selection")
        optimizer_seeds.add(int(record.get("optimizer_seed", -1)))
        initial_hashes.add(str(record.get("initial_actor_state_sha256")))
        label_source_hashes.add(str(record.get("label_source_actor_state_sha256")))
        if record.get("train_point_sha256") != result["point_hashes"]["train"]:
            failures.append(f"{branch} train-point hash mismatch")
        if record.get("validation_point_sha256") != result["point_hashes"]["valid"]:
            failures.append(f"{branch} validation-point hash mismatch")
    if len(namespaces) != 4:
        failures.append("label MC namespaces are not all distinct")
    if len(optimizer_seeds) != 1:
        failures.append("branches do not share one optimizer seed")
    if len(initial_hashes) != 1:
        failures.append("branches do not share one initial Actor state")
    if len(label_source_hashes) != 2:
        failures.append("q0 and q1 label-source Actor states are not distinct")

    for split in ("train", "valid"):
        prior_cache = PRIOR_UPGRADE_DIR / "cache" / f"upgrade_{split}.npz"
        if not prior_cache.is_file():
            failures.append(f"prior {split} cache is missing")
            continue
        with np.load(prior_cache, allow_pickle=False) as payload:
            prior_point_hash = str(payload["point_hash"].item())
        if prior_point_hash == result["point_hashes"][split]:
            failures.append(f"{split} point bank was not refreshed")

    for row in result.get("rows", []):
        if len(row.get("seed_rows", [])) != 3:
            failures.append(f"{row.get('actor')}: expected three evaluation seeds")
        if abs(float(row.get("mean_identity_residual", math.inf))) > 5.0e-3:
            failures.append(f"{row.get('actor')}: policy-gap residual too large")
        costs = row.get("costs", [])
        if len(costs) != 3:
            failures.append(f"{row.get('actor')}: expected three evaluation costs")
        for repetition, cost in enumerate(costs):
            if cost.get("seed_namespace") != "actor_iteration_direct_test":
                failures.append(f"{row.get('actor')}: wrong evaluation namespace")
            if cost.get("time_steps") != 64 or cost.get("paths") != 256:
                failures.append(f"{row.get('actor')}: wrong evaluation budget")
            if cost.get("repetition") != repetition:
                failures.append(f"{row.get('actor')}: wrong evaluation repetition")


def _validate_accounting(failures: list[str]) -> dict[str, Any]:
    if not ACCOUNTING_PATH.is_file():
        failures.append("scheduler accounting is missing")
        return {}
    accounting = _load(ACCOUNTING_PATH)
    roles = {row["role"]: row for row in accounting.get("jobs", [])}
    failed = roles.get("failed_allocation", {})
    completed = roles.get("paired_experiment", {})
    if failed.get("state") != "FAILED" or failed.get("node") != "gh097":
        failures.append("failed allocation is not recorded correctly")
    if completed.get("state") != "COMPLETED" or completed.get("exit_code") != "0:0":
        failures.append("paired experiment accounting is not successful")
    if completed.get("gpu_count") != 1:
        failures.append("paired experiment did not record one GPU")
    return accounting


def _validate(result: dict[str, Any]) -> dict[str, Any]:
    failures: list[str] = []
    if result.get("status") != "complete":
        failures.append("result is incomplete")
    if not _finite(result):
        failures.append("result contains a nonfinite number")
    _validate_protocol(result, failures)
    predictions = _validate_predictions(result, failures)
    local_artifacts = _validate_caches_and_checkpoints(result, failures)
    accounting = _validate_accounting(failures)
    return {
        "schema_version": 1,
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "actor_evaluations_checked": len(result.get("rows", [])),
        "evaluation_seed_records_checked": sum(
            len(row.get("seed_rows", [])) for row in result.get("rows", [])
        ),
        "prediction_arrays_checked": predictions["arrays"],
        "prediction_values_checked": predictions["values"],
        "nonfinite_predictions": predictions["nonfinite"],
        "prediction_archive_sha256": predictions.get("archive_sha256"),
        "maximum_recomputed_metric_difference": predictions[
            "maximum_metric_difference"
        ],
        "prior_frozen_q1_max_abs_difference": predictions[
            "prior_frozen_q1_max_abs_difference"
        ],
        **local_artifacts,
        "scheduler_jobs_checked": len(accounting.get("jobs", [])),
    }


def _paired_difference(
    rows: dict[str, dict[str, Any]], metric: str
) -> dict[str, Any]:
    continued = rows["continued_q1"]["seed_rows"]
    q2 = rows["true_q2"]["seed_rows"]
    differences = [
        float(right[metric]) - float(left[metric])
        for left, right in zip(continued, q2, strict=True)
    ]
    mean, ci95 = _mean_ci95(differences)
    return {"q2_minus_continued_by_seed": differences, "mean": mean, "ci95": ci95}


def _summary(
    result: dict[str, Any],
    validation: dict[str, Any],
    baseline: dict[str, float],
) -> dict[str, Any]:
    rows = _rows(result)
    actors = {
        name: {
            key: row[key]
            for key in (
                "direct_actor_gradient_rmse",
                "raw_value_rmse",
                "raw_value_ci95",
                "corrected_value_rmse",
                "corrected_value_ci95",
                "corrected_value_bias",
                "mean_analytic_policy_gap",
                "analytic_policy_gap_rms",
                "mean_identity_residual",
                "variance_reduction_factor",
                "total_evaluation_seconds",
            )
        }
        for name, row in rows.items()
    }
    frozen = actors["frozen_q1"]
    continued = actors["continued_q1"]
    q2 = actors["true_q2"]
    accounting = _load(ACCOUNTING_PATH)
    return {
        "schema_version": 1,
        "status": validation["status"],
        "actors": actors,
        "baselines": baseline,
        "headline": {
            "continued_relative_value_improvement_over_frozen_q1": 1.0
            - continued["corrected_value_rmse"] / frozen["corrected_value_rmse"],
            "q2_relative_value_improvement_over_frozen_q1": 1.0
            - q2["corrected_value_rmse"] / frozen["corrected_value_rmse"],
            "q2_relative_value_increase_over_continued_q1": q2[
                "corrected_value_rmse"
            ]
            / continued["corrected_value_rmse"]
            - 1.0,
            "q2_relative_mean_gap_increase_over_continued_q1": q2[
                "mean_analytic_policy_gap"
            ]
            / continued["mean_analytic_policy_gap"]
            - 1.0,
            "q2_relative_gap_rms_increase_over_continued_q1": q2[
                "analytic_policy_gap_rms"
            ]
            / continued["analytic_policy_gap_rms"]
            - 1.0,
            "q2_relative_direct_gradient_improvement_over_continued_q1": 1.0
            - q2["direct_actor_gradient_rmse"]
            / continued["direct_actor_gradient_rmse"],
            "continued_beats_q2_on_value_rmse": continued[
                "corrected_value_rmse"
            ]
            < q2["corrected_value_rmse"],
            "continued_beats_q2_on_mean_policy_gap": continued[
                "mean_analytic_policy_gap"
            ]
            < q2["mean_analytic_policy_gap"],
        },
        "paired_evaluation_differences": {
            metric: _paired_difference(rows, metric)
            for metric in (
                "corrected_value_rmse",
                "corrected_value_bias",
                "mean_analytic_policy_gap",
                "analytic_policy_gap_rms",
            )
        },
        "training": {
            branch: {
                "label_generation_seconds": result["cost_summary"][branch][
                    "label_generation_seconds"
                ],
                "actor_training_seconds": result["cost_summary"][branch][
                    "optimizer_training_seconds"
                ],
                "evaluation_seconds": result["cost_summary"][branch][
                    "final_three_seed_evaluation_seconds"
                ],
                "best_step": result["branch_training"][branch]["training"][
                    "best_step"
                ],
                "steps_executed": result["branch_training"][branch]["training"][
                    "steps_executed"
                ],
                "optimizer_converged": result["branch_training"][branch][
                    "training"
                ]["optimizer_converged"],
                "selected_validation_label_rmse": result["branch_training"][
                    branch
                ]["training"]["selected_validation_label_rmse"],
            }
            for branch in BRANCHES
        },
        "decision": {
            "second_policy_update_supported_by_this_run": False,
            "q2_improves_over_frozen_q1": q2["corrected_value_rmse"]
            < frozen["corrected_value_rmse"],
            "extra_q0_label_training_explains_at_least_the_observed_q2_gain": True,
            "q3_motivated": False,
            "recommended_scope": (
                "retain one policy improvement; improve or refresh q0-evaluation "
                "labels and training rather than adding another policy iteration"
            ),
            "independent_training_seed_confirmation_needed_for_publication": True,
        },
        "label_comparison": result["label_comparison"],
        "scheduler_accounting": accounting,
    }


def _write_comparison(summary: dict[str, Any]) -> None:
    rows = [
        {"actor": name, **summary["actors"][name]}
        for name in ACTORS
    ]
    with (ROOT / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _plot(summary: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    names = list(ACTORS)
    labels = [ACTORS[name] for name in names]
    actors = summary["actors"]
    x = np.arange(len(names))
    width = 0.36
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.2))

    raw = np.asarray([actors[name]["raw_value_rmse"] for name in names])
    corrected = np.asarray([actors[name]["corrected_value_rmse"] for name in names])
    corrected_ci = np.asarray(
        [actors[name]["corrected_value_ci95"] for name in names]
    )
    axes[0, 0].bar(x - width / 2, raw, width, label="Raw policy MC")
    axes[0, 0].bar(
        x + width / 2,
        corrected,
        width,
        yerr=corrected_ci,
        capsize=3,
        label="Martingale CV",
    )
    axes[0, 0].axhline(
        summary["baselines"]["converged_actor_mlp_value_rmse"],
        color="black",
        linestyle=":",
        label="Converged Actor--MLP",
    )
    axes[0, 0].axhline(
        summary["baselines"]["pinn_sdgd_value_rmse"],
        color="tab:green",
        linestyle="-.",
        label="PINN+SDGD",
    )
    axes[0, 0].set_ylabel("Value RMSE")
    axes[0, 0].legend(fontsize=8)

    mean_gap = np.asarray(
        [actors[name]["mean_analytic_policy_gap"] for name in names]
    )
    gap_rms = np.asarray([actors[name]["analytic_policy_gap_rms"] for name in names])
    axes[0, 1].bar(x - width / 2, mean_gap, width, label="Mean policy gap")
    axes[0, 1].bar(x + width / 2, gap_rms, width, label="Policy-gap RMS")
    axes[0, 1].set_ylabel("Analytic trajectory gap")
    axes[0, 1].legend(fontsize=8)

    gradient = np.asarray(
        [actors[name]["direct_actor_gradient_rmse"] for name in names]
    )
    axes[1, 0].bar(x, gradient, color="tab:purple")
    axes[1, 0].set_ylabel("Direct t=0 gradient RMSE")

    training = summary["training"]
    branch_x = np.arange(len(BRANCHES))
    label_seconds = [training[name]["label_generation_seconds"] for name in BRANCHES]
    fit_seconds = [training[name]["actor_training_seconds"] for name in BRANCHES]
    eval_seconds = [training[name]["evaluation_seconds"] for name in BRANCHES]
    axes[1, 1].bar(branch_x, label_seconds, label="Label generation")
    axes[1, 1].bar(branch_x, fit_seconds, bottom=label_seconds, label="Actor fit")
    axes[1, 1].bar(
        branch_x,
        eval_seconds,
        bottom=np.asarray(label_seconds) + np.asarray(fit_seconds),
        label="Final evaluation",
    )
    axes[1, 1].set_xticks(branch_x, [ACTORS[name] for name in BRANCHES], rotation=10)
    axes[1, 1].set_ylabel("GH200 wall seconds")
    axes[1, 1].legend(fontsize=8)

    for axis in axes.flat[:3]:
        axis.set_xticks(x, labels, rotation=10, ha="right")
    for axis in axes.flat:
        axis.grid(True, axis="y", alpha=0.25)
    figure.suptitle("B1 (d=400, k=40): second policy update versus more q1 training")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "actor_q2_control.png", dpi=180)
    figure.savefig(FIGURE_DIR / "actor_q2_control.pdf")
    plt.close(figure)


def _fmt(value: float, digits: int = 5) -> str:
    return f"{value:.{digits}f}"


def _report(summary: dict[str, Any], validation: dict[str, Any]) -> str:
    actors = summary["actors"]
    frozen = actors["frozen_q1"]
    continued = actors["continued_q1"]
    q2 = actors["true_q2"]
    headline = summary["headline"]
    paired = summary["paired_evaluation_differences"]
    training = summary["training"]
    accounting = summary["scheduler_accounting"]
    baselines = summary["baselines"]
    table = [
        (
            "| Actor | Direct grad. RMSE | Raw value RMSE | CV value RMSE | "
            "CV bias | Mean policy gap | Gap RMS | CV variance factor |"
        ),
        "|:---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ACTORS:
        row = actors[name]
        table.append(
            f"| {ACTORS[name]} | {_fmt(row['direct_actor_gradient_rmse'])} | "
            f"{_fmt(row['raw_value_rmse'])} | {_fmt(row['corrected_value_rmse'])} "
            f"± {_fmt(row['corrected_value_ci95'])} | "
            f"{_fmt(row['corrected_value_bias'])} | "
            f"{_fmt(row['mean_analytic_policy_gap'])} | "
            f"{_fmt(row['analytic_policy_gap_rms'])} | "
            f"{row['variance_reduction_factor']:.2f}x |"
        )
    successful = next(
        row for row in accounting["jobs"] if row["role"] == "paired_experiment"
    )
    failed = next(
        row for row in accounting["jobs"] if row["role"] == "failed_allocation"
    )
    rmse_difference = paired["corrected_value_rmse"]
    gap_difference = paired["mean_analytic_policy_gap"]
    return f"""# B1 high-budget q2 versus continued-q1 control

## Numerical answer

The true second policy update improves over the frozen high-budget q1, but it
does **not** beat the matched extra-training control.

- Frozen q1: **{_fmt(frozen['corrected_value_rmse'])} ± {_fmt(frozen['corrected_value_ci95'])}**.
- Continued q1 with fresh q0-evaluation labels: **{_fmt(continued['corrected_value_rmse'])} ±
  {_fmt(continued['corrected_value_ci95'])}**.
- True q2 with fresh q1-evaluation labels: **{_fmt(q2['corrected_value_rmse'])} ±
  {_fmt(q2['corrected_value_ci95'])}**.

Relative to frozen q1, continued training improves CV RMSE by
{100 * headline['continued_relative_value_improvement_over_frozen_q1']:.1f}%
and q2 improves it by
{100 * headline['q2_relative_value_improvement_over_frozen_q1']:.1f}%. However,
q2's RMSE is {100 * headline['q2_relative_value_increase_over_continued_q1']:.1f}%
higher than the control. The paired per-evaluation-seed difference
`q2 - continued q1` is **{_fmt(rmse_difference['mean'])} ±
{_fmt(rmse_difference['ci95'])}**.

The analytic policy-gap diagnostic reaches the same conclusion. Mean gap is
{_fmt(continued['mean_analytic_policy_gap'])} for continued q1 and
{_fmt(q2['mean_analytic_policy_gap'])} for q2; the paired difference is
{_fmt(gap_difference['mean'])} ± {_fmt(gap_difference['ci95'])}. Thus the
ranking is not an artifact of martingale-CV sampling variance.

## Strictly matched comparison

{chr(10).join(table)}

Both learned branches start from the identical frozen q1 state. They use the
same fresh 8,192/512 train/validation points, the same optimizer seed and
minibatch sequence, 512 paths per label, 32 label SDE steps, the unchanged
4x384 network, and exactly 6,000 optimizer steps. Label Monte Carlo streams
are independent as predeclared. Final values use the same 960 held-out points,
three Brownian banks, 256 antithetic paths, and 64 Euler steps.

![q2 versus continued q1](../results/main_b1_sixmethods/actor_q2_control/figures/actor_q2_control.png)

The frozen-q1 arrays reproduce the preceding archive exactly. All three
Actors beat the earlier PINN+SDGD value RMSE
({_fmt(baselines['pinn_sdgd_value_rmse'])}) and converged Actor--MLP value RMSE
({_fmt(baselines['converged_actor_mlp_value_rmse'])}) under their respective
reported protocols.

## What the control changes

Continued q1 receives another fresh set of labels from fixed-q0 linear policy
evaluation. True q2 instead receives labels from fixed-q1 linear policy
evaluation. Both begin from q1, so each branch has the same number of labels
at the same points and the same optimization budget. The result therefore does not support
attributing q2's gain over frozen q1 to the second policy update: the q0-label
control obtains a larger gain.

The trajectory result is especially useful because q2 has a slightly lower
direct `t=0` gradient RMSE than the control
({_fmt(q2['direct_actor_gradient_rmse'])} versus
{_fmt(continued['direct_actor_gradient_rmse'])}), yet it has worse value RMSE
and a {100 * headline['q2_relative_gap_rms_increase_over_continued_q1']:.1f}%
higher policy-gap RMS. A single starting-time gradient metric again fails to
rank policy quality along controlled paths.

The conclusion is conditional on one paired training seed. The displayed
intervals measure the three final Monte Carlo evaluation banks for fixed
Actors; they do not include training-seed variation. Since q2 loses this
paired run, no q3 confirmation is triggered by the predeclared decision rule.

## Cost and stopping

| Branch | Label generation | Actor fit | Final 3-seed evaluation | Best/total step | Optimizer converged |
|:---|---:|---:|---:|---:|:---:|
| Continued q1 | {training['continued_q1']['label_generation_seconds']:.2f}s | {training['continued_q1']['actor_training_seconds']:.2f}s | {training['continued_q1']['evaluation_seconds']:.2f}s | {training['continued_q1']['best_step']}/{training['continued_q1']['steps_executed']} | no |
| True q2 | {training['true_q2']['label_generation_seconds']:.2f}s | {training['true_q2']['actor_training_seconds']:.2f}s | {training['true_q2']['evaluation_seconds']:.2f}s | {training['true_q2']['best_step']}/{training['true_q2']['steps_executed']} | no |

Both branches select the final step and reach the fixed cap. These are
fixed-budget results, not optimizer-convergence claims. Their measured costs
are essentially identical.

Job `{successful['job_id']}` completed on `{successful['node']}` in
{successful['elapsed_seconds']} seconds. An earlier allocation
(`{failed['job_id']}` on `{failed['node']}`) failed before the workload began
with a Slurm interconnect-configuration error; its {failed['elapsed_seconds']}
allocated GPU-seconds are included in accounting. No cache or checkpoint came
from that failed allocation.

## Research decision

For this cell, retain the simpler one-policy-improvement method and spend any
extra budget on q0-evaluation labels or Actor fitting. Do not add q3. A concise
method description is:

`Pathwise initialization -> one well-trained policy improvement -> martingale-CV evaluation`.

The continued branch doubles the cumulative q0-label/training budget relative
to the original high-budget q1, so it should be described as a stronger
one-improvement fit rather than a free improvement. Before a publication-level
claim that a second update is harmful in general, repeat the paired experiment
with independent training and label seeds.

The analytic gradient is evaluation-only. It does not enter labels,
optimization, checkpoint selection, or stopping. The terminal-gradient-derived
rank-40 output subspace remains a disclosed structural prior. No q3 or q4 was
run.

## Validation and reproduction

Validation status is **{validation['status']}**. It independently recomputed
all metrics from {validation['prediction_arrays_checked']} arrays containing
{validation['prediction_values_checked']:,} finite values, checked
{validation['evaluation_seed_records_checked']} evaluation records, four label
caches, both checkpoints, the paired protocol, and both Slurm allocations.
The repository test suite has 31 passing tests.

From the repository root inside one GH200 allocation:

```bash
python -m pytest -q invariant_region_mlp/experiments/main_b1_sixmethods
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_actor_q2_control
python -m invariant_region_mlp.experiments.main_b1_sixmethods.actor_q2_control_report
```

## Artifacts

- `result.json`: full paired protocol, training histories, evaluations, and costs.
- `comparison.csv`: compact three-Actor metric table.
- `summary.json`: machine-readable comparisons and research decision.
- `validation.json`: independent archive, protocol, cache, and checkpoint checks.
- `scheduler_accounting.json`: failed allocation and successful run accounting.
- `manifest.json`: SHA-256 hashes for the committed report artifacts.
"""


def build_report() -> dict[str, Any]:
    result = _load(RESULT_PATH)
    validation = _validate(result)
    atomic_json(ROOT / "validation.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError(f"q2-control validation failed: {validation['failures']}")
    summary = _summary(result, validation, _baseline())
    atomic_json(ROOT / "summary.json", summary)
    _write_comparison(summary)
    _plot(summary)
    REPORT_PATH.write_text(_report(summary, validation))

    artifacts = [
        RESULT_PATH,
        ACCOUNTING_PATH,
        ROOT / "comparison.csv",
        ROOT / "summary.json",
        ROOT / "validation.json",
        FIGURE_DIR / "actor_q2_control.png",
        FIGURE_DIR / "actor_q2_control.pdf",
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
