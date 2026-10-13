"""Run the minimal teacher-free Actor Iteration pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time
import traceback
from typing import Any

import numpy as np
import torch

from .actor_convergence import (
    ACTOR_DEPTH,
    ACTOR_WIDTH,
    BudgetLevel,
    checkpoint_payload,
    learn_label_subspace,
    load_converged_actor,
    make_actor,
    predict_actor,
    project_labels,
    train_budget_level,
)
from .actor_iteration import evaluate_policy_mc, terminal_pathwise_labels
from .actor_mlp import evaluate_actor_mlp
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
    seed_everything,
    stable_seed,
    training_points,
)
from .equation import B1Equation
from .metrics import prediction_metrics
from .neural_common import count_parameters


PILOT_DIR = RESULTS_DIR / "actor_iteration_pilot"
CACHE_DIR = PILOT_DIR / "cache"
CHECKPOINT_DIR = PILOT_DIR / "checkpoints"
PREDICTION_DIR = PILOT_DIR / "predictions"

Q0_TRAIN_POINTS = 4096
Q0_VALID_POINTS = 1024
Q0_PATHS = 256
IMPROVEMENT_TRAIN_POINTS = 2048
IMPROVEMENT_VALID_POINTS = 512
IMPROVEMENT_PATHS = 128
POLICY_TIME_STEPS = 16
POLICY_PARENT_BATCH = 16
POLICY_PATH_CHUNK = 32
VALUE_VALID_POINTS = 256
VALUE_VALID_PATHS = 256
VALUE_VALID_SEEDS = 3
DIRECT_TEST_PATHS = 256
DIRECT_TEST_SEEDS = 3
MAX_IMPROVEMENTS = 3
VALUE_ABSOLUTE_TOLERANCE = 1.0e-3
PILOT_CELLS = {(100, 10), (400, 40), (400, 1)}

Q0_BUDGET = BudgetLevel(
    "actor_iteration_q0",
    1,
    Q0_PATHS,
    1,
    Q0_TRAIN_POINTS,
    Q0_VALID_POINTS,
    6000,
    1500,
    5.0e-4,
)
IMPROVEMENT_BUDGET = BudgetLevel(
    "actor_iteration_improvement",
    POLICY_TIME_STEPS,
    IMPROVEMENT_PATHS,
    1,
    IMPROVEMENT_TRAIN_POINTS,
    IMPROVEMENT_VALID_POINTS,
    6000,
    1500,
    1.0e-4,
)


def protocol(d: int, k: int) -> dict[str, Any]:
    return {
        "cell": {"d": d, "k": k},
        "architecture": {
            "hidden_width": ACTOR_WIDTH,
            "hidden_layers": ACTOR_DEPTH,
            "activation": "silu",
            "output_projection": "terminal-label-derived rank k",
        },
        "q0": {
            "training_points": Q0_TRAIN_POINTS,
            "validation_points": Q0_VALID_POINTS,
            "terminal_paths": Q0_PATHS,
            "budget": Q0_BUDGET.__dict__,
        },
        "improvement": {
            "training_points": IMPROVEMENT_TRAIN_POINTS,
            "validation_points": IMPROVEMENT_VALID_POINTS,
            "time_steps": POLICY_TIME_STEPS,
            "paths": IMPROVEMENT_PATHS,
            "maximum_improvements": MAX_IMPROVEMENTS,
            "budget": IMPROVEMENT_BUDGET.__dict__,
        },
        "value_stopping_audit": {
            "points": VALUE_VALID_POINTS,
            "paths": VALUE_VALID_PATHS,
            "seeds": VALUE_VALID_SEEDS,
            "time_steps": POLICY_TIME_STEPS,
            "absolute_tolerance": VALUE_ABSOLUTE_TOLERANCE,
            "consecutive_stable_transitions": 2,
        },
        "final_direct_policy": {
            "test_points": 960,
            "paths": DIRECT_TEST_PATHS,
            "seeds": DIRECT_TEST_SEEDS,
            "time_steps": POLICY_TIME_STEPS,
        },
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_stopping": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }


def _atomic_torch_save(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    torch.save(payload, temporary)
    temporary.replace(path)


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _atomic_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty CSV")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _state_sha256(model: torch.nn.Module) -> str:
    digest = hashlib.sha256()
    for key, value in sorted(model.state_dict().items()):
        array = np.ascontiguousarray(value.detach().cpu().numpy())
        digest.update(key.encode())
        digest.update(str(array.shape).encode())
        digest.update(str(array.dtype).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    mean = float(np.mean(array))
    ci95 = (
        1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array))
        if len(array) > 1
        else 0.0
    )
    return mean, ci95


def _terminal_label_cache(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    split: str,
    device: torch.device,
) -> tuple[np.ndarray, dict[str, Any]]:
    cache = CACHE_DIR / f"d{equation.d}_k{equation.k}" / f"q0_{split}.npz"
    metadata_path = cache.with_suffix(".json")
    point_hash = array_sha256(t, x)
    if cache.is_file() and metadata_path.is_file():
        with np.load(cache, allow_pickle=False) as payload:
            observed_hash = str(payload["point_hash"].item())
            labels = np.asarray(payload["labels"], dtype=np.float64)
        if observed_hash != point_hash:
            raise RuntimeError(f"q0 label point hash mismatch: {cache}")
        return labels, json.loads(metadata_path.read_text())
    labels, cost = terminal_pathwise_labels(
        equation,
        t,
        x,
        paths=Q0_PATHS,
        repetition=0,
        seed_namespace=f"actor_iteration_q0_{split}",
        device=device,
        parent_batch_size=64,
    )
    metadata = {
        "status": "complete",
        "split": split,
        "point_sha256": point_hash,
        "labels_projected": False,
        "cost": cost,
    }
    _atomic_npz(cache, labels=labels, point_hash=np.asarray(point_hash))
    atomic_json(metadata_path, metadata)
    return labels, metadata


def _policy_label_cache(
    equation: B1Equation,
    actor: torch.nn.Module,
    iteration: int,
    t: np.ndarray,
    x: np.ndarray,
    *,
    split: str,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    cache = (
        CACHE_DIR
        / f"d{equation.d}_k{equation.k}"
        / f"iteration_{iteration}_{split}.npz"
    )
    metadata_path = cache.with_suffix(".json")
    point_hash = array_sha256(t, x)
    actor_hash = _state_sha256(actor)
    if cache.is_file() and metadata_path.is_file():
        with np.load(cache, allow_pickle=False) as payload:
            observed_point_hash = str(payload["point_hash"].item())
            observed_actor_hash = str(payload["actor_hash"].item())
            values = np.asarray(payload["values"], dtype=np.float64)
            labels = np.asarray(payload["labels"], dtype=np.float64)
        if observed_point_hash != point_hash or observed_actor_hash != actor_hash:
            raise RuntimeError(f"policy-label cache provenance mismatch: {cache}")
        return values, labels, json.loads(metadata_path.read_text())
    values, labels, cost = evaluate_policy_mc(
        equation,
        actor,
        t,
        x,
        time_steps=POLICY_TIME_STEPS,
        paths=IMPROVEMENT_PATHS,
        repetition=0,
        seed_namespace=f"actor_iteration_{iteration}_{split}_labels",
        device=device,
        parent_batch_size=POLICY_PARENT_BATCH,
        path_chunk_size=POLICY_PATH_CHUNK,
        compute_z=True,
    )
    if labels is None:
        raise RuntimeError("policy improvement did not return pathwise labels")
    metadata = {
        "status": "complete",
        "iteration": iteration,
        "source_actor_iteration": iteration - 1,
        "split": split,
        "point_sha256": point_hash,
        "source_actor_sha256": actor_hash,
        "labels_projected": False,
        "cost": cost,
        "nonterminal_exact_reference_used": False,
    }
    _atomic_npz(
        cache,
        values=values,
        labels=labels,
        point_hash=np.asarray(point_hash),
        actor_hash=np.asarray(actor_hash),
    )
    atomic_json(metadata_path, metadata)
    return values, labels, metadata


def _policy_validation(
    equation: B1Equation,
    actor: torch.nn.Module,
    iteration: int,
    t: np.ndarray,
    x: np.ndarray,
    *,
    device: torch.device,
) -> dict[str, Any]:
    cache = (
        CACHE_DIR
        / f"d{equation.d}_k{equation.k}"
        / f"value_validation_iteration_{iteration}.npz"
    )
    metadata_path = cache.with_suffix(".json")
    point_hash = array_sha256(t, x)
    actor_hash = _state_sha256(actor)
    if cache.is_file() and metadata_path.is_file():
        with np.load(cache, allow_pickle=False) as payload:
            observed_point_hash = str(payload["point_hash"].item())
            observed_actor_hash = str(payload["actor_hash"].item())
            values = np.asarray(payload["values"], dtype=np.float64)
        if observed_point_hash != point_hash or observed_actor_hash != actor_hash:
            raise RuntimeError(f"value-validation cache mismatch: {cache}")
        metadata = json.loads(metadata_path.read_text())
        return {**metadata, "values": values}
    seed_values: list[np.ndarray] = []
    costs: list[dict[str, Any]] = []
    for repetition in range(VALUE_VALID_SEEDS):
        value, _, cost = evaluate_policy_mc(
            equation,
            actor,
            t,
            x,
            time_steps=POLICY_TIME_STEPS,
            paths=VALUE_VALID_PATHS,
            repetition=repetition,
            seed_namespace="actor_iteration_value_validation",
            device=device,
            parent_batch_size=32,
            path_chunk_size=64,
            compute_z=False,
        )
        seed_values.append(value)
        costs.append(cost)
    values = np.stack(seed_values)
    metadata = {
        "status": "complete",
        "iteration": iteration,
        "point_sha256": point_hash,
        "actor_sha256": actor_hash,
        "seed_mean_values": np.mean(values, axis=1).tolist(),
        "mean_value": float(np.mean(values)),
        "costs": costs,
        "analytic_reference_used": False,
    }
    _atomic_npz(
        cache,
        values=values,
        point_hash=np.asarray(point_hash),
        actor_hash=np.asarray(actor_hash),
    )
    atomic_json(metadata_path, metadata)
    return {**metadata, "values": values}


def _value_improvement(
    previous: dict[str, Any], current: dict[str, Any]
) -> dict[str, Any]:
    previous_values = np.asarray(previous["values"], dtype=np.float64)
    current_values = np.asarray(current["values"], dtype=np.float64)
    differences = np.mean(previous_values - current_values, axis=1)
    mean, ci95 = _mean_ci95(differences.tolist())
    stable = abs(mean) <= ci95 + VALUE_ABSOLUTE_TOLERANCE
    return {
        "seed_mean_previous_minus_current": differences.tolist(),
        "mean_previous_minus_current": mean,
        "ci95": ci95,
        "confident_decrease": mean - ci95 > 0.0,
        "confident_increase": mean + ci95 < 0.0,
        "small_relative_to_uncertainty": stable,
    }


def _projected_label_metadata(
    raw: np.ndarray, projected: np.ndarray
) -> dict[str, float]:
    residual = math.sqrt(float(np.mean(np.sum((raw - projected) ** 2, axis=1))))
    scale = math.sqrt(float(np.mean(np.sum(raw**2, axis=1))))
    return {
        "raw_vector_rmse": scale,
        "projection_residual_rmse": residual,
        "relative_projection_residual": residual
        / max(scale, np.finfo(np.float64).tiny),
    }


def _checkpoint_path(d: int, k: int, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"actor_iteration_d{d}_k{k}_q{iteration}.pt"


def _baseline_references(d: int, k: int) -> dict[str, Any]:
    path = (
        RESULTS_DIR
        / "actor_convergence"
        / "protocol_v3"
        / "baseline_comparison"
        / "comparison.csv"
    )
    with path.open(newline="") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["d"]) == d and int(row["k"]) == k
        ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one baseline row for d={d}, k={k}")
    row = matches[0]
    return {
        "vanilla_mlp_value_rmse": float(row["vanilla_mlp_value_rmse"]),
        "double_mlp_value_rmse": float(row["double_mlp_value_rmse"]),
        "dpi_value_rmse": float(row["dpi_value_rmse"]),
        "pinn_sdgd_value_rmse": float(row["pinn_sdgd_value_rmse"]),
        "original_actor_mlp_value_rmse": float(row["actor_mlp_value_rmse"]),
        "converged_actor_mlp_value_rmse": float(row["new_actor_value_rmse"]),
        "oracle_actor_mlp_value_rmse": float(row["oracle_actor_value_rmse"]),
        "best_previous_value_method": row["best_previous_value_method"],
        "best_previous_value_rmse": float(row["best_previous_value_rmse"]),
    }


def _final_evaluation(
    equation: B1Equation,
    actors: list[torch.nn.Module],
    construction_records: list[dict[str, Any]],
    validation_records: list[dict[str, Any]],
    *,
    device: torch.device,
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray], str]:
    config = load_config()["stage_b"]
    points = fixed_test_points(equation.d, equation.k, int(config["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    test_hash = array_sha256(points)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    rows: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {}
    for iteration, (actor, construction, validation) in enumerate(
        zip(actors, construction_records, validation_records, strict=True)
    ):
        direct_seed_rows: list[dict[str, Any]] = []
        direct_predictions: list[np.ndarray] = []
        direct_costs: list[dict[str, Any]] = []
        for repetition in range(DIRECT_TEST_SEEDS):
            value, _, cost = evaluate_policy_mc(
                equation,
                actor,
                times,
                points,
                time_steps=POLICY_TIME_STEPS,
                paths=DIRECT_TEST_PATHS,
                repetition=repetition,
                seed_namespace="actor_iteration_direct_test",
                device=device,
                parent_batch_size=32,
                path_chunk_size=64,
                compute_z=False,
            )
            error = value - exact_u
            direct_seed_rows.append(
                {
                    "repetition": repetition,
                    "value_rmse": float(np.sqrt(np.mean(error**2))),
                    "value_bias": float(np.mean(error)),
                }
            )
            direct_predictions.append(value)
            direct_costs.append(cost)
        direct_mean, direct_ci95 = _mean_ci95(
            [row["value_rmse"] for row in direct_seed_rows]
        )

        actor_mlp_rows: list[dict[str, Any]] = []
        actor_mlp_predictions: list[np.ndarray] = []
        actor_mlp_costs: list[dict[str, Any]] = []
        for repetition in range(int(config["mlp_seeds"])):
            state, cost = evaluate_actor_mlp(
                equation,
                actor,
                times,
                points,
                n=int(config["mlp_n"]),
                M=int(config["mlp_M"]),
                repetition=repetition,
                device=device,
            )
            metrics = prediction_metrics(
                state[:, 0], state[:, 1:], exact_u, exact_z
            )
            actor_mlp_rows.append(
                {"repetition": repetition, **metrics}
            )
            actor_mlp_predictions.append(state[:, 0])
            actor_mlp_costs.append(cost)
        actor_mlp_mean, actor_mlp_ci95 = _mean_ci95(
            [row["value_rmse"] for row in actor_mlp_rows]
        )

        actor_prediction = predict_actor(actor, times, points, device=device)
        actor_error = actor_prediction - exact_z
        actor_gradient_rmse = math.sqrt(
            float(np.mean(np.sum(actor_error**2, axis=1)))
        )
        exact_gradient_scale = math.sqrt(
            float(np.mean(np.sum(exact_z**2, axis=1)))
        )
        predictions[f"q{iteration}_direct_policy_value"] = np.stack(
            direct_predictions
        )
        predictions[f"q{iteration}_actor_mlp_value"] = np.stack(
            actor_mlp_predictions
        )
        predictions[f"q{iteration}_direct_actor"] = actor_prediction
        rows.append(
            {
                "iteration": iteration,
                "direct_policy_value_rmse": direct_mean,
                "direct_policy_value_ci95": direct_ci95,
                "actor_mlp_value_rmse": actor_mlp_mean,
                "actor_mlp_value_ci95": actor_mlp_ci95,
                "actor_gradient_rmse": actor_gradient_rmse,
                "actor_relative_gradient_error": actor_gradient_rmse
                / max(exact_gradient_scale, np.finfo(np.float64).tiny),
                "validation_value_mean": validation["mean_value"],
                "validation_value_improvement": (
                    None
                    if iteration == 0
                    else construction["value_improvement"][
                        "mean_previous_minus_current"
                    ]
                ),
                "validation_value_improvement_ci95": (
                    None
                    if iteration == 0
                    else construction["value_improvement"]["ci95"]
                ),
                "label_generation_seconds": construction[
                    "label_generation_seconds"
                ],
                "actor_training_seconds": construction["training"][
                    "training_seconds"
                ],
                "value_validation_seconds": sum(
                    cost["wall_seconds"] for cost in validation["costs"]
                ),
                "direct_policy_inference_seconds": sum(
                    cost["wall_seconds"] for cost in direct_costs
                ),
                "actor_mlp_inference_seconds": sum(
                    cost["wall_seconds"] for cost in actor_mlp_costs
                ),
                "direct_policy_seed_rows": direct_seed_rows,
                "actor_mlp_seed_rows": actor_mlp_rows,
                "direct_policy_costs": direct_costs,
                "actor_mlp_costs": actor_mlp_costs,
            }
        )
    return rows, predictions, test_hash


def run(d: int, k: int, *, force: bool = False) -> dict[str, Any]:
    result_path = PILOT_DIR / f"result_d{d}_k{k}.json"
    if result_path.is_file() and not force:
        existing = json.loads(result_path.read_text())
        if existing.get("status") == "complete":
            print(json.dumps(existing["headline"], indent=2), flush=True)
            return existing
    if (d, k) not in PILOT_CELLS:
        allowed = ", ".join(
            f"({cell_d},{cell_k})" for cell_d, cell_k in sorted(PILOT_CELLS)
        )
        raise ValueError(f"Actor Iteration pilot supports only {allowed}")

    config = load_config()
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Actor Iteration pilot requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    started = time.perf_counter()

    cell_cache = CACHE_DIR / f"d{d}_k{k}"
    label_basis, subspace_metadata = learn_label_subspace(
        equation, device=device, cache_dir=cell_cache / "label_subspace"
    )
    point_seed = stable_seed("B1", "actor_iteration_points", d, k)
    q0_train_t, q0_train_x = training_points(
        d, k, "actor_iteration_q0_train", point_seed, Q0_TRAIN_POINTS
    )
    q0_valid_t, q0_valid_x = training_points(
        d, k, "actor_iteration_q0_valid", point_seed, Q0_VALID_POINTS
    )
    improve_train_t, improve_train_x = training_points(
        d,
        k,
        "actor_iteration_improvement_train",
        point_seed,
        IMPROVEMENT_TRAIN_POINTS,
    )
    improve_valid_t, improve_valid_x = training_points(
        d,
        k,
        "actor_iteration_improvement_valid",
        point_seed,
        IMPROVEMENT_VALID_POINTS,
    )
    value_valid_t, value_valid_x = training_points(
        d,
        k,
        "actor_iteration_value_validation",
        point_seed,
        VALUE_VALID_POINTS,
    )

    actors: list[torch.nn.Module] = []
    construction_records: list[dict[str, Any]] = []
    validation_records: list[dict[str, Any]] = []

    q0_path = _checkpoint_path(d, k, 0)
    if q0_path.is_file() and not force:
        model, q0_record = load_converged_actor(q0_path, device)
    else:
        q0_train_raw, q0_train_metadata = _terminal_label_cache(
            equation, q0_train_t, q0_train_x, split="train", device=device
        )
        q0_valid_raw, q0_valid_metadata = _terminal_label_cache(
            equation, q0_valid_t, q0_valid_x, split="valid", device=device
        )
        q0_train_y = project_labels(q0_train_raw, label_basis)
        q0_valid_y = project_labels(q0_valid_raw, label_basis)
        initialization_seed = stable_seed(
            "B1", "actor_iteration_model", d, k
        )
        seed_everything(initialization_seed)
        model = make_actor(equation, device, label_basis)
        training = train_budget_level(
            model,
            equation,
            Q0_BUDGET,
            q0_train_t,
            q0_train_x,
            q0_train_y,
            q0_valid_t,
            q0_valid_x,
            q0_valid_y,
            device=device,
            seed=stable_seed("B1", "actor_iteration_optimizer", d, k, 0),
        )
        q0_record = {
            "status": "complete",
            "iteration": 0,
            "label_source": "terminal pathwise gradient only",
            "train_labels": q0_train_metadata,
            "validation_labels": q0_valid_metadata,
            "train_projection": _projected_label_metadata(
                q0_train_raw, q0_train_y
            ),
            "validation_projection": _projected_label_metadata(
                q0_valid_raw, q0_valid_y
            ),
            "label_generation_seconds": q0_train_metadata["cost"][
                "wall_seconds"
            ]
            + q0_valid_metadata["cost"]["wall_seconds"],
            "training": training,
            "analytic_reference_used": False,
        }
        _atomic_torch_save(q0_path, checkpoint_payload(model, q0_record))
    actors.append(model)
    construction_records.append(q0_record)
    q0_validation = _policy_validation(
        equation,
        model,
        0,
        value_valid_t,
        value_valid_x,
        device=device,
    )
    validation_records.append(q0_validation)
    print(
        json.dumps(
            {
                "iteration": 0,
                "validation_value": q0_validation["mean_value"],
                "training_seconds": q0_record["training"]["training_seconds"],
            }
        ),
        flush=True,
    )

    stable_transitions = 0
    stop_reason = "iteration_cap"
    for iteration in range(1, MAX_IMPROVEMENTS + 1):
        checkpoint_path = _checkpoint_path(d, k, iteration)
        source_actor = actors[-1]
        if checkpoint_path.is_file() and not force:
            model, record = load_converged_actor(checkpoint_path, device)
        else:
            _, train_raw, train_label_metadata = _policy_label_cache(
                equation,
                source_actor,
                iteration,
                improve_train_t,
                improve_train_x,
                split="train",
                device=device,
            )
            _, valid_raw, valid_label_metadata = _policy_label_cache(
                equation,
                source_actor,
                iteration,
                improve_valid_t,
                improve_valid_x,
                split="valid",
                device=device,
            )
            train_y = project_labels(train_raw, label_basis)
            valid_y = project_labels(valid_raw, label_basis)
            source_prediction = predict_actor(
                source_actor,
                improve_valid_t,
                improve_valid_x,
                device=device,
            )
            model = make_actor(equation, device, label_basis)
            model.load_state_dict(source_actor.state_dict())
            training = train_budget_level(
                model,
                equation,
                IMPROVEMENT_BUDGET,
                improve_train_t,
                improve_train_x,
                train_y,
                improve_valid_t,
                improve_valid_x,
                valid_y,
                device=device,
                seed=stable_seed(
                    "B1", "actor_iteration_optimizer", d, k, iteration
                ),
            )
            current_prediction = predict_actor(
                model,
                improve_valid_t,
                improve_valid_x,
                device=device,
            )
            record = {
                "status": "complete",
                "iteration": iteration,
                "source_actor_iteration": iteration - 1,
                "source_actor_sha256": train_label_metadata[
                    "source_actor_sha256"
                ],
                "label_source": (
                    "pathwise gradient through controlled-SDE policy evaluation"
                ),
                "train_labels": train_label_metadata,
                "validation_labels": valid_label_metadata,
                "train_projection": _projected_label_metadata(train_raw, train_y),
                "validation_projection": _projected_label_metadata(
                    valid_raw, valid_y
                ),
                "label_generation_seconds": train_label_metadata["cost"][
                    "wall_seconds"
                ]
                + valid_label_metadata["cost"]["wall_seconds"],
                "training": training,
                "actor_validation_prediction_relative_change": math.sqrt(
                    float(
                        np.mean(
                            np.sum(
                                (current_prediction - source_prediction) ** 2,
                                axis=1,
                            )
                        )
                    )
                    / max(
                        float(np.mean(np.sum(current_prediction**2, axis=1))),
                        np.finfo(np.float64).tiny,
                    )
                ),
                "noisy_bellman_gap": 0.5
                * float(np.mean(np.sum((source_prediction - valid_y) ** 2, axis=1))),
                "analytic_reference_used": False,
            }
        actors.append(model)
        current_validation = _policy_validation(
            equation,
            model,
            iteration,
            value_valid_t,
            value_valid_x,
            device=device,
        )
        improvement = _value_improvement(
            validation_records[-1], current_validation
        )
        record["value_improvement"] = improvement
        _atomic_torch_save(checkpoint_path, checkpoint_payload(model, record))
        construction_records.append(record)
        validation_records.append(current_validation)
        if improvement["small_relative_to_uncertainty"]:
            stable_transitions += 1
        else:
            stable_transitions = 0
        print(
            json.dumps(
                {
                    "iteration": iteration,
                    "validation_value": current_validation["mean_value"],
                    "value_improvement": improvement,
                    "stable_transitions": stable_transitions,
                    "label_seconds": record["label_generation_seconds"],
                    "training_seconds": record["training"]["training_seconds"],
                }
            ),
            flush=True,
        )
        if stable_transitions >= 2:
            stop_reason = "two_consecutive_value_plateaus"
            break

    stopping_decision = {
        "selected_iteration": len(actors) - 1,
        "reason": stop_reason,
        "stable_transitions": stable_transitions,
        "analytic_reference_used": False,
    }
    atomic_json(
        PILOT_DIR / f"stopping_decision_d{d}_k{k}.json",
        stopping_decision,
    )

    # Every checkpoint and the label-only stopping decision are frozen above.
    # Analytic references first become visible inside final evaluation below.
    final_rows, predictions, test_hash = _final_evaluation(
        equation,
        actors,
        construction_records,
        validation_records,
        device=device,
    )
    _atomic_npz(
        PREDICTION_DIR / f"actor_iteration_d{d}_k{k}.npz", **predictions
    )
    csv_rows = [
        {
            key: value
            for key, value in row.items()
            if key
            not in {
                "direct_policy_seed_rows",
                "actor_mlp_seed_rows",
                "direct_policy_costs",
                "actor_mlp_costs",
            }
        }
        for row in final_rows
    ]
    _atomic_csv(PILOT_DIR / f"summary_d{d}_k{k}.csv", csv_rows)
    baselines = _baseline_references(d, k)
    selected = final_rows[stopping_decision["selected_iteration"]]
    best_actor_mlp = min(final_rows, key=lambda row: row["actor_mlp_value_rmse"])
    best_direct = min(
        final_rows, key=lambda row: row["direct_policy_value_rmse"]
    )
    headline = {
        "selected_iteration": stopping_decision["selected_iteration"],
        "selected_actor_mlp_value_rmse": selected["actor_mlp_value_rmse"],
        "selected_direct_policy_value_rmse": selected[
            "direct_policy_value_rmse"
        ],
        "best_posthoc_actor_mlp_iteration": best_actor_mlp["iteration"],
        "best_posthoc_actor_mlp_value_rmse": best_actor_mlp[
            "actor_mlp_value_rmse"
        ],
        "best_posthoc_direct_policy_iteration": best_direct["iteration"],
        "best_posthoc_direct_policy_value_rmse": best_direct[
            "direct_policy_value_rmse"
        ],
        "converged_actor_mlp_value_rmse": baselines[
            "converged_actor_mlp_value_rmse"
        ],
        "selected_actor_mlp_beats_converged_actor_mlp": selected[
            "actor_mlp_value_rmse"
        ]
        < baselines["converged_actor_mlp_value_rmse"],
        "selected_direct_policy_beats_converged_actor_mlp": selected[
            "direct_policy_value_rmse"
        ]
        < baselines["converged_actor_mlp_value_rmse"],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "d": d,
        "k": k,
        "device": str(device),
        "protocol": protocol(d, k),
        "label_subspace": subspace_metadata,
        "stopping_decision": stopping_decision,
        "test_point_sha256": test_hash,
        "iterations": final_rows,
        "baselines": baselines,
        "headline": headline,
        "elapsed_seconds": time.perf_counter() - started,
        "analytic_reference_first_used_after_stopping_decision": True,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_stopping": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    atomic_json(result_path, payload)
    print(json.dumps(headline, indent=2), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        run(args.d, args.k, force=args.force)
    except Exception as error:
        PILOT_DIR.mkdir(parents=True, exist_ok=True)
        atomic_json(
            PILOT_DIR / "failure.json",
            {
                "status": "failed",
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()
