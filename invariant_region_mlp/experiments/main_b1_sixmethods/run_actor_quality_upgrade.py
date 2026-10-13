"""Run one higher-quality teacher-free policy improvement for B1 (400,40)."""

from __future__ import annotations

import argparse
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
    load_converged_actor,
    make_actor,
    predict_actor,
    project_labels,
    train_budget_level,
)
from .actor_iteration import evaluate_policy_mc
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
    stable_seed,
    training_points,
)
from .equation import B1Equation
from .run_actor_quality_ceiling import _evaluate_actor


D = 400
K = 40
TRAIN_POINTS = 8192
VALID_POINTS = 512
LABEL_PATHS = 512
LABEL_TIME_STEPS = 32
LABEL_PARENT_BATCH_SIZE = 16
LABEL_PATH_CHUNK_SIZE = 128

OUTPUT_DIR = RESULTS_DIR / "actor_quality_upgrade" / "upgrade"
CACHE_DIR = OUTPUT_DIR / "cache"
CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints"
PREDICTION_DIR = OUTPUT_DIR / "predictions"
CEILING_RESULT = RESULTS_DIR / "actor_quality_upgrade" / "ceiling" / "result.json"
Q0_CHECKPOINT = (
    RESULTS_DIR
    / "actor_iteration_pilot"
    / "checkpoints"
    / "actor_iteration_d400_k40_q0.pt"
)
UPGRADED_CHECKPOINT = CHECKPOINT_DIR / "actor_quality_upgrade_d400_k40_q1.pt"

UPGRADE_BUDGET = BudgetLevel(
    "actor_quality_upgrade_q1",
    LABEL_TIME_STEPS,
    LABEL_PATHS,
    1,
    TRAIN_POINTS,
    VALID_POINTS,
    6000,
    1500,
    1.0e-4,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _state_sha256(model: torch.nn.Module) -> str:
    digest = hashlib.sha256()
    for key, value in sorted(model.state_dict().items()):
        array = np.ascontiguousarray(value.detach().cpu().numpy())
        digest.update(key.encode())
        digest.update(str(array.shape).encode())
        digest.update(str(array.dtype).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _atomic_torch_save(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    torch.save(payload, temporary)
    temporary.replace(path)


def _projection_metadata(
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


def _label_cache(
    equation: B1Equation,
    source_actor: torch.nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    split: str,
    device: torch.device,
    force: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    cache = CACHE_DIR / f"upgrade_{split}.npz"
    metadata_path = cache.with_suffix(".json")
    point_hash = array_sha256(t, x)
    actor_hash = _state_sha256(source_actor)
    if cache.is_file() and metadata_path.is_file() and not force:
        with np.load(cache, allow_pickle=False) as payload:
            if str(payload["point_hash"].item()) != point_hash:
                raise RuntimeError(f"point hash mismatch in {cache}")
            if str(payload["actor_hash"].item()) != actor_hash:
                raise RuntimeError(f"source actor hash mismatch in {cache}")
            values = np.asarray(payload["values"], dtype=np.float64)
            labels = np.asarray(payload["labels"], dtype=np.float64)
        return values, labels, json.loads(metadata_path.read_text())

    values, labels, cost = evaluate_policy_mc(
        equation,
        source_actor,
        t,
        x,
        time_steps=LABEL_TIME_STEPS,
        paths=LABEL_PATHS,
        repetition=0,
        seed_namespace=f"actor_quality_upgrade_q1_{split}_labels",
        device=device,
        parent_batch_size=LABEL_PARENT_BATCH_SIZE,
        path_chunk_size=LABEL_PATH_CHUNK_SIZE,
        compute_z=True,
    )
    if labels is None:
        raise RuntimeError("policy evaluation did not return gradient labels")
    metadata = {
        "status": "complete",
        "split": split,
        "points": len(t),
        "paths_per_label": LABEL_PATHS,
        "sde_steps_per_label": LABEL_TIME_STEPS,
        "point_sha256": point_hash,
        "source_actor_sha256": actor_hash,
        "source_actor": "frozen pathwise q0",
        "label_source": "pathwise gradient through linear policy evaluation",
        "cost": cost,
        "analytic_reference_used": False,
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


def _load_gate() -> dict[str, Any]:
    if not CEILING_RESULT.is_file():
        raise FileNotFoundError(f"Gate A result is missing: {CEILING_RESULT}")
    result = json.loads(CEILING_RESULT.read_text())
    if result.get("status") != "complete":
        raise RuntimeError("Gate A is incomplete")
    if result.get("headline", {}).get("upgrade_gate_passed") is not True:
        raise RuntimeError("Gate A did not show material Actor headroom")
    return result


def run(*, force: bool = False) -> dict[str, Any]:
    result_path = OUTPUT_DIR / "result.json"
    if result_path.is_file() and not force:
        result = json.loads(result_path.read_text())
        if result.get("status") == "complete":
            print(json.dumps(result["headline"], indent=2), flush=True)
            return result
    gate = _load_gate()
    if not Q0_CHECKPOINT.is_file():
        raise FileNotFoundError(f"frozen q0 checkpoint is missing: {Q0_CHECKPOINT}")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Actor-quality upgrade requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    config = load_config()
    equation = B1Equation(
        d=D,
        k=K,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    started = time.perf_counter()
    source_actor, source_metadata = load_converged_actor(Q0_CHECKPOINT, device)
    source_hash = _state_sha256(source_actor)

    point_seed = stable_seed("B1", "actor_iteration_points", D, K)
    train_t, train_x = training_points(
        D, K, "actor_iteration_improvement_train", point_seed, TRAIN_POINTS
    )
    valid_t, valid_x = training_points(
        D, K, "actor_iteration_improvement_valid", point_seed, VALID_POINTS
    )

    if UPGRADED_CHECKPOINT.is_file() and not force:
        model, training_record = load_converged_actor(
            UPGRADED_CHECKPOINT, device
        )
    else:
        _, train_raw, train_label_metadata = _label_cache(
            equation,
            source_actor,
            train_t,
            train_x,
            split="train",
            device=device,
            force=force,
        )
        _, valid_raw, valid_label_metadata = _label_cache(
            equation,
            source_actor,
            valid_t,
            valid_x,
            split="valid",
            device=device,
            force=force,
        )
        label_basis = source_actor.label_basis.detach().cpu().numpy()
        train_y = project_labels(train_raw, label_basis)
        valid_y = project_labels(valid_raw, label_basis)
        source_validation_prediction = predict_actor(
            source_actor, valid_t, valid_x, device=device
        )
        model = make_actor(equation, device, label_basis)
        model.load_state_dict(source_actor.state_dict())
        training = train_budget_level(
            model,
            equation,
            UPGRADE_BUDGET,
            train_t,
            train_x,
            train_y,
            valid_t,
            valid_x,
            valid_y,
            device=device,
            seed=stable_seed("B1", "actor_iteration_optimizer", D, K, 1),
        )
        current_validation_prediction = predict_actor(
            model, valid_t, valid_x, device=device
        )
        training_record = {
            "status": "complete",
            "iteration": 1,
            "source_actor_iteration": 0,
            "source_actor_checkpoint": str(
                Q0_CHECKPOINT.relative_to(RESULTS_DIR.parent)
            ),
            "source_actor_checkpoint_sha256": _sha256(Q0_CHECKPOINT),
            "source_actor_state_sha256": source_hash,
            "source_actor_metadata": source_metadata,
            "label_source": (
                "pathwise gradient through controlled-SDE linear policy evaluation"
            ),
            "train_labels": train_label_metadata,
            "validation_labels": valid_label_metadata,
            "train_projection": _projection_metadata(train_raw, train_y),
            "validation_projection": _projection_metadata(valid_raw, valid_y),
            "label_generation_seconds": (
                train_label_metadata["cost"]["wall_seconds"]
                + valid_label_metadata["cost"]["wall_seconds"]
            ),
            "training": training,
            "actor_validation_prediction_relative_change": math.sqrt(
                float(
                    np.mean(
                        np.sum(
                            (
                                current_validation_prediction
                                - source_validation_prediction
                            )
                            ** 2,
                            axis=1,
                        )
                    )
                )
                / max(
                    float(
                        np.mean(
                            np.sum(current_validation_prediction**2, axis=1)
                        )
                    ),
                    np.finfo(np.float64).tiny,
                )
            ),
            "analytic_reference_used_for_training": False,
            "analytic_reference_used_for_checkpoint_selection": False,
            "policy_improvements": 1,
        }
        _atomic_torch_save(
            UPGRADED_CHECKPOINT, checkpoint_payload(model, training_record)
        )

    # The checkpoint is frozen before exact references are constructed below.
    points = fixed_test_points(D, K, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    evaluation, predictions = _evaluate_actor(
        "upgraded_q1",
        UPGRADED_CHECKPOINT,
        equation,
        times,
        points,
        exact_u,
        exact_z,
        device=device,
    )
    _atomic_npz(PREDICTION_DIR / "upgrade_predictions.npz", **predictions)

    ceiling_rows = {row["actor"]: row for row in gate["rows"]}
    old_q1 = ceiling_rows["teacher_free_q1"]
    converged = ceiling_rows["converged_actor"]
    upgraded_rmse = evaluation["corrected_value_rmse"]
    headline = {
        "old_q1_corrected_value_rmse": old_q1["corrected_value_rmse"],
        "upgraded_q1_corrected_value_rmse": upgraded_rmse,
        "converged_actor_corrected_value_rmse": converged[
            "corrected_value_rmse"
        ],
        "upgraded_relative_improvement_over_old_q1": 1.0
        - upgraded_rmse / old_q1["corrected_value_rmse"],
        "upgraded_beats_old_q1": upgraded_rmse
        < old_q1["corrected_value_rmse"],
        "upgraded_beats_converged_actor": upgraded_rmse
        < converged["corrected_value_rmse"],
        "upgraded_direct_actor_gradient_rmse": evaluation[
            "direct_actor_gradient_rmse"
        ],
        "upgraded_corrected_value_bias": evaluation["corrected_value_bias"],
        "upgraded_mean_analytic_policy_gap": evaluation[
            "mean_analytic_policy_gap"
        ],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "device": str(device),
        "gate_a_headline": gate["headline"],
        "protocol": {
            "d": D,
            "k": K,
            "improvement_training_points": TRAIN_POINTS,
            "improvement_validation_points": VALID_POINTS,
            "paths_per_label": LABEL_PATHS,
            "sde_steps_per_label": LABEL_TIME_STEPS,
            "architecture": {
                "hidden_width": ACTOR_WIDTH,
                "hidden_layers": ACTOR_DEPTH,
                "activation": "silu",
            },
            "policy_improvements": 1,
            "initial_actor": "frozen pathwise q0",
            "final_evaluation": "64-step martingale CV, 256 paths, 3 seeds",
            "analytic_reference_first_used_after_checkpoint_frozen": True,
            "analytic_reference_used_for_training": False,
            "analytic_reference_used_for_checkpoint_selection": False,
            "actor_structure_prior": (
                "terminal-gradient-derived rank-40 output subspace with supplied k"
            ),
            "budget": UPGRADE_BUDGET.__dict__,
        },
        "training": training_record,
        "evaluation": evaluation,
        "headline": headline,
        "elapsed_seconds": time.perf_counter() - started,
        "failed_runs": 0,
        "nan_failures": 0,
        "memory_failures": 0,
    }
    atomic_json(result_path, payload)
    print(json.dumps(headline, indent=2), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        run(force=args.force)
    except Exception as error:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        atomic_json(
            OUTPUT_DIR / "failure.json",
            {
                "status": "failed",
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()
