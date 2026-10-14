"""End-to-end rank-k output-projection ablation for one-step Actor Iteration."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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

from .actor import ActorNet
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
from .actor_iteration import evaluate_policy_mc, evaluate_policy_mc_with_analytic_gap
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
from .neural_common import count_parameters
from .run_actor_no_low_rank_ablation import (
    _load_plain_actor,
    _plain_checkpoint_payload,
)


CELLS = ((100, 10), (400, 40), (400, 1))
Q0_TRAIN_POINTS = 4096
Q0_VALID_POINTS = 1024
Q0_PATHS = 256
TRAIN_POINTS = 8192
VALID_POINTS = 512
LABEL_PATHS = 512
LABEL_TIME_STEPS = 32
LABEL_PARENT_BATCH_SIZE = 16
LABEL_PATH_CHUNK_SIZE = 128
OPTIMIZER_STEPS = 6000
LEARNING_RATE = 1.0e-4
EVALUATION_TIME_STEPS = 64
EVALUATION_PATHS = 256
EVALUATION_SEEDS = 3
EVALUATION_PARENT_BATCH_SIZE = 32
EVALUATION_PATH_CHUNK_SIZE = 64

OUTPUT_DIR = RESULTS_DIR / "actor_iteration_projection_ablation"
CACHE_DIR = OUTPUT_DIR / "cache"
CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints"
PREDICTION_DIR = OUTPUT_DIR / "predictions"
PILOT_DIR = RESULTS_DIR / "actor_iteration_pilot"

Q0_BUDGET = BudgetLevel(
    "actor_iteration_projection_ablation_q0",
    1,
    Q0_PATHS,
    1,
    Q0_TRAIN_POINTS,
    Q0_VALID_POINTS,
    6000,
    1500,
    5.0e-4,
)


@dataclass(frozen=True)
class Phase:
    name: str
    train_point_method: str
    valid_point_method: str
    point_seed_namespace: str
    label_seed_prefix: str
    optimizer_seed_namespace: str


PHASES = (
    Phase(
        "improvement",
        "actor_iteration_improvement_train",
        "actor_iteration_improvement_valid",
        "actor_iteration_points",
        "actor_quality_upgrade_q1",
        "actor_iteration_optimizer",
    ),
    Phase(
        "label_refresh",
        "actor_q2_control_train",
        "actor_q2_control_valid",
        "actor_q2_control_points",
        "actor_q2_control_continued_q1",
        "actor_q2_control_optimizer",
    ),
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


def _parameter_sha256(model: torch.nn.Module) -> str:
    digest = hashlib.sha256()
    for key, value in sorted(model.named_parameters()):
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


def _relative_path(path: Path) -> str:
    return str(path.relative_to(RESULTS_DIR.parent))


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    return (
        float(np.mean(array)),
        1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array)),
    )


def _projection_metadata(raw: np.ndarray, target: np.ndarray) -> dict[str, float]:
    residual = math.sqrt(float(np.mean(np.sum((raw - target) ** 2, axis=1))))
    scale = math.sqrt(float(np.mean(np.sum(raw**2, axis=1))))
    return {
        "raw_vector_rmse": scale,
        "projection_residual_rmse": residual,
        "relative_projection_residual": residual
        / max(scale, np.finfo(np.float64).tiny),
    }


def _phase_budget(structure: str, phase: Phase) -> BudgetLevel:
    return BudgetLevel(
        f"actor_iteration_{structure}_{phase.name}",
        LABEL_TIME_STEPS,
        LABEL_PATHS,
        1,
        TRAIN_POINTS,
        VALID_POINTS,
        OPTIMIZER_STEPS,
        OPTIMIZER_STEPS,
        LEARNING_RATE,
    )


def _optimizer_seed(d: int, k: int, phase: Phase) -> int:
    if phase.name == "improvement":
        return stable_seed("B1", phase.optimizer_seed_namespace, d, k, 1)
    return stable_seed("B1", phase.optimizer_seed_namespace, d, k)


def protocol_record(d: int, k: int) -> dict[str, Any]:
    return {
        "cell": {"d": d, "k": k},
        "structures": {
            "projected": {
                "actor_output": "terminal-label-derived rank-k projection",
                "training_targets": "raw labels projected into the same subspace",
            },
            "unprojected": {
                "actor_output": "plain ambient d-dimensional output",
                "training_targets": "raw ambient labels",
                "label_basis_constructed_or_read": False,
            },
        },
        "architecture": {
            "hidden_width": ACTOR_WIDTH,
            "hidden_layers": ACTOR_DEPTH,
            "activation": "silu",
            "trainable_parameter_count_is_matched": True,
        },
        "q0": {
            "training_points": Q0_TRAIN_POINTS,
            "validation_points": Q0_VALID_POINTS,
            "terminal_paths": Q0_PATHS,
            "budget": asdict(Q0_BUDGET),
            "shared_raw_labels": True,
            "shared_model_initialization": True,
            "shared_optimizer_minibatches": True,
        },
        "phases": [
            {
                **asdict(phase),
                "training_points": TRAIN_POINTS,
                "validation_points": VALID_POINTS,
                "paths_per_label": LABEL_PATHS,
                "sde_steps_per_label": LABEL_TIME_STEPS,
                "budget": asdict(_phase_budget("STRUCTURE", phase)),
                "source_policy": "the arm's frozen q0",
                "shared_points_across_structures": True,
                "common_brownian_numbers_across_structures": True,
                "shared_optimizer_minibatches_across_structures": True,
            }
            for phase in PHASES
        ],
        "phase_2_semantics": (
            "continue q1 using a fresh q0-policy-evaluation label bank; not q2"
        ),
        "final_evaluation": {
            "test_points": 960,
            "paths": EVALUATION_PATHS,
            "seeds": EVALUATION_SEEDS,
            "time_steps": EVALUATION_TIME_STEPS,
            "seed_namespace": "actor_iteration_direct_test",
            "martingale_control_variate": True,
            "common_brownian_numbers_across_structures": True,
        },
        "analytic_reference_first_used_after_both_final_checkpoints_frozen": True,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
        "training_seed_repetitions": 1,
        "policy_improvements": 1,
        "q2_q3_or_q4_permitted": False,
    }


def _q0_checkpoint(d: int, k: int, structure: str) -> Path:
    if structure == "projected":
        return (
            PILOT_DIR
            / "checkpoints"
            / f"actor_iteration_d{d}_k{k}_q0.pt"
        )
    return (
        CHECKPOINT_DIR
        / f"d{d}_k{k}"
        / "unprojected_q0.pt"
    )


def _phase_checkpoint(d: int, k: int, structure: str, phase: Phase) -> Path:
    if structure == "projected" and (d, k) == (400, 40):
        if phase.name == "improvement":
            return (
                RESULTS_DIR
                / "actor_quality_upgrade"
                / "upgrade"
                / "checkpoints"
                / "actor_quality_upgrade_d400_k40_q1.pt"
            )
        return (
            RESULTS_DIR
            / "actor_q2_control"
            / "checkpoints"
            / "actor_q2_control_continued_q1.pt"
        )
    return (
        CHECKPOINT_DIR
        / f"d{d}_k{k}"
        / f"{structure}_{phase.name}.pt"
    )


def _load_actor(
    path: Path,
    equation: B1Equation,
    structure: str,
    device: torch.device,
) -> tuple[torch.nn.Module, dict[str, Any]]:
    if structure == "projected":
        return load_converged_actor(path, device)
    return _load_plain_actor(path, equation, device)


def _checkpoint_payload(
    model: torch.nn.Module,
    equation: B1Equation,
    structure: str,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    if structure == "projected":
        return checkpoint_payload(model, metadata)
    return _plain_checkpoint_payload(model, equation, metadata)


def _load_q0_raw_labels(
    equation: B1Equation,
    split: str,
    t: np.ndarray,
    x: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    path = (
        PILOT_DIR
        / "cache"
        / f"d{equation.d}_k{equation.k}"
        / f"q0_{split}.npz"
    )
    metadata_path = path.with_suffix(".json")
    if not path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(f"shared q0 raw-label cache is missing: {path}")
    with np.load(path, allow_pickle=False) as payload:
        labels = np.asarray(payload["labels"], dtype=np.float64)
        point_hash = str(payload["point_hash"].item())
    expected_hash = array_sha256(t, x)
    if point_hash != expected_hash:
        raise RuntimeError(f"q0 point hash mismatch in {path}")
    if labels.shape != (len(t), equation.d) or not np.all(np.isfinite(labels)):
        raise RuntimeError(f"invalid q0 labels in {path}")
    return labels, {
        **json.loads(metadata_path.read_text()),
        "cache": _relative_path(path),
        "cache_sha256": _sha256(path),
        "raw_label_sha256": array_sha256(labels),
    }


def _paired_initialization_record(
    equation: B1Equation,
    projected_q0: torch.nn.Module,
    device: torch.device,
) -> dict[str, Any]:
    initialization_seed = stable_seed(
        "B1", "actor_iteration_model", equation.d, equation.k
    )
    basis = projected_q0.label_basis.detach().cpu().numpy()
    seed_everything(initialization_seed)
    projected = make_actor(equation, device, basis)
    projected_parameter_hash = _parameter_sha256(projected)
    seed_everything(initialization_seed)
    unprojected = ActorNet(
        equation.d, width=ACTOR_WIDTH, depth=ACTOR_DEPTH
    ).to(device)
    unprojected_parameter_hash = _parameter_sha256(unprojected)
    record = {
        "seed": initialization_seed,
        "projected_parameter_sha256": projected_parameter_hash,
        "unprojected_parameter_sha256": unprojected_parameter_hash,
        "parameter_states_identical": (
            projected_parameter_hash == unprojected_parameter_hash
        ),
        "projected_trainable_parameters": count_parameters(projected),
        "unprojected_trainable_parameters": count_parameters(unprojected),
    }
    del projected, unprojected
    if device.type == "cuda":
        torch.cuda.empty_cache()
    if not record["parameter_states_identical"]:
        raise RuntimeError("q0 parameter initializations are not identical")
    return record


def _train_unprojected_q0(
    equation: B1Equation,
    projected_q0: torch.nn.Module,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_y: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_y: np.ndarray,
    label_metadata: dict[str, Any],
    *,
    device: torch.device,
    force: bool,
) -> tuple[torch.nn.Module, dict[str, Any], dict[str, Any]]:
    path = _q0_checkpoint(equation.d, equation.k, "unprojected")
    initialization = _paired_initialization_record(
        equation, projected_q0, device
    )
    if path.is_file() and not force:
        model, record = _load_actor(path, equation, "unprojected", device)
        return model, record, initialization

    seed_everything(int(initialization["seed"]))
    model = ActorNet(
        equation.d, width=ACTOR_WIDTH, depth=ACTOR_DEPTH
    ).to(device)
    training = train_budget_level(
        model,
        equation,
        Q0_BUDGET,
        train_t,
        train_x,
        train_y,
        valid_t,
        valid_x,
        valid_y,
        device=device,
        seed=stable_seed(
            "B1", "actor_iteration_optimizer", equation.d, equation.k, 0
        ),
    )
    record = {
        "status": "complete",
        "experiment": "actor_iteration_projection_ablation",
        "structure": "unprojected",
        "stage": "q0",
        "label_source": "shared raw terminal pathwise gradients",
        "train_labels": label_metadata["train"],
        "validation_labels": label_metadata["valid"],
        "target_projection": None,
        "initialization": initialization,
        "training": training,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    _atomic_torch_save(
        path,
        _checkpoint_payload(model, equation, "unprojected", record),
    )
    print(
        json.dumps(
            {
                "cell": [equation.d, equation.k],
                "structure": "unprojected",
                "stage": "q0",
                "best_step": training["best_step"],
                "steps_executed": training["steps_executed"],
                "validation_label_rmse": training[
                    "selected_validation_label_rmse"
                ],
            }
        ),
        flush=True,
    )
    return model, record, initialization


def _phase_points(
    d: int, k: int, phase: Phase
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    point_seed = stable_seed("B1", phase.point_seed_namespace, d, k)
    train_t, train_x = training_points(
        d, k, phase.train_point_method, point_seed, TRAIN_POINTS
    )
    valid_t, valid_x = training_points(
        d, k, phase.valid_point_method, point_seed, VALID_POINTS
    )
    return train_t, train_x, valid_t, valid_x


def _policy_label_cache(
    equation: B1Equation,
    source_actor: torch.nn.Module,
    structure: str,
    phase: Phase,
    split: str,
    t: np.ndarray,
    x: np.ndarray,
    *,
    device: torch.device,
    force: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    path = (
        CACHE_DIR
        / f"d{equation.d}_k{equation.k}"
        / f"{structure}_{phase.name}_{split}.npz"
    )
    metadata_path = path.with_suffix(".json")
    point_hash = array_sha256(t, x)
    actor_hash = _state_sha256(source_actor)
    namespace = f"{phase.label_seed_prefix}_{split}_labels"
    if path.is_file() and metadata_path.is_file() and not force:
        with np.load(path, allow_pickle=False) as payload:
            values = np.asarray(payload["values"], dtype=np.float64)
            labels = np.asarray(payload["labels"], dtype=np.float64)
            observed_point_hash = str(payload["point_hash"].item())
            observed_actor_hash = str(payload["actor_hash"].item())
        metadata = json.loads(metadata_path.read_text())
        if observed_point_hash != point_hash or observed_actor_hash != actor_hash:
            raise RuntimeError(f"policy-label cache provenance mismatch: {path}")
        if metadata.get("seed_namespace") != namespace:
            raise RuntimeError(f"policy-label seed mismatch: {path}")
        return values, labels, metadata

    values, labels, cost = evaluate_policy_mc(
        equation,
        source_actor,
        t,
        x,
        time_steps=LABEL_TIME_STEPS,
        paths=LABEL_PATHS,
        repetition=0,
        seed_namespace=namespace,
        device=device,
        parent_batch_size=LABEL_PARENT_BATCH_SIZE,
        path_chunk_size=LABEL_PATH_CHUNK_SIZE,
        compute_z=True,
    )
    if labels is None:
        raise RuntimeError("policy evaluation did not return gradient labels")
    metadata = {
        "status": "complete",
        "structure": structure,
        "phase": phase.name,
        "split": split,
        "points": len(t),
        "paths_per_label": LABEL_PATHS,
        "sde_steps_per_label": LABEL_TIME_STEPS,
        "point_sha256": point_hash,
        "source_actor_state_sha256": actor_hash,
        "seed_namespace": namespace,
        "common_brownian_numbers_across_structures": True,
        "cost": cost,
        "analytic_reference_used": False,
    }
    _atomic_npz(
        path,
        values=values,
        labels=labels,
        point_hash=np.asarray(point_hash),
        actor_hash=np.asarray(actor_hash),
    )
    atomic_json(metadata_path, metadata)
    return values, labels, metadata


def _make_phase_actor(
    equation: B1Equation,
    structure: str,
    initial_actor: torch.nn.Module,
    device: torch.device,
) -> torch.nn.Module:
    if structure == "projected":
        basis = initial_actor.label_basis.detach().cpu().numpy()
        model = make_actor(equation, device, basis)
    else:
        model = ActorNet(
            equation.d, width=ACTOR_WIDTH, depth=ACTOR_DEPTH
        ).to(device)
    model.load_state_dict(initial_actor.state_dict())
    return model


def _external_projected_record(
    path: Path,
    actor: torch.nn.Module,
    phase: Phase,
) -> dict[str, Any]:
    return {
        "status": "complete",
        "experiment": "actor_iteration_projection_ablation",
        "structure": "projected",
        "phase": phase.name,
        "reused_validated_current_checkpoint": True,
        "checkpoint": _relative_path(path),
        "checkpoint_sha256": _sha256(path),
        "actor_state_sha256": _state_sha256(actor),
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }


def _train_phase(
    equation: B1Equation,
    structure: str,
    phase: Phase,
    source_q0: torch.nn.Module,
    initial_actor: torch.nn.Module,
    *,
    device: torch.device,
    force: bool,
) -> tuple[torch.nn.Module, dict[str, Any]]:
    path = _phase_checkpoint(equation.d, equation.k, structure, phase)
    if (
        structure == "projected"
        and (equation.d, equation.k) == (400, 40)
    ):
        if not path.is_file():
            raise FileNotFoundError(f"current projected checkpoint is missing: {path}")
        actor, _ = _load_actor(path, equation, structure, device)
        return actor, _external_projected_record(path, actor, phase)
    if path.is_file() and not force:
        return _load_actor(path, equation, structure, device)

    train_t, train_x, valid_t, valid_x = _phase_points(
        equation.d, equation.k, phase
    )
    _, train_raw, train_label_metadata = _policy_label_cache(
        equation,
        source_q0,
        structure,
        phase,
        "train",
        train_t,
        train_x,
        device=device,
        force=force,
    )
    _, valid_raw, valid_label_metadata = _policy_label_cache(
        equation,
        source_q0,
        structure,
        phase,
        "valid",
        valid_t,
        valid_x,
        device=device,
        force=force,
    )
    if structure == "projected":
        basis = source_q0.label_basis.detach().cpu().numpy()
        train_y = project_labels(train_raw, basis)
        valid_y = project_labels(valid_raw, basis)
    else:
        train_y = train_raw
        valid_y = valid_raw
    model = _make_phase_actor(
        equation, structure, initial_actor, device
    )
    training = train_budget_level(
        model,
        equation,
        _phase_budget(structure, phase),
        train_t,
        train_x,
        train_y,
        valid_t,
        valid_x,
        valid_y,
        device=device,
        seed=_optimizer_seed(equation.d, equation.k, phase),
    )
    if training["steps_executed"] != OPTIMIZER_STEPS:
        raise RuntimeError(
            f"{structure}/{phase.name} did not execute {OPTIMIZER_STEPS} steps"
        )
    record = {
        "status": "complete",
        "experiment": "actor_iteration_projection_ablation",
        "structure": structure,
        "phase": phase.name,
        "source_policy": "the arm's frozen q0",
        "source_q0_checkpoint": _relative_path(
            _q0_checkpoint(equation.d, equation.k, structure)
        ),
        "source_q0_checkpoint_sha256": _sha256(
            _q0_checkpoint(equation.d, equation.k, structure)
        ),
        "source_q0_actor_state_sha256": _state_sha256(source_q0),
        "initial_actor_state_sha256": _state_sha256(initial_actor),
        "train_point_sha256": array_sha256(train_t, train_x),
        "validation_point_sha256": array_sha256(valid_t, valid_x),
        "train_labels": train_label_metadata,
        "validation_labels": valid_label_metadata,
        "train_target_projection": (
            _projection_metadata(train_raw, train_y)
            if structure == "projected"
            else None
        ),
        "validation_target_projection": (
            _projection_metadata(valid_raw, valid_y)
            if structure == "projected"
            else None
        ),
        "labels_projected": structure == "projected",
        "label_generation_seconds": (
            train_label_metadata["cost"]["wall_seconds"]
            + valid_label_metadata["cost"]["wall_seconds"]
        ),
        "optimizer_seed": _optimizer_seed(equation.d, equation.k, phase),
        "training": training,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
    }
    _atomic_torch_save(
        path,
        _checkpoint_payload(model, equation, structure, record),
    )
    print(
        json.dumps(
            {
                "cell": [equation.d, equation.k],
                "structure": structure,
                "phase": phase.name,
                "best_step": training["best_step"],
                "steps_executed": training["steps_executed"],
                "label_seconds": record["label_generation_seconds"],
                "training_seconds": training["training_seconds"],
                "validation_label_rmse": training[
                    "selected_validation_label_rmse"
                ],
            }
        ),
        flush=True,
    )
    return model, record


def _evaluate_actor(
    name: str,
    actor: torch.nn.Module,
    checkpoint: Path,
    equation: B1Equation,
    times: np.ndarray,
    points: np.ndarray,
    exact_u: np.ndarray,
    exact_z: np.ndarray,
    *,
    device: torch.device,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    direct_actor = predict_actor(actor, times, points, device=device)
    direct_error = direct_actor - exact_z
    direct_gradient_rmse = math.sqrt(
        float(np.mean(np.sum(direct_error * direct_error, axis=1)))
    )
    arrays: dict[str, list[np.ndarray]] = {
        "raw": [],
        "corrected": [],
        "policy_gap": [],
        "raw_estimator_variance": [],
        "corrected_estimator_variance": [],
        "policy_gap_estimator_variance": [],
    }
    seed_rows: list[dict[str, Any]] = []
    costs: list[dict[str, Any]] = []
    for repetition in range(EVALUATION_SEEDS):
        (
            raw,
            corrected,
            policy_gap,
            raw_variance,
            corrected_variance,
            policy_gap_variance,
            cost,
        ) = evaluate_policy_mc_with_analytic_gap(
            equation,
            actor,
            times,
            points,
            time_steps=EVALUATION_TIME_STEPS,
            paths=EVALUATION_PATHS,
            repetition=repetition,
            seed_namespace="actor_iteration_direct_test",
            device=device,
            parent_batch_size=EVALUATION_PARENT_BATCH_SIZE,
            path_chunk_size=EVALUATION_PATH_CHUNK_SIZE,
            coarse_time_steps=16,
            maximum_time_steps=64,
        )
        raw_error = raw - exact_u
        corrected_error = corrected - exact_u
        identity_residual = corrected_error - policy_gap
        seed_rows.append(
            {
                "repetition": repetition,
                "raw_value_rmse": float(np.sqrt(np.mean(raw_error**2))),
                "corrected_value_rmse": float(
                    np.sqrt(np.mean(corrected_error**2))
                ),
                "raw_value_bias": float(np.mean(raw_error)),
                "corrected_value_bias": float(np.mean(corrected_error)),
                "mean_analytic_policy_gap": float(np.mean(policy_gap)),
                "analytic_policy_gap_rms": float(
                    np.sqrt(np.mean(policy_gap**2))
                ),
                "mean_identity_residual": float(np.mean(identity_residual)),
                "identity_residual_rms": float(
                    np.sqrt(np.mean(identity_residual**2))
                ),
            }
        )
        for key, value in (
            ("raw", raw),
            ("corrected", corrected),
            ("policy_gap", policy_gap),
            ("raw_estimator_variance", raw_variance),
            ("corrected_estimator_variance", corrected_variance),
            ("policy_gap_estimator_variance", policy_gap_variance),
        ):
            arrays[key].append(value)
        costs.append(cost)
    stacked = {key: np.stack(value) for key, value in arrays.items()}
    raw_rmse, raw_ci95 = _mean_ci95(
        [row["raw_value_rmse"] for row in seed_rows]
    )
    corrected_rmse, corrected_ci95 = _mean_ci95(
        [row["corrected_value_rmse"] for row in seed_rows]
    )
    corrected_bias, corrected_bias_ci95 = _mean_ci95(
        [row["corrected_value_bias"] for row in seed_rows]
    )
    mean_gap, mean_gap_ci95 = _mean_ci95(
        [row["mean_analytic_policy_gap"] for row in seed_rows]
    )
    mean_residual, mean_residual_ci95 = _mean_ci95(
        [row["mean_identity_residual"] for row in seed_rows]
    )
    gap_rms, gap_rms_ci95 = _mean_ci95(
        [row["analytic_policy_gap_rms"] for row in seed_rows]
    )
    raw_variance = float(np.mean(stacked["raw_estimator_variance"]))
    corrected_variance = float(
        np.mean(stacked["corrected_estimator_variance"])
    )
    row = {
        "actor": name,
        "checkpoint": _relative_path(checkpoint),
        "checkpoint_sha256": _sha256(checkpoint),
        "actor_state_sha256": _state_sha256(actor),
        "direct_actor_gradient_rmse": direct_gradient_rmse,
        "raw_value_rmse": raw_rmse,
        "raw_value_ci95": raw_ci95,
        "corrected_value_rmse": corrected_rmse,
        "corrected_value_ci95": corrected_ci95,
        "corrected_value_bias": corrected_bias,
        "corrected_value_bias_ci95": corrected_bias_ci95,
        "mean_analytic_policy_gap": mean_gap,
        "mean_analytic_policy_gap_ci95": mean_gap_ci95,
        "analytic_policy_gap_rms": gap_rms,
        "analytic_policy_gap_rms_ci95": gap_rms_ci95,
        "mean_identity_residual": mean_residual,
        "mean_identity_residual_ci95": mean_residual_ci95,
        "raw_estimator_variance": raw_variance,
        "corrected_estimator_variance": corrected_variance,
        "variance_reduction_factor": raw_variance
        / max(corrected_variance, np.finfo(np.float64).tiny),
        "seed_rows": seed_rows,
        "costs": costs,
        "total_evaluation_seconds": float(
            sum(cost["wall_seconds"] for cost in costs)
        ),
    }
    prediction_arrays = {
        f"{name}_{key}": value for key, value in stacked.items()
    }
    prediction_arrays[f"{name}_direct_actor"] = direct_actor
    return row, prediction_arrays


def _paired_difference(
    projected: dict[str, Any],
    unprojected: dict[str, Any],
    metric: str,
) -> dict[str, Any]:
    values = [
        float(right[metric]) - float(left[metric])
        for left, right in zip(
            projected["seed_rows"],
            unprojected["seed_rows"],
            strict=True,
        )
    ]
    mean, ci95 = _mean_ci95(values)
    return {
        "unprojected_minus_projected_by_seed": values,
        "mean": mean,
        "ci95": ci95,
    }


def _comparison(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    projected = rows["projected"]
    unprojected = rows["unprojected"]
    return {
        "unprojected_to_projected_value_rmse_ratio": unprojected[
            "corrected_value_rmse"
        ]
        / projected["corrected_value_rmse"],
        "unprojected_relative_value_rmse_change": unprojected[
            "corrected_value_rmse"
        ]
        / projected["corrected_value_rmse"]
        - 1.0,
        "unprojected_to_projected_mean_policy_gap_ratio": unprojected[
            "mean_analytic_policy_gap"
        ]
        / projected["mean_analytic_policy_gap"],
        "unprojected_to_projected_gap_rms_ratio": unprojected[
            "analytic_policy_gap_rms"
        ]
        / projected["analytic_policy_gap_rms"],
        "unprojected_beats_projected_on_value_rmse": unprojected[
            "corrected_value_rmse"
        ]
        < projected["corrected_value_rmse"],
        "paired_differences": {
            metric: _paired_difference(projected, unprojected, metric)
            for metric in (
                "corrected_value_rmse",
                "corrected_value_bias",
                "mean_analytic_policy_gap",
                "analytic_policy_gap_rms",
            )
        },
        "interpretation_requires_independent_training_seed_confirmation": True,
    }


def run_cell(d: int, k: int, *, force: bool = False) -> dict[str, Any]:
    if (d, k) not in CELLS:
        raise ValueError(f"unsupported cell {(d, k)}")
    result_path = OUTPUT_DIR / f"result_d{d}_k{k}.json"
    if result_path.is_file() and not force:
        result = json.loads(result_path.read_text())
        if result.get("status") == "complete":
            print(json.dumps(result["headline"], indent=2), flush=True)
            return result
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("projection ablation requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    config = load_config()
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    started = time.perf_counter()
    protocol = protocol_record(d, k)

    point_seed = stable_seed("B1", "actor_iteration_points", d, k)
    q0_train_t, q0_train_x = training_points(
        d, k, "actor_iteration_q0_train", point_seed, Q0_TRAIN_POINTS
    )
    q0_valid_t, q0_valid_x = training_points(
        d, k, "actor_iteration_q0_valid", point_seed, Q0_VALID_POINTS
    )
    q0_train_raw, q0_train_metadata = _load_q0_raw_labels(
        equation, "train", q0_train_t, q0_train_x
    )
    q0_valid_raw, q0_valid_metadata = _load_q0_raw_labels(
        equation, "valid", q0_valid_t, q0_valid_x
    )
    shared_q0_labels = {
        "train": q0_train_metadata,
        "valid": q0_valid_metadata,
        "train_point_sha256": array_sha256(q0_train_t, q0_train_x),
        "validation_point_sha256": array_sha256(q0_valid_t, q0_valid_x),
    }

    projected_q0_path = _q0_checkpoint(d, k, "projected")
    if not projected_q0_path.is_file():
        raise FileNotFoundError(f"projected q0 is missing: {projected_q0_path}")
    projected_q0, projected_q0_metadata = _load_actor(
        projected_q0_path, equation, "projected", device
    )
    unprojected_q0, unprojected_q0_record, initialization = (
        _train_unprojected_q0(
            equation,
            projected_q0,
            q0_train_t,
            q0_train_x,
            q0_train_raw,
            q0_valid_t,
            q0_valid_x,
            q0_valid_raw,
            {"train": q0_train_metadata, "valid": q0_valid_metadata},
            device=device,
            force=force,
        )
    )
    actors = {"projected": projected_q0, "unprojected": unprojected_q0}
    construction: dict[str, dict[str, Any]] = {
        "projected": {
            "q0": {
                "reused_validated_current_checkpoint": True,
                "checkpoint": _relative_path(projected_q0_path),
                "checkpoint_sha256": _sha256(projected_q0_path),
                "actor_state_sha256": _state_sha256(projected_q0),
                "training_summary": {
                    key: projected_q0_metadata["training"][key]
                    for key in (
                        "best_step",
                        "steps_executed",
                        "selected_validation_label_rmse",
                        "optimizer_converged",
                        "training_seconds",
                    )
                },
            }
        },
        "unprojected": {"q0": unprojected_q0_record},
    }

    for phase in PHASES:
        for structure in ("projected", "unprojected"):
            initial_actor = actors[structure]
            actor, record = _train_phase(
                equation,
                structure,
                phase,
                projected_q0 if structure == "projected" else unprojected_q0,
                initial_actor,
                device=device,
                force=force,
            )
            actors[structure] = actor
            construction[structure][phase.name] = record

    # Exact references are first constructed after both final checkpoints exist.
    final_checkpoints = {
        structure: _phase_checkpoint(d, k, structure, PHASES[-1])
        for structure in ("projected", "unprojected")
    }
    for structure, path in final_checkpoints.items():
        if not path.is_file():
            raise RuntimeError(f"{structure} final checkpoint was not frozen")
    points = fixed_test_points(d, k, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    rows: list[dict[str, Any]] = []
    prediction_arrays: dict[str, np.ndarray] = {}
    for structure in ("projected", "unprojected"):
        row, arrays = _evaluate_actor(
            structure,
            actors[structure],
            final_checkpoints[structure],
            equation,
            times,
            points,
            exact_u,
            exact_z,
            device=device,
        )
        rows.append(row)
        prediction_arrays.update(arrays)
        print(
            json.dumps(
                {
                    "cell": [d, k],
                    "actor": structure,
                    "corrected_value_rmse": row["corrected_value_rmse"],
                    "corrected_value_ci95": row["corrected_value_ci95"],
                    "mean_analytic_policy_gap": row[
                        "mean_analytic_policy_gap"
                    ],
                    "direct_actor_gradient_rmse": row[
                        "direct_actor_gradient_rmse"
                    ],
                }
            ),
            flush=True,
        )
    _atomic_npz(
        PREDICTION_DIR / f"paired_predictions_d{d}_k{k}.npz",
        **prediction_arrays,
    )
    row_lookup = {row["actor"]: row for row in rows}
    comparison = _comparison(row_lookup)
    headline = {
        "cell": {"d": d, "k": k},
        "projected_corrected_value_rmse": row_lookup["projected"][
            "corrected_value_rmse"
        ],
        "unprojected_corrected_value_rmse": row_lookup["unprojected"][
            "corrected_value_rmse"
        ],
        "unprojected_relative_value_rmse_change": comparison[
            "unprojected_relative_value_rmse_change"
        ],
        "unprojected_beats_projected": comparison[
            "unprojected_beats_projected_on_value_rmse"
        ],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "device": str(device),
        "protocol": protocol,
        "shared_q0_labels": shared_q0_labels,
        "paired_q0_initialization": initialization,
        "construction": construction,
        "rows": rows,
        "comparison": comparison,
        "headline": headline,
        "test_point_sha256": array_sha256(points),
        "elapsed_seconds": time.perf_counter() - started,
        "failed_runs": 0,
        "nan_failures": 0,
        "memory_failures": 0,
        "q2_q3_or_q4_run": False,
        "limitations": [
            "Each structure uses one training seed.",
            (
                "Evaluation confidence intervals condition on the frozen "
                "trained Actors and exclude training-seed variability."
            ),
            (
                "The unprojected arm removes the Actor output prior, while the "
                "benchmark PDE still has its declared effective dimension k."
            ),
        ],
    }
    atomic_json(result_path, payload)
    print(json.dumps(headline, indent=2), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d", type=int)
    parser.add_argument("--k", type=int)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--describe-protocol", action="store_true")
    args = parser.parse_args()
    if args.describe_protocol:
        print(
            json.dumps(
                {f"d{d}_k{k}": protocol_record(d, k) for d, k in CELLS},
                indent=2,
                sort_keys=True,
            )
        )
        return
    if args.all:
        selected = CELLS
    elif args.d is not None and args.k is not None:
        selected = ((args.d, args.k),)
    else:
        parser.error("pass --all or both --d and --k")
    for d, k in selected:
        try:
            run_cell(d, k, force=args.force)
        except Exception as error:
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            atomic_json(
                OUTPUT_DIR / f"failure_d{d}_k{k}.json",
                {
                    "status": "failed",
                    "d": d,
                    "k": k,
                    "error": repr(error),
                    "traceback": traceback.format_exc(),
                },
            )
            raise


if __name__ == "__main__":
    main()
