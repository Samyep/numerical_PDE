"""Progressive, non-oracle actor-label and optimization convergence study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np
import torch

from .actor import ActorNet
from .common import array_sha256, atomic_json, seed_everything, stable_seed
from .equation import B1Equation
from .mlp_adapter import evaluate_mlp
from .neural_common import clone_state_dict, count_parameters


ACTOR_WIDTH = 384
ACTOR_DEPTH = 4
ACTOR_BATCH_SIZE = 1024
ACTOR_LEARNING_RATE = 5.0e-4
ACTOR_MIN_LEARNING_RATE = 1.0e-6
ACTOR_WEIGHT_DECAY = 1.0e-5
VALIDATION_EVERY = 100
SCHEDULER_PATIENCE = 6
STOPPING_PATIENCE = 8
PARAMETER_UPDATE_TOLERANCE = 5.0e-3

AUDIT_POINTS = 2048
AUDIT_N = 3
AUDIT_M = 16
AUDIT_BANK_REPLICAS = 4

AUDIT_STABILITY_TOLERANCE = 0.05
PREDICTION_STABILITY_TOLERANCE = 0.05
NOISE_FLOOR_FACTOR = 1.10


@dataclass(frozen=True)
class BudgetLevel:
    name: str
    n: int
    M: int
    replicas: int
    training_points: int
    validation_points: int
    max_steps: int
    minimum_steps: int


BUDGET_LEVELS = (
    BudgetLevel("level_1", 2, 16, 4, 8192, 2048, 6000, 1500),
    BudgetLevel("level_2", 3, 8, 4, 12288, 3072, 10000, 2000),
    BudgetLevel("level_3", 3, 12, 6, 16384, 4096, 16000, 2500),
    BudgetLevel("level_4", 3, 16, 8, 24576, 6144, 24000, 3000),
    BudgetLevel("level_5", 3, 20, 10, 32768, 8192, 32000, 4000),
)


def budget_protocol() -> dict[str, Any]:
    return {
        "architecture": {
            "hidden_width": ACTOR_WIDTH,
            "hidden_layers": ACTOR_DEPTH,
            "activation": "silu",
            "input": "ambient (t,x)",
            "output": "ambient z",
        },
        "optimizer": {
            "name": "AdamW",
            "batch_size": ACTOR_BATCH_SIZE,
            "learning_rate": ACTOR_LEARNING_RATE,
            "minimum_learning_rate": ACTOR_MIN_LEARNING_RATE,
            "weight_decay": ACTOR_WEIGHT_DECAY,
            "validation_every": VALIDATION_EVERY,
            "scheduler_patience": SCHEDULER_PATIENCE,
            "stopping_patience": STOPPING_PATIENCE,
            "parameter_update_tolerance": PARAMETER_UPDATE_TOLERANCE,
        },
        "audit": {
            "points": AUDIT_POINTS,
            "n": AUDIT_N,
            "M": AUDIT_M,
            "replicas_per_bank": AUDIT_BANK_REPLICAS,
            "banks": 2,
        },
        "convergence": {
            "minimum_budget_levels": 2,
            "consecutive_stable_transitions": 2,
            "audit_excess_rmse_relative_change_tolerance": AUDIT_STABILITY_TOLERANCE,
            "prediction_relative_change_tolerance": PREDICTION_STABILITY_TOLERANCE,
            "noise_floor_factor": NOISE_FLOOR_FACTOR,
        },
        "levels": [asdict(level) for level in BUDGET_LEVELS],
        "nonterminal_oracle_reference_used_for_training": False,
        "nonterminal_oracle_reference_used_for_selection": False,
    }


def make_actor(equation: B1Equation, device: torch.device) -> ActorNet:
    return ActorNet(equation.d, width=ACTOR_WIDTH, depth=ACTOR_DEPTH).to(device)


def label_chunk_size(d: int, n: int, M: int) -> int:
    """Bound the largest root normal array while retaining useful GPU batches."""

    coordinate_budget = 500_000_000
    coordinates_per_point = max(1, int(d) * int(M) ** int(n))
    return max(8, min(512, coordinate_budget // coordinates_per_point))


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez(temporary, **arrays)
    temporary.replace(path)


def generate_label_ensemble(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    replicas: int,
    repetition_offset: int,
    seed_namespace: str,
    device: torch.device,
    cache_path: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Average independent pathwise-baseline labels and cache the result."""

    point_hash = array_sha256(np.asarray(t), np.asarray(x))
    metadata_path = cache_path.with_suffix(".json")
    if cache_path.is_file() and metadata_path.is_file():
        with np.load(cache_path, allow_pickle=False) as payload:
            observed_hash = str(payload["point_hash"].item())
            labels = np.asarray(payload["labels"], dtype=np.float64)
        import json

        metadata = json.loads(metadata_path.read_text())
        if observed_hash != point_hash:
            raise RuntimeError(f"label-cache point hash mismatch: {cache_path}")
        expected = (len(t), equation.d)
        if labels.shape != expected or not np.all(np.isfinite(labels)):
            raise RuntimeError(f"invalid label cache {cache_path}: {labels.shape}")
        return labels, metadata

    if replicas < 1:
        raise ValueError("replicas must be positive")
    total = np.zeros((len(t), equation.d), dtype=np.float64)
    total_squared_norm = 0.0
    aggregate_counts = {
        "terminal_samples": 0,
        "transition_samples": 0,
        "generator_calls": 0,
        "recursively_evaluated_states": 0,
        "recursive_calls": 0,
        "gaussian_coordinate_work": 0,
    }
    replicate_seconds: list[float] = []
    peak_gpu_bytes = 0
    chunk = label_chunk_size(equation.d, n, M)
    started = time.perf_counter()
    for replica in range(replicas):
        state, cost = evaluate_mlp(
            equation,
            t,
            x,
            n=n,
            M=M,
            method="actor_label",
            repetition=repetition_offset + replica,
            device=device,
            chunk_size=chunk,
            seed_namespace=seed_namespace,
        )
        labels = state[:, 1:]
        if not np.all(np.isfinite(labels)):
            raise FloatingPointError(
                f"nonfinite labels for {seed_namespace}, replica {replica}"
            )
        total += labels
        total_squared_norm += float(np.sum(labels * labels, dtype=np.float64))
        replicate_seconds.append(float(cost["wall_seconds"]))
        peak_gpu_bytes = max(peak_gpu_bytes, int(cost["peak_gpu_bytes"] or 0))
        for key in aggregate_counts:
            aggregate_counts[key] += int(cost[key])

    mean = total / replicas
    if replicas > 1:
        centered_sum_squares = max(
            0.0,
            total_squared_norm
            - replicas * float(np.sum(mean * mean, dtype=np.float64)),
        )
        sample_variance_sum = centered_sum_squares / (replicas - 1)
        spread_rmse = math.sqrt(sample_variance_sum / len(t))
        standard_error_rmse = math.sqrt(sample_variance_sum / (len(t) * replicas))
    else:
        spread_rmse = float("nan")
        standard_error_rmse = float("nan")

    metadata: dict[str, Any] = {
        "status": "complete",
        "d": equation.d,
        "k": equation.k,
        "points": len(t),
        "point_sha256": point_hash,
        "n": n,
        "M": M,
        "replicas": replicas,
        "repetition_offset": repetition_offset,
        "seed_namespace": seed_namespace,
        "chunk_size": chunk,
        "label_method": "pathwise terminal + level-matched baseline",
        "ensemble_spread_vector_rmse": spread_rmse,
        "ensemble_standard_error_vector_rmse": standard_error_rmse,
        "replicate_wall_seconds": replicate_seconds,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak_gpu_bytes,
        **aggregate_counts,
        "nonterminal_oracle_reference_used": False,
    }
    _atomic_npz(
        cache_path,
        labels=mean,
        point_hash=np.asarray(point_hash),
    )
    atomic_json(metadata_path, metadata)
    return mean, metadata


def predict_actor(
    model: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    device: torch.device,
    batch_size: int = 4096,
) -> np.ndarray:
    pieces: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(t), batch_size):
            stop = min(start + batch_size, len(t))
            prediction = model.predict(
                torch.as_tensor(t[start:stop], dtype=torch.float32, device=device),
                torch.as_tensor(x[start:stop], dtype=torch.float32, device=device),
            )
            pieces.append(prediction.cpu().numpy().astype(np.float64))
    return np.concatenate(pieces, axis=0)


def _validation_loss(
    model: ActorNet,
    t: torch.Tensor,
    x: torch.Tensor,
    y: torch.Tensor,
    *,
    batch_size: int = 4096,
) -> float:
    total = 0.0
    model.eval()
    with torch.no_grad():
        for start in range(0, len(t), batch_size):
            stop = min(start + batch_size, len(t))
            error = model.predict(t[start:stop], x[start:stop]) - y[start:stop]
            total += float(torch.sum(error.square()).item())
    return total / len(t)


def _relative_parameter_distance(
    first: dict[str, torch.Tensor], second: dict[str, torch.Tensor]
) -> float:
    numerator = 0.0
    denominator = 0.0
    for key in first:
        left = first[key].to(dtype=torch.float64)
        right = second[key].to(dtype=torch.float64)
        numerator += float(torch.sum((left - right).square()).item())
        denominator += float(torch.sum(right.square()).item())
    return math.sqrt(numerator / max(denominator, np.finfo(np.float64).tiny))


def train_budget_level(
    model: ActorNet,
    equation: B1Equation,
    level: BudgetLevel,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_y: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_y: np.ndarray,
    *,
    device: torch.device,
    seed: int,
) -> dict[str, Any]:
    """Fit one level until its stochastic validation optimization plateaus."""

    seed_everything(seed)
    start_state = clone_state_dict(model)
    tensors = {
        "train_t": torch.as_tensor(train_t, dtype=torch.float32, device=device),
        "train_x": torch.as_tensor(train_x, dtype=torch.float32, device=device),
        "train_y": torch.as_tensor(train_y, dtype=torch.float32, device=device),
        "valid_t": torch.as_tensor(valid_t, dtype=torch.float32, device=device),
        "valid_x": torch.as_tensor(valid_x, dtype=torch.float32, device=device),
        "valid_y": torch.as_tensor(valid_y, dtype=torch.float32, device=device),
    }
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=ACTOR_LEARNING_RATE,
        weight_decay=ACTOR_WEIGHT_DECAY,
    )
    generator = torch.Generator(device=device).manual_seed(seed)
    best_loss = float("inf")
    best_step = 0
    best_state: dict[str, torch.Tensor] | None = None
    bad_validations = 0
    reductions = 0
    history: list[dict[str, float | int]] = []
    last_eval = [parameter.detach().clone() for parameter in model.parameters()]
    stopping_reason = "maximum_steps"
    optimizer_converged = False

    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    model.train()
    for step in range(1, level.max_steps + 1):
        index = torch.randint(
            len(train_t),
            (min(ACTOR_BATCH_SIZE, len(train_t)),),
            generator=generator,
            device=device,
        )
        prediction = model.predict(
            tensors["train_t"][index], tensors["train_x"][index]
        )
        vector_loss = torch.sum(
            (prediction - tensors["train_y"][index]).square(), dim=-1
        ).mean()
        loss = vector_loss / equation.d
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
        optimizer.step()

        if step % VALIDATION_EVERY != 0 and step != level.max_steps:
            continue

        valid_loss = _validation_loss(
            model,
            tensors["valid_t"],
            tensors["valid_x"],
            tensors["valid_y"],
        )
        numerator = torch.zeros((), dtype=torch.float64, device=device)
        denominator = torch.zeros((), dtype=torch.float64, device=device)
        with torch.no_grad():
            for parameter, previous in zip(model.parameters(), last_eval, strict=True):
                numerator += torch.sum(
                    (parameter.detach().to(torch.float64) - previous.to(torch.float64)).square()
                )
                denominator += torch.sum(previous.to(torch.float64).square())
                previous.copy_(parameter.detach())
        relative_update = math.sqrt(
            float(numerator.item())
            / max(float(denominator.item()), np.finfo(np.float64).tiny)
        )
        learning_rate = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "step": step,
                "training_vector_mse": float(vector_loss.item()),
                "validation_vector_mse": valid_loss,
                "learning_rate": learning_rate,
                "relative_parameter_update": relative_update,
            }
        )
        if valid_loss < best_loss * (1.0 - 1.0e-4):
            best_loss = valid_loss
            best_step = step
            best_state = clone_state_dict(model)
            bad_validations = 0
        else:
            bad_validations += 1

        if bad_validations >= SCHEDULER_PATIENCE:
            old_rate = float(optimizer.param_groups[0]["lr"])
            new_rate = max(ACTOR_MIN_LEARNING_RATE, old_rate * 0.25)
            if new_rate < old_rate:
                for group in optimizer.param_groups:
                    group["lr"] = new_rate
                reductions += 1
                bad_validations = 0

        at_minimum_rate = (
            float(optimizer.param_groups[0]["lr"])
            <= ACTOR_MIN_LEARNING_RATE * (1.0 + 1.0e-12)
        )
        if (
            step >= level.minimum_steps
            and at_minimum_rate
            and bad_validations >= STOPPING_PATIENCE
            and relative_update <= PARAMETER_UPDATE_TOLERANCE
        ):
            stopping_reason = "validation_and_parameter_plateau"
            optimizer_converged = True
            break
        model.train()

    if best_state is None:
        raise RuntimeError(f"{level.name} produced no actor checkpoint")
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak_gpu_bytes: int | None = int(torch.cuda.max_memory_allocated(device))
    else:
        peak_gpu_bytes = None
    training_seconds = time.perf_counter() - started
    model.load_state_dict(best_state)
    model.eval()

    final_train_loss = _validation_loss(
        model, tensors["train_t"], tensors["train_x"], tensors["train_y"]
    )
    final_valid_loss = _validation_loss(
        model, tensors["valid_t"], tensors["valid_x"], tensors["valid_y"]
    )
    return {
        "budget_level": level.name,
        "d": equation.d,
        "k": equation.k,
        "seed": seed,
        "parameter_count": count_parameters(model),
        "best_step": best_step,
        "steps_executed": step,
        "best_validation_vector_mse": best_loss,
        "selected_training_vector_mse": final_train_loss,
        "selected_validation_vector_mse": final_valid_loss,
        "selected_validation_label_rmse": math.sqrt(final_valid_loss),
        "learning_rate_reductions": reductions,
        "final_learning_rate": float(optimizer.param_groups[0]["lr"]),
        "optimizer_converged": optimizer_converged,
        "stopping_reason": stopping_reason,
        "relative_parameter_change_from_level_start": _relative_parameter_distance(
            best_state, start_state
        ),
        "training_seconds": training_seconds,
        "peak_gpu_bytes": peak_gpu_bytes,
        "history": history,
        "nonterminal_oracle_reference_used_for_training": False,
        "nonterminal_oracle_reference_used_for_checkpoint_selection": False,
    }


def audit_record(
    prediction: np.ndarray,
    audit_bank_a: np.ndarray,
    audit_bank_b: np.ndarray,
    *,
    previous_prediction: np.ndarray | None,
) -> dict[str, float | None]:
    consensus = 0.5 * (audit_bank_a + audit_bank_b)
    target_mse = float(np.mean(np.sum((prediction - consensus) ** 2, axis=1)))
    # If A and B are independent averages with equal variance, the consensus
    # variance is one quarter of E||A-B||^2.
    noise_floor_mse = 0.25 * float(
        np.mean(np.sum((audit_bank_a - audit_bank_b) ** 2, axis=1))
    )
    excess_mse = max(0.0, target_mse - noise_floor_mse)
    if previous_prediction is None:
        prediction_change = None
    else:
        change = math.sqrt(
            float(np.mean(np.sum((prediction - previous_prediction) ** 2, axis=1)))
        )
        scale = math.sqrt(float(np.mean(np.sum(prediction**2, axis=1))))
        prediction_change = change / max(scale, np.finfo(np.float64).tiny)
    return {
        "audit_target_rmse": math.sqrt(target_mse),
        "audit_consensus_noise_floor_rmse": math.sqrt(noise_floor_mse),
        "audit_noise_corrected_rmse": math.sqrt(excess_mse),
        "prediction_relative_change": prediction_change,
    }


def convergence_decision(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Apply the frozen label-only rule; oracle diagnostic keys are ignored."""

    items = list(records)
    if len(items) < 2:
        return {
            "converged": False,
            "reason": "minimum_budget_levels_not_reached",
            "stable_transitions": 0,
        }
    latest = items[-1]
    if (
        bool(latest["optimizer_converged"])
        and float(latest["audit_target_rmse"])
        <= NOISE_FLOOR_FACTOR * float(latest["audit_consensus_noise_floor_rmse"])
    ):
        return {
            "converged": True,
            "reason": "stochastic_audit_noise_floor_reached",
            "stable_transitions": 0,
        }

    stable: list[bool] = []
    changes: list[float] = []
    for previous, current in zip(items[:-1], items[1:], strict=True):
        old = float(previous["audit_noise_corrected_rmse"])
        new = float(current["audit_noise_corrected_rmse"])
        relative_change = abs(new - old) / max(old, np.finfo(np.float64).tiny)
        changes.append(relative_change)
        prediction_change = current.get("prediction_relative_change")
        stable.append(
            bool(current["optimizer_converged"])
            and relative_change <= AUDIT_STABILITY_TOLERANCE
            and prediction_change is not None
            and float(prediction_change) <= PREDICTION_STABILITY_TOLERANCE
        )
    consecutive = 0
    for value in reversed(stable):
        if not value:
            break
        consecutive += 1
    return {
        "converged": consecutive >= 2,
        "reason": (
            "two_consecutive_budget_plateaus"
            if consecutive >= 2
            else "additional_budget_required"
        ),
        "stable_transitions": consecutive,
        "audit_relative_changes": changes,
    }


def checkpoint_payload(
    model: ActorNet, metadata: dict[str, Any]
) -> dict[str, Any]:
    return {
        "state_dict": clone_state_dict(model),
        "architecture": {
            "d": model.d,
            "width": ACTOR_WIDTH,
            "depth": ACTOR_DEPTH,
            "activation": "silu",
        },
        "metadata": metadata,
    }


def load_converged_actor(path: Path, device: torch.device) -> tuple[ActorNet, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = ActorNet(
        int(architecture["d"]),
        width=int(architecture["width"]),
        depth=int(architecture["depth"]),
    ).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload["metadata"]

