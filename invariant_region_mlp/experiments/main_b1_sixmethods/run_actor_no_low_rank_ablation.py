"""Focused value-only ablation of the actor's rank-k output projection.

This runner continues the protocol-v2 ambient-output actor from its level-5
checkpoint through the confirmation schedule selected by protocol v3.  It
uses the cached, unprojected level-5 labels and evaluates only the Stage-B
value error.  No analytic gradient reference is constructed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .actor import ActorNet
from .actor_convergence import (
    ACTOR_DEPTH,
    ACTOR_WIDTH,
    BUDGET_LEVELS,
    generate_label_ensemble,
    train_budget_level,
)
from .actor_mlp import evaluate_actor_mlp
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
from .neural_common import clone_state_dict, count_parameters


CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence"
CACHE_DIR = CONVERGENCE_DIR / "cache"
PROJECTED_DIR = CONVERGENCE_DIR / "protocol_v3"
ABLATION_DIR = CONVERGENCE_DIR / "no_low_rank_ablation"
CHECKPOINT_DIR = ABLATION_DIR / "checkpoints"
PREDICTION_DIR = ABLATION_DIR / "predictions"

START_LEVEL = "level_5"
CONFIRMATION_LEVELS = (
    "level_6_confirmation",
    "level_7_confirmation",
    "level_8_confirmation",
)


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


def _load_plain_actor(
    checkpoint: Path, equation: B1Equation, device: torch.device
) -> tuple[ActorNet, dict[str, Any]]:
    payload = torch.load(checkpoint, map_location=device, weights_only=False)
    state_dict = payload["state_dict"]
    if "label_basis" in state_dict:
        raise RuntimeError(
            f"checkpoint unexpectedly contains a label basis: {checkpoint}"
        )
    architecture = payload.get("architecture", {})
    if int(architecture.get("d", -1)) != equation.d:
        raise RuntimeError(f"checkpoint dimension mismatch: {checkpoint}")
    if int(architecture.get("width", -1)) != ACTOR_WIDTH:
        raise RuntimeError(f"checkpoint width mismatch: {checkpoint}")
    if int(architecture.get("depth", -1)) != ACTOR_DEPTH:
        raise RuntimeError(f"checkpoint depth mismatch: {checkpoint}")
    model = ActorNet(equation.d, width=ACTOR_WIDTH, depth=ACTOR_DEPTH).to(device)
    model.load_state_dict(state_dict)
    model.eval()
    return model, payload["metadata"]


def _plain_checkpoint_payload(
    model: ActorNet, equation: B1Equation, metadata: dict[str, Any]
) -> dict[str, Any]:
    return {
        "state_dict": clone_state_dict(model),
        "architecture": {
            "kind": "ambient_actor_without_output_projection",
            "d": equation.d,
            "width": ACTOR_WIDTH,
            "depth": ACTOR_DEPTH,
            "activation": "silu",
            "output_projection_rank": None,
        },
        "metadata": metadata,
    }


def _raw_level_five_data(
    equation: B1Equation, *, device: torch.device
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[str, Any],
]:
    level_index = next(
        index for index, level in enumerate(BUDGET_LEVELS) if level.name == START_LEVEL
    )
    level = BUDGET_LEVELS[level_index]
    seed = stable_seed("B1", "actor_convergence", equation.d, equation.k)
    train_t, train_x = training_points(
        equation.d,
        equation.k,
        f"actor_convergence_{START_LEVEL}_train",
        seed,
        level.training_points,
    )
    valid_t, valid_x = training_points(
        equation.d,
        equation.k,
        f"actor_convergence_{START_LEVEL}_valid",
        seed,
        level.validation_points,
    )
    cell_cache = CACHE_DIR / f"d{equation.d}_k{equation.k}"
    train_cache = cell_cache / f"{START_LEVEL}_train.npz"
    valid_cache = cell_cache / f"{START_LEVEL}_valid.npz"
    if not train_cache.is_file() or not valid_cache.is_file():
        raise FileNotFoundError(
            "focused ablation requires the existing level-5 raw-label caches"
        )
    train_y, train_metadata = generate_label_ensemble(
        equation,
        train_t,
        train_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=100_000 + 1000 * level_index,
        seed_namespace=f"actor_convergence_{START_LEVEL}_train",
        device=device,
        cache_path=train_cache,
    )
    valid_y, valid_metadata = generate_label_ensemble(
        equation,
        valid_t,
        valid_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=200_000 + 1000 * level_index,
        seed_namespace=f"actor_convergence_{START_LEVEL}_valid",
        device=device,
        cache_path=valid_cache,
    )
    metadata = {
        "source_level": START_LEVEL,
        "labels_projected": False,
        "train": train_metadata,
        "validation": valid_metadata,
    }
    return train_t, train_x, train_y, valid_t, valid_x, valid_y, metadata


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    mean = float(np.mean(array))
    radius = (
        1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array))
        if len(array) > 1
        else 0.0
    )
    return mean, radius


def _comparison_references(d: int, k: int) -> dict[str, Any]:
    projected_path = PROJECTED_DIR / "raw" / f"actor_converged_d{d}_k{k}.json"
    comparison_path = PROJECTED_DIR / "baseline_comparison" / "comparison.csv"
    projected = json.loads(projected_path.read_text())
    projected_rows = {
        int(row["repetition"]): float(row["value_rmse"])
        for row in projected["rows"]
    }
    with comparison_path.open(newline="") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["d"]) == d and int(row["k"]) == k
        ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one baseline-comparison row for d={d}, k={k}")
    baseline = matches[0]
    return {
        "projected_selected_level": projected["selected_level"],
        "projected_value_rmse_by_repetition": projected_rows,
        "projected_value_rmse": float(baseline["new_actor_value_rmse"]),
        "projected_value_ci95": float(baseline["new_actor_value_ci95"]),
        "best_previous_method": baseline["best_previous_value_method"],
        "best_previous_value_rmse": float(baseline["best_previous_value_rmse"]),
    }


def _evaluate_value(
    equation: B1Equation,
    model: ActorNet,
    *,
    device: torch.device,
) -> tuple[list[dict[str, Any]], np.ndarray, str]:
    stage = load_config()["stage_b"]
    points = fixed_test_points(equation.d, equation.k, int(stage["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    reference_u = equation.exact_u(times, points)
    values = np.empty((int(stage["mlp_seeds"]), len(points)), dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for repetition in range(int(stage["mlp_seeds"])):
        state, cost = evaluate_actor_mlp(
            equation,
            model,
            times,
            points,
            n=int(stage["mlp_n"]),
            M=int(stage["mlp_M"]),
            repetition=repetition,
            device=device,
        )
        prediction = state[:, 0]
        values[repetition] = prediction
        error = prediction - reference_u
        rows.append(
            {
                "repetition": repetition,
                "value_rmse": float(np.sqrt(np.mean(error**2))),
                "value_bias": float(np.mean(error)),
                "inference_seconds": float(cost["wall_seconds"]),
                "peak_gpu_bytes": cost["peak_gpu_bytes"],
                "terminal_samples": int(cost["terminal_samples"]),
                "transition_samples": int(cost["transition_samples"]),
                "actor_calls": int(cost["actor_calls"]),
                "actor_points": int(cost["actor_points"]),
            }
        )
    return rows, values, array_sha256(points)


def run(d: int, k: int, *, force: bool = False) -> dict[str, Any]:
    result_path = ABLATION_DIR / f"result_d{d}_k{k}.json"
    if result_path.is_file() and not force:
        existing = json.loads(result_path.read_text())
        if existing.get("status") == "complete":
            print(json.dumps(existing["value_summary"], indent=2), flush=True)
            return existing

    config = load_config()
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("no-low-rank actor ablation requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    started = time.perf_counter()

    start_checkpoint = (
        CONVERGENCE_DIR / "checkpoints" / f"actor_d{d}_k{k}_{START_LEVEL}.pt"
    )
    if not start_checkpoint.is_file():
        raise FileNotFoundError(
            f"missing unprojected starting checkpoint: {start_checkpoint}"
        )
    model, start_metadata = _load_plain_actor(start_checkpoint, equation, device)
    (
        train_t,
        train_x,
        train_y,
        valid_t,
        valid_x,
        valid_y,
        label_metadata,
    ) = _raw_level_five_data(equation, device=device)

    records: list[dict[str, Any]] = []
    for level_name in CONFIRMATION_LEVELS:
        level = next(level for level in BUDGET_LEVELS if level.name == level_name)
        training = train_budget_level(
            model,
            equation,
            level,
            train_t,
            train_x,
            train_y,
            valid_t,
            valid_x,
            valid_y,
            device=device,
            seed=stable_seed(
                "B1", "actor_convergence_optimizer", d, k, level.name
            ),
        )
        record = {
            "status": "complete",
            "d": d,
            "k": k,
            "budget": {
                "name": level.name,
                "n": level.n,
                "M": level.M,
                "replicas": level.replicas,
                "training_points": level.training_points,
                "validation_points": level.validation_points,
                "max_steps": level.max_steps,
                "minimum_steps": level.minimum_steps,
                "learning_rate": level.learning_rate,
                "reuse_labels_from": level.reuse_labels_from,
            },
            "training": training,
            "labels_projected": False,
        }
        records.append(record)
        checkpoint_path = CHECKPOINT_DIR / f"actor_d{d}_k{k}_{level.name}.pt"
        _atomic_torch_save(
            checkpoint_path,
            _plain_checkpoint_payload(model, equation, record),
        )
        print(
            json.dumps(
                {
                    "level": level.name,
                    "best_step": training["best_step"],
                    "steps_executed": training["steps_executed"],
                    "validation_label_rmse": training[
                        "selected_validation_label_rmse"
                    ],
                    "training_seconds": training["training_seconds"],
                }
            ),
            flush=True,
        )

    rows, values, test_hash = _evaluate_value(equation, model, device=device)
    _atomic_npz(
        PREDICTION_DIR / f"value_predictions_d{d}_k{k}.npz",
        value=values,
        repetition=np.arange(len(rows), dtype=np.int64),
    )
    references = _comparison_references(d, k)
    if references["projected_selected_level"] != CONFIRMATION_LEVELS[-1]:
        raise RuntimeError(
            "focused ablation target does not match the projected actor's "
            "selected level"
        )
    rmse = [float(row["value_rmse"]) for row in rows]
    mean_rmse, ci95 = _mean_ci95(rmse)
    projected_by_seed = references["projected_value_rmse_by_repetition"]
    paired_differences = [
        row["value_rmse"] - projected_by_seed[int(row["repetition"])]
        for row in rows
    ]
    paired_mean, paired_ci95 = _mean_ci95(paired_differences)
    summary = {
        "no_low_rank_value_rmse": mean_rmse,
        "no_low_rank_value_ci95": ci95,
        "projected_actor_value_rmse": references["projected_value_rmse"],
        "projected_actor_value_ci95": references["projected_value_ci95"],
        "no_low_rank_over_projected": (
            mean_rmse / references["projected_value_rmse"]
        ),
        "paired_no_low_rank_minus_projected": paired_mean,
        "paired_difference_ci95": paired_ci95,
        "best_previous_method": references["best_previous_method"],
        "best_previous_value_rmse": references["best_previous_value_rmse"],
        "no_low_rank_over_best_previous": (
            mean_rmse / references["best_previous_value_rmse"]
        ),
        "beats_best_previous": mean_rmse < references["best_previous_value_rmse"],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "actor output-projection ablation",
        "d": d,
        "k": k,
        "device": str(device),
        "architecture": {
            "input": "ambient (t,x)",
            "hidden_width": ACTOR_WIDTH,
            "hidden_layers": ACTOR_DEPTH,
            "activation": "silu",
            "output": "ambient z",
            "output_projection": None,
            "parameter_count": count_parameters(model),
        },
        "start_checkpoint": str(start_checkpoint),
        "start_level_metadata": start_metadata,
        "confirmation_levels": records,
        "labels": label_metadata,
        "selected_level": CONFIRMATION_LEVELS[-1],
        "evaluation": {
            "metric_scope": "value only",
            "test_points": int(config["stage_b"]["test_points"]),
            "mlp_seeds": int(config["stage_b"]["mlp_seeds"]),
            "mlp_n": int(config["stage_b"]["mlp_n"]),
            "mlp_M": int(config["stage_b"]["mlp_M"]),
            "paired_seed_namespace": "B1|actor_mlp",
            "test_point_sha256": test_hash,
            "gradient_reference_evaluated": False,
        },
        "rows": rows,
        "value_summary": summary,
        "elapsed_seconds": time.perf_counter() - started,
        "nonterminal_oracle_reference_used_for_training": False,
        "nonterminal_oracle_reference_used_for_selection": False,
    }
    atomic_json(result_path, payload)
    print(json.dumps(summary, indent=2), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d", type=int, default=400)
    parser.add_argument("--k", type=int, default=40)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    run(args.d, args.k, force=args.force)


if __name__ == "__main__":
    main()
