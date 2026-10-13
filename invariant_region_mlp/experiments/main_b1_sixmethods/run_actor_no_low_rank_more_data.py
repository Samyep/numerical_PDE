"""Value-only no-projection actor ablation with twice as many labeled states."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .actor_convergence import (
    ACTOR_DEPTH,
    ACTOR_WIDTH,
    BUDGET_LEVELS,
    BudgetLevel,
    generate_label_ensemble,
    train_budget_level,
)
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    load_config,
    stable_seed,
    training_points,
)
from .equation import B1Equation
from .neural_common import count_parameters
from .run_actor_no_low_rank_ablation import (
    _atomic_npz,
    _atomic_torch_save,
    _comparison_references,
    _evaluate_value,
    _load_plain_actor,
    _mean_ci95,
    _plain_checkpoint_payload,
    _raw_level_five_data,
)


CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence"
MORE_DATA_DIR = CONVERGENCE_DIR / "no_low_rank_more_data"
CACHE_DIR = MORE_DATA_DIR / "cache"
CHECKPOINT_DIR = MORE_DATA_DIR / "checkpoints"
PREDICTION_DIR = MORE_DATA_DIR / "predictions"
START_LEVEL = "level_5"


MORE_DATA_LEVELS = (
    (
        "level_6_confirmation",
        BudgetLevel(
            "more_data_level_6",
            3,
            20,
            10,
            65_536,
            16_384,
            12_000,
            4_000,
            1.0e-4,
            "level_5_plus_extra",
        ),
    ),
    (
        "level_7_confirmation",
        BudgetLevel(
            "more_data_level_7",
            3,
            20,
            10,
            65_536,
            16_384,
            8_000,
            3_000,
            2.5e-5,
            "level_5_plus_extra",
        ),
    ),
    (
        "level_8_confirmation",
        BudgetLevel(
            "more_data_level_8",
            3,
            20,
            10,
            65_536,
            16_384,
            6_000,
            2_000,
            5.0e-6,
            "level_5_plus_extra",
        ),
    ),
)


def _extra_raw_data(
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
    seed = stable_seed("B1", "actor_no_low_rank_more_data", equation.d, equation.k)
    train_t, train_x = training_points(
        equation.d,
        equation.k,
        "actor_no_low_rank_more_data_extra_train",
        seed,
        level.training_points,
    )
    valid_t, valid_x = training_points(
        equation.d,
        equation.k,
        "actor_no_low_rank_more_data_extra_valid",
        seed,
        level.validation_points,
    )
    cell_cache = CACHE_DIR / f"d{equation.d}_k{equation.k}"
    train_y, train_metadata = generate_label_ensemble(
        equation,
        train_t,
        train_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=600_000,
        seed_namespace="actor_no_low_rank_more_data_extra_train",
        device=device,
        cache_path=cell_cache / "extra_train.npz",
    )
    valid_y, valid_metadata = generate_label_ensemble(
        equation,
        valid_t,
        valid_x,
        n=level.n,
        M=level.M,
        replicas=level.replicas,
        repetition_offset=700_000,
        seed_namespace="actor_no_low_rank_more_data_extra_valid",
        device=device,
        cache_path=cell_cache / "extra_valid.npz",
    )
    metadata = {
        "train": train_metadata,
        "validation": valid_metadata,
        "train_point_sha256": array_sha256(train_t, train_x),
        "validation_point_sha256": array_sha256(valid_t, valid_x),
        "labels_projected": False,
    }
    return train_t, train_x, train_y, valid_t, valid_x, valid_y, metadata


def _paired_summary(
    rows: list[dict[str, Any]], reference_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    reference = {
        int(row["repetition"]): float(row["value_rmse"])
        for row in reference_rows
    }
    differences = [
        float(row["value_rmse"]) - reference[int(row["repetition"])]
        for row in rows
    ]
    mean, ci95 = _mean_ci95(differences)
    return {
        "mean_difference": mean,
        "difference_ci95": ci95,
        "worse_seeds": sum(value > 0.0 for value in differences),
        "total_seeds": len(differences),
    }


def _reference_result(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    rows = payload["rows"]
    test_hash = payload.get("test_point_sha256")
    if test_hash is None:
        test_hash = payload["evaluation"]["test_point_sha256"]
    value_rmse, value_ci95 = _mean_ci95(
        [float(row["value_rmse"]) for row in rows]
    )
    return {
        "value_rmse": value_rmse,
        "value_ci95": value_ci95,
        "test_point_sha256": test_hash,
        "rows": rows,
    }


def run(d: int, k: int, *, force: bool = False) -> dict[str, Any]:
    result_path = MORE_DATA_DIR / f"result_d{d}_k{k}.json"
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
        raise RuntimeError("more-data actor ablation requires an allocated GPU")
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
        base_train_t,
        base_train_x,
        base_train_y,
        base_valid_t,
        base_valid_x,
        base_valid_y,
        base_metadata,
    ) = _raw_level_five_data(equation, device=device)
    (
        extra_train_t,
        extra_train_x,
        extra_train_y,
        extra_valid_t,
        extra_valid_x,
        extra_valid_y,
        extra_metadata,
    ) = _extra_raw_data(equation, device=device)

    train_t = np.concatenate((base_train_t, extra_train_t))
    train_x = np.concatenate((base_train_x, extra_train_x))
    train_y = np.concatenate((base_train_y, extra_train_y))
    valid_t = np.concatenate((base_valid_t, extra_valid_t))
    valid_x = np.concatenate((base_valid_x, extra_valid_x))
    valid_y = np.concatenate((base_valid_y, extra_valid_y))
    expected_train = MORE_DATA_LEVELS[0][1].training_points
    expected_valid = MORE_DATA_LEVELS[0][1].validation_points
    if len(train_t) != expected_train or len(valid_t) != expected_valid:
        raise RuntimeError("combined more-data arrays have unexpected sizes")
    if array_sha256(base_train_t, base_train_x) == array_sha256(
        extra_train_t, extra_train_x
    ):
        raise RuntimeError("extra training points duplicate the base point array")

    records: list[dict[str, Any]] = []
    for source_name, level in MORE_DATA_LEVELS:
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
                "B1", "actor_convergence_optimizer", d, k, source_name
            ),
        )
        record = {
            "status": "complete",
            "d": d,
            "k": k,
            "source_confirmation_level": source_name,
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

    rows, predictions, test_hash = _evaluate_value(equation, model, device=device)
    _atomic_npz(
        PREDICTION_DIR / f"value_predictions_d{d}_k{k}.npz",
        value=predictions,
        repetition=np.arange(len(rows), dtype=np.int64),
    )
    comparison = _comparison_references(d, k)
    previous = _reference_result(
        CONVERGENCE_DIR / "no_low_rank_ablation" / f"result_d{d}_k{k}.json"
    )
    level_five = _reference_result(
        CONVERGENCE_DIR / "raw" / f"actor_converged_d{d}_k{k}.json"
    )
    if test_hash != previous["test_point_sha256"]:
        raise RuntimeError("more-data and previous ablation test hashes differ")
    if test_hash != level_five["test_point_sha256"]:
        raise RuntimeError("more-data and level-5 test hashes differ")

    rmse = [float(row["value_rmse"]) for row in rows]
    mean_rmse, ci95 = _mean_ci95(rmse)
    projected_rows = [
        {"repetition": repetition, "value_rmse": value}
        for repetition, value in comparison[
            "projected_value_rmse_by_repetition"
        ].items()
    ]
    summary = {
        "more_data_value_rmse": mean_rmse,
        "more_data_value_ci95": ci95,
        "projected_actor_value_rmse": comparison["projected_value_rmse"],
        "more_data_over_projected": mean_rmse
        / comparison["projected_value_rmse"],
        "paired_vs_projected": _paired_summary(rows, projected_rows),
        "previous_no_projection_value_rmse": previous["value_rmse"],
        "more_data_over_previous_no_projection": mean_rmse
        / previous["value_rmse"],
        "paired_vs_previous_no_projection": _paired_summary(
            rows, previous["rows"]
        ),
        "level_five_no_projection_value_rmse": level_five["value_rmse"],
        "more_data_over_level_five": mean_rmse / level_five["value_rmse"],
        "paired_vs_level_five": _paired_summary(rows, level_five["rows"]),
        "best_previous_method": comparison["best_previous_method"],
        "best_previous_value_rmse": comparison["best_previous_value_rmse"],
        "more_data_over_best_previous": mean_rmse
        / comparison["best_previous_value_rmse"],
        "beats_best_previous": mean_rmse
        < comparison["best_previous_value_rmse"],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "no-projection actor with twice as many labeled states",
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
        "data": {
            "base_training_points": len(base_train_t),
            "extra_training_points": len(extra_train_t),
            "total_training_points": len(train_t),
            "base_validation_points": len(base_valid_t),
            "extra_validation_points": len(extra_valid_t),
            "total_validation_points": len(valid_t),
            "n": MORE_DATA_LEVELS[0][1].n,
            "M": MORE_DATA_LEVELS[0][1].M,
            "replicas": MORE_DATA_LEVELS[0][1].replicas,
            "labels_projected": False,
            "base": base_metadata,
            "extra": extra_metadata,
            "combined_train_point_sha256": array_sha256(train_t, train_x),
            "combined_validation_point_sha256": array_sha256(valid_t, valid_x),
        },
        "optimization_levels": records,
        "selected_level": MORE_DATA_LEVELS[-1][1].name,
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
