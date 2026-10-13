"""Run the frozen-actor martingale control-variate evaluation pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import time
import traceback
from typing import Any

import numpy as np
import torch

from .actor_convergence import load_converged_actor
from .actor_iteration import evaluate_policy_mc_with_martingale_cv
from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    fixed_test_points,
    load_config,
)
from .equation import B1Equation


CELLS = ((100, 10), (400, 40), (400, 1))
ACTOR_ITERATIONS = (0, 1)
TIME_STEPS = (16, 32, 64)
PATHS = 256
SEEDS = 3
PARENT_BATCH_SIZE = 32
PATH_CHUNK_SIZE = 64
OUTPUT_DIR = RESULTS_DIR / "policy_evaluation_cv_pilot"
PREDICTION_DIR = OUTPUT_DIR / "predictions"
ACTOR_ITERATION_DIR = RESULTS_DIR / "actor_iteration_pilot"
CHECKPOINT_DIR = ACTOR_ITERATION_DIR / "checkpoints"


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    mean = float(np.mean(array))
    ci95 = 1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array))
    return mean, ci95


def _atomic_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _checkpoint_path(d: int, k: int, iteration: int) -> Path:
    return CHECKPOINT_DIR / f"actor_iteration_d{d}_k{k}_q{iteration}.pt"


def _legacy_predictions(d: int, k: int) -> dict[str, np.ndarray]:
    path = (
        ACTOR_ITERATION_DIR
        / "predictions"
        / f"actor_iteration_d{d}_k{k}.npz"
    )
    with np.load(path, allow_pickle=False) as payload:
        return {
            f"q{iteration}": np.asarray(
                payload[f"q{iteration}_direct_policy_value"], dtype=np.float64
            )
            for iteration in ACTOR_ITERATIONS
        }


def _cell_evaluation(
    d: int,
    k: int,
    *,
    device: torch.device,
) -> tuple[list[dict[str, Any]], dict[str, np.ndarray], dict[str, Any]]:
    config = load_config()
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    points = fixed_test_points(d, k, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact = equation.exact_u(times, points)
    legacy = _legacy_predictions(d, k)
    rows: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {}
    checkpoint_metadata: dict[str, Any] = {}

    for iteration in ACTOR_ITERATIONS:
        checkpoint = _checkpoint_path(d, k, iteration)
        if not checkpoint.is_file():
            raise FileNotFoundError(f"missing frozen actor checkpoint: {checkpoint}")
        actor, metadata = load_converged_actor(checkpoint, device)
        checkpoint_metadata[f"q{iteration}"] = {
            "path": str(checkpoint.relative_to(RESULTS_DIR.parent)),
            "label_source": metadata["label_source"],
            "analytic_reference_used": metadata["analytic_reference_used"],
        }
        for time_steps in TIME_STEPS:
            raw_by_seed: list[np.ndarray] = []
            corrected_by_seed: list[np.ndarray] = []
            raw_variance_by_seed: list[np.ndarray] = []
            corrected_variance_by_seed: list[np.ndarray] = []
            seed_rows: list[dict[str, Any]] = []
            costs: list[dict[str, Any]] = []
            for repetition in range(SEEDS):
                (
                    raw,
                    corrected,
                    raw_variance,
                    corrected_variance,
                    cost,
                ) = evaluate_policy_mc_with_martingale_cv(
                    equation,
                    actor,
                    times,
                    points,
                    time_steps=time_steps,
                    paths=PATHS,
                    repetition=repetition,
                    seed_namespace="actor_iteration_direct_test",
                    device=device,
                    parent_batch_size=PARENT_BATCH_SIZE,
                    path_chunk_size=PATH_CHUNK_SIZE,
                    coarse_time_steps=TIME_STEPS[0],
                    maximum_time_steps=TIME_STEPS[-1],
                )
                raw_error = raw - exact
                corrected_error = corrected - exact
                raw_mean_variance = float(np.mean(raw_variance))
                corrected_mean_variance = float(np.mean(corrected_variance))
                seed_rows.append(
                    {
                        "repetition": repetition,
                        "raw_value_rmse": float(
                            np.sqrt(np.mean(np.square(raw_error)))
                        ),
                        "corrected_value_rmse": float(
                            np.sqrt(np.mean(np.square(corrected_error)))
                        ),
                        "raw_value_bias": float(np.mean(raw_error)),
                        "corrected_value_bias": float(
                            np.mean(corrected_error)
                        ),
                        "raw_estimator_variance": raw_mean_variance,
                        "corrected_estimator_variance": corrected_mean_variance,
                        "variance_reduction_factor": raw_mean_variance
                        / max(
                            corrected_mean_variance,
                            np.finfo(np.float64).tiny,
                        ),
                        "mean_control_variate_shift": float(
                            np.mean(corrected - raw)
                        ),
                    }
                )
                raw_by_seed.append(raw)
                corrected_by_seed.append(corrected)
                raw_variance_by_seed.append(raw_variance)
                corrected_variance_by_seed.append(corrected_variance)
                costs.append(cost)

            raw_stack = np.stack(raw_by_seed)
            corrected_stack = np.stack(corrected_by_seed)
            raw_rmse, raw_rmse_ci95 = _mean_ci95(
                [row["raw_value_rmse"] for row in seed_rows]
            )
            corrected_rmse, corrected_rmse_ci95 = _mean_ci95(
                [row["corrected_value_rmse"] for row in seed_rows]
            )
            paired_improvement, paired_improvement_ci95 = _mean_ci95(
                [
                    row["raw_value_rmse"] - row["corrected_value_rmse"]
                    for row in seed_rows
                ]
            )
            raw_estimator_variance = float(
                np.mean(np.stack(raw_variance_by_seed))
            )
            corrected_estimator_variance = float(
                np.mean(np.stack(corrected_variance_by_seed))
            )
            raw_cross_seed_variance = float(
                np.mean(np.var(raw_stack, axis=0, ddof=1))
            )
            corrected_cross_seed_variance = float(
                np.mean(np.var(corrected_stack, axis=0, ddof=1))
            )
            legacy_difference = (
                float(np.max(np.abs(raw_stack - legacy[f"q{iteration}"])))
                if time_steps == TIME_STEPS[0]
                else None
            )
            key = f"d{d}_k{k}_q{iteration}_n{time_steps}"
            predictions[f"{key}_raw"] = raw_stack
            predictions[f"{key}_corrected"] = corrected_stack
            rows.append(
                {
                    "d": d,
                    "k": k,
                    "actor_iteration": iteration,
                    "time_steps": time_steps,
                    "paths": PATHS,
                    "seeds": SEEDS,
                    "raw_value_rmse": raw_rmse,
                    "raw_value_ci95": raw_rmse_ci95,
                    "corrected_value_rmse": corrected_rmse,
                    "corrected_value_ci95": corrected_rmse_ci95,
                    "paired_rmse_improvement": paired_improvement,
                    "paired_rmse_improvement_ci95": paired_improvement_ci95,
                    "raw_value_bias": float(
                        np.mean([row["raw_value_bias"] for row in seed_rows])
                    ),
                    "corrected_value_bias": float(
                        np.mean(
                            [row["corrected_value_bias"] for row in seed_rows]
                        )
                    ),
                    "raw_estimator_variance": raw_estimator_variance,
                    "corrected_estimator_variance": corrected_estimator_variance,
                    "antithetic_variance_reduction_factor": (
                        raw_estimator_variance
                        / max(
                            corrected_estimator_variance,
                            np.finfo(np.float64).tiny,
                        )
                    ),
                    "raw_cross_seed_variance": raw_cross_seed_variance,
                    "corrected_cross_seed_variance": corrected_cross_seed_variance,
                    "cross_seed_variance_reduction_factor": (
                        raw_cross_seed_variance
                        / max(
                            corrected_cross_seed_variance,
                            np.finfo(np.float64).tiny,
                        )
                    ),
                    "legacy_n16_max_abs_difference": legacy_difference,
                    "wall_seconds": float(
                        sum(cost["wall_seconds"] for cost in costs)
                    ),
                    "peak_gpu_bytes": int(
                        max(cost["peak_gpu_bytes"] or 0 for cost in costs)
                    ),
                    "seed_rows": seed_rows,
                    "costs": costs,
                }
            )

    discretization: list[dict[str, Any]] = []
    for iteration in ACTOR_ITERATIONS:
        for estimator in ("raw", "corrected"):
            coarse = predictions[f"d{d}_k{k}_q{iteration}_n16_{estimator}"]
            fine = predictions[f"d{d}_k{k}_q{iteration}_n64_{estimator}"]
            difference = fine - coarse
            discretization.append(
                {
                    "actor_iteration": iteration,
                    "estimator": estimator,
                    "paired_n64_minus_n16_bias": float(np.mean(difference)),
                    "paired_n64_minus_n16_rms": float(
                        np.sqrt(np.mean(np.square(difference)))
                    ),
                }
            )
    cell_metadata = {
        "d": d,
        "k": k,
        "test_point_sha256": array_sha256(points),
        "checkpoint_metadata": checkpoint_metadata,
        "discretization": discretization,
    }
    return rows, predictions, cell_metadata


def run(*, force: bool = False) -> dict[str, Any]:
    result_path = OUTPUT_DIR / "result.json"
    if result_path.is_file() and not force:
        existing = json.loads(result_path.read_text())
        if existing.get("status") == "complete":
            print(json.dumps(existing["headline"], indent=2), flush=True)
            return existing
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("policy-evaluation CV pilot requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    started = time.perf_counter()
    rows: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {}
    cells: list[dict[str, Any]] = []
    for d, k in CELLS:
        cell_rows, cell_predictions, metadata = _cell_evaluation(
            d, k, device=device
        )
        rows.extend(cell_rows)
        predictions.update(cell_predictions)
        cells.append(metadata)
        q1_rows = {
            row["time_steps"]: row
            for row in cell_rows
            if row["actor_iteration"] == 1
        }
        print(
            json.dumps(
                {
                    "d": d,
                    "k": k,
                    "q1": {
                        steps: {
                            "raw_rmse": q1_rows[steps]["raw_value_rmse"],
                            "corrected_rmse": q1_rows[steps][
                                "corrected_value_rmse"
                            ],
                            "variance_reduction": q1_rows[steps][
                                "antithetic_variance_reduction_factor"
                            ],
                        }
                        for steps in TIME_STEPS
                    },
                }
            ),
            flush=True,
        )

    _atomic_npz(PREDICTION_DIR / "paired_predictions.npz", **predictions)
    csv_rows = [
        {
            key: value
            for key, value in row.items()
            if key not in {"seed_rows", "costs"}
        }
        for row in rows
    ]
    _atomic_csv(OUTPUT_DIR / "comparison.csv", csv_rows)
    q1_n64 = [
        row
        for row in rows
        if row["actor_iteration"] == 1 and row["time_steps"] == 64
    ]
    headline = {
        "q1_n64_corrected_beats_raw_cells": sum(
            row["corrected_value_rmse"] < row["raw_value_rmse"]
            for row in q1_n64
        ),
        "q1_n64_mean_antithetic_variance_reduction_factor": float(
            np.mean(
                [
                    row["antithetic_variance_reduction_factor"]
                    for row in q1_n64
                ]
            )
        ),
        "maximum_legacy_n16_absolute_difference": float(
            max(
                row["legacy_n16_max_abs_difference"]
                for row in rows
                if row["time_steps"] == 16
            )
        ),
        "q1_n64": [
            {
                "d": row["d"],
                "k": row["k"],
                "raw_value_rmse": row["raw_value_rmse"],
                "corrected_value_rmse": row["corrected_value_rmse"],
                "raw_value_bias": row["raw_value_bias"],
                "corrected_value_bias": row["corrected_value_bias"],
                "antithetic_variance_reduction_factor": row[
                    "antithetic_variance_reduction_factor"
                ],
            }
            for row in q1_n64
        ],
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "device": str(device),
        "protocol": {
            "cells": [list(cell) for cell in CELLS],
            "actor_iterations": list(ACTOR_ITERATIONS),
            "time_steps": list(TIME_STEPS),
            "paths": PATHS,
            "seeds": SEEDS,
            "test_points": 960,
            "same_brownian_for_raw_and_corrected": True,
            "nested_brownian_across_time_steps": True,
            "additional_actor_calls_for_control_variate": 0,
            "actor_retraining": False,
            "analytic_reference_use": "final evaluation only",
            "actor_structure_prior": (
                "terminal-gradient-derived rank-k output subspace with supplied k"
            ),
            "nonlinear_mlp_teacher_used": False,
        },
        "cells": cells,
        "rows": rows,
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
