"""Compare frozen actors under identical martingale-CV policy evaluation."""

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

from .actor_convergence import load_converged_actor, predict_actor
from .actor_iteration import evaluate_policy_mc_with_analytic_gap
from .common import RESULTS_DIR, array_sha256, atomic_json, fixed_test_points, load_config
from .equation import B1Equation


D = 400
K = 40
TIME_STEPS = 64
PATHS = 256
SEEDS = 3
PARENT_BATCH_SIZE = 32
PATH_CHUNK_SIZE = 64

OUTPUT_DIR = RESULTS_DIR / "actor_quality_upgrade" / "ceiling"
PREDICTION_DIR = OUTPUT_DIR / "predictions"
ACTOR_ITERATION_DIR = RESULTS_DIR / "actor_iteration_pilot"
CONVERGENCE_DIR = RESULTS_DIR / "actor_convergence" / "protocol_v3"
ACTORS = {
    "teacher_free_q1": (
        ACTOR_ITERATION_DIR
        / "checkpoints"
        / "actor_iteration_d400_k40_q1.pt"
    ),
    "converged_actor": (
        CONVERGENCE_DIR
        / "checkpoints"
        / "actor_converged_d400_k40.pt"
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mean_ci95(values: list[float]) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    return (
        float(np.mean(array)),
        1.96 * float(np.std(array, ddof=1)) / math.sqrt(len(array)),
    )


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def _evaluate_actor(
    name: str,
    checkpoint: Path,
    equation: B1Equation,
    times: np.ndarray,
    points: np.ndarray,
    exact_u: np.ndarray,
    exact_z: np.ndarray,
    *,
    device: torch.device,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    actor, checkpoint_metadata = load_converged_actor(checkpoint, device)
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
    for repetition in range(SEEDS):
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
            time_steps=TIME_STEPS,
            paths=PATHS,
            repetition=repetition,
            seed_namespace="actor_iteration_direct_test",
            device=device,
            parent_batch_size=PARENT_BATCH_SIZE,
            path_chunk_size=PATH_CHUNK_SIZE,
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
                "raw_estimator_variance": float(np.mean(raw_variance)),
                "corrected_estimator_variance": float(
                    np.mean(corrected_variance)
                ),
                "policy_gap_estimator_variance": float(
                    np.mean(policy_gap_variance)
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
    mean_identity_residual, mean_identity_residual_ci95 = _mean_ci95(
        [row["mean_identity_residual"] for row in seed_rows]
    )
    raw_variance = float(np.mean(stacked["raw_estimator_variance"]))
    corrected_variance = float(
        np.mean(stacked["corrected_estimator_variance"])
    )
    row = {
        "actor": name,
        "checkpoint": str(checkpoint.relative_to(RESULTS_DIR.parent)),
        "checkpoint_sha256": _sha256(checkpoint),
        "checkpoint_metadata": checkpoint_metadata,
        "direct_actor_gradient_rmse": direct_gradient_rmse,
        "raw_value_rmse": raw_rmse,
        "raw_value_ci95": raw_ci95,
        "corrected_value_rmse": corrected_rmse,
        "corrected_value_ci95": corrected_ci95,
        "corrected_value_bias": corrected_bias,
        "corrected_value_bias_ci95": corrected_bias_ci95,
        "mean_analytic_policy_gap": mean_gap,
        "mean_analytic_policy_gap_ci95": mean_gap_ci95,
        "policy_gap_fraction_of_corrected_bias": mean_gap
        / max(abs(corrected_bias), np.finfo(np.float64).tiny),
        "mean_identity_residual": mean_identity_residual,
        "mean_identity_residual_ci95": mean_identity_residual_ci95,
        "raw_estimator_variance": raw_variance,
        "corrected_estimator_variance": corrected_variance,
        "variance_reduction_factor": raw_variance
        / max(corrected_variance, np.finfo(np.float64).tiny),
        "seed_rows": seed_rows,
        "costs": costs,
    }
    prediction_arrays = {
        f"{name}_{key}": value for key, value in stacked.items()
    }
    prediction_arrays[f"{name}_direct_actor"] = direct_actor
    return row, prediction_arrays


def run(*, force: bool = False) -> dict[str, Any]:
    result_path = OUTPUT_DIR / "result.json"
    if result_path.is_file() and not force:
        result = json.loads(result_path.read_text())
        if result.get("status") == "complete":
            print(json.dumps(result["headline"], indent=2), flush=True)
            return result
    missing = [str(path) for path in ACTORS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing actor checkpoints: {missing}")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Actor-quality ceiling requires an allocated GPU")
    torch.set_float32_matmul_precision("high")
    config = load_config()
    equation = B1Equation(
        d=D,
        k=K,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    points = fixed_test_points(D, K, int(config["stage_b"]["test_points"]))
    times = np.zeros(len(points), dtype=np.float64)
    exact_u = equation.exact_u(times, points)
    exact_z = equation.exact_z(times, points)
    started = time.perf_counter()
    rows: list[dict[str, Any]] = []
    prediction_arrays: dict[str, np.ndarray] = {}
    for name, checkpoint in ACTORS.items():
        row, actor_arrays = _evaluate_actor(
            name,
            checkpoint,
            equation,
            times,
            points,
            exact_u,
            exact_z,
            device=device,
        )
        rows.append(row)
        prediction_arrays.update(actor_arrays)
        print(
            json.dumps(
                {
                    "actor": name,
                    "direct_actor_gradient_rmse": row[
                        "direct_actor_gradient_rmse"
                    ],
                    "corrected_value_rmse": row["corrected_value_rmse"],
                    "corrected_value_bias": row["corrected_value_bias"],
                    "mean_analytic_policy_gap": row[
                        "mean_analytic_policy_gap"
                    ],
                }
            ),
            flush=True,
        )

    row_lookup = {row["actor"]: row for row in rows}
    teacher_free = row_lookup["teacher_free_q1"]
    converged = row_lookup["converged_actor"]
    legacy_path = (
        RESULTS_DIR
        / "policy_evaluation_cv_pilot"
        / "predictions"
        / "paired_predictions.npz"
    )
    legacy_differences: dict[str, float] = {}
    with np.load(legacy_path, allow_pickle=False) as legacy:
        for estimator in ("raw", "corrected"):
            key = f"d400_k40_q1_n64_{estimator}"
            observed = prediction_arrays[f"teacher_free_q1_{estimator}"]
            legacy_differences[estimator] = float(
                np.max(np.abs(observed - legacy[key]))
            )
    relative_rmse_improvement = 1.0 - (
        converged["corrected_value_rmse"]
        / teacher_free["corrected_value_rmse"]
    )
    material_headroom = (
        relative_rmse_improvement >= 0.10
        and converged["corrected_value_rmse"]
        + converged["corrected_value_ci95"]
        < teacher_free["corrected_value_rmse"]
        - teacher_free["corrected_value_ci95"]
    )
    headline = {
        "teacher_free_q1_corrected_value_rmse": teacher_free[
            "corrected_value_rmse"
        ],
        "converged_actor_corrected_value_rmse": converged[
            "corrected_value_rmse"
        ],
        "converged_actor_relative_rmse_improvement": relative_rmse_improvement,
        "material_actor_headroom": material_headroom,
        "upgrade_gate_passed": material_headroom,
        "legacy_q1_n64_max_abs_difference": max(legacy_differences.values()),
    }
    payload = {
        "schema_version": 1,
        "status": "complete",
        "device": str(device),
        "protocol": {
            "d": D,
            "k": K,
            "test_points": len(points),
            "test_point_sha256": array_sha256(points),
            "time_steps": TIME_STEPS,
            "paths": PATHS,
            "seeds": SEEDS,
            "same_brownian_across_actors": True,
            "same_evaluator_as_policy_evaluation_cv_pilot": True,
            "analytic_reference_use": "post-training diagnostic only",
            "analytic_reference_used_for_training": False,
            "analytic_reference_used_for_checkpoint_selection": False,
            "actor_structure_prior": (
                "terminal-gradient-derived rank-k output subspace with supplied k"
            ),
        },
        "rows": rows,
        "legacy_differences": legacy_differences,
        "headline": headline,
        "elapsed_seconds": time.perf_counter() - started,
        "failed_runs": 0,
        "nan_failures": 0,
        "memory_failures": 0,
    }
    _atomic_npz(PREDICTION_DIR / "ceiling_predictions.npz", **prediction_arrays)
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
