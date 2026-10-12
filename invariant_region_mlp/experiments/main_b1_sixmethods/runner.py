"""Resumable execution for one B1 stage/cell/method."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import time
import traceback
from typing import Any, Callable

import numpy as np
import torch

from .actor import actor_checkpoint_path, train_actor
from .actor_mlp import (
    actor_linearization_identity,
    evaluate_actor_mlp,
)
from .common import (
    CHECKPOINT_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    ensure_result_dirs,
    fixed_test_points,
    load_config,
)
from .deep_bsde_adapter import (
    deep_bsde_checkpoint_path,
    predict_deep_bsde,
    train_deep_bsde,
)
from .dpi_adapter import dpi_checkpoint_path, predict_dpi, train_dpi
from .equation import B1Equation
from .metrics import prediction_metrics
from .mlp_adapter import evaluate_mlp
from .sdgd_pinn_adapter import (
    predict_sdgd_pinn,
    sdgd_checkpoint_path,
    train_sdgd_pinn,
)


METHODS = (
    "vanilla_mlp",
    "double_mlp",
    "actor_mlp",
    "dpi",
    "pinn_sdgd",
    "deep_bsde",
)


def _paths(stage: str, method: str, d: int, k: int) -> tuple[Path, Path]:
    raw = RESULTS_DIR / "raw" / stage
    raw.mkdir(parents=True, exist_ok=True)
    return (
        raw / f"{method}_d{d}_k{k}.json",
        raw / f"{method}_d{d}_k{k}_predictions.npz",
    )


def _save_predictions(path: Path, predictions: list[dict[str, Any]]) -> None:
    payload: dict[str, np.ndarray] = {}
    for item in predictions:
        seed = item["seed"]
        payload[f"value_seed{seed}"] = item["value"]
        if item["z"] is not None:
            payload[f"z_seed{seed}"] = item["z"]
    np.savez_compressed(path, **payload)


def _base_row(
    stage: str,
    method: str,
    equation: B1Equation,
    repetition: int,
    points: np.ndarray,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "method": method,
        "d": equation.d,
        "k": equation.k,
        "K": equation.K,
        "coefficient_scale": None if equation.k == 1 else equation.coefficient_scale,
        "repetition": repetition,
        "test_points": len(points),
        "test_point_sha256": array_sha256(points),
        "test_time": 0.0,
        "test_distribution": "x~N(0,I_d)",
    }


def run_method(
    stage: str,
    method: str,
    d: int,
    k: int,
    *,
    force: bool = False,
) -> dict[str, Any]:
    if stage not in {"stage_a", "stage_b"}:
        raise ValueError(stage)
    if method not in METHODS:
        raise ValueError(method)
    ensure_result_dirs()
    output_path, prediction_path = _paths(stage, method, d, k)
    if output_path.is_file() and not force:
        with output_path.open() as handle:
            return json.load(handle)

    config = load_config()
    stage_config = config[stage]
    scale = float(config["pde"]["frozen_coefficient_scale"])
    equation = B1Equation(d=d, k=k, coefficient_scale=scale)
    point_count = int(stage_config["test_points"])
    points = fixed_test_points(d, k, point_count)
    t = np.zeros(point_count, dtype=np.float64)
    reference_u = equation.exact_u(t, points)
    reference_z = equation.exact_z(t, points)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("Stage A/B methods require an allocated GPU")

    rows: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    shared: dict[str, Any] = {}
    started = time.perf_counter()

    if method in {"vanilla_mlp", "double_mlp"}:
        for repetition in range(int(stage_config["mlp_seeds"])):
            state, cost = evaluate_mlp(
                equation,
                t,
                points,
                n=int(stage_config["mlp_n"]),
                M=int(stage_config["mlp_M"]),
                method=method,
                repetition=repetition,
                device=device,
            )
            row = _base_row(stage, method, equation, repetition, points)
            row.update(prediction_metrics(state[:, 0], state[:, 1:], reference_u, reference_z))
            row.update({f"cost_{key}": value for key, value in cost.items()})
            row.update({"training_seconds": 0.0, "inference_seconds": cost["wall_seconds"], "parameter_count": 0})
            rows.append(row)
            predictions.append({"seed": repetition, "value": state[:, 0], "z": state[:, 1:]})

    elif method == "actor_mlp":
        actor_path = actor_checkpoint_path(CHECKPOINT_DIR, equation)
        actor, actor_metadata = train_actor(
            equation, device=device, checkpoint_path=actor_path, force=force
        )
        identity_error = actor_linearization_identity(actor, equation, device)
        if identity_error > 2e-5:
            raise RuntimeError(f"actor linearization identity failed: {identity_error}")
        shared["actor"] = actor_metadata
        shared["actor_identity_max_abs_error"] = identity_error
        for repetition in range(int(stage_config["mlp_seeds"])):
            state, cost = evaluate_actor_mlp(
                equation,
                actor,
                t,
                points,
                n=int(stage_config["mlp_n"]),
                M=int(stage_config["mlp_M"]),
                repetition=repetition,
                device=device,
            )
            row = _base_row(stage, method, equation, repetition, points)
            row.update(prediction_metrics(state[:, 0], state[:, 1:], reference_u, reference_z))
            row.update({f"cost_{key}": value for key, value in cost.items()})
            row.update(
                {
                    "training_seconds": actor_metadata["total_actor_seconds"],
                    "inference_seconds": cost["wall_seconds"],
                    "parameter_count": actor_metadata["parameter_count"],
                    "actor_identity_max_abs_error": identity_error,
                }
            )
            rows.append(row)
            predictions.append({"seed": repetition, "value": state[:, 0], "z": state[:, 1:]})

    else:
        neural_seeds = int(stage_config["neural_seeds"])
        for repetition in range(neural_seeds):
            if method == "dpi":
                checkpoint = dpi_checkpoint_path(CHECKPOINT_DIR, equation, repetition, stage)
                model, metadata = train_dpi(
                    equation,
                    seed=repetition,
                    stage=stage,
                    device=device,
                    checkpoint_path=checkpoint,
                    force=force,
                )
                value, z, inference = predict_dpi(model, equation, points, device)
            elif method == "deep_bsde":
                checkpoint = deep_bsde_checkpoint_path(CHECKPOINT_DIR, equation, repetition, stage)
                model, metadata = train_deep_bsde(
                    equation,
                    seed=repetition,
                    stage=stage,
                    device=device,
                    checkpoint_path=checkpoint,
                    force=force,
                )
                value, z, inference = predict_deep_bsde(model, equation, points, device)
            else:
                checkpoint = sdgd_checkpoint_path(CHECKPOINT_DIR, equation, repetition, stage)
                model, metadata = train_sdgd_pinn(
                    equation,
                    seed=repetition,
                    stage=stage,
                    device=device,
                    checkpoint_path=checkpoint,
                    force=force,
                )
                value, z, inference = predict_sdgd_pinn(model, equation, points, device)
            row = _base_row(stage, method, equation, repetition, points)
            row.update(prediction_metrics(value, z, reference_u, reference_z))
            row.update(
                {
                    "training_seconds": metadata["training_seconds"],
                    "inference_seconds": inference,
                    "parameter_count": metadata["parameter_count"],
                    "peak_gpu_bytes": metadata["peak_gpu_bytes"],
                    "optimizer_steps": metadata["optimizer_steps"],
                }
            )
            rows.append(row)
            predictions.append({"seed": repetition, "value": value, "z": z})
            shared[f"seed_{repetition}"] = metadata
            del model
            torch.cuda.empty_cache()

    payload = {
        "schema_version": 1,
        "status": "complete",
        "stage": stage,
        "method": method,
        "d": d,
        "k": k,
        "device": str(device),
        "elapsed_seconds": time.perf_counter() - started,
        "reference_only_used_after_training": method in {"actor_mlp", "dpi", "deep_bsde", "pinn_sdgd"},
        "rows": rows,
        "shared": shared,
    }
    _save_predictions(prediction_path, predictions)
    atomic_json(output_path, payload)
    print(json.dumps({key: payload[key] for key in ("status", "stage", "method", "d", "k", "elapsed_seconds")}, indent=2))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("stage_a", "stage_b"))
    parser.add_argument("--method", required=True, choices=METHODS)
    parser.add_argument("--d", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        run_method(args.stage, args.method, args.d, args.k, force=args.force)
    except Exception as error:
        failure_path, _ = _paths(args.stage, args.method, args.d, args.k)
        atomic_json(
            failure_path.with_name(failure_path.stem + "_failure.json"),
            {
                "status": "failed",
                "stage": args.stage,
                "method": args.method,
                "d": args.d,
                "k": args.k,
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()
