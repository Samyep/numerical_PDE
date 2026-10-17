"""Frozen post-selection evaluation of Tangent-Actor checkpoints."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .benchmarks import BENCHMARK_NAMES, DPIOUBenchmark, make_benchmark
from .common import (
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    stable_seed,
)
from .diagnostics import quadratic_bregman_diagnostic
from .models import load_actor
from .references import dpi_exact_reference, hjb_log_reference
from .tangent_solver import evaluate_policy


def _native_reference(name: str) -> float:
    if name == "hjb_log_d100":
        return float(
            hjb_log_reference(np.asarray([0.0]), np.zeros((1, 100)))[0]
        )
    if name == "allen_cahn_d100":
        audit = json.loads((RESULTS_DIR / "reference_accuracy_audit.json").read_text())
        return float(audit["allen_cahn"]["native_u_0_0_extrapolated"])
    benchmark = DPIOUBenchmark()
    return float(
        dpi_exact_reference(benchmark, np.asarray([0.0]), np.zeros((1, 100)))[0]
    )


def _metrics(prediction: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    error = np.asarray(prediction) - np.asarray(reference)
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
    }


def evaluate(
    benchmark_name: str,
    trial: str,
    seed: int,
    iteration: int,
    device: torch.device,
    *,
    paths: int,
    time_steps: int,
    replicates: int,
    native_paths: int,
) -> dict[str, Any]:
    ensure_dirs()
    evaluation_process_started = time.perf_counter()
    benchmark = make_benchmark(benchmark_name)
    training_path = RAW_DIR / f"actor_{benchmark_name}_{trial}_s{seed}.json"
    training = json.loads(training_path.read_text())
    phase = training["phases"][iteration - 1]
    actor, payload = load_actor(phase["checkpoint"], benchmark, device)
    with np.load(RESULTS_DIR / f"evaluation_bank_{benchmark_name}.npz") as bank:
        indices = bank["test_indices"]
        t = bank["t"][indices]
        x = bank["x"][indices]
        reference = bank["value"][indices]
        reference_z = bank["z"][indices]

    predictions: list[np.ndarray] = []
    plain_predictions: list[np.ndarray] = []
    replicate_metadata: list[dict[str, Any]] = []
    evaluation_seconds = 0.0
    for replicate in range(replicates):
        result, metadata = evaluate_policy(
            benchmark,
            actor,
            t,
            x,
            paths=paths,
            time_steps=time_steps,
            seed=stable_seed(
                "original_benchmark_final_evaluation",
                benchmark_name,
                seed,
                iteration,
                time_steps,
                replicate,
            ),
            device=device,
            parent_batch=8,
            path_chunk=64,
        )
        key = "mcv" if benchmark.actor_kind == "gradient" else "plain"
        predictions.append(result[key])
        plain_predictions.append(result["plain"])
        replicate_metadata.append(metadata)
        evaluation_seconds += float(metadata["wall_seconds"])
    matrix = np.asarray(predictions)
    plain_matrix = np.asarray(plain_predictions)
    mean_prediction = matrix.mean(axis=0)
    mean_plain = plain_matrix.mean(axis=0)
    corrected_metrics = _metrics(mean_prediction, reference)
    plain_metrics = _metrics(mean_plain, reference)
    per_replicate_rmse = [
        _metrics(row, reference)["rmse"] for row in matrix
    ]
    per_replicate_plain_rmse = [
        _metrics(row, reference)["rmse"] for row in plain_matrix
    ]

    tt = torch.as_tensor(t, dtype=torch.float32, device=device)
    xx = torch.as_tensor(x, dtype=torch.float32, device=device)
    actor_started = time.perf_counter()
    with torch.no_grad():
        r, q = actor.predict(tt, xx)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    actor_inference_seconds = time.perf_counter() - actor_started
    if q is not None:
        direct = q.detach().double().cpu().numpy()
        difference = direct - reference_z
        actor_diagnostics: dict[str, Any] = {
            "target": "z=sigma*grad(u)",
            "component_rmse": float(np.sqrt(np.mean(np.square(difference)))),
            "vector_rmse": float(
                np.sqrt(np.mean(np.square(difference).sum(axis=1)))
            ),
            "relative_l2": float(
                np.linalg.norm(difference) / np.linalg.norm(reference_z)
            ),
        }
    else:
        if r is None:
            raise RuntimeError("Value Actor returned no value")
        direct = r.detach().double().cpu().numpy()
        actor_diagnostics = {
            "target": "value",
            **_metrics(direct, reference),
        }

    native_t = np.asarray([0.0])
    native_x = np.zeros((1, benchmark.d), dtype=np.float64)
    native_replicates: list[float] = []
    native_plain_replicates: list[float] = []
    native_metadata: list[dict[str, Any]] = []
    for replicate in range(max(replicates, 4)):
        result, metadata = evaluate_policy(
            benchmark,
            actor,
            native_t,
            native_x,
            paths=native_paths,
            time_steps=time_steps,
            seed=stable_seed(
                "original_benchmark_native_final",
                benchmark_name,
                seed,
                iteration,
                replicate,
            ),
            device=device,
            parent_batch=1,
            path_chunk=64,
        )
        key = "mcv" if benchmark.actor_kind == "gradient" else "plain"
        native_replicates.append(float(result[key][0]))
        native_plain_replicates.append(float(result["plain"][0]))
        native_metadata.append(metadata)
        evaluation_seconds += float(metadata["wall_seconds"])
    native_reference = _native_reference(benchmark_name)
    native_mean = float(np.mean(native_replicates))
    native_plain_mean = float(np.mean(native_plain_replicates))

    training_peak_candidates: list[int] = []
    training_nonfinite_count = 0
    for completed_phase in training["phases"][:iteration]:
        fit_peak = completed_phase["training"].get("peak_gpu_bytes")
        if fit_peak is not None:
            training_peak_candidates.append(int(fit_peak))
        for label_split in ("train", "validation", "replication"):
            label_metadata = completed_phase["labels"].get(label_split, {})
            label_peak = label_metadata.get("peak_gpu_bytes")
            if label_peak is not None:
                training_peak_candidates.append(int(label_peak))
            training_nonfinite_count += int(
                label_metadata.get("nonfinite_count", 0)
            )
    evaluation_peak_candidates = [
        int(metadata["peak_gpu_bytes"])
        for metadata in replicate_metadata + native_metadata
        if metadata.get("peak_gpu_bytes") is not None
    ]
    evaluation_nonfinite_count = sum(
        int(metadata.get("nonfinite_count", 0))
        for metadata in replicate_metadata + native_metadata
    )

    if benchmark.actor_kind == "gradient":
        diagnostic_indices = np.linspace(
            0, len(reference) - 1, min(96, len(reference)), dtype=int
        )
        bregman = quadratic_bregman_diagnostic(
            benchmark,
            actor,
            t[diagnostic_indices],
            x[diagnostic_indices],
            paths=128,
            time_steps=time_steps,
            seed=stable_seed(
                "original_benchmark_bregman",
                benchmark_name,
                seed,
                iteration,
                time_steps,
            ),
            device=device,
        )
        observed_gap = float(
            np.mean(mean_prediction[diagnostic_indices] - reference[diagnostic_indices])
        )
        bregman["observed_mcv_policy_minus_exact_mean"] = observed_gap
        bregman["identity_residual"] = (
            observed_gap - float(bregman["mean_integrated_bregman_gap"])
        )
        bregman["observed_gap_uses_frozen_final_evaluation"] = True
    else:
        bregman = {
            "available": False,
            "reason": (
                "Allen-Cahn has a value-dependent cubic Bregman term; the "
                "quadratic-gradient identity does not apply"
            ),
        }

    prediction_path = PREDICTION_DIR / (
        f"actor_{benchmark_name}_{trial}_s{seed}_a{iteration}_n{time_steps}.npz"
    )
    np.savez_compressed(
        prediction_path,
        corrected=matrix,
        plain=plain_matrix,
        test_indices=indices,
        native=np.asarray(native_replicates),
        native_plain=np.asarray(native_plain_replicates),
        direct_actor=direct,
    )
    evaluation_process_seconds = time.perf_counter() - evaluation_process_started
    result = {
        "method": "terminal_aware_tangent_actor_mcv",
        "benchmark": benchmark_name,
        "trial": trial,
        "seed": seed,
        "iteration": iteration,
        "checkpoint": phase["checkpoint"],
        "configuration": training["configuration"],
        "parameter_count": sum(parameter.numel() for parameter in actor.parameters()),
        "states": len(reference),
        "paths_per_replicate": paths,
        "replicates": replicates,
        "time_steps": time_steps,
        "corrected_metrics_of_replicate_mean": corrected_metrics,
        "plain_metrics_of_replicate_mean": plain_metrics,
        "corrected_replicate_rmse_mean": float(np.mean(per_replicate_rmse)),
        "corrected_replicate_rmse_std": float(
            np.std(per_replicate_rmse, ddof=1)
        )
        if len(per_replicate_rmse) > 1
        else 0.0,
        "plain_replicate_rmse_mean": float(np.mean(per_replicate_plain_rmse)),
        "plain_replicate_rmse_std": float(
            np.std(per_replicate_plain_rmse, ddof=1)
        )
        if len(per_replicate_plain_rmse) > 1
        else 0.0,
        "sampling_variance_component_mse": float(
            np.mean(np.var(matrix, axis=0, ddof=1)) / replicates
        )
        if replicates > 1
        else None,
        "plain_sampling_variance_component_mse": float(
            np.mean(np.var(plain_matrix, axis=0, ddof=1)) / replicates
        )
        if replicates > 1
        else None,
        "mean_variance_reduction": float(
            np.mean([item["variance_reduction"] for item in replicate_metadata])
        ),
        "mean_weight_ess_fraction": float(
            np.mean([item["weight_ess_fraction"] for item in replicate_metadata])
        ),
        "actor_diagnostics": actor_diagnostics,
        "bregman_diagnostics": bregman,
        "actor_inference_seconds": actor_inference_seconds,
        "peak_training_gpu_bytes": max(training_peak_candidates, default=None),
        "peak_inference_gpu_bytes": max(evaluation_peak_candidates, default=None),
        "peak_gpu_bytes": max(
            training_peak_candidates + evaluation_peak_candidates,
            default=None,
        ),
        "training_nonfinite_count": training_nonfinite_count,
        "evaluation_nonfinite_count": evaluation_nonfinite_count,
        "native": {
            "point": "t=0,x=0",
            "is_published_native_point": benchmark_name
            in {"hjb_log_d100", "allen_cahn_d100"},
            "reference": native_reference,
            "corrected_mean": native_mean,
            "corrected_absolute_error": abs(native_mean - native_reference),
            "corrected_replicate_standard_error": float(
                np.std(native_replicates, ddof=1) / math.sqrt(len(native_replicates))
            ),
            "plain_mean": native_plain_mean,
            "plain_absolute_error": abs(native_plain_mean - native_reference),
            "paths_per_replicate": native_paths,
            "replicates": len(native_replicates),
            "mean_variance_reduction": float(
                np.mean([item["variance_reduction"] for item in native_metadata])
            ),
        },
        "training_seconds": phase["cumulative_training_seconds"],
        "evaluation_seconds": evaluation_seconds,
        "evaluation_process_seconds_including_diagnostics": evaluation_process_seconds,
        "diagnostic_and_io_overhead_seconds": max(
            evaluation_process_seconds - evaluation_seconds, 0.0
        ),
        "total_seconds": float(phase["cumulative_training_seconds"])
        + evaluation_seconds,
        "total_seconds_including_diagnostics": float(
            phase["cumulative_training_seconds"]
        )
        + evaluation_process_seconds,
        "prediction_file": str(prediction_path),
        "ground_truth_used_for_training_or_selection": False,
        "frozen_before_test_reference": True,
    }
    output = RAW_DIR / (
        f"actor_evaluation_{benchmark_name}_{trial}_s{seed}_a{iteration}_n{time_steps}.json"
    )
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARK_NAMES, required=True)
    parser.add_argument("--trial", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--iteration", type=int, required=True)
    parser.add_argument("--paths", type=int, default=256)
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--native-paths", type=int, default=8192)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = evaluate(
        args.benchmark,
        args.trial,
        args.seed,
        args.iteration,
        torch.device(args.device),
        paths=args.paths,
        time_steps=args.steps,
        replicates=args.replicates,
        native_paths=args.native_paths,
    )
    print(
        json.dumps(
            {
                "status": "complete",
                "rmse": result["corrected_metrics_of_replicate_mean"]["rmse"],
                "native_absolute_error": result["native"]["corrected_absolute_error"],
            }
        )
    )


if __name__ == "__main__":
    main()

