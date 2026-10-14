"""High-accuracy evaluation of frozen Actor candidates.

This module never trains or changes an Actor.  It re-evaluates the five
predeclared final-profile candidates with common random numbers at larger MC
or time-discretization budgets, keeping those two effects separate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import torch

from .actor_method import evaluate_policy_cv, load_actor
from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    file_sha256,
    model_metrics,
    require_cuda,
    stable_seed,
)
from .equations import Equation, equation_key, make_equation
from .run_experiment import _banks


EVALUATION_BUDGETS: tuple[dict[str, int | str], ...] = (
    {
        "name": "paths1024_steps64",
        "paths": 1024,
        "time_steps": 64,
    },
    {
        "name": "paths256_steps128",
        "paths": 256,
        "time_steps": 128,
    },
)


def _git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unavailable"


def _evaluate_split(
    equation: Equation,
    actor: torch.nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Evaluate 2,048 fixed points while retaining N=960 timing."""

    if len(t) != 2048 or len(x) != 2048:
        raise ValueError("sensitivity evaluation expects the frozen 2,048-point bank")
    common = {
        "paths": paths,
        "time_steps": time_steps,
        "device": device,
        "path_chunk": min(64, paths),
    }
    first_value, first_variance, first_meta = evaluate_policy_cv(
        equation,
        actor,
        t[:960],
        x[:960],
        seed=seed,
        **common,
    )
    second_value, second_variance, second_meta = evaluate_policy_cv(
        equation,
        actor,
        t[960:],
        x[960:],
        seed=stable_seed(seed, "tail_1088"),
        **common,
    )
    return (
        np.concatenate((first_value, second_value)),
        np.concatenate((first_variance, second_variance)),
        {
            **first_meta,
            "points": 2048,
            "wall_seconds_960": first_meta["wall_seconds"],
            "wall_seconds_2048": first_meta["wall_seconds"]
            + second_meta["wall_seconds"],
            "peak_gpu_bytes": max(
                int(first_meta["peak_gpu_bytes"] or 0),
                int(second_meta["peak_gpu_bytes"] or 0),
            ),
            "split_evaluation": [first_meta, second_meta],
        },
    )


def run(seed_index: int) -> None:
    equation = make_equation("hjb_log")
    key = equation_key(equation)
    device = require_cuda()
    run_path = RESULTS_DIR / "runs" / key / "ours" / f"seed{seed_index}_final.json"
    if not run_path.is_file():
        raise FileNotFoundError(run_path)
    source_run = json.loads(run_path.read_text())
    if source_run.get("profile") != "final":
        raise ValueError(f"expected final-profile source, got {source_run.get('profile')}")
    actor_names = [entry["actor"] for entry in source_run["evaluations"]]
    if actor_names != ["q0", "q1", "q1_refresh", "q2", "q3"]:
        raise ValueError(f"unexpected Actor candidates: {actor_names}")

    banks = _banks(equation)
    checkpoint_root = CHECKPOINT_DIR / key / "ours" / f"seed{seed_index}" / "final"
    prediction_root = (
        PREDICTION_DIR / key / "ours" / f"seed{seed_index}" / "evaluation_sensitivity"
    )
    evaluations: list[dict[str, Any]] = []
    started = time.perf_counter()
    for budget in EVALUATION_BUDGETS:
        budget_name = str(budget["name"])
        paths = int(budget["paths"])
        time_steps = int(budget["time_steps"])
        # Every candidate in this seed/budget sees the same Brownian stream.
        evaluation_seed = stable_seed(
            "actor_sensitivity_common_random_numbers",
            key,
            seed_index,
            budget_name,
        )
        for actor_name in actor_names:
            checkpoint_path = checkpoint_root / f"{actor_name}.pt"
            if not checkpoint_path.is_file():
                raise FileNotFoundError(checkpoint_path)
            actor, checkpoint_metadata = load_actor(checkpoint_path, device)
            value, estimator_variance, metadata = _evaluate_split(
                equation,
                actor,
                banks["test_t"],
                banks["test_x"],
                paths=paths,
                time_steps=time_steps,
                seed=evaluation_seed,
                device=device,
            )
            metrics = model_metrics(value, banks["test_u"])
            rms_mc_standard_error = float(
                np.sqrt(np.mean(np.maximum(estimator_variance, 0.0)))
            )
            mean_estimator_variance = float(
                np.mean(np.maximum(estimator_variance, 0.0))
            )
            observed_mse = float(metrics["rmse"] ** 2)
            prediction_path = prediction_root / f"{actor_name}_{budget_name}.npz"
            atomic_npz(
                prediction_path,
                value=value,
                estimator_variance=estimator_variance,
                truth_u=banks["test_u"],
                test_points_sha256=np.asarray(array_sha256(banks["test_x"])),
            )
            evaluations.append(
                {
                    "actor": actor_name,
                    "budget": budget_name,
                    "paths": paths,
                    "time_steps": time_steps,
                    "common_random_numbers_seed": evaluation_seed,
                    "checkpoint_path": str(checkpoint_path),
                    "checkpoint_sha256": file_sha256(checkpoint_path),
                    "path_training_seconds": checkpoint_metadata[
                        "path_training_seconds"
                    ],
                    "inference_seconds_960": metadata["wall_seconds_960"],
                    "inference_seconds_2048": metadata["wall_seconds_2048"],
                    "peak_gpu_bytes": metadata["peak_gpu_bytes"],
                    "mean_mc_standard_error": float(
                        np.mean(np.sqrt(np.maximum(estimator_variance, 0.0)))
                    ),
                    "rms_mc_standard_error": rms_mc_standard_error,
                    "mean_estimator_variance": mean_estimator_variance,
                    "mc_variance_fraction_of_observed_mse": (
                        mean_estimator_variance / observed_mse
                        if observed_mse > 0.0
                        else 0.0
                    ),
                    "noise_corrected_rmse_diagnostic": float(
                        np.sqrt(max(observed_mse - mean_estimator_variance, 0.0))
                    ),
                    "prediction_path": str(prediction_path),
                    **metrics,
                }
            )
            del actor
            torch.cuda.empty_cache()

    selected_by_budget: dict[str, str] = {}
    for budget in EVALUATION_BUDGETS:
        name = str(budget["name"])
        eligible = [row for row in evaluations if row["budget"] == name]
        selected_by_budget[name] = min(eligible, key=lambda row: row["rmse"])["actor"]
    payload = {
        "method": "ours_actor_iteration_cv_evaluation_sensitivity",
        "equation": key,
        "seed_index": seed_index,
        "source_profile": "final",
        "source_run_path": str(run_path),
        "source_run_sha256": file_sha256(run_path),
        "source_selected_actor_at_256_paths_64_steps": source_run["selected_actor"],
        "gpu_name": torch.cuda.get_device_name(device),
        "git_commit": _git_head(),
        "frozen_actor_weights": True,
        "training_performed": False,
        "test_used_for_actor_ranking": True,
        "candidate_common_random_numbers": True,
        "selected_actor_by_budget": selected_by_budget,
        "evaluations": evaluations,
        "job_wall_seconds": time.perf_counter() - started,
    }
    output_path = (
        RESULTS_DIR / "evaluation_sensitivity" / f"hjb_log_seed{seed_index}.json"
    )
    atomic_json(output_path, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, choices=(0, 1, 2), required=True)
    arguments = parser.parse_args()
    run(arguments.seed)


if __name__ == "__main__":
    main()
