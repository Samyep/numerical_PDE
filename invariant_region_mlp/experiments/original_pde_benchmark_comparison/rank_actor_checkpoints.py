"""High-precision reference-free ranking of frozen Actor checkpoints."""

from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .benchmarks import BENCHMARK_NAMES, make_benchmark
from .common import RAW_DIR, RESULTS_DIR, atomic_json, stable_seed
from .models import load_actor
from .run_actor import _policy_bank
from .tangent_solver import evaluate_policy


def rank(
    benchmark_name: str,
    seed: int,
    device: torch.device,
    *,
    paths: int,
    time_steps: int,
    replicates: int,
) -> dict[str, Any]:
    benchmark = make_benchmark(benchmark_name)
    t, x = _policy_bank(benchmark)
    candidates: list[dict[str, Any]] = []
    pattern = str(RAW_DIR / f"actor_{benchmark_name}_*_s{seed}.json")
    for filename in sorted(glob.glob(pattern)):
        record = json.loads(Path(filename).read_text())
        label_estimator = str(record["configuration"]["label_estimator"])
        for phase in record["phases"]:
            candidates.append(
                {
                    "trial": record["trial"],
                    "iteration": int(phase["iteration"]),
                    "checkpoint": phase["checkpoint"],
                    "training_seconds": float(phase["cumulative_training_seconds"]),
                    "heldout_label_rmse": float(
                        phase["training"]["selected_validation_rmse"]
                    ),
                    "label_replication_rmse": float(
                        phase["labels"]["validation_label_replication_rmse"]
                    ),
                    "label_estimator": label_estimator,
                    "eligible_for_primary_selection": label_estimator == "mcv",
                }
            )
    if not candidates:
        raise RuntimeError(f"no Actor candidates found for {benchmark_name}")

    started = time.perf_counter()
    for candidate in candidates:
        actor, _ = load_actor(candidate["checkpoint"], benchmark, device)
        predictions: list[np.ndarray] = []
        metadata: list[dict[str, Any]] = []
        for replicate in range(replicates):
            result, meta = evaluate_policy(
                benchmark,
                actor,
                t,
                x,
                paths=paths,
                time_steps=time_steps,
                seed=stable_seed(
                    "original_benchmark_final_reference_free_ranking",
                    benchmark_name,
                    replicate,
                ),
                device=device,
                parent_batch=8,
                path_chunk=64,
            )
            key = "mcv" if benchmark.actor_kind == "gradient" else "plain"
            predictions.append(result[key])
            metadata.append(meta)
        matrix = np.asarray(predictions)
        state_prediction = matrix.mean(axis=0)
        candidate["policy_value"] = float(np.mean(state_prediction))
        candidate["state_standard_error"] = float(
            np.std(state_prediction, ddof=1) / math.sqrt(len(state_prediction))
        )
        replicate_means = matrix.mean(axis=1)
        candidate["mc_standard_error"] = float(
            np.std(replicate_means, ddof=1) / math.sqrt(replicates)
        ) if replicates > 1 else None
        candidate["mean_variance_reduction"] = float(
            np.mean([item["variance_reduction"] for item in metadata])
        )
        candidate["mean_weight_ess_fraction"] = float(
            np.mean([item["weight_ess_fraction"] for item in metadata])
        )
        candidate["nonfinite_count"] = int(
            sum(item["nonfinite_count"] for item in metadata)
        )
        candidate["_state_prediction"] = state_prediction

    candidates.sort(key=lambda row: (row["policy_value"], row["heldout_label_rmse"]))
    eligible = [row for row in candidates if row["eligible_for_primary_selection"]]
    if not eligible:
        raise RuntimeError("no MCV-training candidate is eligible for primary selection")
    selected_candidate = min(
        eligible,
        key=lambda row: (row["policy_value"], row["heldout_label_rmse"]),
    )
    best = selected_candidate["_state_prediction"]
    for candidate in candidates:
        difference = candidate.pop("_state_prediction") - best
        candidate["candidate_minus_best_policy_value"] = float(np.mean(difference))
        candidate["paired_state_standard_error_vs_best"] = float(
            np.std(difference, ddof=1) / math.sqrt(len(difference))
        )
        candidate["ground_truth_used"] = False
        if not candidate["eligible_for_primary_selection"]:
            candidate["selection_role"] = "plain_FK_training_ablation_only"
        elif candidate is selected_candidate:
            candidate["selection_role"] = "primary_MCV_training_selected"
        else:
            candidate["selection_role"] = "primary_MCV_training_candidate"

    result = {
        "benchmark": benchmark_name,
        "seed": seed,
        "stage": "post-training frozen reference-free checkpoint ranking",
        "criterion": (
            "minimum replicated common-bank tangent-policy value among "
            "MCV-training candidates; plain-FK training remains an ablation"
        ),
        "states": len(t),
        "paths_per_replicate": paths,
        "time_steps": time_steps,
        "replicates": replicates,
        "common_random_numbers_across_candidates": True,
        "ranking": candidates,
        "selected": {
            key: selected_candidate[key]
            for key in (
                "trial",
                "iteration",
                "checkpoint",
                "training_seconds",
                "label_estimator",
            )
        },
        "ranking_seconds": time.perf_counter() - started,
        "ground_truth_used": False,
        "final_test_opened": False,
    }
    atomic_json(RESULTS_DIR / f"final_actor_selection_{benchmark_name}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARK_NAMES, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--paths", type=int, default=512)
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = rank(
        args.benchmark,
        args.seed,
        torch.device(args.device),
        paths=args.paths,
        time_steps=args.steps,
        replicates=args.replicates,
    )
    print(json.dumps(result["selected"], sort_keys=True))


if __name__ == "__main__":
    main()
