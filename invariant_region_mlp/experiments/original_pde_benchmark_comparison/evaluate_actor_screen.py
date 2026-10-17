"""Descriptive final-test screen run only after Actor selection is frozen."""

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
from .evaluate_actor import _select_reference_z
from .models import load_actor
from .tangent_solver import evaluate_policy


def _metrics(prediction: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    error = np.asarray(prediction) - np.asarray(reference)
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
        "bias": float(np.mean(error)),
    }


def evaluate_screen(
    benchmark_name: str,
    device: torch.device,
    *,
    paths: int,
    time_steps: int,
) -> dict[str, Any]:
    selection_path = RESULTS_DIR / f"final_actor_selection_{benchmark_name}.json"
    if not selection_path.exists():
        raise RuntimeError("freeze final_actor_selection before opening test references")
    selection = json.loads(selection_path.read_text())
    if not selection.get("final_test_opened") is False:
        raise RuntimeError("selection artifact lacks the pre-test freeze declaration")
    benchmark = make_benchmark(benchmark_name)
    with np.load(RESULTS_DIR / f"evaluation_bank_{benchmark_name}.npz") as bank:
        indices = bank["test_indices"]
        t = bank["t"][indices]
        x = bank["x"][indices]
        reference = bank["value"][indices]
        reference_z = _select_reference_z(bank, indices, benchmark.actor_kind)

    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    pattern = str(RAW_DIR / f"actor_{benchmark_name}_*_s0.json")
    for filename in sorted(glob.glob(pattern)):
        training = json.loads(Path(filename).read_text())
        for phase in training["phases"]:
            actor, _ = load_actor(phase["checkpoint"], benchmark, device)
            result, metadata = evaluate_policy(
                benchmark,
                actor,
                t,
                x,
                paths=paths,
                time_steps=time_steps,
                seed=stable_seed(
                    "original_benchmark_postselection_screen",
                    benchmark_name,
                    time_steps,
                ),
                device=device,
                parent_batch=8,
                path_chunk=64,
            )
            key = "mcv" if benchmark.actor_kind == "gradient" else "plain"
            tt = torch.as_tensor(t, dtype=torch.float32, device=device)
            xx = torch.as_tensor(x, dtype=torch.float32, device=device)
            with torch.no_grad():
                r, q = actor.predict(tt, xx)
            direct = q if q is not None else r
            if direct is None:
                raise RuntimeError("Actor returned no direct output")
            if q is not None and reference_z is None:
                raise RuntimeError("Gradient Actor screen has no gradient reference")
            target = reference_z if q is not None else reference
            direct_np = direct.detach().double().cpu().numpy()
            direct_error = direct_np - target
            rows.append(
                {
                    "trial": training["trial"],
                    "iteration": int(phase["iteration"]),
                    "checkpoint": phase["checkpoint"],
                    "training_seconds": float(phase["cumulative_training_seconds"]),
                    "configuration": training["configuration"],
                    "value_metrics": _metrics(result[key], reference),
                    "plain_value_metrics": _metrics(result["plain"], reference),
                    "direct_actor_target": benchmark.actor_kind,
                    "direct_actor_component_rmse": float(
                        np.sqrt(np.mean(np.square(direct_error)))
                    ),
                    "direct_actor_vector_rmse": float(
                        np.sqrt(np.mean(np.square(direct_error).sum(axis=-1)))
                    ) if direct_error.ndim == 2 else None,
                    "mc_standard_error_mean": float(
                        metadata["mean_mcv_mc_se"]
                        if benchmark.actor_kind == "gradient"
                        else metadata["mean_plain_mc_se"]
                    ),
                    "variance_reduction": float(metadata["variance_reduction"]),
                    "weight_ess_fraction": float(metadata["weight_ess_fraction"]),
                    "evaluation_seconds": float(metadata["wall_seconds"]),
                    "nonfinite_count": int(metadata["nonfinite_count"]),
                }
            )
    result = {
        "benchmark": benchmark_name,
        "purpose": "descriptive test-error envelope after frozen reference-free selection",
        "selected_before_test": selection["selected"],
        "states": len(reference),
        "paths": paths,
        "time_steps": time_steps,
        "common_random_numbers_across_candidates": True,
        "rows": rows,
        "wall_seconds": time.perf_counter() - started,
        "used_for_model_selection": False,
        "final_test_used_descriptively": True,
    }
    atomic_json(RESULTS_DIR / f"actor_postselection_screen_{benchmark_name}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARK_NAMES, required=True)
    parser.add_argument("--paths", type=int, default=128)
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = evaluate_screen(
        args.benchmark,
        torch.device(args.device),
        paths=args.paths,
        time_steps=args.steps,
    )
    best = min(result["rows"], key=lambda row: row["value_metrics"]["rmse"])
    print(
        json.dumps(
            {
                "status": "complete",
                "descriptive_best_trial": best["trial"],
                "descriptive_best_iteration": best["iteration"],
                "descriptive_best_rmse": best["value_metrics"]["rmse"],
                "selection_unchanged": True,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
