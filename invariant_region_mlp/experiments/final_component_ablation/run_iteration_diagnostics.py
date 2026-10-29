"""Evaluate every frozen Full-TA actor iterate on the formal test protocol."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    evaluate_policy as evaluate_original,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy as evaluate_joint,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor import (
    make_driver,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor import (
    make_benchmark,
)

from .protocol import (
    DIMENSION,
    EVALUATION_PATHS,
    EVALUATION_STEPS,
    PDES,
    PREDICTION_DIR,
    RAW_DIR,
    TRAINING_SEEDS,
    UPDATES,
    ensure_dirs,
    stable_seed,
)
from .run_ablation import (
    _array_hash,
    _file_hash,
    _initial_actor,
    _is_hjb,
    _load_local_actor,
    _load_reused_actor,
    _local_checkpoint,
    _test_bank,
    _truth,
    _write_json,
)


def _training_record(pde: str, seed: int) -> dict[str, Any]:
    path = RAW_DIR / f"training_{pde}_d100_a0_full_s{seed}.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "complete":
        raise RuntimeError(f"incomplete Full-TA training record: {path}")
    return payload


def _diagnostic_path(pde: str, seed: int, iteration: int) -> Path:
    return RAW_DIR / (
        f"iteration_diagnostic_{pde}_d100_a0_full_s{seed}_a{iteration}_rep0.json"
    )


def _prediction_path(pde: str, seed: int, iteration: int) -> Path:
    return PREDICTION_DIR / (
        f"iteration_diagnostic_{pde}_d100_a0_full_s{seed}_a{iteration}_rep0.npz"
    )


def _checkpoint_path(pde: str, seed: int, iteration: int) -> Path | None:
    if iteration == 0:
        return None
    local = _local_checkpoint(pde, "a0_full", seed, iteration)
    if local.exists():
        return local
    # Reused source checkpoints are resolved and validated by the shared loader.
    from .run_ablation import _reuse_checkpoint

    reused = _reuse_checkpoint(pde, seed, iteration)
    if reused is None or not reused.exists():
        raise FileNotFoundError(
            f"missing Full-TA checkpoint for {pde} seed={seed} a{iteration}"
        )
    return reused


def _load_actor(
    pde: str,
    seed: int,
    iteration: int,
    problem: Any,
    device: torch.device,
) -> tuple[torch.nn.Module, Path | None]:
    checkpoint = _checkpoint_path(pde, seed, iteration)
    if iteration == 0:
        return _initial_actor(pde, "a0_full", seed, problem, device), None
    local = _local_checkpoint(pde, "a0_full", seed, iteration)
    if local.exists():
        return (
            _load_local_actor(local, pde, "a0_full", problem, device),
            checkpoint,
        )
    return _load_reused_actor(pde, seed, iteration, problem, device), checkpoint


def _reuse_selected_result(
    pde: str,
    seed: int,
    iteration: int,
    training: dict[str, Any],
) -> dict[str, Any] | None:
    if iteration != int(training["selected_iteration"]):
        return None
    source = RAW_DIR / f"evaluation_{pde}_d100_a0_full_s{seed}_rep0.json"
    if not source.exists():
        return None
    payload = json.loads(source.read_text(encoding="utf-8"))
    checkpoint = _checkpoint_path(pde, seed, iteration)
    result = {
        "pde": pde,
        "seed": seed,
        "iteration": iteration,
        "selected_iteration": int(training["selected_iteration"]),
        "is_selected_full_actor": True,
        "estimator": "mcv",
        "states": int(payload["states"]),
        "steps": int(payload["steps"]),
        "paths": int(payload["paths"]),
        "value_rmse": float(payload["mcv_rmse_same_paths"]),
        "rmae": float(payload["rmae"]),
        "bias": float(payload["bias"]),
        "mc_rmse_contribution": float(payload["mc_rmse_contribution"]),
        "weight_ess_fraction": float(payload["weight_ess_fraction"]),
        "nonfinite_count": int(payload["nonfinite_count"]),
        "inference_seconds": float(payload["inference_seconds"]),
        "test_bank": payload["test_bank"],
        "test_bank_array_sha256": payload["test_bank_array_sha256"],
        "reference_array_sha256": payload["reference_array_sha256"],
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": _file_hash(checkpoint),
        "reused_formal_primary_evaluation": True,
        "source_evaluation": str(source),
        "prediction_path": payload["prediction_path"],
        "common_brownian_seed": stable_seed(
            "component_ablation_evaluation", pde, seed, 0
        ),
        "status": "complete",
    }
    _write_json(_diagnostic_path(pde, seed, iteration), result)
    return result


def evaluate_iteration(
    pde: str,
    seed: int,
    iteration: int,
    device: torch.device,
    *,
    force: bool = False,
) -> dict[str, Any]:
    ensure_dirs()
    destination = _diagnostic_path(pde, seed, iteration)
    if destination.exists() and not force:
        payload = json.loads(destination.read_text(encoding="utf-8"))
        if payload.get("status") == "complete":
            print(f"ITERATION_SKIP pde={pde} seed={seed} a{iteration}", flush=True)
            return payload

    training = _training_record(pde, seed)
    reused = _reuse_selected_result(pde, seed, iteration, training)
    if reused is not None and not force:
        print(
            f"ITERATION_REUSE pde={pde} seed={seed} a{iteration} "
            f"rmse={reused['value_rmse']:.6e}",
            flush=True,
        )
        return reused

    problem = (
        make_benchmark(pde, DIMENSION)
        if _is_hjb(pde)
        else make_driver(pde, DIMENSION)
    )
    actor, checkpoint = _load_actor(pde, seed, iteration, problem, device)
    t, x, bank_path = _test_bank(pde)
    if len(t) != 960:
        raise RuntimeError(f"{pde} test bank has {len(t)} states, expected 960")
    truth, _, _ = _truth(pde, problem, t, x)
    evaluation_seed = stable_seed(
        "component_ablation_evaluation", pde, seed, 0
    )
    started = time.perf_counter()
    if _is_hjb(pde):
        output, metadata = evaluate_original(
            problem,
            actor,
            t,
            x,
            paths=EVALUATION_PATHS,
            time_steps=EVALUATION_STEPS,
            seed=evaluation_seed,
            device=device,
            parent_batch=8,
            path_chunk=64,
        )
        prediction = output["mcv"]
    else:
        output, metadata = evaluate_joint(
            problem,
            actor,
            t,
            x,
            paths=EVALUATION_PATHS,
            time_steps=EVALUATION_STEPS,
            seed=evaluation_seed,
            device=device,
            parent_batch=8,
            path_chunk=64,
        )
        prediction = problem.working_to_original(output["mcv"])
    error = prediction - truth
    prediction_path = _prediction_path(pde, seed, iteration)
    np.savez_compressed(
        prediction_path,
        truth=truth,
        prediction=prediction,
        pair_variance=output["mcv_pair_variance"],
    )
    result = {
        "pde": pde,
        "seed": seed,
        "iteration": iteration,
        "selected_iteration": int(training["selected_iteration"]),
        "is_selected_full_actor": iteration == int(training["selected_iteration"]),
        "estimator": "mcv",
        "states": len(t),
        "steps": EVALUATION_STEPS,
        "paths": EVALUATION_PATHS,
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "rmae": float(np.sum(np.abs(error)) / np.sum(np.abs(truth))),
        "bias": float(np.mean(error)),
        "mc_rmse_contribution": float(
            np.sqrt(
                np.mean(
                    output["mcv_pair_variance"] / (EVALUATION_PATHS // 2)
                )
            )
        ),
        "weight_ess_fraction": float(metadata["weight_ess_fraction"]),
        "nonfinite_count": int(metadata["nonfinite_count"]),
        "inference_seconds": float(metadata["wall_seconds"]),
        "invocation_seconds": time.perf_counter() - started,
        "test_bank": str(bank_path),
        "test_bank_file_sha256": _file_hash(bank_path),
        "test_bank_array_sha256": _array_hash(t, x),
        "reference_array_sha256": _array_hash(truth),
        "checkpoint": str(checkpoint) if checkpoint is not None else None,
        "checkpoint_sha256": (
            _file_hash(checkpoint) if checkpoint is not None else None
        ),
        "reused_formal_primary_evaluation": False,
        "prediction_path": str(prediction_path),
        "prediction_file_sha256": _file_hash(prediction_path),
        "common_brownian_seed": evaluation_seed,
        "status": "complete",
    }
    _write_json(destination, result)
    print(
        f"ITERATION_DONE pde={pde} seed={seed} a{iteration} "
        f"rmse={result['value_rmse']:.6e} seconds={result['inference_seconds']:.1f}",
        flush=True,
    )
    del actor
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return result


def audit(pdes: list[str], seeds: list[int], iterations: list[int]) -> None:
    rows = []
    for pde in pdes:
        for seed in seeds:
            training = _training_record(pde, seed)
            for iteration in iterations:
                checkpoint = _checkpoint_path(pde, seed, iteration)
                rows.append(
                    {
                        "pde": pde,
                        "seed": seed,
                        "iteration": iteration,
                        "selected_iteration": int(training["selected_iteration"]),
                        "checkpoint": str(checkpoint) if checkpoint else None,
                        "checkpoint_exists": checkpoint is None or checkpoint.exists(),
                        "output_exists": _diagnostic_path(
                            pde, seed, iteration
                        ).exists(),
                    }
                )
    missing = [row for row in rows if not row["checkpoint_exists"]]
    print(json.dumps({"rows": rows, "missing": missing}, indent=2))
    if missing:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pde", choices=PDES, action="append")
    parser.add_argument("--seed", type=int, choices=TRAINING_SEEDS, action="append")
    parser.add_argument("--iterations", default="0,1,2,3,4")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    pdes = args.pde or list(PDES)
    seeds = args.seed or list(TRAINING_SEEDS)
    iterations = [int(token) for token in args.iterations.split(",")]
    if any(iteration < 0 or iteration > UPDATES for iteration in iterations):
        parser.error(f"iterations must lie in [0, {UPDATES}]")
    if args.audit_only:
        audit(pdes, seeds, iterations)
        return
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("iteration diagnostics require a CUDA allocation")
    for pde in pdes:
        for seed in seeds:
            for iteration in iterations:
                evaluate_iteration(
                    pde, seed, iteration, device, force=args.force
                )


if __name__ == "__main__":
    main()
