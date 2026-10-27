"""Adapter for the frozen mature joint Tangent-Actor implementation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import torch

from invariant_region_mlp.experiments.main_highdim_expansion import (
    run_dense_actor as _frozen,
)

from .benchmark import make_benchmark, sample_states
from .protocol import (
    ACTOR_UPDATES,
    CACHE_DIR,
    CHECKPOINT_DIR,
    JOINT_ACTOR_CONFIG,
    PDE,
    PREDICTION_DIR,
    PREFLIGHT_AUDIT,
    PRIMARY_TA_TRAINING_SEEDS,
    RAW_DIR,
    REVISION,
    TA_EVALUATION_PATHS,
    TA_INFERENCE_REPLICATES,
    TA_PRIMARY_EULER_STEPS,
    ensure_directories,
    file_sha256,
    load_query_bank,
    query_path,
    stable_seed,
    write_json_atomic,
)


def _run_id(dimension: int, seed: int) -> str:
    return f"{PDE}_d{int(dimension)}_mature_s{int(seed)}"


def _states(
    dimension: int, count: int, *, split: str, seed: int
):
    return sample_states(
        dimension,
        count,
        seed=stable_seed("ta_states", dimension, split, seed),
    )


def _load_bank(_pde: str, dimension: int):
    t, x, truth, truth_z = load_query_bank(dimension)
    return t, x, truth, truth_z, query_path(dimension) / "bank.npz"


def _configure() -> None:
    _frozen.make_dense_fourier = make_benchmark
    _frozen._states = _states
    _frozen._run_id = _run_id
    _frozen.load_test_bank = _load_bank
    _frozen.file_sha256 = file_sha256
    _frozen.ACTORS_UPDATES = ACTOR_UPDATES
    _frozen.ACTOR_UPDATES = ACTOR_UPDATES
    _frozen.CACHE_DIR = CACHE_DIR
    _frozen.CHECKPOINT_DIR = CHECKPOINT_DIR
    _frozen.EVALUATION_PATHS = TA_EVALUATION_PATHS
    _frozen.INFERENCE_REPLICATES = TA_INFERENCE_REPLICATES
    _frozen.JOINT_ACTOR_CONFIG = JOINT_ACTOR_CONFIG
    _frozen.PREDICTION_DIR = PREDICTION_DIR
    _frozen.PRIMARY_EULER_STEPS = TA_PRIMARY_EULER_STEPS
    _frozen.RAW_DIR = RAW_DIR
    _frozen.REVISION = REVISION
    _frozen.TRAINING_SEED = PRIMARY_TA_TRAINING_SEEDS[0]
    _frozen.VALIDATION_PATH = PREFLIGHT_AUDIT
    _frozen.ensure_directories = ensure_directories
    _frozen.stable_seed = stable_seed
    _frozen.write_json_atomic = write_json_atomic


def _normalize_training(result: dict) -> dict:
    result["study_revision"] = REVISION
    result["pde"] = PDE
    result["method"] = "terminal_aware_joint_tangent_actor_mcv"
    result["source_solver_module"] = (
        "invariant_region_mlp.experiments.main_highdim_expansion.run_dense_actor"
    )
    result["adapter_only"] = True
    path = RAW_DIR / "actor" / f"training_{result['run_id']}.json"
    write_json_atomic(path, result)
    return result


def train(dimension: int, seed: int, device: torch.device) -> dict:
    _configure()
    return _normalize_training(_frozen.train(int(dimension), int(seed), device))


def evaluate(
    dimension: int,
    seed: int,
    device: torch.device,
    *,
    actor_names: tuple[str, ...] = ("selected", "a0"),
    steps: tuple[int, ...] = (TA_PRIMARY_EULER_STEPS,),
) -> list[dict]:
    _configure()
    rows = _frozen.evaluate(
        int(dimension),
        int(seed),
        device,
        actor_names=actor_names,
        steps=steps,
    )
    normalized: list[dict] = []
    run_id = _run_id(dimension, seed)
    for row in rows:
        row["study_revision"] = REVISION
        row["pde"] = PDE
        row["source_solver_module"] = (
            "invariant_region_mlp.experiments.main_highdim_expansion.run_dense_actor"
        )
        output = (
            RAW_DIR
            / "actor_evaluations"
            / run_id
            / row["actor_request"]
            / f"steps{row['steps']}_rep{row['inference_replicate']:02d}.json"
        )
        write_json_atomic(output, row)
        normalized.append(row)
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--train-only", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--actors", default="selected,a0")
    parser.add_argument("--steps", default=str(TA_PRIMARY_EULER_STEPS))
    args = parser.parse_args()
    if args.train_only and args.evaluate_only:
        raise ValueError("choose at most one of --train-only and --evaluate-only")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    if not args.evaluate_only:
        result = train(args.dimension, args.seed, device)
        print(json.dumps({
            "training": result["status"],
            "selected": result["selection"]["selected_actor"],
        }), flush=True)
    if not args.train_only:
        evaluate(
            args.dimension,
            args.seed,
            device,
            actor_names=tuple(value for value in args.actors.split(",") if value),
            steps=tuple(int(value) for value in args.steps.split(",") if value),
        )


if __name__ == "__main__":
    main()


__all__ = ["evaluate", "train"]
