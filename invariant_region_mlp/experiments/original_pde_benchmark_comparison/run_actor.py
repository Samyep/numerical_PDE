"""Staged training CLI for our Terminal-Aware Tangent-Actor method."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .benchmarks import BENCHMARK_NAMES, OriginalBenchmark, make_benchmark
from .common import (
    CHECKPOINT_DIR,
    RAW_DIR,
    atomic_json,
    atomic_torch_save,
    ensure_dirs,
    stable_seed,
    synchronize,
)
from .models import TerminalAwareActor, load_actor
from .tangent_solver import ActorConfig, evaluate_policy, fit_actor, generate_labels


TRIALS: dict[str, ActorConfig] = {
    "pilot192": ActorConfig(
        training_points=2048, validation_points=256, paths_per_label=64,
        label_time_steps=16, optimizer_steps=2000, width=192, depth=3,
        activation="silu", learning_rate=2e-4, batch_size=512,
    ),
    "wide384": ActorConfig(
        training_points=4096, validation_points=256, paths_per_label=64,
        label_time_steps=24, optimizer_steps=3000, width=384, depth=4,
        activation="silu", learning_rate=2e-4, batch_size=512,
    ),
    "deep_elu": ActorConfig(
        training_points=4096, validation_points=256, paths_per_label=64,
        label_time_steps=24, optimizer_steps=3000, width=256, depth=5,
        activation="elu", learning_rate=2e-4, batch_size=512,
    ),
    "data8192": ActorConfig(
        training_points=8192, validation_points=384, paths_per_label=64,
        label_time_steps=32, optimizer_steps=4000, width=256, depth=4,
        activation="silu", learning_rate=2e-4, batch_size=512,
    ),
    "quality128": ActorConfig(
        training_points=4096, validation_points=384, paths_per_label=128,
        label_time_steps=32, optimizer_steps=4000, width=256, depth=4,
        activation="silu", learning_rate=2e-4, batch_size=512,
    ),
    "focused": ActorConfig(
        training_points=8192, validation_points=512, paths_per_label=128,
        label_time_steps=32, optimizer_steps=5000, width=384, depth=4,
        activation="silu", learning_rate=1.5e-4, batch_size=512,
    ),
    "focused_plain": ActorConfig(
        training_points=8192, validation_points=512, paths_per_label=128,
        label_time_steps=32, optimizer_steps=5000, width=384, depth=4,
        activation="silu", learning_rate=1.5e-4, batch_size=512,
        label_estimator="plain",
    ),
}


def _distribution(benchmark: OriginalBenchmark) -> str:
    if benchmark.name == "hjb_log_d100":
        return "mixed"
    if benchmark.name == "dpi_ou_gmm_d100":
        return "official"
    return "native"


def _state_bank(
    benchmark: OriginalBenchmark,
    trial: str,
    seed: int,
    iteration: int,
    config: ActorConfig,
) -> dict[str, np.ndarray]:
    distribution = config.training_distribution
    train_t, train_x = benchmark.sample_states(
        config.training_points,
        seed=stable_seed("original_benchmark_train_states", benchmark.name, trial, seed, iteration),
        purpose=distribution,
    )
    # Fix validation locations across Actor rounds.  Their labels still use an
    # independent stream and change with the frozen source policy.
    valid_t, valid_x = benchmark.sample_states(
        config.validation_points,
        seed=stable_seed("original_benchmark_validation_states", benchmark.name, trial, seed),
        purpose=distribution,
    )
    return {
        "train_t": train_t,
        "train_x": train_x,
        "valid_t": valid_t,
        "valid_x": valid_x,
    }


def _phase_labels(
    benchmark: OriginalBenchmark,
    trial: str,
    seed: int,
    iteration: int,
    config: ActorConfig,
    source: TerminalAwareActor,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    bank = _state_bank(benchmark, trial, seed, iteration, config)
    common = {
        "paths": config.paths_per_label,
        "time_steps": config.label_time_steps,
        "device": device,
        "parent_batch": config.parent_batch,
        "path_chunk": config.path_chunk,
        "estimator": config.label_estimator,
    }
    train, train_meta = generate_labels(
        benchmark,
        source,
        bank["train_t"],
        bank["train_x"],
        seed=stable_seed("original_benchmark_labels", benchmark.name, trial, seed, iteration, "train"),
        **common,
    )
    valid, valid_meta = generate_labels(
        benchmark,
        source,
        bank["valid_t"],
        bank["valid_x"],
        seed=stable_seed("original_benchmark_labels", benchmark.name, trial, seed, iteration, "valid"),
        **common,
    )
    replicate, replicate_meta = generate_labels(
        benchmark,
        source,
        bank["valid_t"],
        bank["valid_x"],
        seed=stable_seed("original_benchmark_labels", benchmark.name, trial, seed, iteration, "replicate"),
        **common,
    )
    difference = valid - replicate
    arrays = {**bank, "train_labels": train, "valid_labels": valid}
    return arrays, {
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "validation_label_replication_rmse": float(
            np.sqrt(np.mean(np.square(difference)))
        ),
        "validation_label_replication_vector_rmse": float(
            np.sqrt(np.mean(np.square(difference).sum(axis=1)))
        )
        if difference.ndim == 2
        else None,
        "wall_seconds": float(
            train_meta["wall_seconds"]
            + valid_meta["wall_seconds"]
            + replicate_meta["wall_seconds"]
        ),
        "ground_truth_used": False,
    }


def _policy_bank(benchmark: OriginalBenchmark) -> tuple[np.ndarray, np.ndarray]:
    purpose = "function" if benchmark.name == "hjb_log_d100" else (
        "official" if benchmark.name == "dpi_ou_gmm_d100" else "native"
    )
    return benchmark.sample_states(
        96,
        seed=stable_seed("original_benchmark_common_policy_bank", benchmark.name),
        purpose=purpose,
    )


def _policy_validation(
    benchmark: OriginalBenchmark,
    trial: str,
    seed: int,
    iteration: int,
    source: TerminalAwareActor,
    actor: TerminalAwareActor,
    device: torch.device,
) -> dict[str, Any]:
    t, x = _policy_bank(benchmark)
    evaluation_seed = stable_seed(
        "original_benchmark_common_policy_paths", benchmark.name, seed, iteration
    )
    common = {
        "paths": 256,
        "time_steps": 32,
        "seed": evaluation_seed,
        "device": device,
        "parent_batch": 8,
        "path_chunk": 64,
    }
    old, old_meta = evaluate_policy(benchmark, source, t, x, **common)
    new, new_meta = evaluate_policy(benchmark, actor, t, x, **common)
    key = "mcv" if benchmark.actor_kind == "gradient" else "plain"
    paired = old[key] - new[key]
    return {
        "interpretation": "positive source_minus_new means lower tangent-policy value",
        "bank": "fixed reference-free common bank shared by all trials",
        "states": len(t),
        "source_mean": float(np.mean(old[key])),
        "new_mean": float(np.mean(new[key])),
        "source_minus_new_mean": float(np.mean(paired)),
        "source_minus_new_state_se": float(np.std(paired, ddof=1) / np.sqrt(len(paired))),
        "source_metadata": old_meta,
        "new_metadata": new_meta,
        "ground_truth_used": False,
    }


def train_trial(
    benchmark_name: str,
    trial: str,
    seed: int,
    iterations: int,
    device: torch.device,
) -> dict[str, Any]:
    ensure_dirs()
    benchmark = make_benchmark(benchmark_name)
    config = replace(TRIALS[trial], training_distribution=_distribution(benchmark))
    output_path = RAW_DIR / f"actor_{benchmark_name}_{trial}_s{seed}.json"
    phases: list[dict[str, Any]] = []
    cumulative = 0.0
    source = TerminalAwareActor(
        benchmark,
        width=config.width,
        depth=config.depth,
        activation=config.activation,
    ).to(device)
    source.eval()
    if output_path.exists():
        existing = json.loads(output_path.read_text())
        if existing["configuration"] != asdict(config):
            raise RuntimeError("existing trial configuration differs")
        phases = existing.get("phases", [])
        if phases:
            source, _ = load_actor(phases[-1]["checkpoint"], benchmark, device)
            cumulative = float(phases[-1]["cumulative_training_seconds"])
    if len(phases) >= iterations:
        return json.loads(output_path.read_text())

    synchronize(device)
    resource_started = time.perf_counter()
    for iteration in range(len(phases) + 1, iterations + 1):
        synchronize(device)
        phase_started = time.perf_counter()
        arrays, label_meta = _phase_labels(
            benchmark, trial, seed, iteration, config, source, device
        )
        actor, training, trace = fit_actor(
            benchmark,
            config,
            source,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_labels"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_labels"],
            seed=stable_seed("original_benchmark_fit", benchmark.name, trial, seed, iteration),
            device=device,
        )
        policy = _policy_validation(
            benchmark, trial, seed, iteration, source, actor, device
        )
        checkpoint = CHECKPOINT_DIR / (
            f"actor_{benchmark_name}_{trial}_s{seed}_a{iteration}.pt"
        )
        atomic_torch_save(
            checkpoint,
            {
                "state_dict": {
                    key: value.detach().cpu().clone()
                    for key, value in actor.state_dict().items()
                },
                "benchmark": benchmark_name,
                "trial": trial,
                "seed": seed,
                "iteration": iteration,
                "width": config.width,
                "depth": config.depth,
                "activation": config.activation,
                "training_label_estimator": config.label_estimator,
            },
        )
        synchronize(device)
        phase_seconds = time.perf_counter() - phase_started
        cumulative += phase_seconds
        phase = {
            "iteration": iteration,
            "checkpoint": str(checkpoint),
            "phase_seconds": phase_seconds,
            "cumulative_training_seconds": cumulative,
            "label_seconds_reported": label_meta["wall_seconds"],
            "fit_seconds_reported": training["wall_seconds"],
            "labels": label_meta,
            "training": training,
            "trace": trace,
            "reference_free_policy_validation": policy,
            "selection_criterion": "held-out tangent labels and common reference-free policy value",
            "ground_truth_used_for_selection": False,
        }
        phases.append(phase)
        source = actor
        result = {
            "method": "terminal_aware_tangent_actor_mcv",
            "benchmark": benchmark_name,
            "trial": trial,
            "seed": seed,
            "configuration": asdict(config),
            "benchmark_metadata": benchmark.metadata(),
            "phases": phases,
            "cumulative_training_seconds": cumulative,
            "resource_seconds_this_invocation": time.perf_counter() - resource_started,
            "ground_truth_used_for_training_or_selection": False,
            "oracle_subspace_projection": False,
        }
        atomic_json(output_path, result)
        print(
            f"{benchmark_name} {trial} seed={seed} a{iteration} "
            f"time={cumulative:.1f}s policy={policy['new_mean']:.8g} "
            f"improvement={policy['source_minus_new_mean']:.6g}",
            flush=True,
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARK_NAMES, required=True)
    parser.add_argument("--trial", choices=tuple(TRIALS), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("a CUDA allocation is required")
    result = train_trial(
        args.benchmark, args.trial, args.seed, args.iterations, device
    )
    print(
        json.dumps(
            {
                "status": "complete",
                "benchmark": args.benchmark,
                "trial": args.trial,
                "seed": args.seed,
                "phases": len(result["phases"]),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

