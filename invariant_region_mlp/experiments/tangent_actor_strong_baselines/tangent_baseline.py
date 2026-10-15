"""Predeclared Tangent-Actor methods with end-to-end phase accounting."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import time
from typing import Any, Literal

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    evaluate_tangent_mc,
    make_targets,
    tangent_labels,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.models import (
    ConstantTangentActor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    make_actor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    train_actor,
)

from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    atomic_json,
    atomic_torch_save,
    driver_for,
    parameter_count,
    stable_seed,
    synchronize,
)


TangentMethod = Literal["tangent_actor", "standard_tangent_actor"]


def _points(
    benchmark: str,
    seed: int,
    iteration: int,
    split: str,
    count: int,
    d: int,
    T: float,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(
        stable_seed("strong_tangent_points", benchmark, seed, iteration, split)
    )
    t = rng.uniform(0.0, np.nextafter(T, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, d)).astype(np.float64)
    return t, x


def _generate_phase_labels(
    benchmark: str,
    training_seed: int,
    iteration: int,
    source: nn.Module,
    config: FitConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    driver = driver_for(benchmark)
    train_t, train_x = _points(
        benchmark,
        training_seed,
        iteration,
        "train",
        config.training_points,
        driver.d,
        driver.T,
    )
    valid_t, valid_x = _points(
        benchmark,
        training_seed,
        iteration,
        "validation",
        config.validation_points,
        driver.d,
        driver.T,
    )
    common = dict(
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    train_value, train_q, train_meta = tangent_labels(
        driver,
        source,
        train_t,
        train_x,
        seed=stable_seed("strong_tangent_labels", benchmark, training_seed, iteration, "train"),
        **common,
    )
    valid_value, valid_q, valid_meta = tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed("strong_tangent_labels", benchmark, training_seed, iteration, "valid"),
        **common,
    )
    replicate_value, replicate_q, replicate_meta = tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed("strong_tangent_labels", benchmark, training_seed, iteration, "replicate"),
        **common,
    )
    train_r, train_q, train_projection = make_targets(driver, train_value, train_q)
    valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_q)
    value_replication_rmse = float(
        np.sqrt(np.mean(np.square(valid_value - replicate_value)))
    )
    if valid_q is None:
        gradient_replication_component_rmse = None
        gradient_replication_vector_rmse = None
    else:
        assert replicate_q is not None
        difference = valid_q - replicate_q
        gradient_replication_component_rmse = float(np.sqrt(np.mean(difference**2)))
        gradient_replication_vector_rmse = float(
            np.sqrt(np.mean(np.square(difference).sum(axis=1)))
        )
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_r": train_r,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_r": valid_r,
    }
    if train_q is not None and valid_q is not None:
        arrays.update(train_q=train_q, valid_q=valid_q)
    metadata = {
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "train_projection": train_projection,
        "validation_projection": valid_projection,
        "value_label_replication_rmse": value_replication_rmse,
        "gradient_label_replication_component_rmse": gradient_replication_component_rmse,
        "gradient_label_replication_vector_rmse": gradient_replication_vector_rmse,
        "wall_seconds": float(
            train_meta["wall_seconds"]
            + valid_meta["wall_seconds"]
            + replicate_meta["wall_seconds"]
        ),
        "ground_truth_used": False,
    }
    return arrays, metadata


def train_tangent_method(
    method: TangentMethod,
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    iterations: int = 3,
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    actor_kind = "terminal_aware" if method == "tangent_actor" else "standard"
    config = FitConfig(
        training_points=2048,
        validation_points=256,
        paths_per_label=128,
        label_time_steps=16,
        optimizer_steps=2500,
        width=192,
        depth=3,
        learning_rate=2e-4,
        batch_size=512,
        parent_batch=8,
        path_chunk=32,
        checkpoint_interval=125,
        gradient_clip=10.0,
    )
    actors: list[nn.Module] = [ConstantTangentActor(driver).to(device)]
    phases: list[dict[str, Any]] = []
    cumulative = 0.0
    synchronize(device)
    overall_started = time.perf_counter()
    for iteration in range(1, iterations + 1):
        synchronize(device)
        phase_started = time.perf_counter()
        arrays, label_meta = _generate_phase_labels(
            benchmark, seed, iteration, actors[-1], config, device
        )
        initial = (
            actors[-1]
            if getattr(actors[-1], "actor_kind", None) == actor_kind
            else None
        )
        actor, training, trace = train_actor(
            driver,
            config,
            actor_kind,
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_r"],
            arrays.get("train_q"),
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_r"],
            arrays.get("valid_q"),
            seed=stable_seed("strong_tangent_fit", method, benchmark, seed, iteration),
            device=device,
            initial=initial,
        )
        checkpoint_path = CHECKPOINT_DIR / f"{method}_{benchmark}_s{seed}_a{iteration}.pt"
        atomic_torch_save(
            checkpoint_path,
            {
                "state_dict": {
                    name: value.detach().cpu().clone()
                    for name, value in actor.state_dict().items()
                },
                "actor_kind": actor_kind,
                "width": config.width,
                "depth": config.depth,
                "benchmark": benchmark,
                "seed": seed,
                "iteration": iteration,
            },
        )
        synchronize(device)
        actual_phase_seconds = time.perf_counter() - phase_started
        cumulative += actual_phase_seconds
        phase = {
            "method": method,
            "benchmark": benchmark,
            "seed": seed,
            "iteration": iteration,
            "checkpoint": str(checkpoint_path),
            "phase_seconds": actual_phase_seconds,
            "cumulative_training_seconds": cumulative,
            "label_seconds_reported": label_meta["wall_seconds"],
            "fit_seconds_reported": training["wall_seconds"],
            "checkpoint_io_and_other_seconds": actual_phase_seconds
            - label_meta["wall_seconds"]
            - training["wall_seconds"],
            "labels": label_meta,
            "training": training,
            "trace": trace,
            "parameter_count": parameter_count(actor),
            "selection_criterion": "minimum normalized held-out noisy tangent-label loss",
            "final_test_used_for_selection": False,
        }
        phases.append(phase)
        actors.append(actor)
        print(
            f"{method} {benchmark} seed={seed} a{iteration} "
            f"cumulative={cumulative:.2f}s val={training['best_validation_normalized_loss']:.6g}",
            flush=True,
        )
    result = {
        "method": method,
        "benchmark": benchmark,
        "seed": seed,
        "actor_kind": actor_kind,
        "configuration": asdict(config),
        "phases": phases,
        "total_training_seconds": time.perf_counter() - overall_started,
        "ground_truth_used_for_training_or_selection": False,
        "analytic_nonterminal_training_labels": False,
        "oracle_subspace_projection": False,
    }
    atomic_json(RAW_DIR / f"{method}_{benchmark}_s{seed}_main.json", result)
    return result


def load_tangent_actor(
    path: Path, device: torch.device
) -> tuple[nn.Module, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    driver = driver_for(payload["benchmark"])
    actor = make_actor(
        payload["actor_kind"], driver, width=payload["width"], depth=payload["depth"]
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor, payload


def terminal_tangent_actor(benchmark: str, device: torch.device) -> nn.Module:
    """The zero-residual terminal-aware network A_g=(g,sqrt(2) grad g)."""
    driver = driver_for(benchmark)
    actor = make_actor("terminal_aware", driver, width=192, depth=3).to(device)
    actor.eval()
    return actor


def evaluate_tangent_actor(
    benchmark: str,
    actor: nn.Module,
    x: np.ndarray,
    device: torch.device,
    *,
    paths: int,
    time_steps: int = 32,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    driver = driver_for(benchmark)
    t = np.zeros(len(x), dtype=np.float64)
    output, metadata = evaluate_tangent_mc(
        driver,
        actor,
        t,
        x,
        paths=paths,
        time_steps=time_steps,
        seed=seed,
        device=device,
        parent_batch=8,
        path_chunk=64 if paths >= 64 else paths,
    )
    return (
        np.asarray(driver.working_to_original(output["working_value"]), dtype=np.float64),
        np.asarray(output["mc_variance_of_mean"], dtype=np.float64),
        metadata,
    )


__all__ = [
    "evaluate_tangent_actor", "load_tangent_actor", "terminal_tangent_actor",
    "train_tangent_method",
]
