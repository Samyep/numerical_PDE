"""Deep Picard Iteration baseline adapted from the authors' public formulation.

Each round regresses value and z labels from the Picard map.  The gradient
label uses a terminal pathwise derivative and a centered Bismut estimator for
the driver integral, avoiding the uncentered estimator's infinite variance at
sample times close to the parent time.
"""

from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    make_actor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    train_actor,
)

from .common import (
    BUDGETS,
    CHECKPOINT_DIR,
    RAW_DIR,
    atomic_json,
    atomic_torch_save,
    driver_for,
    parameter_count,
    stable_seed,
    synchronize,
)


class ZeroPicardActor(nn.Module):
    actor_kind = "zero_picard"

    def __init__(self, d: int, has_gradient: bool) -> None:
        super().__init__()
        self.d, self.has_gradient = d, has_gradient

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        del t, T
        value = torch.zeros(x.shape[:-1], dtype=x.dtype, device=x.device)
        return value, torch.zeros_like(x) if self.has_gradient else None


@torch.no_grad()
def picard_labels(
    benchmark: str,
    source: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 4,
    path_chunk: int = 32,
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    if paths % path_chunk:
        raise ValueError("paths must be divisible by path_chunk")
    driver = driver_for(benchmark)
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, driver.d)
    values = np.empty(len(times), dtype=np.float64)
    gradients = np.empty((len(times), driver.d), dtype=np.float64) if driver.has_gradient else None
    generator = torch.Generator(device=device).manual_seed(seed)
    f_sum = f_square_sum = centered_square_sum = 0.0
    count = nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for begin in range(0, len(times), parent_batch):
        end = min(begin + parent_batch, len(times))
        bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
        bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
        tau = torch.clamp(driver.T - bt, min=1e-6)
        value_sum = torch.zeros(len(bt), dtype=torch.float64, device=device)
        q_sum = (
            torch.zeros((len(bt), driver.d), dtype=torch.float64, device=device)
            if driver.has_gradient
            else None
        )
        for _ in range(0, paths, path_chunk):
            terminal_normal = torch.randn(
                (len(bt), path_chunk, driver.d), device=device, generator=generator
            )
            terminal_x = bx[:, None, :] + driver.sigma * torch.sqrt(tau)[:, None, None] * terminal_normal
            terminal_value = driver.terminal(terminal_x.reshape(-1, driver.d)).reshape(len(bt), path_chunk)
            value_sum += terminal_value.double().sum(dim=1)
            if driver.has_gradient:
                assert q_sum is not None
                terminal_q = driver.sigma * driver.terminal_gradient_torch(
                    terminal_x.reshape(-1, driver.d)
                ).reshape(len(bt), path_chunk, driver.d)
                q_sum += terminal_q.double().sum(dim=1)

            uniform = torch.clamp(
                torch.rand((len(bt), path_chunk), device=device, generator=generator),
                min=1e-5,
            )
            delta = tau[:, None] * uniform
            normal = torch.randn(
                (len(bt), path_chunk, driver.d), device=device, generator=generator
            )
            sample_t = bt[:, None] + delta
            sample_x = bx[:, None, :] + driver.sigma * torch.sqrt(delta)[:, :, None] * normal
            flat_t = sample_t.reshape(-1)
            flat_x = sample_x.reshape(-1, driver.d)
            source_r, source_q = source.predict(flat_t, flat_x, driver.T)
            sample_f = driver.driver(source_r, source_q).reshape(len(bt), path_chunk)
            value_sum += (tau[:, None] * sample_f).double().sum(dim=1)
            f_sum += float(sample_f.double().sum().item())
            f_square_sum += float(sample_f.double().square().sum().item())
            if driver.has_gradient:
                assert q_sum is not None
                parent_x = bx[:, None, :].expand(-1, path_chunk, -1).reshape(-1, driver.d)
                parent_r, parent_q = source.predict(flat_t, parent_x, driver.T)
                parent_f = driver.driver(parent_r, parent_q).reshape(len(bt), path_chunk)
                centered = sample_f - parent_f
                centered_square_sum += float(centered.double().square().sum().item())
                q_integral = (
                    tau[:, None, None]
                    * centered[:, :, None]
                    * normal
                    / torch.sqrt(delta)[:, :, None]
                )
                q_sum += q_integral.double().sum(dim=1)
            nonfinite += int((~torch.isfinite(sample_f)).sum().item())
            count += len(bt) * path_chunk
        values[begin:end] = (value_sum / paths).cpu().numpy()
        if gradients is not None:
            assert q_sum is not None
            gradients[begin:end] = (q_sum / paths).cpu().numpy()
    synchronize(device)
    if nonfinite or not np.all(np.isfinite(values)) or (
        gradients is not None and not np.all(np.isfinite(gradients))
    ):
        raise FloatingPointError("nonfinite Deep Picard labels")
    f_mean = f_sum / count
    f_variance = max(f_square_sum / count - f_mean**2, 0.0)
    metadata = {
        "kind": "one-time-randomized Picard Feynman-Kac with centered Bismut gradient",
        "points": len(times),
        "paths": paths,
        "seed": seed,
        "driver_mean": f_mean,
        "driver_std": math.sqrt(f_variance),
        "centered_driver_rms": math.sqrt(centered_square_sum / count)
        if driver.has_gradient
        else None,
        "nonfinite_count": nonfinite,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "ground_truth_used": False,
    }
    return values, gradients, metadata


def _points(
    benchmark: str, seed: int, round_index: int, split: str, count: int, d: int, T: float
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(
        stable_seed("deep_picard_points", benchmark, seed, round_index, split)
    )
    return (
        rng.uniform(0.0, np.nextafter(T, 0.0), count).astype(np.float64),
        rng.standard_normal((count, d)).astype(np.float64),
    )


def train_deep_picard(
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    width: int = 256,
    depth: int = 3,
    training_points: int = 2048,
    validation_points: int = 256,
    paths: int = 128,
    epochs: int = 20,
    batch_size: int = 512,
    learning_rate: float = 5e-4,
    maximum_rounds: int = 40,
    budgets: tuple[float, ...] = BUDGETS,
    calibration_rounds: int | None = None,
    trial_tag: str = "",
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    source: nn.Module = ZeroPicardActor(driver.d, driver.has_gradient).to(device)
    optimizer_steps = epochs * math.ceil(training_points / batch_size)
    config = FitConfig(
        training_points=training_points,
        validation_points=validation_points,
        paths_per_label=paths,
        label_time_steps=1,
        optimizer_steps=optimizer_steps,
        width=width,
        depth=depth,
        learning_rate=learning_rate,
        batch_size=batch_size,
        parent_batch=4,
        path_chunk=32,
        checkpoint_interval=max(1, optimizer_steps // 8),
        gradient_clip=10.0,
    )
    rounds = calibration_rounds or maximum_rounds
    cumulative = 0.0
    next_budget = 0
    history: list[dict[str, Any]] = []
    checkpoints: list[dict[str, Any]] = []
    synchronize(device)
    overall_started = time.perf_counter()
    for round_index in range(1, rounds + 1):
        synchronize(device)
        round_started = time.perf_counter()
        train_t, train_x = _points(
            benchmark, seed, round_index, "train", training_points, driver.d, driver.T
        )
        valid_t, valid_x = _points(
            benchmark, seed, round_index, "valid", validation_points, driver.d, driver.T
        )
        train_value, train_q, train_meta = picard_labels(
            benchmark,
            source,
            train_t,
            train_x,
            paths=paths,
            seed=stable_seed("deep_picard_labels", benchmark, seed, round_index, "train"),
            device=device,
        )
        valid_value, valid_q, valid_meta = picard_labels(
            benchmark,
            source,
            valid_t,
            valid_x,
            paths=paths,
            seed=stable_seed("deep_picard_labels", benchmark, seed, round_index, "valid"),
            device=device,
        )
        train_r = np.clip(train_value, driver.r_low, driver.r_high)
        valid_r = np.clip(valid_value, driver.r_low, driver.r_high)
        initial = source if getattr(source, "actor_kind", None) == "terminal_aware" else None
        actor, training, trace = train_actor(
            driver,
            config,
            "terminal_aware",
            train_t,
            train_x,
            train_r,
            train_q,
            valid_t,
            valid_x,
            valid_r,
            valid_q,
            seed=stable_seed("deep_picard_fit", benchmark, seed, round_index),
            device=device,
            initial=initial,
        )
        synchronize(device)
        round_seconds_without_save = time.perf_counter() - round_started
        cumulative += round_seconds_without_save
        crossed: list[float] = []
        while next_budget < len(budgets) and cumulative >= budgets[next_budget]:
            crossed.append(budgets[next_budget])
            next_budget += 1
        is_final = round_index == rounds
        should_save = bool(crossed or is_final)
        checkpoint_path: Path | None = None
        if should_save:
            label = f"b{int(crossed[-1])}" if crossed else f"round{round_index}"
            checkpoint_path = CHECKPOINT_DIR / f"deep_picard_{benchmark}_s{seed}{trial_tag}_{label}.pt"
            atomic_torch_save(
                checkpoint_path,
                {
                    "state_dict": {
                        name: value.detach().cpu().clone()
                        for name, value in actor.state_dict().items()
                    },
                    "benchmark": benchmark,
                    "seed": seed,
                    "round": round_index,
                    "width": width,
                    "depth": depth,
                },
            )
            synchronize(device)
            io_seconds = time.perf_counter() - round_started - round_seconds_without_save
            cumulative += io_seconds
        row = {
            "round": round_index,
            "round_seconds": time.perf_counter() - round_started,
            "cumulative_training_seconds": cumulative,
            "train_label_seconds": train_meta["wall_seconds"],
            "validation_label_seconds": valid_meta["wall_seconds"],
            "fit_seconds": training["wall_seconds"],
            "validation_normalized_loss": training["best_validation_normalized_loss"],
            "validation_value_rmse": training["selected_validation_value_rmse"],
            "validation_gradient_vector_rmse": training[
                "selected_validation_gradient_vector_rmse"
            ],
            "selected_optimizer_step": training["selected_step"],
            "checkpoint": None if checkpoint_path is None else str(checkpoint_path),
            "crossed_budgets": crossed,
            "nonfinite_count": train_meta["nonfinite_count"] + valid_meta["nonfinite_count"],
            "trace": trace,
        }
        history.append(row)
        for target in crossed:
            checkpoints.append(
                {
                    "method": "deep_picard",
                    "benchmark": benchmark,
                    "seed": seed,
                    "target_budget_seconds": target,
                    "training_seconds": cumulative,
                    "round": round_index,
                    "validation_objective": training["best_validation_normalized_loss"],
                    "checkpoint": str(checkpoint_path),
                    "parameter_count": parameter_count(actor),
                    "selection_criterion": "first complete Picard round crossing predeclared time budget",
                    "final_test_used_for_selection": False,
                }
            )
        source = actor
        print(
            f"deep_picard {benchmark} seed={seed} round={round_index} "
            f"cumulative={cumulative:.2f}s val={training['best_validation_normalized_loss']:.5g}",
            flush=True,
        )
        # The common ladder ends at B5.  Continuing after the first complete
        # round that crosses it spends compute with no predeclared comparison
        # point and can make later, network-sourced rounds dominate wall time.
        if next_budget == len(budgets):
            break
    result = {
        "method": "deep_picard",
        "benchmark": benchmark,
        "seed": seed,
        "configuration": asdict(config),
        "maximum_rounds": rounds,
        "history": history,
        "checkpoints": checkpoints,
        "total_training_seconds": time.perf_counter() - overall_started,
        "ground_truth_used_for_training_or_selection": False,
        "official_configuration_delta": "authors use 4x512, 4096 points, 4096 paths and 40 rounds; the calibrated common-budget configuration uses 3x256, 2048 points and 128 paths, and stops after the complete round crossing B5",
    }
    suffix = f"cal{rounds}{trial_tag}" if calibration_rounds else "main"
    atomic_json(RAW_DIR / f"deep_picard_{benchmark}_s{seed}_{suffix}.json", result)
    return result


def load_deep_picard(path: Path, device: torch.device) -> tuple[nn.Module, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    driver = driver_for(payload["benchmark"])
    actor = make_actor(
        "terminal_aware", driver, width=payload["width"], depth=payload["depth"]
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor, payload


__all__ = ["load_deep_picard", "picard_labels", "train_deep_picard"]
