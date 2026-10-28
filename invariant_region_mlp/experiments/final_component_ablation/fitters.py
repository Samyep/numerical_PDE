"""Source-faithful fitting with a configurable Actor parameterization."""

from __future__ import annotations

from dataclasses import asdict
import math
import time
from typing import Any, Callable

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything as seed_original,
    synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    ActorConfig,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.common import (
    seed_everything as seed_joint,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    _fit_metrics,
)


def fit_original(
    benchmark: Any,
    config: ActorConfig,
    source: nn.Module,
    factory: Callable[[], nn.Module],
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_labels: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_labels: np.ndarray,
    *,
    seed: int,
    device: torch.device,
) -> tuple[nn.Module, dict[str, Any], list[dict[str, float]]]:
    """Exact original-HJB optimizer/checkpoint rule, generalized by factory."""
    seed_original(seed)
    actor = factory().to(device)
    actor.load_state_dict(source.state_dict())
    tt = torch.as_tensor(train_t, dtype=torch.float32, device=device)
    tx = torch.as_tensor(train_x, dtype=torch.float32, device=device)
    ty = torch.as_tensor(train_labels, dtype=torch.float32, device=device)
    vt = torch.as_tensor(valid_t, dtype=torch.float32, device=device)
    vx = torch.as_tensor(valid_x, dtype=torch.float32, device=device)
    vy = torch.as_tensor(valid_labels, dtype=torch.float32, device=device)
    scale = max(float(np.sqrt(np.mean(np.square(train_labels)))), 1e-3)
    optimizer = torch.optim.AdamW(
        actor.parameters(), lr=config.learning_rate, weight_decay=1e-6
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.optimizer_steps
    )
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = set(
        range(
            config.checkpoint_interval,
            config.optimizer_steps + 1,
            config.checkpoint_interval,
        )
    )
    checkpoints.update((0, 1, config.optimizer_steps))
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    trace: list[dict[str, float]] = []

    def prediction(t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        r, q = actor.predict(t, x)
        output = q if benchmark.actor_kind == "gradient" else r
        if output is None:
            raise RuntimeError("Actor returned the wrong target")
        return output

    def record(step: int) -> None:
        nonlocal best_loss, best_state
        actor.eval()
        with torch.no_grad():
            train_mse = float((prediction(tt, tx) - ty).double().square().mean())
            valid_mse = float((prediction(vt, vx) - vy).double().square().mean())
        normalized = valid_mse / scale**2
        trace.append({
            "step": float(step),
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "train_rmse": math.sqrt(train_mse),
            "validation_rmse": math.sqrt(valid_mse),
            "validation_normalized_mse": normalized,
        })
        if normalized < best_loss:
            best_loss = normalized
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in actor.state_dict().items()
            }

    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    record(0)
    for step in range(1, config.optimizer_steps + 1):
        actor.train()
        index = torch.randint(
            len(tt),
            (min(config.batch_size, len(tt)),),
            generator=generator,
            device=device,
        )
        loss = (prediction(tt[index], tx[index]) - ty[index]).square().mean() / scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)
    if best_state is None:
        raise RuntimeError("no original Actor checkpoint selected")
    actor.load_state_dict(best_state)
    actor.eval()
    synchronize(device)
    selected = min(trace, key=lambda row: row["validation_normalized_mse"])
    return actor, {
        "config": asdict(config),
        "warm_started": True,
        "target_scale": scale,
        "selected_step": int(selected["step"]),
        "selected_validation_rmse": selected["validation_rmse"],
        "selected_validation_normalized_mse": selected["validation_normalized_mse"],
        "selected_by": "minimum held-out noisy tangent-label error",
        "ground_truth_used": False,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "parameters": sum(parameter.numel() for parameter in actor.parameters()),
    }, trace


def fit_joint(
    driver: Any,
    config: FitConfig,
    source: nn.Module,
    factory: Callable[[], nn.Module],
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_r: np.ndarray,
    train_q: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_r: np.ndarray,
    valid_q: np.ndarray,
    *,
    seed: int,
    device: torch.device,
) -> tuple[nn.Module, dict[str, Any], list[dict[str, Any]]]:
    """Exact mature joint optimizer/checkpoint rule, generalized by factory."""
    seed_joint(seed)
    actor = factory().to(device)
    actor.load_state_dict(source.state_dict())
    tensor = lambda value: torch.as_tensor(value, dtype=torch.float32, device=device)
    tt, tx, tr, tq = map(tensor, (train_t, train_x, train_r, train_q))
    vt, vx, vr, vq = map(tensor, (valid_t, valid_x, valid_r, valid_q))
    r_scale = max(float(np.std(train_r)), 1e-2)
    q_scale = max(float(np.sqrt(np.mean(np.square(train_q)))), 1e-3)
    optimizer = torch.optim.AdamW(
        actor.parameters(), lr=config.learning_rate, weight_decay=1e-6
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.optimizer_steps
    )
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = set(
        range(
            config.checkpoint_interval,
            config.optimizer_steps + 1,
            config.checkpoint_interval,
        )
    )
    checkpoints.update((0, 1, config.optimizer_steps))
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    trace: list[dict[str, Any]] = []
    gradient_norms: list[float] = []
    clipped = 0

    def record(step: int) -> None:
        nonlocal best_loss, best_state, gradient_norms, clipped
        train_metrics = _fit_metrics(
            actor, driver, tt, tx, tr, tq, r_scale=r_scale, q_scale=q_scale
        )
        valid_metrics = _fit_metrics(
            actor, driver, vt, vx, vr, vq, r_scale=r_scale, q_scale=q_scale
        )
        row: dict[str, Any] = {
            "step": step,
            "learning_rate": optimizer.param_groups[0]["lr"],
            **{f"train_{key}": value for key, value in train_metrics.items()},
            **{f"validation_{key}": value for key, value in valid_metrics.items()},
            "gradient_norm_mean_since_checkpoint": float(np.mean(gradient_norms))
            if gradient_norms
            else None,
            "gradient_norm_max_since_checkpoint": float(np.max(gradient_norms))
            if gradient_norms
            else None,
            "gradient_clip_fraction_since_checkpoint": clipped / len(gradient_norms)
            if gradient_norms
            else None,
        }
        trace.append(row)
        score = float(valid_metrics["normalized_loss"])
        if score < best_loss:
            best_loss = score
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in actor.state_dict().items()
            }
        gradient_norms, clipped = [], 0

    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    record(0)
    for step in range(1, config.optimizer_steps + 1):
        actor.train()
        index = torch.randint(
            len(tt),
            (min(config.batch_size, len(tt)),),
            generator=generator,
            device=device,
        )
        predicted_r, predicted_q = actor.predict(tt[index], tx[index], driver.T)
        if predicted_q is None:
            raise RuntimeError("joint Actor omitted q")
        loss = (predicted_r - tr[index]).square().mean() / r_scale**2
        loss = loss + (predicted_q - tq[index]).square().mean() / q_scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        norm = float(
            torch.nn.utils.clip_grad_norm_(actor.parameters(), config.gradient_clip)
        )
        gradient_norms.append(norm)
        clipped += int(norm > config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)
    if best_state is None:
        raise RuntimeError("no joint Actor checkpoint selected")
    actor.load_state_dict(best_state)
    actor.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    selected = min(trace, key=lambda row: row["validation_normalized_loss"])
    return actor, {
        "config": asdict(config),
        "actor_kind": getattr(actor, "actor_kind", "terminal_aware"),
        "parameters": sum(parameter.numel() for parameter in actor.parameters()),
        "warm_started": True,
        "value_label_scale": r_scale,
        "gradient_label_scale": q_scale,
        "selected_step": int(selected["step"]),
        "selected_by": "minimum normalized held-out noisy tangent-label loss",
        "best_validation_normalized_loss": best_loss,
        "selected_validation_value_rmse": selected["validation_value_rmse"],
        "selected_validation_gradient_component_rmse": selected[
            "validation_gradient_component_rmse"
        ],
        "selected_validation_gradient_vector_rmse": selected[
            "validation_gradient_vector_rmse"
        ],
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "used_exact_reference": False,
    }, trace


__all__ = ["fit_joint", "fit_original"]
