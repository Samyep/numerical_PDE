"""Fitting for the preregistered unbounded joint terminal-aware actor."""

from __future__ import annotations

from dataclasses import asdict
import math
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
    synchronize,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)

from .joint_actor import JointTerminalActor
from .benchmarks import DenseFourierDriver


def _metrics(
    actor: JointTerminalActor,
    driver: DenseFourierDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    r: torch.Tensor,
    q: torch.Tensor,
    *,
    r_scale: float,
    q_scale: float,
    chunk: int = 1024,
) -> dict[str, float]:
    r_square = q_square = q_vector_square = 0.0
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), chunk):
            end = min(begin + chunk, len(t))
            pr, pq = actor.predict(t[begin:end], x[begin:end], driver.T)
            re = (pr - r[begin:end]).double()
            qe = (pq - q[begin:end]).double()
            r_square += float(re.square().sum().item())
            q_square += float(qe.square().sum().item())
            q_vector_square += float(qe.square().sum(dim=1).sum().item())
    r_mse = r_square / len(t)
    q_mse = q_square / (len(t) * driver.d)
    return {
        "normalized_loss": r_mse / r_scale**2 + q_mse / q_scale**2,
        "value_rmse": math.sqrt(r_mse),
        "gradient_component_rmse": math.sqrt(q_mse),
        "gradient_vector_rmse": math.sqrt(q_vector_square / len(t)),
    }


def fit_joint_actor(
    driver: DenseFourierDriver,
    config: FitConfig,
    source: JointTerminalActor,
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
) -> tuple[JointTerminalActor, dict[str, Any], list[dict[str, Any]]]:
    seed_everything(seed)
    actor = JointTerminalActor(
        driver, width=config.width, depth=config.depth
    ).to(device)
    actor.load_state_dict(source.state_dict())
    tensors = [
        torch.as_tensor(value, dtype=torch.float32, device=device)
        for value in (
            train_t,
            train_x,
            train_r,
            train_q,
            valid_t,
            valid_x,
            valid_r,
            valid_q,
        )
    ]
    tt, tx, tr, tq, vt, vx, vr, vq = tensors
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
        train = _metrics(
            actor, driver, tt, tx, tr, tq, r_scale=r_scale, q_scale=q_scale
        )
        valid = _metrics(
            actor, driver, vt, vx, vr, vq, r_scale=r_scale, q_scale=q_scale
        )
        row: dict[str, Any] = {
            "step": step,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            **{f"train_{key}": value for key, value in train.items()},
            **{f"validation_{key}": value for key, value in valid.items()},
            "gradient_norm_mean_since_checkpoint": (
                float(np.mean(gradient_norms)) if gradient_norms else None
            ),
            "gradient_norm_max_since_checkpoint": (
                float(np.max(gradient_norms)) if gradient_norms else None
            ),
            "gradient_clip_fraction_since_checkpoint": (
                clipped / len(gradient_norms) if gradient_norms else None
            ),
        }
        trace.append(row)
        if valid["normalized_loss"] < best_loss:
            best_loss = valid["normalized_loss"]
            best_state = {
                name: value.detach().cpu().clone()
                for name, value in actor.state_dict().items()
            }
        gradient_norms, clipped = [], 0

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
        pr, pq = actor.predict(tt[index], tx[index], driver.T)
        loss = (pr - tr[index]).square().mean() / r_scale**2
        loss = loss + (pq - tq[index]).square().mean() / q_scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        norm = float(
            torch.nn.utils.clip_grad_norm_(
                actor.parameters(), config.gradient_clip
            ).item()
        )
        gradient_norms.append(norm)
        clipped += int(norm > config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)
    if best_state is None:
        raise RuntimeError("no joint actor checkpoint recorded")
    actor.load_state_dict(best_state)
    actor.eval()
    synchronize(device)
    selected = min(trace, key=lambda row: row["validation_normalized_loss"])
    return actor, {
        "config": asdict(config),
        "warm_started": True,
        "value_label_scale": r_scale,
        "gradient_label_component_scale": q_scale,
        "selected_step": int(selected["step"]),
        "selected_by": "minimum normalized held-out noisy joint tangent-label loss",
        "best_validation_normalized_loss": float(best_loss),
        "selected_validation_value_rmse": selected["validation_value_rmse"],
        "selected_validation_gradient_component_rmse": selected[
            "validation_gradient_component_rmse"
        ],
        "selected_validation_gradient_vector_rmse": selected[
            "validation_gradient_vector_rmse"
        ],
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": (
            int(torch.cuda.max_memory_allocated(device))
            if device.type == "cuda"
            else None
        ),
        "used_exact_reference": False,
        **actor.specification(),
    }, trace


__all__ = ["fit_joint_actor"]
