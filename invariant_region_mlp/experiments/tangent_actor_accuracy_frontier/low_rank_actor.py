"""Non-oracle learned low-rank Gradient-Actor ablation."""

from __future__ import annotations

from dataclasses import asdict
import math
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.models import (
    count_parameters,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    _fit_metrics,
)

from .common import atomic_torch_save, seed_everything, synchronize


def _network(input_dim: int, output_dim: int, width: int, depth: int) -> nn.Sequential:
    layers: list[nn.Module] = []
    previous = input_dim
    for _ in range(depth):
        layer = nn.Linear(previous, width)
        nn.init.xavier_uniform_(layer.weight)
        nn.init.zeros_(layer.bias)
        layers.extend((layer, nn.SiLU()))
        previous = width
    final = nn.Linear(previous, output_dim)
    nn.init.zeros_(final.weight)
    nn.init.zeros_(final.bias)
    layers.append(final)
    return nn.Sequential(*layers)


class LearnedLowRankTerminalAwareActor(nn.Module):
    """Full-ambient network with a learned d-by-r Gradient output factor."""

    actor_kind = "terminal_aware_learned_low_rank"

    def __init__(self, driver, *, width: int, depth: int, rank: int, seed: int) -> None:
        super().__init__()
        if rank <= 0 or rank > driver.d:
            raise ValueError(rank)
        self.driver, self.width, self.depth, self.rank = driver, width, depth, rank
        self.network = _network(driver.d + 1, 1 + rank, width, depth)
        generator = torch.Generator().manual_seed(seed)
        random = torch.randn(driver.d, rank, generator=generator)
        basis, _ = torch.linalg.qr(random, mode="reduced")
        self.gradient_basis = nn.Parameter(basis)

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float):
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        output = self.network(inputs)
        tau = torch.clamp(T - t.reshape(-1), min=0.0)
        terminal = self.driver.terminal(x)
        span = self.driver.r_high - self.driver.r_low
        unit = torch.clamp(
            (terminal - self.driver.r_low) / span,
            min=torch.finfo(x.dtype).eps,
            max=1.0 - torch.finfo(x.dtype).eps,
        )
        terminal_logit = torch.log(unit) - torch.log1p(-unit)
        r = self.driver.r_low + span * torch.sigmoid(
            terminal_logit + tau * output[:, 0]
        )
        terminal_q = self.driver.sigma * self.driver.terminal_gradient_torch(x)
        q = terminal_q + tau[:, None] * (output[:, 1:] @ self.gradient_basis.T)
        return r, q


def make_low_rank_actor(driver, config: FitConfig, rank: int, *, seed: int, device):
    seed_everything(seed)
    actor = LearnedLowRankTerminalAwareActor(
        driver, width=config.width, depth=config.depth, rank=rank, seed=seed + 11
    ).to(device)
    actor.eval()
    return actor


def save_low_rank_actor(path, actor, *, driver, config, metadata) -> None:
    atomic_torch_save(path, {
        "state_dict": {key: value.detach().cpu().clone() for key, value in actor.state_dict().items()},
        "driver": driver.metadata(),
        "fit_config": asdict(config),
        "rank": actor.rank,
        "metadata": metadata,
    })


def load_low_rank_actor(path, driver, config: FitConfig, device):
    payload = torch.load(path, map_location=device, weights_only=True)
    actor = LearnedLowRankTerminalAwareActor(
        driver, width=config.width, depth=config.depth,
        rank=int(payload["rank"]), seed=0,
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor


def _tensor(array: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(array, dtype=torch.float32, device=device)


def train_low_rank_actor(
    driver,
    config: FitConfig,
    rank: int,
    arrays: dict[str, np.ndarray],
    *,
    seed: int,
    device: torch.device,
    initial: LearnedLowRankTerminalAwareActor,
) -> tuple[LearnedLowRankTerminalAwareActor, dict[str, Any], list[dict[str, Any]]]:
    seed_everything(seed)
    actor = make_low_rank_actor(driver, config, rank, seed=seed, device=device)
    actor.load_state_dict(initial.state_dict())
    tt, tx, tr, tq = (
        _tensor(arrays["train_t"], device), _tensor(arrays["train_x"], device),
        _tensor(arrays["train_r"], device), _tensor(arrays["train_q"], device),
    )
    vt, vx, vr, vq = (
        _tensor(arrays["valid_t"], device), _tensor(arrays["valid_x"], device),
        _tensor(arrays["valid_r"], device), _tensor(arrays["valid_q"], device),
    )
    r_scale = max(float(np.std(arrays["train_r"])), 1e-2)
    q_scale = max(float(np.sqrt(np.mean(np.square(arrays["train_q"])))), 1e-3)
    optimizer = torch.optim.AdamW(actor.parameters(), lr=config.learning_rate, weight_decay=1e-6)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.optimizer_steps)
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = set(range(config.checkpoint_interval, config.optimizer_steps + 1, config.checkpoint_interval))
    checkpoints.update((0, 1, config.optimizer_steps))
    best_loss, best_state = float("inf"), None
    trace: list[dict[str, Any]] = []

    def record(step: int) -> None:
        nonlocal best_loss, best_state
        train = _fit_metrics(actor, driver, tt, tx, tr, tq, r_scale=r_scale, q_scale=q_scale)
        valid = _fit_metrics(actor, driver, vt, vx, vr, vq, r_scale=r_scale, q_scale=q_scale)
        row: dict[str, Any] = {"step": step, "learning_rate": optimizer.param_groups[0]["lr"]}
        row.update({f"train_{key}": value for key, value in train.items()})
        row.update({f"validation_{key}": value for key, value in valid.items()})
        row["generalization_gap_value_rmse"] = float(valid["value_rmse"]) - float(train["value_rmse"])
        row["generalization_gap_gradient_vector_rmse"] = float(valid["gradient_vector_rmse"]) - float(train["gradient_vector_rmse"])
        trace.append(row)
        loss = float(valid["normalized_loss"])
        if loss < best_loss:
            best_loss = loss
            best_state = {name: value.detach().cpu().clone() for name, value in actor.state_dict().items()}

    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    record(0)
    for step in range(1, config.optimizer_steps + 1):
        actor.train()
        index = torch.randint(len(tt), (min(config.batch_size, len(tt)),), generator=generator, device=device)
        predicted_r, predicted_q = actor.predict(tt[index], tx[index], driver.T)
        loss = (predicted_r - tr[index]).square().mean() / r_scale**2
        loss = loss + (predicted_q - tq[index]).square().mean() / q_scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)
    if best_state is None:
        raise RuntimeError("no low-rank checkpoint recorded")
    actor.load_state_dict(best_state)
    actor.eval()
    synchronize(device)
    selected = min(trace, key=lambda row: row["validation_normalized_loss"])
    training = {
        "config": asdict(config),
        "actor_kind": actor.actor_kind,
        "rank": rank,
        "parameters": count_parameters(actor),
        "selected_step": int(selected["step"]),
        "selected_by": "minimum normalized held-out noisy tangent-label loss",
        "best_validation_normalized_loss": best_loss,
        "selected_train_value_rmse": selected["train_value_rmse"],
        "selected_validation_value_rmse": selected["validation_value_rmse"],
        "selected_train_gradient_vector_rmse": selected["train_gradient_vector_rmse"],
        "selected_validation_gradient_vector_rmse": selected["validation_gradient_vector_rmse"],
        "selected_gradient_generalization_gap_vector_rmse": selected["generalization_gap_gradient_vector_rmse"],
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "learned_basis_initialized_without_oracle_information": True,
        "used_exact_reference": False,
    }
    return actor, training, trace


__all__ = [
    "LearnedLowRankTerminalAwareActor", "load_low_rank_actor", "make_low_rank_actor",
    "save_low_rank_actor", "train_low_rank_actor",
]
