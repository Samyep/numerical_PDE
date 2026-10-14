"""No-projection teacher-free Actor Iteration and martingale-CV evaluation."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import math
from pathlib import Path
import time
from typing import Any, Iterator

import numpy as np
import torch
from torch import nn

from .common import (
    array_sha256,
    atomic_json,
    atomic_npz,
    atomic_torch_save,
    model_metrics,
    seed_everything,
    stable_seed,
    vector_relative_error,
)
from .equations import Equation, equation_key, torch_terminal, torch_terminal_z
from .models import ActorNet, clone_state_dict, count_parameters


@dataclass(frozen=True)
class ActorPhase:
    name: str
    source: str | None
    training_points: int
    validation_points: int
    paths_per_label: int
    label_time_steps: int
    optimizer_steps: int
    width: int = 384
    depth: int = 4
    learning_rate: float = 1.0e-4
    batch_size: int = 1024


@contextmanager
def frozen(module: nn.Module) -> Iterator[None]:
    flags = [parameter.requires_grad for parameter in module.parameters()]
    training = module.training
    for parameter in module.parameters():
        parameter.requires_grad_(False)
        parameter.grad = None
    module.eval()
    try:
        yield
    finally:
        for parameter, flag in zip(module.parameters(), flags, strict=True):
            parameter.requires_grad_(flag)
        module.train(training)


def antithetic_normals(
    shape: tuple[int, ...],
    *,
    generator: torch.Generator,
    device: torch.device,
    dtype: torch.dtype = torch.float32,
) -> torch.Tensor:
    if len(shape) < 2 or shape[1] % 2:
        raise ValueError("path count must be even")
    half = torch.randn(
        (shape[0], shape[1] // 2, *shape[2:]),
        generator=generator,
        device=device,
        dtype=dtype,
    )
    return torch.cat((half, -half), dim=1)


def terminal_pathwise_labels(
    equation: Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 16,
) -> tuple[np.ndarray, dict[str, Any]]:
    if paths % 2:
        raise ValueError("paths must be even")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, equation.d)
    labels = np.empty_like(states)
    generator = torch.Generator(device=device).manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for start in range(0, len(times), parent_batch):
        stop = min(start + parent_batch, len(times))
        batch_t = torch.as_tensor(times[start:stop], dtype=torch.float32, device=device)
        batch_x = torch.as_tensor(states[start:stop], dtype=torch.float32, device=device)
        normal = antithetic_normals(
            (len(batch_t), paths, equation.d), generator=generator, device=device
        )
        tau = torch.clamp(equation.T - batch_t, min=0.0)
        terminal_x = batch_x[:, None, :] + equation.sigma * torch.sqrt(tau)[:, None, None] * normal
        with torch.no_grad():
            labels[start:stop] = torch_terminal_z(equation, terminal_x).mean(dim=1).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    return labels, {
        "kind": "terminal_pathwise_gradient",
        "points": len(times),
        "paths": paths,
        "seed": seed,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "analytic_nonterminal_reference_used": False,
    }


def _pathwise_policy_chunk(
    equation: Equation,
    source: ActorNet,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    batch, paths, steps, dimension = normals.shape
    if dimension != equation.d:
        raise ValueError("normal dimension mismatch")
    initial = x.detach().clone().requires_grad_(True)
    times = t.detach()
    dt = torch.clamp(equation.T - times, min=0.0) / steps
    root_dt = torch.sqrt(dt)
    state = initial[:, None, :].expand(batch, paths, equation.d)
    running = torch.zeros((batch, paths), dtype=x.dtype, device=x.device)
    with frozen(source), torch.enable_grad():
        for step in range(steps):
            step_t = times + step * dt
            policy = source.predict(
                step_t[:, None].expand(batch, paths).reshape(-1),
                state.reshape(-1, equation.d),
                equation.T,
            ).reshape(batch, paths, equation.d)
            dw = root_dt[:, None, None] * normals[:, :, step, :]
            running = running + 0.5 * torch.sum(policy.square(), dim=-1) * dt[:, None]
            state = state - equation.sigma * policy * dt[:, None, None] + equation.sigma * dw
        value = (torch_terminal(equation, state) + running).mean(dim=1)
        z = equation.sigma * torch.autograd.grad(value.sum(), initial)[0]
    return value.detach(), z.detach()


def policy_gradient_labels(
    equation: Equation,
    source: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 8,
    path_chunk: int = 64,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, equation.d)
    values = np.empty(len(times), dtype=np.float64)
    labels = np.empty_like(states)
    generator = torch.Generator(device=device).manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for start in range(0, len(times), parent_batch):
        stop = min(start + parent_batch, len(times))
        batch_t = torch.as_tensor(times[start:stop], dtype=torch.float32, device=device)
        batch_x = torch.as_tensor(states[start:stop], dtype=torch.float32, device=device)
        value_sum = torch.zeros(len(batch_t), dtype=torch.float64, device=device)
        z_sum = torch.zeros((len(batch_t), equation.d), dtype=torch.float64, device=device)
        for _ in range(0, paths, path_chunk):
            normal = antithetic_normals(
                (len(batch_t), path_chunk, time_steps, equation.d),
                generator=generator,
                device=device,
            )
            chunk_value, chunk_z = _pathwise_policy_chunk(
                equation, source, batch_t, batch_x, normal
            )
            value_sum += path_chunk * chunk_value.to(torch.float64)
            z_sum += path_chunk * chunk_z.to(torch.float64)
        values[start:stop] = (value_sum / paths).cpu().numpy()
        labels[start:stop] = (z_sum / paths).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    if not np.all(np.isfinite(labels)):
        raise FloatingPointError("nonfinite policy-gradient labels")
    return values, labels, {
        "kind": "full_controlled_trajectory_pathwise_gradient",
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "source_actor_frozen": True,
        "state_dependence_preserved": True,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "analytic_nonterminal_reference_used": False,
    }


def train_actor_supervised(
    equation: Equation,
    phase: ActorPhase,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_y: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_y: np.ndarray,
    *,
    seed: int,
    device: torch.device,
    initial: ActorNet | None,
) -> tuple[ActorNet, dict[str, Any], list[dict[str, Any]]]:
    seed_everything(seed)
    actor = ActorNet(
        equation.d, width=phase.width, depth=phase.depth
    ).to(device)
    if initial is not None:
        if initial.width != phase.width or initial.depth != phase.depth:
            raise ValueError("warm-start architecture mismatch")
        actor.load_state_dict(initial.state_dict())
    tensors = {
        "train_t": torch.as_tensor(train_t, dtype=torch.float32, device=device),
        "train_x": torch.as_tensor(train_x, dtype=torch.float32, device=device),
        "train_y": torch.as_tensor(train_y, dtype=torch.float32, device=device),
        "valid_t": torch.as_tensor(valid_t, dtype=torch.float32, device=device),
        "valid_x": torch.as_tensor(valid_x, dtype=torch.float32, device=device),
        "valid_y": torch.as_tensor(valid_y, dtype=torch.float32, device=device),
    }
    optimizer = torch.optim.AdamW(actor.parameters(), lr=phase.learning_rate, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer,
        milestones=[max(1, phase.optimizer_steps // 2), max(2, 4 * phase.optimizer_steps // 5)],
        gamma=0.2,
    )
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = sorted(set([phase.optimizer_steps // 2, phase.optimizer_steps]))
    best_loss = float("inf")
    best_step = 0
    best_state: dict[str, torch.Tensor] | None = None
    history: list[dict[str, Any]] = []
    trajectory: list[dict[str, Any]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for step in range(1, phase.optimizer_steps + 1):
        actor.train()
        index = torch.randint(
            len(train_t),
            (min(phase.batch_size, len(train_t)),),
            generator=generator,
            device=device,
        )
        error = actor.predict(
            tensors["train_t"][index], tensors["train_x"][index], equation.T
        ) - tensors["train_y"][index]
        vector_mse = torch.mean(torch.sum(error.square(), dim=-1))
        loss = vector_mse / equation.d
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step % 250 == 0 or step in checkpoints:
            actor.eval()
            total = 0.0
            with torch.no_grad():
                for start in range(0, len(valid_t), 1024):
                    stop = min(start + 1024, len(valid_t))
                    valid_error = actor.predict(
                        tensors["valid_t"][start:stop],
                        tensors["valid_x"][start:stop],
                        equation.T,
                    ) - tensors["valid_y"][start:stop]
                    total += float(torch.sum(valid_error.square()).item())
            valid_loss = total / len(valid_t)
            record = {
                "step": step,
                "training_vector_mse": float(vector_mse.item()),
                "validation_vector_mse": valid_loss,
                "learning_rate": float(optimizer.param_groups[0]["lr"]),
                "elapsed_seconds": time.perf_counter() - started,
            }
            history.append(record)
            if valid_loss < best_loss:
                best_loss = valid_loss
                best_step = step
                best_state = clone_state_dict(actor)
            if step in checkpoints:
                trajectory.append({**record, "state_dict": clone_state_dict(actor)})
    if best_state is None:
        raise RuntimeError("actor training produced no checkpoint")
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    actor.load_state_dict(best_state)
    actor.eval()
    metadata = {
        "phase": asdict(phase),
        "seed": seed,
        "parameter_count": count_parameters(actor),
        "warm_started": initial is not None,
        "best_step": best_step,
        "best_validation_vector_mse": best_loss,
        "optimizer_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "history": history,
        "selection": "held-out pathwise-label vector MSE within phase",
        "analytic_reference_used_for_training": False,
    }
    return actor, metadata, trajectory


def predict_actor(
    actor: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    equation: Equation,
    device: torch.device,
) -> np.ndarray:
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, equation.d)
    output: list[np.ndarray] = []
    actor.eval()
    with torch.no_grad():
        for start in range(0, len(times), 4096):
            stop = min(start + 4096, len(times))
            output.append(
                actor.predict(
                    torch.as_tensor(times[start:stop], dtype=torch.float32, device=device),
                    torch.as_tensor(states[start:stop], dtype=torch.float32, device=device),
                    equation.T,
                ).cpu().numpy()
            )
    return np.concatenate(output).astype(np.float64)


def evaluate_policy_cv(
    equation: Equation,
    actor: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 16,
    path_chunk: int = 64,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Return CV value estimates and per-point estimator variances."""

    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, equation.d)
    output = np.empty(len(times), dtype=np.float64)
    variances = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    pair_count = paths // 2
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor), torch.no_grad():
        for start in range(0, len(times), parent_batch):
            stop = min(start + parent_batch, len(times))
            batch_t = torch.as_tensor(times[start:stop], dtype=torch.float32, device=device)
            initial = torch.as_tensor(states[start:stop], dtype=torch.float32, device=device)
            corrected_sum = torch.zeros(len(batch_t), dtype=torch.float64, device=device)
            paired_chunks: list[torch.Tensor] = []
            dt = torch.clamp(equation.T - batch_t, min=0.0) / time_steps
            root_dt = torch.sqrt(dt)
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (len(batch_t), path_chunk, time_steps, equation.d),
                    generator=generator,
                    device=device,
                )
                state = initial[:, None, :].expand(len(batch_t), path_chunk, equation.d)
                running = torch.zeros((len(batch_t), path_chunk), dtype=torch.float32, device=device)
                martingale = torch.zeros_like(running)
                for step in range(time_steps):
                    step_t = batch_t + step * dt
                    policy = actor.predict(
                        step_t[:, None].expand(len(batch_t), path_chunk).reshape(-1),
                        state.reshape(-1, equation.d),
                        equation.T,
                    ).reshape(len(batch_t), path_chunk, equation.d)
                    dw = root_dt[:, None, None] * normals[:, :, step, :]
                    running += 0.5 * torch.sum(policy.square(), dim=-1) * dt[:, None]
                    martingale += torch.sum(policy * dw, dim=-1)
                    state = state - equation.sigma * policy * dt[:, None, None] + equation.sigma * dw
                corrected = torch_terminal(equation, state) + running - martingale
                corrected_sum += corrected.sum(dim=1).to(torch.float64)
                half = path_chunk // 2
                paired = 0.5 * (corrected[:, :half] + corrected[:, half:])
                paired_chunks.append(paired.to(torch.float64))
            output[start:stop] = (corrected_sum / paths).cpu().numpy()
            all_pairs = torch.cat(paired_chunks, dim=1)
            sample_var = torch.var(all_pairs, dim=1, correction=1 if pair_count > 1 else 0)
            variances[start:stop] = (torch.clamp(sample_var, min=0.0) / pair_count).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    return output, variances, {
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "antithetic": True,
        "martingale_control_variate": "-integral q^T dW",
        "seed": seed,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
    }


def save_actor(path: Path, actor: ActorNet, metadata: dict[str, Any]) -> None:
    atomic_torch_save(
        path,
        {
            "state_dict": clone_state_dict(actor),
            "architecture": {
                "d": actor.d,
                "width": actor.width,
                "depth": actor.depth,
                "output_projection": False,
            },
            "metadata": metadata,
        },
    )


def load_actor(path: Path, device: torch.device) -> tuple[ActorNet, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    actor = ActorNet(
        int(architecture["d"]),
        width=int(architecture["width"]),
        depth=int(architecture["depth"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor, payload["metadata"]
