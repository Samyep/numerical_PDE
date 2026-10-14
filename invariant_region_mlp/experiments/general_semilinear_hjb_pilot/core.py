"""Shared controlled simulation, pathwise VJP labels, training, and plain MC."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import time
from typing import Any, Iterator

import numpy as np
import torch
from torch import nn

from .common import seed_everything
from .models import ActorNet
from .problems import ControlledHJB


@dataclass(frozen=True)
class ActorConfig:
    training_points: int = 2048
    validation_points: int = 256
    paths_per_label: int = 128
    label_time_steps: int = 16
    optimizer_steps: int = 2500
    width: int = 192
    depth: int = 3
    learning_rate: float = 2.0e-4
    batch_size: int = 512
    parent_batch: int = 8
    path_chunk: int = 32


@contextmanager
def frozen(module: nn.Module) -> Iterator[None]:
    flags = [parameter.requires_grad for parameter in module.parameters()]
    was_training = module.training
    for parameter in module.parameters():
        parameter.requires_grad_(False)
        parameter.grad = None
    module.eval()
    try:
        yield
    finally:
        for parameter, flag in zip(module.parameters(), flags, strict=True):
            parameter.requires_grad_(flag)
        module.train(was_training)


def antithetic_normals(
    shape: tuple[int, ...], *, generator: torch.Generator, device: torch.device
) -> torch.Tensor:
    if len(shape) < 2 or shape[1] % 2:
        raise ValueError("path count must be even")
    half = torch.randn(
        (shape[0], shape[1] // 2, *shape[2:]),
        generator=generator,
        dtype=torch.float32,
        device=device,
    )
    return torch.cat((half, -half), dim=1)


def _pathwise_chunk(
    problem: ControlledHJB,
    source: nn.Module,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    batch, paths, steps, dimension = normals.shape
    if dimension != problem.d:
        raise ValueError("Brownian dimension mismatch")
    initial = x.detach().clone().requires_grad_(True)
    times = t.detach()
    dt = torch.clamp(problem.T - times, min=0.0) / steps
    state = initial[:, None, :].expand(batch, paths, dimension)
    running = torch.zeros((batch, paths), dtype=x.dtype, device=x.device)
    with torch.enable_grad():
        for step in range(steps):
            step_t = times + step * dt
            flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
            flat_x = state.reshape(-1, dimension)
            action = source.predict(flat_t, flat_x, problem.T).reshape(batch, paths, problem.m)
            cost = problem.running_cost(
                flat_t, flat_x, action.reshape(-1, problem.m)
            ).reshape(batch, paths)
            drift = problem.drift(flat_t, flat_x, action.reshape(-1, problem.m)).reshape(
                batch, paths, dimension
            )
            running = running + cost * dt[:, None]
            state = (
                state
                + drift * dt[:, None, None]
                + problem.sigma * torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
            )
        value = (problem.terminal(state) + running).mean(dim=1)
        gradient = torch.autograd.grad(value.sum(), initial, create_graph=False)[0]
    return value.detach(), gradient.detach()


def policy_gradient_labels(
    problem: ControlledHJB,
    source: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int,
    path_chunk: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, problem.d)
    values = np.empty(len(times), dtype=np.float64)
    gradients = np.empty_like(states)
    generator = torch.Generator(device=device).manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(source):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            batch_t = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            batch_x = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            value_sum = torch.zeros(end - begin, dtype=torch.float64, device=device)
            gradient_sum = torch.zeros((end - begin, problem.d), dtype=torch.float64, device=device)
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, problem.d),
                    generator=generator,
                    device=device,
                )
                chunk_value, chunk_gradient = _pathwise_chunk(
                    problem, source, batch_t, batch_x, normals
                )
                value_sum += path_chunk * chunk_value.double()
                gradient_sum += path_chunk * chunk_gradient.double()
            values[begin:end] = (value_sum / paths).cpu().numpy()
            gradients[begin:end] = (gradient_sum / paths).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(gradients)):
        raise FloatingPointError("nonfinite pathwise labels")
    return values, gradients, {
        "kind": "full_discretized_controlled_trajectory_vjp",
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "parent_batch": parent_batch,
        "path_chunk": path_chunk,
        "seed": seed,
        "source_actor_frozen": True,
        "state_dependence_preserved": True,
        "analytic_nonterminal_reference_used": False,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
    }


def greedy_targets(
    problem: ControlledHJB, t: np.ndarray, x: np.ndarray, gradient: np.ndarray, device: torch.device
) -> np.ndarray:
    with torch.no_grad():
        target = problem.greedy_action(
            torch.as_tensor(t, dtype=torch.float32, device=device),
            torch.as_tensor(x, dtype=torch.float32, device=device),
            torch.as_tensor(gradient, dtype=torch.float32, device=device),
        )
    result = target.cpu().numpy().astype(np.float64)
    if result.shape != (len(t), problem.m) or not np.all(np.isfinite(result)):
        raise ValueError("invalid greedy targets")
    return result


def train_actor(
    problem: ControlledHJB,
    config: ActorConfig,
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
) -> tuple[ActorNet, dict[str, Any], list[dict[str, float]]]:
    seed_everything(seed)
    actor = ActorNet(problem, width=config.width, depth=config.depth).to(device)
    if initial is not None:
        actor.load_state_dict(initial.state_dict())
    tensors = [
        torch.as_tensor(value, dtype=torch.float32, device=device)
        for value in (train_t, train_x, train_y, valid_t, valid_x, valid_y)
    ]
    tt, tx, ty, vt, vx, vy = tensors
    optimizer = torch.optim.AdamW(actor.parameters(), lr=config.learning_rate, weight_decay=1e-6)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.optimizer_steps)
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = sorted(set([1, config.optimizer_steps // 4, config.optimizer_steps // 2,
                              3 * config.optimizer_steps // 4, config.optimizer_steps]))
    best_state: dict[str, torch.Tensor] | None = None
    best_loss = float("inf")
    trace: list[dict[str, float]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for step in range(1, config.optimizer_steps + 1):
        index = torch.randint(
            len(tt), (min(config.batch_size, len(tt)),), generator=generator, device=device
        )
        prediction = actor.predict(tt[index], tx[index], problem.T)
        loss = (prediction - ty[index]).square().mean()
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            actor.eval()
            with torch.no_grad():
                valid_prediction = actor.predict(vt, vx, problem.T)
                valid_loss = float((valid_prediction - vy).square().mean().item())
                train_loss = float(loss.item())
            actor.train()
            trace.append({"step": float(step), "train_mse": train_loss, "validation_mse": valid_loss})
            if valid_loss < best_loss:
                best_loss = valid_loss
                best_state = {key: value.detach().cpu().clone() for key, value in actor.state_dict().items()}
    if best_state is None:
        raise RuntimeError("no actor checkpoint")
    actor.load_state_dict(best_state)
    actor.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    return actor, {
        "config": asdict(config),
        "best_validation_action_mse": best_loss,
        "selected_by": "held_out_greedy_label_mse_only",
        "used_exact_reference": False,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "parameters": sum(parameter.numel() for parameter in actor.parameters()),
    }, trace


@torch.no_grad()
def evaluate_policy_plain_mc(
    problem: ControlledHJB,
    policy: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 8,
    path_chunk: int = 64,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Ordinary payoff MC. There is deliberately no martingale control variate."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, problem.d)
    generator = torch.Generator(device=device).manual_seed(seed)
    outputs = {
        key: np.empty(len(times), dtype=np.float64)
        for key in ("value", "mc_variance_of_mean", "running_cost", "mean_abs_action", "saturation_rate")
    }
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(policy):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            sums = {key: torch.zeros(end - begin, dtype=torch.float64, device=device)
                    for key in ("payoff", "pair", "pair2", "running", "abs_action", "saturated")}
            action_observations = 0
            pair_observations = 0
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, problem.d),
                    generator=generator,
                    device=device,
                )
                dt = torch.clamp(problem.T - bt, min=0.0) / time_steps
                state = bx[:, None, :].expand(end - begin, path_chunk, problem.d)
                running = torch.zeros((end - begin, path_chunk), device=device)
                abs_action = torch.zeros_like(running)
                saturated = torch.zeros_like(running)
                for step in range(time_steps):
                    step_t = bt + step * dt
                    flat_t = step_t[:, None].expand(end - begin, path_chunk).reshape(-1)
                    flat_x = state.reshape(-1, problem.d)
                    action = policy.predict(flat_t, flat_x, problem.T).reshape(
                        end - begin, path_chunk, problem.m
                    )
                    cost = problem.running_cost(flat_t, flat_x, action.reshape(-1, problem.m)).reshape(
                        end - begin, path_chunk
                    )
                    drift = problem.drift(flat_t, flat_x, action.reshape(-1, problem.m)).reshape(
                        end - begin, path_chunk, problem.d
                    )
                    running += cost * dt[:, None]
                    abs_action += action.abs().mean(dim=-1)
                    if hasattr(problem, "A"):
                        saturated += (action[..., 0].abs() >= 0.99 * problem.A).float()
                    state = state + drift * dt[:, None, None] + problem.sigma * torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
                payoff = problem.terminal(state) + running
                half = path_chunk // 2
                pair_payoff = 0.5 * (payoff[:, :half] + payoff[:, half:])
                sums["payoff"] += payoff.double().sum(dim=1)
                sums["pair"] += pair_payoff.double().sum(dim=1)
                sums["pair2"] += pair_payoff.double().square().sum(dim=1)
                sums["running"] += running.double().sum(dim=1)
                sums["abs_action"] += abs_action.double().sum(dim=1)
                sums["saturated"] += saturated.double().sum(dim=1)
                action_observations += path_chunk * time_steps
                pair_observations += half
            mean = sums["payoff"] / paths
            pair_sample_var = torch.clamp(
                (sums["pair2"] - sums["pair"].square() / pair_observations)
                / max(pair_observations - 1, 1),
                min=0.0,
            )
            outputs["value"][begin:end] = mean.cpu().numpy()
            outputs["mc_variance_of_mean"][begin:end] = (
                pair_sample_var / pair_observations
            ).cpu().numpy()
            outputs["running_cost"][begin:end] = (sums["running"] / paths).cpu().numpy()
            outputs["mean_abs_action"][begin:end] = (sums["abs_action"] / action_observations).cpu().numpy()
            outputs["saturation_rate"][begin:end] = (sums["saturated"] / action_observations).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    return outputs, {
        "estimator": "ordinary_controlled_policy_payoff_mc",
        "martingale_control_variate": False,
        "antithetic_sampling": True,
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "mean_value_mc_se": float(np.sqrt(np.sum(outputs["mc_variance_of_mean"])) / len(times)),
    }


@torch.no_grad()
def actor_actions(
    problem: ControlledHJB, policy: nn.Module, t: np.ndarray, x: np.ndarray, device: torch.device
) -> np.ndarray:
    prediction = policy.predict(
        torch.as_tensor(t, dtype=torch.float32, device=device),
        torch.as_tensor(x, dtype=torch.float32, device=device),
        problem.T,
    )
    return prediction.cpu().numpy().astype(np.float64)

