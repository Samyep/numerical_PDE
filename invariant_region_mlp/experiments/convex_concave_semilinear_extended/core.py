"""Linear weighted Feynman--Kac evaluation and supervised tangent updates."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import time
from typing import Any, Iterator

import numpy as np
import torch
from torch import nn

from .common import seed_everything
from .drivers import ConcaveSemilinearDriver
from .models import TangentActor


@dataclass(frozen=True)
class ActorConfig:
    training_points: int = 2048
    validation_points: int = 256
    paths_per_label: int = 128
    label_time_steps: int = 16
    optimizer_steps: int = 2500
    width: int = 192
    depth: int = 3
    learning_rate: float = 2e-4
    batch_size: int = 512
    parent_batch: int = 8
    path_chunk: int = 32


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
    shape: tuple[int, ...], *, generator: torch.Generator, device: torch.device
) -> torch.Tensor:
    if shape[1] % 2:
        raise ValueError("path count must be even")
    half = torch.randn(
        (shape[0], shape[1] // 2, *shape[2:]),
        generator=generator,
        dtype=torch.float32,
        device=device,
    )
    return torch.cat((half, -half), dim=1)


def _weighted_rollout_chunk(
    driver: ConcaveSemilinearDriver,
    actor: nn.Module,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    need_z: bool,
) -> tuple[torch.Tensor, torch.Tensor | None, torch.Tensor, dict[str, float]]:
    batch, paths, steps, dimension = normals.shape
    if dimension != driver.d:
        raise ValueError("Brownian dimension mismatch")
    initial = x.detach().clone().requires_grad_(need_z)
    times = t.detach()
    dt = torch.clamp(driver.T - times, min=0.0) / steps
    state = initial[:, None, :].expand(batch, paths, dimension)
    weight = torch.ones((batch, paths), dtype=x.dtype, device=x.device)
    source_integral = torch.zeros_like(weight)
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = 0
    context = torch.enable_grad() if need_z else torch.no_grad()
    with context:
        for step in range(steps):
            step_t = times + step * dt
            flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
            flat_x = state.reshape(-1, dimension)
            r, q = actor.predict(flat_t, flat_x, driver.T)
            c, b, h = driver.tangent_coefficients(r, q)
            c = c.reshape(batch, paths)
            b = b.reshape(batch, paths, dimension)
            h = h.reshape(batch, paths)
            c_min = min(c_min, float(c.detach().min().item()))
            c_max = max(c_max, float(c.detach().max().item()))
            c_sum += float(c.detach().double().sum().item())
            c_square_sum += float(c.detach().double().square().sum().item())
            c_count += c.numel()
            source_integral = source_integral + weight * h * dt[:, None]
            weight = weight * torch.exp(c * dt[:, None])
            state = (
                state
                + driver.sigma * b * dt[:, None, None]
                + driver.sigma * torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
            )
        payoff = weight * driver.terminal(state) + source_integral
        if need_z:
            value = payoff.mean(dim=1)
            z = driver.sigma * torch.autograd.grad(value.sum(), initial)[0]
        else:
            z = None
    diagnostics = {
        "c_min": c_min,
        "c_max": c_max,
        "c_sum": c_sum,
        "c_square_sum": c_square_sum,
        "c_count": c_count,
        "nonfinite_count": float(
            (~torch.isfinite(payoff)).sum().item() + (~torch.isfinite(weight)).sum().item()
        ),
    }
    return payoff.detach(), None if z is None else z.detach(), weight.detach(), diagnostics


def tangent_labels(
    driver: ConcaveSemilinearDriver,
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
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("invalid path chunk")
    times = np.asarray(t).reshape(-1)
    states = np.asarray(x).reshape(-1, driver.d)
    values = np.empty(len(times), dtype=np.float64)
    z_labels = np.empty_like(states) if driver.has_gradient else None
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = 0
    nonfinite = 0
    gradient_norms: list[np.ndarray] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(source):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            value_sum = torch.zeros(end - begin, dtype=torch.float64, device=device)
            z_sum = torch.zeros((end - begin, driver.d), dtype=torch.float64, device=device)
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                payoff, z, weight, diagnostics = _weighted_rollout_chunk(
                    driver, source, bt, bx, normals, need_z=driver.has_gradient
                )
                value_sum += payoff.double().sum(dim=1)
                if z is not None:
                    z_sum += path_chunk * z.double()
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                c_sum += diagnostics["c_sum"]
                c_square_sum += diagnostics["c_square_sum"]
                c_count += diagnostics["c_count"]
                nonfinite += int(diagnostics["nonfinite_count"])
            values[begin:end] = (value_sum / paths).cpu().numpy()
            if z_labels is not None:
                batch_z = (z_sum / paths).cpu().numpy()
                z_labels[begin:end] = batch_z
                gradient_norms.append(np.linalg.norm(batch_z, axis=1))
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    count = len(times) * paths
    if nonfinite or not np.all(np.isfinite(values)) or (
        z_labels is not None and not np.all(np.isfinite(z_labels))
    ):
        raise FloatingPointError("nonfinite weighted tangent labels")
    gradient_values = np.concatenate(gradient_norms) if gradient_norms else np.zeros(1)
    c_mean = c_sum / c_count
    c_variance = max(c_square_sum / c_count - c_mean**2, 0.0)
    weight_mean = weight_sum / count
    weight_variance = max(weight_square_sum / count - weight_mean**2, 0.0)
    return values, z_labels, {
        "kind": "weighted_linear_feynman_kac_full_trajectory_vjp",
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "source_actor_frozen": True,
        "state_dependence_preserved": True,
        "exponential_weight_differentiated": bool(driver.has_gradient),
        "analytic_nonterminal_reference_used": False,
        "c_min": c_min,
        "c_max": c_max,
        "c_mean": c_mean,
        "c_std": float(np.sqrt(c_variance)),
        "weight_min": weight_min,
        "weight_max": weight_max,
        "weight_mean": weight_mean,
        "weight_std": float(np.sqrt(weight_variance)),
        "weight_ess_fraction": weight_sum**2 / max(count * weight_square_sum, np.finfo(float).tiny),
        "gradient_norm_mean": float(np.mean(gradient_values)),
        "gradient_norm_max": float(np.max(gradient_values)),
        "nonfinite_count": nonfinite,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
    }


def make_targets(
    driver: ConcaveSemilinearDriver,
    values: np.ndarray,
    z: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray | None, dict[str, float]]:
    r = np.clip(values, driver.r_low, driver.r_high)
    diagnostics = {
        "raw_value_min": float(np.min(values)),
        "raw_value_max": float(np.max(values)),
        "value_projection_rate": float(np.mean((values < driver.r_low) | (values > driver.r_high))),
    }
    return r, z, diagnostics


def train_actor(
    driver: ConcaveSemilinearDriver,
    config: ActorConfig,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_r: np.ndarray,
    train_q: np.ndarray | None,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_r: np.ndarray,
    valid_q: np.ndarray | None,
    *,
    seed: int,
    device: torch.device,
    initial: TangentActor | None,
) -> tuple[TangentActor, dict[str, Any], list[dict[str, float]]]:
    seed_everything(seed)
    actor = TangentActor(driver, width=config.width, depth=config.depth).to(device)
    if initial is not None:
        actor.load_state_dict(initial.state_dict())
    tt = torch.as_tensor(train_t, dtype=torch.float32, device=device)
    tx = torch.as_tensor(train_x, dtype=torch.float32, device=device)
    tr = torch.as_tensor(train_r, dtype=torch.float32, device=device)
    vt = torch.as_tensor(valid_t, dtype=torch.float32, device=device)
    vx = torch.as_tensor(valid_x, dtype=torch.float32, device=device)
    vr = torch.as_tensor(valid_r, dtype=torch.float32, device=device)
    tq = None if train_q is None else torch.as_tensor(train_q, dtype=torch.float32, device=device)
    vq = None if valid_q is None else torch.as_tensor(valid_q, dtype=torch.float32, device=device)
    r_scale = max(float(np.std(train_r)), 1e-2)
    q_scale = (
        None if train_q is None else max(float(np.sqrt(np.mean(np.square(train_q)))), 1e-3)
    )
    optimizer = torch.optim.AdamW(actor.parameters(), lr=config.learning_rate, weight_decay=1e-6)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.optimizer_steps)
    generator = torch.Generator(device=device).manual_seed(seed + 1)
    checkpoints = sorted(set([1, config.optimizer_steps // 4, config.optimizer_steps // 2,
                              3 * config.optimizer_steps // 4, config.optimizer_steps]))
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    trace: list[dict[str, float]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for step in range(1, config.optimizer_steps + 1):
        index = torch.randint(len(tt), (min(config.batch_size, len(tt)),), generator=generator, device=device)
        predicted_r, predicted_q = actor.predict(tt[index], tx[index], driver.T)
        r_mse = (predicted_r - tr[index]).square().mean()
        loss = r_mse / r_scale**2
        if tq is not None:
            assert predicted_q is not None and q_scale is not None
            q_mse = (predicted_q - tq[index]).square().mean()
            loss = loss + q_mse / q_scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            actor.eval()
            with torch.no_grad():
                pr, pq = actor.predict(vt, vx, driver.T)
                valid_r_mse = float((pr - vr).square().mean().item())
                if vq is not None:
                    assert pq is not None and q_scale is not None
                    valid_q_mse = float((pq - vq).square().mean().item())
                else:
                    valid_q_mse = 0.0
                valid_loss = valid_r_mse / r_scale**2 + (
                    0.0 if q_scale is None else valid_q_mse / q_scale**2
                )
            actor.train()
            trace.append({
                "step": float(step),
                "validation_normalized_loss": valid_loss,
                "validation_value_mse": valid_r_mse,
                "validation_gradient_mse": valid_q_mse,
            })
            if valid_loss < best_loss:
                best_loss = valid_loss
                best_state = {key: value.detach().cpu().clone() for key, value in actor.state_dict().items()}
    if best_state is None:
        raise RuntimeError("no checkpoint")
    actor.load_state_dict(best_state)
    actor.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    selected = min(trace, key=lambda item: item["validation_normalized_loss"])
    return actor, {
        "config": asdict(config),
        "value_label_scale": r_scale,
        "gradient_label_scale": q_scale,
        "best_validation_normalized_loss": best_loss,
        "best_validation_value_mse": selected["validation_value_mse"],
        "best_validation_gradient_mse": selected["validation_gradient_mse"],
        "selected_step": int(selected["step"]),
        "selected_by": "normalized held-out tangent-label loss only",
        "used_exact_reference": False,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
        "parameters": sum(parameter.numel() for parameter in actor.parameters()),
    }, trace


@torch.no_grad()
def evaluate_tangent_mc(
    driver: ConcaveSemilinearDriver,
    actor: nn.Module,
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
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("invalid path chunk")
    times = np.asarray(t).reshape(-1)
    states = np.asarray(x).reshape(-1, driver.d)
    values = np.empty(len(times), dtype=np.float64)
    variance = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = 0
    nonfinite = 0
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            payoff_sum = torch.zeros(end - begin, dtype=torch.float64, device=device)
            pair_sum = torch.zeros_like(payoff_sum)
            pair_square_sum = torch.zeros_like(payoff_sum)
            pair_count = 0
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                payoff, _, weight, diagnostics = _weighted_rollout_chunk(
                    driver, actor, bt, bx, normals, need_z=False
                )
                half = path_chunk // 2
                pair = 0.5 * (payoff[:, :half] + payoff[:, half:])
                payoff_sum += payoff.double().sum(dim=1)
                pair_sum += pair.double().sum(dim=1)
                pair_square_sum += pair.double().square().sum(dim=1)
                pair_count += half
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                c_sum += diagnostics["c_sum"]
                c_square_sum += diagnostics["c_square_sum"]
                c_count += diagnostics["c_count"]
                nonfinite += int(diagnostics["nonfinite_count"])
            values[begin:end] = (payoff_sum / paths).cpu().numpy()
            pair_var = torch.clamp(
                (pair_square_sum - pair_sum.square() / pair_count) / max(pair_count - 1, 1),
                min=0.0,
            )
            variance[begin:end] = (pair_var / pair_count).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    if nonfinite or not np.all(np.isfinite(values)):
        raise FloatingPointError("nonfinite weighted MC evaluation")
    count = len(times) * paths
    c_mean = c_sum / c_count
    c_variance = max(c_square_sum / c_count - c_mean**2, 0.0)
    weight_mean = weight_sum / count
    weight_variance = max(weight_square_sum / count - weight_mean**2, 0.0)
    return {"working_value": values, "mc_variance_of_mean": variance}, {
        "estimator": "plain_weighted_linear_feynman_kac_mc",
        "martingale_control_variate": False,
        "antithetic": True,
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "c_min": c_min,
        "c_max": c_max,
        "c_mean": c_mean,
        "c_std": float(np.sqrt(c_variance)),
        "weight_min": weight_min,
        "weight_max": weight_max,
        "weight_mean": weight_mean,
        "weight_std": float(np.sqrt(weight_variance)),
        "weight_ess_fraction": weight_sum**2 / max(count * weight_square_sum, np.finfo(float).tiny),
        "nonfinite_count": nonfinite,
        "mean_value_mc_se": float(np.sqrt(np.sum(variance)) / len(times)),
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
    }


@torch.no_grad()
def actor_tangents(
    driver: ConcaveSemilinearDriver,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray | None]:
    r, q = actor.predict(
        torch.as_tensor(t, dtype=torch.float32, device=device),
        torch.as_tensor(x, dtype=torch.float32, device=device),
        driver.T,
    )
    return r.cpu().numpy().astype(np.float64), (
        None if q is None else q.cpu().numpy().astype(np.float64)
    )
