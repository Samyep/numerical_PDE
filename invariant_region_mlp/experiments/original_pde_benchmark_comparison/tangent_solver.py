"""Linear tangent Feynman--Kac labels, fitting, and MCV evaluation."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass
import math
import time
from typing import Any, Iterator

import numpy as np
import torch
from torch import nn

from .benchmarks import OriginalBenchmark
from .common import seed_everything, synchronize
from .models import Activation, TerminalAwareActor


@dataclass(frozen=True)
class ActorConfig:
    training_points: int
    validation_points: int
    paths_per_label: int
    label_time_steps: int
    optimizer_steps: int
    width: int
    depth: int
    activation: Activation
    learning_rate: float
    batch_size: int
    parent_batch: int = 8
    path_chunk: int = 32
    checkpoint_interval: int = 250
    gradient_clip: float = 10.0
    training_distribution: str = "mixed"
    label_estimator: str = "mcv"


@contextmanager
def frozen(module: nn.Module) -> Iterator[None]:
    flags = [parameter.requires_grad for parameter in module.parameters()]
    mode = module.training
    for parameter in module.parameters():
        parameter.requires_grad_(False)
        parameter.grad = None
    module.eval()
    try:
        yield
    finally:
        for parameter, flag in zip(module.parameters(), flags, strict=True):
            parameter.requires_grad_(flag)
        module.train(mode)


def antithetic_normals(
    batch: int,
    paths: int,
    steps: int,
    dimension: int,
    *,
    generator: torch.Generator,
    device: torch.device,
) -> torch.Tensor:
    if paths % 2:
        raise ValueError("antithetic path count must be even")
    half = torch.randn(
        (batch, paths // 2, steps, dimension),
        generator=generator,
        dtype=torch.float32,
        device=device,
    )
    return torch.cat((half, -half), dim=1)


def _rollout_chunk(
    benchmark: OriginalBenchmark,
    actor: TerminalAwareActor,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    mcv: bool,
    need_gradient: bool,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor | None, dict[str, float]]:
    """Return plain/corrected payoffs and the corrected payoff gradient label."""
    batch, paths, steps, dimension = normals.shape
    if dimension != benchmark.d:
        raise ValueError("Brownian dimension mismatch")
    if mcv and benchmark.actor_kind != "gradient":
        raise ValueError("MCV needs an already trained Gradient Actor")
    initial = x.detach().clone().requires_grad_(need_gradient)
    times = t.detach().reshape(-1)
    dt = torch.clamp(benchmark.T - times, min=0.0) / steps
    state = initial[:, None, :].expand(batch, paths, dimension)
    weight = torch.ones((batch, paths), dtype=x.dtype, device=x.device)
    integral = torch.zeros_like(weight)
    martingale = torch.zeros_like(weight)
    c_min, c_max = float("inf"), -float("inf")
    nonfinite = 0
    context = torch.enable_grad() if need_gradient else torch.no_grad()
    with context:
        for step in range(steps):
            step_t = times + step * dt
            flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
            flat_x = state.reshape(-1, dimension)
            r, q = actor.predict(flat_t, flat_x)
            c, b, h = benchmark.tangent_coefficients(r, q)
            c = c.reshape(batch, paths)
            b = b.reshape(batch, paths, dimension)
            h = h.reshape(batch, paths)
            dw = torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
            if mcv:
                if q is None:
                    raise RuntimeError("missing q for MCV")
                martingale = martingale + weight * (
                    q.reshape(batch, paths, dimension) * dw
                ).sum(dim=-1)
            integral = integral + weight * h * dt[:, None]
            weight = weight * torch.exp(c * dt[:, None])
            drift = benchmark.base_drift(flat_t, flat_x).reshape(
                batch, paths, dimension
            )
            state = state + (drift + benchmark.sigma * b) * dt[:, None, None]
            state = state + benchmark.sigma * dw
            c_min = min(c_min, float(c.detach().min().item()))
            c_max = max(c_max, float(c.detach().max().item()))
        plain = weight * benchmark.terminal(state) + integral
        corrected = plain - martingale if mcv else plain
        if need_gradient:
            value = corrected.mean(dim=1)
            z = benchmark.sigma * torch.autograd.grad(value.sum(), initial)[0]
        else:
            z = None
    nonfinite += int((~torch.isfinite(plain)).sum().item())
    nonfinite += int((~torch.isfinite(corrected)).sum().item())
    nonfinite += int((~torch.isfinite(weight)).sum().item())
    return (
        plain.detach(),
        corrected.detach(),
        weight.detach(),
        None if z is None else z.detach(),
        {
            "c_min": c_min,
            "c_max": c_max,
            "nonfinite_count": float(nonfinite),
            "martingale_sum": float(martingale.detach().double().sum().item()),
            "martingale_square_sum": float(
                martingale.detach().double().square().sum().item()
            ),
        },
    )


def generate_labels(
    benchmark: OriginalBenchmark,
    source: TerminalAwareActor,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int,
    path_chunk: int,
    estimator: str,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Generate Value labels or z labels according to the benchmark Actor."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    use_mcv = estimator == "mcv" and benchmark.actor_kind == "gradient"
    if estimator not in {"plain", "mcv"}:
        raise KeyError(estimator)
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, benchmark.d)
    output_shape = (
        (len(times), benchmark.d)
        if benchmark.actor_kind == "gradient"
        else (len(times),)
    )
    labels = np.empty(output_shape, dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    martingale_sum = martingale_square_sum = 0.0
    c_min, c_max = float("inf"), -float("inf")
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(source):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            value_sum = torch.zeros(end - begin, dtype=torch.float64, device=device)
            z_sum = torch.zeros(
                (end - begin, benchmark.d), dtype=torch.float64, device=device
            )
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    end - begin,
                    path_chunk,
                    time_steps,
                    benchmark.d,
                    generator=generator,
                    device=device,
                )
                plain, corrected, weight, z, diagnostics = _rollout_chunk(
                    benchmark,
                    source,
                    bt,
                    bx,
                    normals,
                    mcv=use_mcv,
                    need_gradient=benchmark.actor_kind == "gradient",
                )
                value_sum += corrected.double().sum(dim=1)
                if z is not None:
                    z_sum += path_chunk * z.double()
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                martingale_sum += diagnostics["martingale_sum"]
                martingale_square_sum += diagnostics["martingale_square_sum"]
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                nonfinite += int(diagnostics["nonfinite_count"])
            if benchmark.actor_kind == "gradient":
                labels[begin:end] = (z_sum / paths).cpu().numpy()
            else:
                labels[begin:end] = (value_sum / paths).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    count = len(times) * paths
    if nonfinite or not np.all(np.isfinite(labels)):
        raise FloatingPointError("nonfinite tangent labels")
    weight_ess = weight_sum**2 / max(
        count * weight_square_sum, np.finfo(float).tiny
    )
    martingale_mean = martingale_sum / count
    martingale_variance = max(
        martingale_square_sum / count - martingale_mean**2, 0.0
    )
    return labels, {
        "kind": "differentiable_linear_tangent_feynman_kac",
        "actor_target": benchmark.actor_kind,
        "estimator": "direct_alpha1_mcv" if use_mcv else "plain_fk",
        "alpha": 1.0 if use_mcv else 0.0,
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "source_actor_frozen": True,
        "complete_corrected_payoff_differentiated": use_mcv,
        "left_endpoint_predictable_martingale": use_mcv,
        "analytic_nonterminal_reference_used": False,
        "ground_truth_used": False,
        "weight_ess_fraction": float(weight_ess),
        "c_min": c_min,
        "c_max": c_max,
        "martingale_path_mean": float(martingale_mean),
        "martingale_path_std": float(math.sqrt(martingale_variance)),
        "label_rms": float(np.sqrt(np.mean(np.square(labels)))),
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
    }


def fit_actor(
    benchmark: OriginalBenchmark,
    config: ActorConfig,
    source: TerminalAwareActor,
    train_t: np.ndarray,
    train_x: np.ndarray,
    train_labels: np.ndarray,
    valid_t: np.ndarray,
    valid_x: np.ndarray,
    valid_labels: np.ndarray,
    *,
    seed: int,
    device: torch.device,
) -> tuple[TerminalAwareActor, dict[str, Any], list[dict[str, float]]]:
    seed_everything(seed)
    actor = TerminalAwareActor(
        benchmark,
        width=config.width,
        depth=config.depth,
        activation=config.activation,
    ).to(device)
    if (
        source.width == actor.width
        and source.depth == actor.depth
        and source.activation_name == actor.activation_name
    ):
        actor.load_state_dict(source.state_dict())
        warm_started = True
    else:
        warm_started = False
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
        range(config.checkpoint_interval, config.optimizer_steps + 1, config.checkpoint_interval)
    )
    checkpoints.update((0, 1, config.optimizer_steps))
    trace: list[dict[str, float]] = []
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None

    def prediction(tensor_t: torch.Tensor, tensor_x: torch.Tensor) -> torch.Tensor:
        r, q = actor.predict(tensor_t, tensor_x)
        output = q if benchmark.actor_kind == "gradient" else r
        if output is None:
            raise RuntimeError("Actor returned the wrong target")
        return output

    def record(step: int) -> None:
        nonlocal best_loss, best_state
        actor.eval()
        with torch.no_grad():
            train_error = prediction(tt, tx) - ty
            valid_error = prediction(vt, vx) - vy
            train_mse = float(train_error.double().square().mean().item())
            valid_mse = float(valid_error.double().square().mean().item())
        normalized = valid_mse / scale**2
        trace.append(
            {
                "step": float(step),
                "learning_rate": float(optimizer.param_groups[0]["lr"]),
                "train_rmse": math.sqrt(train_mse),
                "validation_rmse": math.sqrt(valid_mse),
                "validation_normalized_mse": normalized,
            }
        )
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
        error = prediction(tt[index], tx[index]) - ty[index]
        loss = error.square().mean() / scale**2
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(actor.parameters(), config.gradient_clip)
        optimizer.step()
        scheduler.step()
        if step in checkpoints:
            record(step)
    if best_state is None:
        raise RuntimeError("no Actor checkpoint selected")
    actor.load_state_dict(best_state)
    actor.eval()
    synchronize(device)
    elapsed = time.perf_counter() - started
    selected = min(trace, key=lambda row: row["validation_normalized_mse"])
    return actor, {
        "config": asdict(config),
        "warm_started": warm_started,
        "target_scale": scale,
        "selected_step": int(selected["step"]),
        "selected_validation_rmse": selected["validation_rmse"],
        "selected_validation_normalized_mse": selected[
            "validation_normalized_mse"
        ],
        "selected_by": "minimum held-out noisy tangent-label error",
        "ground_truth_used": False,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        **actor.specification(),
    }, trace


@torch.no_grad()
def evaluate_policy(
    benchmark: OriginalBenchmark,
    actor: TerminalAwareActor,
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
    """Evaluate plain FK and, for Gradient Actors, direct-alpha MCV together."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("invalid evaluation path layout")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, benchmark.d)
    plain_mean = np.empty(len(times), dtype=np.float64)
    corrected_mean = np.empty(len(times), dtype=np.float64)
    plain_pair_variance = np.empty(len(times), dtype=np.float64)
    corrected_pair_variance = np.empty(len(times), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            plain_pairs: list[torch.Tensor] = []
            corrected_pairs: list[torch.Tensor] = []
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    end - begin,
                    path_chunk,
                    time_steps,
                    benchmark.d,
                    generator=generator,
                    device=device,
                )
                plain, corrected, weight, _, diagnostics = _rollout_chunk(
                    benchmark,
                    actor,
                    bt,
                    bx,
                    normals,
                    mcv=benchmark.actor_kind == "gradient",
                    need_gradient=False,
                )
                half = path_chunk // 2
                plain_pairs.append(0.5 * (plain[:, :half] + plain[:, half:]))
                corrected_pairs.append(
                    0.5 * (corrected[:, :half] + corrected[:, half:])
                )
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                nonfinite += int(diagnostics["nonfinite_count"])
            local_plain = torch.cat(plain_pairs, dim=1).double()
            local_corrected = torch.cat(corrected_pairs, dim=1).double()
            plain_mean[begin:end] = local_plain.mean(dim=1).cpu().numpy()
            corrected_mean[begin:end] = local_corrected.mean(dim=1).cpu().numpy()
            plain_pair_variance[begin:end] = local_plain.var(
                dim=1, unbiased=True
            ).cpu().numpy()
            corrected_pair_variance[begin:end] = local_corrected.var(
                dim=1, unbiased=True
            ).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    count = len(times) * paths
    return {
        "plain": plain_mean,
        "mcv": corrected_mean,
        "plain_pair_variance": plain_pair_variance,
        "mcv_pair_variance": corrected_pair_variance,
    }, {
        "paths": paths,
        "antithetic_pairs": paths // 2,
        "time_steps": time_steps,
        "seed": seed,
        "mcv": benchmark.actor_kind == "gradient",
        "mcv_alpha": 1.0 if benchmark.actor_kind == "gradient" else None,
        "weight_ess_fraction": float(
            weight_sum**2
            / max(count * weight_square_sum, np.finfo(float).tiny)
        ),
        "mean_plain_pair_variance": float(np.mean(plain_pair_variance)),
        "mean_mcv_pair_variance": float(np.mean(corrected_pair_variance)),
        "variance_reduction": float(
            np.mean(plain_pair_variance) / max(np.mean(corrected_pair_variance), 1e-300)
        ),
        "mean_plain_mc_se": float(
            np.mean(np.sqrt(plain_pair_variance / (paths // 2)))
        ),
        "mean_mcv_mc_se": float(
            np.mean(np.sqrt(corrected_pair_variance / (paths // 2)))
        ),
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
    }


__all__ = [
    "ActorConfig",
    "antithetic_normals",
    "evaluate_policy",
    "fit_actor",
    "frozen",
    "generate_labels",
    "_rollout_chunk",
]
