"""Teacher-free controlled-SDE policy evaluation and actor improvement."""

from __future__ import annotations

from contextlib import contextmanager
import math
import time
from typing import Any, Iterator

import numpy as np
import torch
from torch import nn

from .actor import ActorNet
from .common import stable_seed
from .equation import B1Equation


def torch_terminal(equation: B1Equation, x: torch.Tensor) -> torch.Tensor:
    """Evaluate the specified terminal condition without an analytic solution."""

    directions = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
    log_weights = torch.as_tensor(equation.logc, dtype=x.dtype, device=x.device)
    logits = x @ directions.T + log_weights
    return -torch.logsumexp(logits, dim=-1)


def torch_terminal_z(equation: B1Equation, x: torch.Tensor) -> torch.Tensor:
    """Explicit terminal derivative, which is part of the PDE specification."""

    directions = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
    log_weights = torch.as_tensor(equation.logc, dtype=x.dtype, device=x.device)
    logits = x @ directions.T + log_weights
    probabilities = torch.softmax(logits, dim=-1)
    return -equation.sigma * (probabilities @ directions)


@contextmanager
def frozen_actor(actor: nn.Module) -> Iterator[None]:
    """Freeze parameters while preserving derivatives with respect to states."""

    flags = [parameter.requires_grad for parameter in actor.parameters()]
    was_training = actor.training
    for parameter in actor.parameters():
        parameter.requires_grad_(False)
        parameter.grad = None
    actor.eval()
    try:
        yield
    finally:
        for parameter, flag in zip(actor.parameters(), flags, strict=True):
            parameter.requires_grad_(flag)
        actor.train(was_training)


def antithetic_normals(
    shape: tuple[int, ...],
    *,
    generator: torch.Generator,
    device: torch.device,
    dtype: torch.dtype,
) -> torch.Tensor:
    """Draw normals antithetically along the path axis (dimension one)."""

    if len(shape) < 2 or shape[1] < 2 or shape[1] % 2:
        raise ValueError("the antithetic path count must be a positive even number")
    half_shape = (shape[0], shape[1] // 2, *shape[2:])
    half = torch.randn(
        half_shape,
        generator=generator,
        device=device,
        dtype=dtype,
    )
    return torch.cat((half, -half), dim=1)


def refine_antithetic_normals(
    coarse: torch.Tensor,
    refinement: int,
    *,
    generator: torch.Generator,
) -> torch.Tensor:
    """Brownian-bridge refinement that preserves every coarse increment.

    If ``Z`` is a coarse standard normal and ``r`` is the refinement, the
    refined normals satisfy ``sum_i Z_i / sqrt(r) = Z``.  Centered independent
    bridge innovations give each refined normal unit variance and zero
    cross-covariance.  Antithetic pairing is retained along the path axis.
    """

    if coarse.ndim != 4 or coarse.shape[1] % 2:
        raise ValueError("coarse normals must have an even antithetic path axis")
    if refinement < 1:
        raise ValueError("refinement must be positive")
    if refinement == 1:
        return coarse
    batch, paths, steps, dimension = coarse.shape
    half_paths = paths // 2
    coarse_half = coarse[:, :half_paths]
    innovations = torch.randn(
        (batch, half_paths, steps, refinement, dimension),
        generator=generator,
        device=coarse.device,
        dtype=coarse.dtype,
    )
    bridge = innovations - innovations.mean(dim=3, keepdim=True)
    fine_half = coarse_half.unsqueeze(3) / math.sqrt(refinement) + bridge
    fine_half = fine_half.reshape(
        batch, half_paths, steps * refinement, dimension
    )
    return torch.cat((fine_half, -fine_half), dim=1)


def coarsen_normals(normals: torch.Tensor, time_steps: int) -> torch.Tensor:
    """Aggregate standard-normal Brownian increments to ``time_steps``."""

    if normals.ndim != 4 or time_steps < 1:
        raise ValueError("normals must be rank four and time_steps positive")
    batch, paths, fine_steps, dimension = normals.shape
    if fine_steps % time_steps:
        raise ValueError("fine time steps must be divisible by target time steps")
    block = fine_steps // time_steps
    if block == 1:
        return normals
    reshaped = normals.reshape(
        batch, paths, time_steps, block, dimension
    )
    return reshaped.sum(dim=3) / math.sqrt(block)


def policy_payoffs_with_martingale_from_normals(
    equation: B1Equation,
    actor: ActorNet,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return raw, martingale-corrected, and stochastic-integral payoffs."""

    if normals.ndim != 4:
        raise ValueError("normals must have shape (batch, paths, steps, d)")
    batch, paths, steps, dimension = normals.shape
    if batch != len(x) or dimension != equation.d or len(t) != batch:
        raise ValueError("incompatible policy rollout shapes")
    if paths < 1 or steps < 1:
        raise ValueError("paths and time steps must be positive")

    dtype = next(actor.parameters()).dtype
    device = x.device
    times = t.detach().to(device=device, dtype=dtype).reshape(-1)
    initial = x.detach().to(device=device, dtype=dtype).reshape(
        batch, equation.d
    )
    draws = normals.to(device=device, dtype=dtype)
    dt = torch.clamp(equation.T - times, min=0.0) / steps
    root_dt = torch.sqrt(dt)

    with frozen_actor(actor), torch.no_grad():
        state = initial[:, None, :].expand(batch, paths, equation.d)
        running_cost = torch.zeros((batch, paths), dtype=dtype, device=device)
        martingale = torch.zeros((batch, paths), dtype=dtype, device=device)
        for step in range(steps):
            step_time = times + step * dt
            flat_time = step_time[:, None].expand(batch, paths).reshape(-1)
            flat_state = state.reshape(batch * paths, equation.d)
            policy = actor.predict(flat_time, flat_state, equation.T).reshape(
                batch, paths, equation.d
            )
            brownian_increment = root_dt[:, None, None] * draws[:, :, step, :]
            running_cost = running_cost + (
                0.5 * torch.sum(policy.square(), dim=-1) * dt[:, None]
            )
            martingale = martingale + torch.sum(
                policy * brownian_increment, dim=-1
            )
            state = (
                state
                - equation.sigma * policy * dt[:, None, None]
                + equation.sigma
                * root_dt[:, None, None]
                * draws[:, :, step, :]
            )
        raw = torch_terminal(equation, state) + running_cost
        corrected = raw - martingale
    return raw.detach(), corrected.detach(), martingale.detach()


def evaluate_policy_mc_with_martingale_cv(
    equation: B1Equation,
    actor: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    time_steps: int,
    paths: int,
    repetition: int,
    seed_namespace: str,
    device: torch.device,
    parent_batch_size: int,
    path_chunk_size: int,
    coarse_time_steps: int = 16,
    maximum_time_steps: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Evaluate paired raw/CV estimators and antithetic-pair variances.

    The returned variances estimate the variance of each finite-path sample
    mean.  They use independent antithetic-pair averages, so their scaling
    accounts for covariance within each antithetic pair.
    """

    if paths % path_chunk_size:
        raise ValueError("paths must be divisible by path_chunk_size")
    if path_chunk_size % 2:
        raise ValueError("path_chunk_size must be even")
    if maximum_time_steps % coarse_time_steps:
        raise ValueError("maximum_time_steps must refine coarse_time_steps")
    if time_steps < coarse_time_steps or maximum_time_steps % time_steps:
        raise ValueError("time_steps must be a supported nested resolution")

    times = np.array(t, dtype=np.float64, copy=True).reshape(-1)
    states = np.array(x, dtype=np.float64, copy=True, order="C").reshape(
        len(times), equation.d
    )
    raw_values = np.empty(len(times), dtype=np.float64)
    corrected_values = np.empty(len(times), dtype=np.float64)
    raw_estimator_variances = np.empty(len(times), dtype=np.float64)
    corrected_estimator_variances = np.empty(len(times), dtype=np.float64)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    pair_count = paths // 2
    refinement = maximum_time_steps // coarse_time_steps

    for batch_index, start in enumerate(range(0, len(times), parent_batch_size)):
        stop = min(start + parent_batch_size, len(times))
        batch_t = torch.as_tensor(
            times[start:stop], dtype=torch.float32, device=device
        )
        batch_x = torch.as_tensor(
            states[start:stop], dtype=torch.float32, device=device
        )
        shape = (len(batch_t),)
        raw_sum = torch.zeros(shape, dtype=torch.float64, device=device)
        corrected_sum = torch.zeros(shape, dtype=torch.float64, device=device)
        raw_pair_sum = torch.zeros(shape, dtype=torch.float64, device=device)
        corrected_pair_sum = torch.zeros(
            shape, dtype=torch.float64, device=device
        )
        raw_pair_square_sum = torch.zeros(
            shape, dtype=torch.float64, device=device
        )
        corrected_pair_square_sum = torch.zeros(
            shape, dtype=torch.float64, device=device
        )
        coarse_seed = stable_seed(
            "B1",
            "actor_iteration",
            seed_namespace,
            equation.d,
            equation.k,
            repetition,
            batch_index,
        )
        bridge_seed = stable_seed(
            "B1",
            "policy_evaluation_cv_bridge",
            seed_namespace,
            equation.d,
            equation.k,
            repetition,
            batch_index,
        )
        coarse_generator = torch.Generator(device=device).manual_seed(coarse_seed)
        bridge_generator = torch.Generator(device=device).manual_seed(bridge_seed)
        for _ in range(0, paths, path_chunk_size):
            coarse = antithetic_normals(
                (
                    len(batch_t),
                    path_chunk_size,
                    coarse_time_steps,
                    equation.d,
                ),
                generator=coarse_generator,
                device=device,
                dtype=torch.float32,
            )
            if time_steps == coarse_time_steps:
                normals = coarse
            else:
                fine = refine_antithetic_normals(
                    coarse, refinement, generator=bridge_generator
                )
                normals = coarsen_normals(fine, time_steps)
            raw, corrected, _ = policy_payoffs_with_martingale_from_normals(
                equation, actor, batch_t, batch_x, normals
            )
            raw_sum += raw.mean(dim=1).to(torch.float64) * path_chunk_size
            corrected_sum += (
                corrected.mean(dim=1).to(torch.float64) * path_chunk_size
            )
            half = path_chunk_size // 2
            raw_pairs = 0.5 * (raw[:, :half] + raw[:, half:])
            corrected_pairs = 0.5 * (
                corrected[:, :half] + corrected[:, half:]
            )
            raw_pairs = raw_pairs.to(torch.float64)
            corrected_pairs = corrected_pairs.to(torch.float64)
            raw_pair_sum += raw_pairs.sum(dim=1)
            corrected_pair_sum += corrected_pairs.sum(dim=1)
            raw_pair_square_sum += raw_pairs.square().sum(dim=1)
            corrected_pair_square_sum += corrected_pairs.square().sum(dim=1)

        raw_mean = raw_sum / paths
        corrected_mean = corrected_sum / paths
        raw_pair_variance = (
            raw_pair_square_sum - raw_pair_sum.square() / pair_count
        ) / (pair_count - 1)
        corrected_pair_variance = (
            corrected_pair_square_sum
            - corrected_pair_sum.square() / pair_count
        ) / (pair_count - 1)
        raw_values[start:stop] = raw_mean.cpu().numpy()
        corrected_values[start:stop] = corrected_mean.cpu().numpy()
        raw_estimator_variances[start:stop] = (
            torch.clamp(raw_pair_variance, min=0.0) / pair_count
        ).cpu().numpy()
        corrected_estimator_variances[start:stop] = (
            torch.clamp(corrected_pair_variance, min=0.0) / pair_count
        ).cpu().numpy()

    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak_gpu_bytes: int | None = int(torch.cuda.max_memory_allocated(device))
    else:
        peak_gpu_bytes = None
    metadata = {
        "time_steps": time_steps,
        "coarse_time_steps": coarse_time_steps,
        "maximum_time_steps": maximum_time_steps,
        "paths": paths,
        "antithetic_pairs": pair_count,
        "repetition": repetition,
        "seed_namespace": seed_namespace,
        "parent_batch_size": parent_batch_size,
        "path_chunk_size": path_chunk_size,
        "controlled_drift": "-sqrt(2) * q",
        "martingale_control_variate": "-integral q^T dW",
        "raw_and_corrected_share_paths": True,
        "additional_actor_calls_for_control_variate": 0,
        "nested_brownian_refinement": time_steps > coarse_time_steps,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak_gpu_bytes,
        "path_steps": len(times) * paths * time_steps,
        "actor_parameters_frozen": True,
    }
    arrays = (
        raw_values,
        corrected_values,
        raw_estimator_variances,
        corrected_estimator_variances,
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise FloatingPointError("nonfinite martingale-CV policy result")
    return (*arrays, metadata)


def policy_value_and_z_from_normals(
    equation: B1Equation,
    actor: ActorNet,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    compute_z: bool,
) -> tuple[torch.Tensor, torch.Tensor | None]:
    """Euler policy evaluation, optionally differentiated through every step.

    ``normals`` has shape ``(batch, paths, time_steps, d)``.  Summing the
    batch-level mean values before ``autograd.grad`` is a vector-Jacobian
    product; independent rollouts ensure that it returns one gradient per
    initial state without constructing a dense Jacobian.
    """

    if normals.ndim != 4:
        raise ValueError("normals must have shape (batch, paths, steps, d)")
    batch, paths, steps, dimension = normals.shape
    if batch != len(x) or dimension != equation.d or len(t) != batch:
        raise ValueError("incompatible policy rollout shapes")
    if paths < 1 or steps < 1:
        raise ValueError("paths and time steps must be positive")

    dtype = next(actor.parameters()).dtype
    device = x.device
    times = t.detach().to(device=device, dtype=dtype).reshape(-1)
    initial = x.detach().to(device=device, dtype=dtype).reshape(batch, equation.d)
    initial = initial.clone().requires_grad_(compute_z)
    draws = normals.to(device=device, dtype=dtype)
    dt = torch.clamp(equation.T - times, min=0.0) / steps

    context = torch.enable_grad() if compute_z else torch.no_grad()
    with frozen_actor(actor), context:
        state = initial[:, None, :].expand(batch, paths, equation.d)
        running_cost = torch.zeros((batch, paths), dtype=dtype, device=device)
        root_dt = torch.sqrt(dt)
        for step in range(steps):
            step_time = times + step * dt
            flat_time = step_time[:, None].expand(batch, paths).reshape(-1)
            flat_state = state.reshape(batch * paths, equation.d)
            policy = actor.predict(flat_time, flat_state, equation.T).reshape(
                batch, paths, equation.d
            )
            running_cost = running_cost + (
                0.5 * torch.sum(policy.square(), dim=-1) * dt[:, None]
            )
            state = (
                state
                - equation.sigma * policy * dt[:, None, None]
                + equation.sigma
                * root_dt[:, None, None]
                * draws[:, :, step, :]
            )
        samples = torch_terminal(equation, state) + running_cost
        value = samples.mean(dim=1)
        if compute_z:
            gradient = torch.autograd.grad(value.sum(), initial)[0]
            z = equation.sigma * gradient
        else:
            z = None
    return value.detach(), None if z is None else z.detach()


def evaluate_policy_mc(
    equation: B1Equation,
    actor: ActorNet,
    t: np.ndarray,
    x: np.ndarray,
    *,
    time_steps: int,
    paths: int,
    repetition: int,
    seed_namespace: str,
    device: torch.device,
    parent_batch_size: int,
    path_chunk_size: int,
    compute_z: bool,
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    """Chunked policy Monte Carlo with pathwise input gradients."""

    if paths % path_chunk_size:
        raise ValueError("paths must be divisible by path_chunk_size")
    if path_chunk_size % 2:
        raise ValueError("path_chunk_size must be even for antithetic sampling")
    # ``fixed_test_points`` can return a read-only view.  Own these buffers so
    # that ``torch.as_tensor`` never exposes PyTorch to non-writable storage.
    times = np.array(t, dtype=np.float64, copy=True).reshape(-1)
    states = np.array(x, dtype=np.float64, copy=True, order="C").reshape(
        len(times), equation.d
    )
    values = np.empty(len(times), dtype=np.float64)
    gradients = (
        np.empty((len(times), equation.d), dtype=np.float64)
        if compute_z
        else None
    )
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for batch_index, start in enumerate(range(0, len(times), parent_batch_size)):
        stop = min(start + parent_batch_size, len(times))
        batch_t = torch.as_tensor(times[start:stop], dtype=torch.float32, device=device)
        batch_x = torch.as_tensor(
            states[start:stop], dtype=torch.float32, device=device
        )
        value_sum = torch.zeros(len(batch_t), dtype=torch.float64, device=device)
        z_sum = (
            torch.zeros(
                (len(batch_t), equation.d), dtype=torch.float64, device=device
            )
            if compute_z
            else None
        )
        seed = stable_seed(
            "B1",
            "actor_iteration",
            seed_namespace,
            equation.d,
            equation.k,
            repetition,
            batch_index,
        )
        generator = torch.Generator(device=device).manual_seed(seed)
        for _ in range(0, paths, path_chunk_size):
            normals = antithetic_normals(
                (len(batch_t), path_chunk_size, time_steps, equation.d),
                generator=generator,
                device=device,
                dtype=torch.float32,
            )
            chunk_value, chunk_z = policy_value_and_z_from_normals(
                equation,
                actor,
                batch_t,
                batch_x,
                normals,
                compute_z=compute_z,
            )
            value_sum += chunk_value.to(torch.float64) * path_chunk_size
            if compute_z:
                if chunk_z is None or z_sum is None:
                    raise RuntimeError("pathwise gradient was not returned")
                z_sum += chunk_z.to(torch.float64) * path_chunk_size
        values[start:stop] = (value_sum / paths).cpu().numpy()
        if compute_z:
            if gradients is None or z_sum is None:
                raise RuntimeError("gradient output was not allocated")
            gradients[start:stop] = (z_sum / paths).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak_gpu_bytes: int | None = int(torch.cuda.max_memory_allocated(device))
    else:
        peak_gpu_bytes = None
    elapsed = time.perf_counter() - started
    metadata = {
        "time_steps": time_steps,
        "paths": paths,
        "repetition": repetition,
        "seed_namespace": seed_namespace,
        "parent_batch_size": parent_batch_size,
        "path_chunk_size": path_chunk_size,
        "compute_z": compute_z,
        "controlled_drift": "-sqrt(2) * q",
        "wall_seconds": elapsed,
        "peak_gpu_bytes": peak_gpu_bytes,
        "path_steps": len(times) * paths * time_steps,
        "actor_parameters_frozen": True,
        "state_graph_detached": False,
    }
    if not np.all(np.isfinite(values)):
        raise FloatingPointError("nonfinite policy values")
    if gradients is not None and not np.all(np.isfinite(gradients)):
        raise FloatingPointError("nonfinite pathwise policy gradients")
    return values, gradients, metadata


def terminal_pathwise_labels(
    equation: B1Equation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    repetition: int,
    seed_namespace: str,
    device: torch.device,
    parent_batch_size: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Estimate E[sqrt(2) grad g(X_T)] using terminal derivatives only."""

    if paths % 2:
        raise ValueError("terminal path count must be even")
    times = np.array(t, dtype=np.float64, copy=True).reshape(-1)
    states = np.array(x, dtype=np.float64, copy=True, order="C").reshape(
        len(times), equation.d
    )
    labels = np.empty((len(times), equation.d), dtype=np.float64)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for batch_index, start in enumerate(range(0, len(times), parent_batch_size)):
        stop = min(start + parent_batch_size, len(times))
        batch_t = torch.as_tensor(times[start:stop], dtype=torch.float32, device=device)
        batch_x = torch.as_tensor(
            states[start:stop], dtype=torch.float32, device=device
        )
        seed = stable_seed(
            "B1",
            "actor_iteration",
            seed_namespace,
            equation.d,
            equation.k,
            repetition,
            batch_index,
        )
        generator = torch.Generator(device=device).manual_seed(seed)
        normals = antithetic_normals(
            (len(batch_t), paths, equation.d),
            generator=generator,
            device=device,
            dtype=torch.float32,
        )
        tau = torch.clamp(equation.T - batch_t, min=0.0)
        terminal_state = (
            batch_x[:, None, :]
            + equation.sigma * torch.sqrt(tau)[:, None, None] * normals
        )
        with torch.no_grad():
            terminal_z = torch_terminal_z(equation, terminal_state)
            labels[start:stop] = terminal_z.mean(dim=1).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak_gpu_bytes: int | None = int(torch.cuda.max_memory_allocated(device))
    else:
        peak_gpu_bytes = None
    metadata = {
        "points": len(times),
        "paths": paths,
        "repetition": repetition,
        "seed_namespace": seed_namespace,
        "parent_batch_size": parent_batch_size,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak_gpu_bytes,
        "terminal_samples": len(times) * paths,
        "label_source": "terminal pathwise gradient only",
        "nonterminal_mlp_used": False,
        "nonterminal_exact_gradient_used": False,
    }
    if not np.all(np.isfinite(labels)):
        raise FloatingPointError("nonfinite terminal pathwise labels")
    return labels, metadata


def rms_prediction_change(current: np.ndarray, previous: np.ndarray) -> float:
    difference = np.asarray(current, dtype=np.float64) - np.asarray(
        previous, dtype=np.float64
    )
    scale = math.sqrt(float(np.mean(np.sum(np.asarray(current) ** 2, axis=-1))))
    change = math.sqrt(float(np.mean(np.sum(difference**2, axis=-1))))
    return change / max(scale, np.finfo(np.float64).tiny)
