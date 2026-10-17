"""Post-selection diagnostics that may use analytic benchmark references."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.special import roots_laguerre
import torch

from .benchmarks import DPIOUBenchmark, HJBLogBenchmark, OriginalBenchmark
from .common import synchronize
from .models import TerminalAwareActor
from .tangent_solver import antithetic_normals, frozen


def exact_z_torch(
    benchmark: OriginalBenchmark,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    quadrature_order: int = 64,
) -> torch.Tensor:
    """Return the independent exact ``z=sigma*grad(u)`` on evaluation paths."""
    times = t.reshape(-1)
    states = x.reshape(len(times), benchmark.d)
    if isinstance(benchmark, DPIOUBenchmark):
        return benchmark.sigma * benchmark.exact_gradient(times, states)
    if not isinstance(benchmark, HJBLogBenchmark):
        raise TypeError("an exact gradient is unavailable for this benchmark")

    tau = torch.clamp(benchmark.T - times, min=0.0)
    radius2 = states.square().sum(dim=1)
    nodes_np, weights_np = roots_laguerre(quadrature_order)
    nodes = torch.as_tensor(nodes_np, dtype=states.dtype, device=states.device)
    log_weights = torch.as_tensor(
        np.log(weights_np), dtype=states.dtype, device=states.device
    )
    scale = 1.0 + 2.0 * benchmark.d * tau + radius2
    scaled = nodes[None, :] / scale[:, None]
    denominator = 1.0 + 4.0 * tau[:, None] * scaled
    log_relative = (
        nodes[None, :]
        - scaled
        - 0.5 * benchmark.d * torch.log(denominator)
        - scaled * radius2[:, None] / denominator
    )
    probabilities = torch.softmax(log_weights[None, :] + log_relative, dim=1)
    ratio = (probabilities * scaled / denominator).sum(dim=1)
    gradient = 2.0 * states * ratio[:, None]
    terminal = tau <= 8.0 * torch.finfo(states.dtype).eps
    if terminal.any():
        gradient[terminal] = 2.0 * states[terminal] / (
            1.0 + radius2[terminal, None]
        )
    return benchmark.sigma * gradient


@torch.no_grad()
def quadratic_bregman_diagnostic(
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
) -> dict[str, Any]:
    """Estimate ``0.5 E^q int ||q-z*||^2 ds`` on controlled paths."""
    if not isinstance(benchmark, (HJBLogBenchmark, DPIOUBenchmark)):
        raise TypeError("quadratic Bregman diagnostic requires an exact z")
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, benchmark.d)
    generator = torch.Generator(device=device).manual_seed(seed)
    state_means = np.empty(len(times), dtype=np.float64)
    state_pair_variances = np.empty(len(times), dtype=np.float64)
    integrated_squared_error_sum = 0.0
    horizon_sum = 0.0
    nonfinite = 0

    synchronize(device)
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            dt = torch.clamp(benchmark.T - bt, min=0.0) / time_steps
            pair_chunks: list[torch.Tensor] = []
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    end - begin,
                    path_chunk,
                    time_steps,
                    benchmark.d,
                    generator=generator,
                    device=device,
                )
                state = bx[:, None, :].expand(end - begin, path_chunk, benchmark.d)
                gap = torch.zeros(
                    (end - begin, path_chunk), dtype=bx.dtype, device=device
                )
                for step in range(time_steps):
                    step_t = bt + step * dt
                    flat_t = step_t[:, None].expand(end - begin, path_chunk).reshape(-1)
                    flat_x = state.reshape(-1, benchmark.d)
                    _, q = actor.predict(flat_t, flat_x)
                    if q is None:
                        raise RuntimeError("quadratic benchmark returned no Gradient Actor")
                    exact_z = exact_z_torch(benchmark, flat_t, flat_x)
                    squared = (q - exact_z).square().sum(dim=1).reshape(
                        end - begin, path_chunk
                    )
                    gap = gap + 0.5 * squared * dt[:, None]
                    _, b, _ = benchmark.tangent_coefficients(None, q)
                    dw = torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
                    drift = benchmark.base_drift(flat_t, flat_x).reshape(
                        end - begin, path_chunk, benchmark.d
                    )
                    state = state + (
                        drift + benchmark.sigma * b.reshape_as(state)
                    ) * dt[:, None, None] + benchmark.sigma * dw
                half = path_chunk // 2
                pairs = 0.5 * (gap[:, :half] + gap[:, half:])
                pair_chunks.append(pairs.double())
                integrated_squared_error_sum += float((2.0 * gap).double().sum().item())
                horizon_sum += float(
                    (benchmark.T - bt).double().sum().item() * path_chunk
                )
                nonfinite += int((~torch.isfinite(gap)).sum().item())
            all_pairs = torch.cat(pair_chunks, dim=1)
            state_means[begin:end] = all_pairs.mean(dim=1).cpu().numpy()
            state_pair_variances[begin:end] = all_pairs.var(
                dim=1, unbiased=True
            ).cpu().numpy()
    synchronize(device)
    pair_count = paths // 2
    mean_gap = float(np.mean(state_means))
    path_mc_se = float(
        math.sqrt(np.sum(state_pair_variances / pair_count) / len(times) ** 2)
    )
    return {
        "identity": "v^q-u*=0.5 E^q integral ||q-z*||^2 ds",
        "states": len(times),
        "paths": paths,
        "antithetic_pairs": pair_count,
        "time_steps": time_steps,
        "seed": seed,
        "mean_integrated_bregman_gap": mean_gap,
        "state_standard_error": float(
            np.std(state_means, ddof=1) / math.sqrt(len(state_means))
        ) if len(state_means) > 1 else 0.0,
        "path_sampling_standard_error_of_state_mean": path_mc_se,
        "trajectory_vector_rmse": float(
            math.sqrt(integrated_squared_error_sum / max(horizon_sum, 1e-300))
        ),
        "nonfinite_count": nonfinite,
        "ground_truth_used_for_training_or_selection": False,
        "post_selection_only": True,
    }


__all__ = ["exact_z_torch", "quadratic_bregman_diagnostic"]
