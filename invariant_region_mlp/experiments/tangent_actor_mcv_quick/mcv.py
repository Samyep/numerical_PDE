"""Predictable left-endpoint martingale control variate for weighted FK paths."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    antithetic_normals,
    frozen,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)


def antithetic_pair_average(paths: torch.Tensor) -> torch.Tensor:
    """Pair first-half draws with their negated second-half counterparts."""
    if paths.ndim < 2 or paths.shape[1] % 2:
        raise ValueError("antithetic path axis must have even length")
    half = paths.shape[1] // 2
    return 0.5 * (paths[:, :half] + paths[:, half:])


def variance_of_antithetic_mean(paths: torch.Tensor) -> torch.Tensor:
    """Variance of a path mean using independent antithetic pair averages."""
    pairs = antithetic_pair_average(paths).double()
    if pairs.shape[1] < 2:
        return torch.zeros(pairs.shape[0], dtype=torch.float64, device=pairs.device)
    return pairs.var(dim=1, unbiased=True) / pairs.shape[1]


def calibrated_alpha(
    payoff_pairs: np.ndarray,
    martingale_pairs: np.ndarray,
    *,
    denominator_tolerance: float = 1e-20,
) -> tuple[float, dict[str, float]]:
    """Estimate one scalar coefficient after centering within each state."""
    payoff = np.asarray(payoff_pairs, dtype=np.float64)
    martingale = np.asarray(martingale_pairs, dtype=np.float64)
    if payoff.shape != martingale.shape or payoff.ndim != 2:
        raise ValueError((payoff.shape, martingale.shape))
    centered_payoff = payoff - payoff.mean(axis=1, keepdims=True)
    centered_martingale = martingale - martingale.mean(axis=1, keepdims=True)
    numerator = float(np.sum(centered_payoff * centered_martingale))
    denominator = float(np.sum(centered_martingale**2))
    alpha = 0.0 if denominator <= denominator_tolerance else numerator / denominator
    return alpha, {
        "covariance_numerator": numerator,
        "variance_denominator": denominator,
        "denominator_tolerance": denominator_tolerance,
        "denominator_was_zero": denominator <= denominator_tolerance,
    }


def _weighted_rollout_chunk_mcv(
    driver: ConcaveSemilinearDriver,
    actor: nn.Module,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    collect_martingale: bool,
) -> tuple[torch.Tensor, torch.Tensor | None, torch.Tensor, dict[str, float]]:
    """Run the unchanged tangent FK dynamics and optionally accumulate H."""
    batch, paths, steps, dimension = normals.shape
    if dimension != driver.d:
        raise ValueError("Brownian dimension mismatch")
    times = t.detach()
    dt = torch.clamp(driver.T - times, min=0.0) / steps
    state = x.detach()[:, None, :].expand(batch, paths, dimension)
    weight = torch.ones((batch, paths), dtype=x.dtype, device=x.device)
    source_integral = torch.zeros_like(weight)
    martingale = torch.zeros_like(weight) if collect_martingale else None
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = 0
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
        # dW is the unscaled Brownian increment. q already equals sqrt(2) grad u.
        dw = torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
        if collect_martingale:
            if q is None:
                raise ValueError("MCV requires a Gradient Actor")
            assert martingale is not None
            q_view = q.reshape(batch, paths, dimension)
            # Predictable left endpoint: C_m, q(t_m, X_m), then increment dW_m.
            martingale = martingale + weight * torch.sum(q_view * dw, dim=-1)
        source_integral = source_integral + weight * h * dt[:, None]
        weight = weight * torch.exp(c * dt[:, None])
        state = state + driver.sigma * b * dt[:, None, None] + driver.sigma * dw
    payoff = weight * driver.terminal(state) + source_integral
    diagnostics = {
        "c_min": c_min,
        "c_max": c_max,
        "c_sum": c_sum,
        "c_square_sum": c_square_sum,
        "c_count": float(c_count),
        "nonfinite_count": float(
            (~torch.isfinite(payoff)).sum().item()
            + (~torch.isfinite(weight)).sum().item()
            + (0 if martingale is None else (~torch.isfinite(martingale)).sum().item())
        ),
    }
    return payoff.detach(), None if martingale is None else martingale.detach(), weight.detach(), diagnostics


@torch.no_grad()
def evaluate_pair_samples(
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
    collect_martingale: bool = True,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Return independent antithetic pair samples for P and H."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, driver.d)
    pair_count = paths // 2
    payoff_pairs = np.empty((len(times), pair_count), dtype=np.float64)
    martingale_pairs = (
        np.empty_like(payoff_pairs) if collect_martingale else None
    )
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    raw_payoff_var_sum = pair_payoff_var_sum = 0.0
    raw_h_var_sum = pair_h_var_sum = 0.0
    pair_cov_sum = pair_p_centered_sum = pair_h_centered_sum = 0.0
    pair_moment_count = 0
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = nonfinite = 0
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            pair_offset = 0
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                payoff, martingale, weight, diagnostics = _weighted_rollout_chunk_mcv(
                    driver, actor, bt, bx, normals,
                    collect_martingale=collect_martingale,
                )
                payoff_pair = antithetic_pair_average(payoff)
                local_pairs = path_chunk // 2
                payoff_pairs[begin:end, pair_offset:pair_offset + local_pairs] = (
                    payoff_pair.double().cpu().numpy()
                )
                raw_payoff_var_sum += float(
                    payoff.double().var(dim=1, unbiased=True).sum().item()
                )
                pair_payoff_var_sum += float(
                    payoff_pair.double().var(dim=1, unbiased=True).sum().item()
                )
                if collect_martingale:
                    assert martingale is not None and martingale_pairs is not None
                    martingale_pair = antithetic_pair_average(martingale)
                    martingale_pairs[begin:end, pair_offset:pair_offset + local_pairs] = (
                        martingale_pair.double().cpu().numpy()
                    )
                    raw_h_var_sum += float(
                        martingale.double().var(dim=1, unbiased=True).sum().item()
                    )
                    pair_h_var_sum += float(
                        martingale_pair.double().var(dim=1, unbiased=True).sum().item()
                    )
                    p_centered = payoff_pair.double() - payoff_pair.double().mean(dim=1, keepdim=True)
                    h_centered = martingale_pair.double() - martingale_pair.double().mean(dim=1, keepdim=True)
                    pair_cov_sum += float((p_centered * h_centered).sum().item())
                    pair_p_centered_sum += float(p_centered.square().sum().item())
                    pair_h_centered_sum += float(h_centered.square().sum().item())
                    pair_moment_count += p_centered.numel()
                pair_offset += local_pairs
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                c_sum += diagnostics["c_sum"]
                c_square_sum += diagnostics["c_square_sum"]
                c_count += int(diagnostics["c_count"])
                nonfinite += int(diagnostics["nonfinite_count"])
            if pair_offset != pair_count:
                raise AssertionError((pair_offset, pair_count))
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    elapsed = time.perf_counter() - started
    if nonfinite or not np.all(np.isfinite(payoff_pairs)) or (
        martingale_pairs is not None and not np.all(np.isfinite(martingale_pairs))
    ):
        raise FloatingPointError("nonfinite MCV rollout")
    count = len(times) * paths
    c_mean = c_sum / c_count
    c_variance = max(c_square_sum / c_count - c_mean**2, 0.0)
    weight_mean = weight_sum / count
    weight_variance = max(weight_square_sum / count - weight_mean**2, 0.0)
    payoff_antithetic_ratio = pair_payoff_var_sum / max(
        0.5 * raw_payoff_var_sum, np.finfo(float).tiny
    )
    if pair_h_centered_sum > 0.0 and pair_p_centered_sum > 0.0:
        pair_correlation = pair_cov_sum / np.sqrt(
            pair_p_centered_sum * pair_h_centered_sum
        )
    else:
        pair_correlation = 0.0
    result = {"payoff_pairs": payoff_pairs}
    if martingale_pairs is not None:
        result["martingale_pairs"] = martingale_pairs
    return result, {
        "estimator": "plain_weighted_linear_feynman_kac_plus_optional_mcv",
        "martingale_control_variate": collect_martingale,
        "actor_forward_passes_per_step": 1,
        "q_scaled_again_by_sqrt2": False,
        "left_endpoint_preincrement": True,
        "antithetic": True,
        "independent_units": "antithetic pair averages",
        "paths": paths,
        "pairs": pair_count,
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
        "weight_ess_fraction": weight_sum**2 / max(
            count * weight_square_sum, np.finfo(float).tiny
        ),
        "payoff_antithetic_variance_ratio": payoff_antithetic_ratio,
        "martingale_antithetic_variance_ratio": (
            None if not collect_martingale else
            (pair_h_var_sum / max(0.5 * raw_h_var_sum, np.finfo(float).tiny))
        ),
        "payoff_martingale_pair_correlation": pair_correlation,
        "nonfinite_count": nonfinite,
        "wall_seconds": elapsed,
        "peak_gpu_bytes": peak,
    }


__all__ = [
    "_weighted_rollout_chunk_mcv", "antithetic_pair_average", "calibrated_alpha",
    "evaluate_pair_samples", "variance_of_antithetic_mean",
]
