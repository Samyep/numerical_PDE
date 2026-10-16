"""Differentiable tangent labels with a predictable martingale correction."""

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


def _weighted_rollout_chunk_mcv_labels(
    driver: ConcaveSemilinearDriver,
    source: nn.Module,
    t: torch.Tensor,
    x: torch.Tensor,
    normals: torch.Tensor,
    *,
    need_z: bool,
) -> tuple[torch.Tensor, torch.Tensor | None, torch.Tensor, torch.Tensor, dict[str, float]]:
    """Evaluate P-H and differentiate the complete corrected pathwise payoff."""
    batch, paths, steps, dimension = normals.shape
    if dimension != driver.d:
        raise ValueError("Brownian dimension mismatch")
    initial = x.detach().clone().requires_grad_(need_z)
    times = t.detach()
    dt = torch.clamp(driver.T - times, min=0.0) / steps
    state = initial[:, None, :].expand(batch, paths, dimension)
    weight = torch.ones((batch, paths), dtype=x.dtype, device=x.device)
    source_integral = torch.zeros_like(weight)
    martingale = torch.zeros_like(weight)
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = 0
    context = torch.enable_grad() if need_z else torch.no_grad()
    with context:
        for step in range(steps):
            step_t = times + step * dt
            flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
            flat_x = state.reshape(-1, dimension)
            r, q = source.predict(flat_t, flat_x, driver.T)
            if q is None:
                raise ValueError("MCV training requires a Gradient Actor")
            c, b, h = driver.tangent_coefficients(r, q)
            c = c.reshape(batch, paths)
            b = b.reshape(batch, paths, dimension)
            h = h.reshape(batch, paths)
            q = q.reshape(batch, paths, dimension)
            c_min = min(c_min, float(c.detach().min().item()))
            c_max = max(c_max, float(c.detach().max().item()))
            c_sum += float(c.detach().double().sum().item())
            c_square_sum += float(c.detach().double().square().sum().item())
            c_count += c.numel()
            # q already represents sqrt(2) grad u. Use the predictable left
            # endpoint and unscaled dW, before updating C or X.
            dw = torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
            martingale = martingale + weight * torch.sum(q * dw, dim=-1)
            source_integral = source_integral + weight * h * dt[:, None]
            weight = weight * torch.exp(c * dt[:, None])
            state = state + driver.sigma * b * dt[:, None, None] + driver.sigma * dw
        plain_payoff = weight * driver.terminal(state) + source_integral
        corrected_payoff = plain_payoff - martingale
        if need_z:
            corrected_value = corrected_payoff.mean(dim=1)
            z = driver.sigma * torch.autograd.grad(
                corrected_value.sum(), initial
            )[0]
        else:
            z = None
    diagnostics = {
        "c_min": c_min,
        "c_max": c_max,
        "c_sum": c_sum,
        "c_square_sum": c_square_sum,
        "c_count": float(c_count),
        "nonfinite_count": float(
            (~torch.isfinite(corrected_payoff)).sum().item()
            + (~torch.isfinite(martingale)).sum().item()
            + (~torch.isfinite(weight)).sum().item()
        ),
    }
    return (
        corrected_payoff.detach(),
        None if z is None else z.detach(),
        weight.detach(),
        martingale.detach(),
        diagnostics,
    )


def mcv_tangent_labels(
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
    need_z: bool = True,
) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    """Generate direct-alpha MCV Value and full pathwise Gradient labels."""
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("invalid path chunk")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, driver.d)
    values = np.empty(len(times), dtype=np.float64)
    z_labels = np.empty_like(states) if need_z else None
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    correction_sum = correction_square_sum = 0.0
    c_min, c_max = float("inf"), -float("inf")
    c_sum = c_square_sum = 0.0
    c_count = nonfinite = 0
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
            z_sum = torch.zeros(
                (end - begin, driver.d), dtype=torch.float64, device=device
            )
            for _ in range(0, paths, path_chunk):
                normals = antithetic_normals(
                    (end - begin, path_chunk, time_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                payoff, z, weight, correction, diagnostics = (
                    _weighted_rollout_chunk_mcv_labels(
                        driver, source, bt, bx, normals, need_z=need_z
                    )
                )
                value_sum += payoff.double().sum(dim=1)
                if z is not None:
                    z_sum += path_chunk * z.double()
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                correction_sum += float(correction.double().sum().item())
                correction_square_sum += float(
                    correction.double().square().sum().item()
                )
                c_min = min(c_min, diagnostics["c_min"])
                c_max = max(c_max, diagnostics["c_max"])
                c_sum += diagnostics["c_sum"]
                c_square_sum += diagnostics["c_square_sum"]
                c_count += int(diagnostics["c_count"])
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
    if nonfinite or not np.all(np.isfinite(values)) or (
        z_labels is not None and not np.all(np.isfinite(z_labels))
    ):
        raise FloatingPointError("nonfinite MCV tangent labels")
    count = len(times) * paths
    c_mean = c_sum / c_count
    c_variance = max(c_square_sum / c_count - c_mean**2, 0.0)
    weight_mean = weight_sum / count
    weight_variance = max(weight_square_sum / count - weight_mean**2, 0.0)
    correction_mean = correction_sum / count
    correction_variance = max(
        correction_square_sum / count - correction_mean**2, 0.0
    )
    gradient_values = np.concatenate(gradient_norms) if gradient_norms else np.zeros(1)
    return values, z_labels, {
        "kind": "weighted_linear_feynman_kac_full_trajectory_vjp_mcv",
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "alpha": 1.0,
        "source_actor_frozen": True,
        "state_dependence_preserved": True,
        "complete_corrected_payoff_differentiated": bool(need_z),
        "mcv_q_detached": False,
        "mcv_weight_detached": False,
        "mcv_state_detached": False,
        "left_endpoint_preincrement": True,
        "q_scaled_again_by_sqrt2": False,
        "actor_forward_passes_per_step": 1,
        "analytic_nonterminal_reference_used": False,
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
        "martingale_path_mean": correction_mean,
        "martingale_path_std": float(np.sqrt(correction_variance)),
        "gradient_norm_mean": float(np.mean(gradient_values)),
        "gradient_norm_max": float(np.max(gradient_values)),
        "nonfinite_count": nonfinite,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": peak,
    }


__all__ = ["_weighted_rollout_chunk_mcv_labels", "mcv_tangent_labels"]
