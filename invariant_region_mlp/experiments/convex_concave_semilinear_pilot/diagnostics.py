"""Post-training Bregman-regret diagnostics using the independent reference."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .core import antithetic_normals, frozen
from .drivers import ConcaveSemilinearDriver
from .reference import PeriodicReference, ReferenceTangentActor


@torch.no_grad()
def bregman_identity_audit(
    driver: ConcaveSemilinearDriver,
    actor: nn.Module,
    reference: PeriodicReference,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    time_steps: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 4,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Compare v^A-u* with the weighted integral of D_{-f}(Y*,A)."""
    if paths % 2:
        raise ValueError("paths must be even")
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(-1, driver.d)
    exact_actor = ReferenceTangentActor(driver, reference).to(device)
    generator = torch.Generator(device=device).manual_seed(seed)
    lhs_values = np.empty(len(times), dtype=np.float64)
    rhs_values = np.empty(len(times), dtype=np.float64)
    lhs_variances = np.empty(len(times), dtype=np.float64)
    rhs_variances = np.empty(len(times), dtype=np.float64)
    divergence_min = float("inf")
    divergence_negative = divergence_count = 0
    weight_min, weight_max = float("inf"), -float("inf")
    nonfinite = 0
    started = time.perf_counter()
    with frozen(actor), frozen(exact_actor):
        for begin in range(0, len(times), parent_batch):
            end = min(begin + parent_batch, len(times))
            bt = torch.as_tensor(times[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
            batch = end - begin
            dt = torch.clamp(driver.T - bt, min=0.0) / time_steps
            normals = antithetic_normals(
                (batch, paths, time_steps, driver.d), generator=generator, device=device
            )
            state = bx[:, None, :].expand(batch, paths, driver.d)
            weight = torch.ones((batch, paths), dtype=bx.dtype, device=device)
            source = torch.zeros_like(weight)
            rhs = torch.zeros_like(weight)
            initial_exact, _ = exact_actor.predict(bt, bx, driver.T)
            for step in range(time_steps):
                step_t = bt + step * dt
                flat_t = step_t[:, None].expand(batch, paths).reshape(-1)
                flat_x = state.reshape(-1, driver.d)
                r, q = actor.predict(flat_t, flat_x, driver.T)
                exact_y, exact_z = exact_actor.predict(flat_t, flat_x, driver.T)
                c, b, h = driver.tangent_coefficients(r, q)
                divergence = driver.bregman_divergence(exact_y, exact_z, r, q)
                c = c.reshape(batch, paths)
                b = b.reshape(batch, paths, driver.d)
                h = h.reshape(batch, paths)
                divergence = divergence.reshape(batch, paths)
                divergence_min = min(divergence_min, float(divergence.min().item()))
                divergence_negative += int((divergence < -2e-6).sum().item())
                divergence_count += divergence.numel()
                rhs = rhs + weight * divergence * dt[:, None]
                source = source + weight * h * dt[:, None]
                weight = weight * torch.exp(c * dt[:, None])
                state = (
                    state
                    + driver.sigma * b * dt[:, None, None]
                    + driver.sigma * torch.sqrt(dt)[:, None, None] * normals[:, :, step, :]
                )
            payoff = weight * driver.terminal(state) + source
            lhs = payoff - initial_exact[:, None]
            lhs_values[begin:end] = lhs.double().mean(dim=1).cpu().numpy()
            rhs_values[begin:end] = rhs.double().mean(dim=1).cpu().numpy()
            half = paths // 2
            lhs_pairs = 0.5 * (lhs[:, :half] + lhs[:, half:])
            rhs_pairs = 0.5 * (rhs[:, :half] + rhs[:, half:])
            lhs_variances[begin:end] = (
                lhs_pairs.double().var(dim=1, unbiased=True).cpu().numpy() / half
            )
            rhs_variances[begin:end] = (
                rhs_pairs.double().var(dim=1, unbiased=True).cpu().numpy() / half
            )
            weight_min = min(weight_min, float(weight.min().item()))
            weight_max = max(weight_max, float(weight.max().item()))
            nonfinite += int((~torch.isfinite(payoff)).sum().item())
            nonfinite += int((~torch.isfinite(rhs)).sum().item())
    discrepancy = lhs_values - rhs_values
    scale = max(float(np.sqrt(np.mean(lhs_values**2))), np.finfo(float).tiny)
    summary = {
        "identity": "v^A-u*=E^A integral C_s D_{-f}(Y*,A) ds",
        "posthoc_reference_only": True,
        "points": len(times),
        "paths": paths,
        "time_steps": time_steps,
        "seed": seed,
        "lhs_mean": float(np.mean(lhs_values)),
        "rhs_mean": float(np.mean(rhs_values)),
        "mean_discrepancy": float(np.mean(discrepancy)),
        "discrepancy_rmse_across_states": float(np.sqrt(np.mean(discrepancy**2))),
        "relative_discrepancy_rmse": float(np.sqrt(np.mean(discrepancy**2)) / scale),
        "lhs_mean_mc_se": float(np.sqrt(np.sum(lhs_variances)) / len(times)),
        "rhs_mean_mc_se": float(np.sqrt(np.sum(rhs_variances)) / len(times)),
        "correlation_across_states": float(np.corrcoef(lhs_values, rhs_values)[0, 1]),
        "divergence_min": divergence_min,
        "divergence_negative_fraction_below_minus_2e-6": divergence_negative / divergence_count,
        "terminal_weight_min": weight_min,
        "terminal_weight_max": weight_max,
        "nonfinite_count": nonfinite,
        "wall_seconds": time.perf_counter() - started,
    }
    arrays = {
        "t": times,
        "x": states,
        "lhs_value_gap": lhs_values,
        "rhs_weighted_bregman": rhs_values,
        "lhs_mc_variance_of_mean": lhs_variances,
        "rhs_mc_variance_of_mean": rhs_variances,
    }
    return summary, arrays
