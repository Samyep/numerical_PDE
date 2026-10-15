"""SCaSML level-one correction for an ordinary soft-terminal PINN."""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.scasml import (
    pinn_residual_hutchinson,
)

from .common import driver_for, synchronize
from .models import RawTanhPINN


def predict_scasml_original_level1(
    model: RawTanhPINN,
    benchmark: str,
    x: np.ndarray,
    device: torch.device,
    *,
    paths: int,
    probes: int,
    seed: int,
    parent_batch: int = 8,
    path_chunk: int = 8,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Estimate base + terminal defect + one sampled residual integral.

    The hard-terminal pilot omitted the terminal defect because it was exactly
    zero by construction. It is required for the corrected ordinary PINN.
    """
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("paths must be divisible by an even path_chunk")
    driver = driver_for(benchmark)
    states = np.asarray(x, dtype=np.float32)
    prediction = np.empty(len(states), dtype=np.float64)
    variance_of_mean = np.empty(len(states), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    residual_rms: list[float] = []
    terminal_defect_rms: list[float] = []
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    model.eval()
    for begin in range(0, len(states), parent_batch):
        end = min(begin + parent_batch, len(states))
        bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
        t0 = torch.zeros(len(bx), dtype=bx.dtype, device=device)
        with torch.no_grad():
            base = model.value(driver, t0, bx).double()
        payoff_sum = torch.zeros(len(bx), dtype=torch.float64, device=device)
        pair_sum = torch.zeros_like(payoff_sum)
        pair_square_sum = torch.zeros_like(payoff_sum)
        pair_count = 0
        for _ in range(0, paths, path_chunk):
            half = path_chunk // 2
            uniform_half = torch.clamp(
                torch.rand((len(bx), half), device=device, generator=generator),
                min=1e-5, max=1.0 - 1e-5,
            )
            uniform = torch.cat((uniform_half, uniform_half), dim=1)
            sample_t = driver.T * uniform
            normal_t_half = torch.randn(
                (len(bx), half, driver.d), device=device, generator=generator
            )
            normal_t = torch.cat((normal_t_half, -normal_t_half), dim=1)
            sample_x = (
                bx[:, None, :]
                + driver.sigma * torch.sqrt(sample_t)[:, :, None] * normal_t
            )
            normal_tail_half = torch.randn(
                (len(bx), half, driver.d), device=device, generator=generator
            )
            normal_tail = torch.cat((normal_tail_half, -normal_tail_half), dim=1)
            terminal_x = (
                sample_x
                + driver.sigma
                * torch.sqrt(driver.T - sample_t)[:, :, None]
                * normal_tail
            )
            with torch.enable_grad():
                residual, diagnostics = pinn_residual_hutchinson(
                    model, benchmark, sample_t.reshape(-1),
                    sample_x.reshape(-1, driver.d), probes=probes,
                    generator=generator,
                )
            with torch.no_grad():
                terminal_time = torch.full(
                    (len(bx) * path_chunk,), driver.T,
                    dtype=bx.dtype, device=device,
                )
                terminal_flat = terminal_x.reshape(-1, driver.d)
                terminal_defect = (
                    driver.terminal(terminal_flat)
                    - model.value(driver, terminal_time, terminal_flat)
                ).reshape(len(bx), path_chunk)
            payoff = terminal_defect + driver.T * residual.reshape(len(bx), path_chunk)
            payoff_sum += payoff.double().sum(dim=1)
            paired = 0.5 * (payoff[:, :half] + payoff[:, half:])
            pair_sum += paired.double().sum(dim=1)
            pair_square_sum += paired.double().square().sum(dim=1)
            pair_count += half
            residual_rms.append(diagnostics["residual_rms"])
            terminal_defect_rms.append(
                float(torch.sqrt(torch.mean(terminal_defect.square())).item())
            )
            nonfinite += int((~torch.isfinite(payoff)).sum().item())
        prediction[begin:end] = driver.working_to_original(
            (base + payoff_sum / paths).cpu().numpy()
        )
        pair_variance = torch.clamp(
            (pair_square_sum - pair_sum.square() / pair_count) / max(pair_count - 1, 1),
            min=0.0,
        )
        variance_of_mean[begin:end] = (pair_variance / pair_count).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    if nonfinite or not np.all(np.isfinite(prediction)):
        raise FloatingPointError("nonfinite corrected SCaSML prediction")
    metadata = {
        "method": "scasml_original", "benchmark": benchmark,
        "defect_picard_level": 1, "paths": paths,
        "hutchinson_probes": probes, "antithetic": True, "seed": seed,
        "inference_seconds": elapsed,
        "mean_mc_se_across_workload": float(
            np.sqrt(np.sum(variance_of_mean)) / len(states)
        ),
        "mean_sampled_residual_rms": float(np.mean(residual_rms)),
        "mean_terminal_defect_rms": float(np.mean(terminal_defect_rms)),
        "terminal_defect_included": True,
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda" else None,
        "adaptation": "SCaSML first defect correction with terminal and PDE defects",
    }
    return prediction, variance_of_mean, metadata


__all__ = ["predict_scasml_original_level1"]
