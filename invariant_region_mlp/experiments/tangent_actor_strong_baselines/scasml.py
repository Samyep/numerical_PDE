"""One-level SCaSML defect correction for a hard-terminal SDGD-PINN.

For a surrogate u_hat with exact terminal data, the error e=u-u_hat obeys a
semilinear PDE with source equal to the surrogate PDE residual.  The first
SCaSML/Picard correction is E[int R_hat(s,X_s) ds].  We estimate the trace in
R_hat with Hutchinson probes, a standard high-dimensional PINN adaptation.
"""

from __future__ import annotations

import math
import time
from typing import Any

import numpy as np
import torch

from .common import driver_for, stable_seed, synchronize
from .models import HardTerminalValue


def pinn_residual_hutchinson(
    model: HardTerminalValue,
    benchmark: str,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    probes: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, dict[str, float]]:
    driver = driver_for(benchmark)
    local_t = t.detach().requires_grad_(True)
    local_x = x.detach().requires_grad_(True)
    value = model.value(driver, local_t, local_x)
    u_t = torch.autograd.grad(
        value.sum(), local_t, create_graph=False, retain_graph=True
    )[0]
    gradient = torch.autograd.grad(
        value.sum(), local_x, create_graph=True, retain_graph=True
    )[0]
    trace = torch.zeros(len(local_x), dtype=local_x.dtype, device=local_x.device)
    for index in range(probes):
        direction = 2.0 * torch.randint(
            0,
            2,
            local_x.shape,
            device=local_x.device,
            generator=generator,
            dtype=torch.int64,
        ).to(local_x.dtype) - 1.0
        hessian_vector = torch.autograd.grad(
            torch.sum(gradient * direction),
            local_x,
            create_graph=False,
            retain_graph=index + 1 < probes,
        )[0]
        trace = trace + torch.sum(hessian_vector * direction, dim=-1)
    trace = trace / probes
    z = driver.sigma * gradient.detach()
    residual = u_t.detach() + trace.detach() + driver.driver(value.detach(), z)
    return residual, {
        "value_min": float(value.detach().min().item()),
        "value_max": float(value.detach().max().item()),
        "residual_rms": float(torch.sqrt(torch.mean(residual.square())).item()),
        "residual_max_abs": float(residual.abs().max().item()),
    }


def predict_scasml_level1(
    model: HardTerminalValue,
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
    if paths % path_chunk:
        raise ValueError("paths must be divisible by path_chunk")
    if path_chunk % 2:
        raise ValueError("path_chunk must be even for antithetic sampling")
    driver = driver_for(benchmark)
    states = np.asarray(x, dtype=np.float32)
    prediction = np.empty(len(states), dtype=np.float64)
    variance_of_mean = np.empty(len(states), dtype=np.float64)
    generator = torch.Generator(device=device).manual_seed(seed)
    residual_rms: list[float] = []
    residual_max = 0.0
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    model.eval()
    for begin in range(0, len(states), parent_batch):
        end = min(begin + parent_batch, len(states))
        bx = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
        bt = torch.zeros(len(bx), dtype=torch.float32, device=device)
        with torch.no_grad():
            base_value = model.value(driver, bt, bx).double()
        payoff_sum = torch.zeros(len(bx), dtype=torch.float64, device=device)
        pair_sum = torch.zeros_like(payoff_sum)
        pair_square_sum = torch.zeros_like(payoff_sum)
        pair_count = 0
        for _ in range(0, paths, path_chunk):
            half = path_chunk // 2
            uniform_half = torch.clamp(
                torch.rand((len(bx), half), device=device, generator=generator),
                min=1e-5,
            )
            uniform = torch.cat((uniform_half, uniform_half), dim=1)
            delta = driver.T * uniform
            normal_half = torch.randn(
                (len(bx), half, driver.d), device=device, generator=generator
            )
            normal = torch.cat((normal_half, -normal_half), dim=1)
            sample_x = bx[:, None, :] + driver.sigma * torch.sqrt(delta)[:, :, None] * normal
            sample_t = delta
            with torch.enable_grad():
                residual, diagnostics = pinn_residual_hutchinson(
                    model,
                    benchmark,
                    sample_t.reshape(-1),
                    sample_x.reshape(-1, driver.d),
                    probes=probes,
                    generator=generator,
                )
            residual = residual.reshape(len(bx), path_chunk)
            payoff = driver.T * residual
            payoff_sum += payoff.double().sum(dim=1)
            paired = 0.5 * (payoff[:, :half] + payoff[:, half:])
            pair_sum += paired.double().sum(dim=1)
            pair_square_sum += paired.double().square().sum(dim=1)
            pair_count += half
            residual_rms.append(diagnostics["residual_rms"])
            residual_max = max(residual_max, diagnostics["residual_max_abs"])
            nonfinite += int((~torch.isfinite(payoff)).sum().item())
        correction = payoff_sum / paths
        prediction[begin:end] = driver.working_to_original(
            (base_value + correction).cpu().numpy()
        )
        pair_variance = torch.clamp(
            (pair_square_sum - pair_sum.square() / pair_count) / max(pair_count - 1, 1),
            min=0.0,
        )
        variance_of_mean[begin:end] = (pair_variance / pair_count).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    if nonfinite or not np.all(np.isfinite(prediction)):
        raise FloatingPointError("nonfinite SCaSML prediction")
    metadata = {
        "method": "scasml_l1",
        "benchmark": benchmark,
        "defect_picard_level": 1,
        "paths": paths,
        "hutchinson_probes": probes,
        "antithetic": True,
        "seed": seed,
        "inference_seconds": elapsed,
        "mean_mc_se_across_workload": float(np.sqrt(np.sum(variance_of_mean)) / len(states)),
        "mean_sampled_residual_rms": float(np.mean(residual_rms)),
        "max_sampled_residual_abs": residual_max,
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "adaptation": "official SCaSML defect equation, first correction, Hutchinson Laplacian",
    }
    return prediction, variance_of_mean, metadata


__all__ = ["pinn_residual_hutchinson", "predict_scasml_level1"]
