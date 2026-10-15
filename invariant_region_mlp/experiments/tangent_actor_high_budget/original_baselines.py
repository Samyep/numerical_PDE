"""Compute-matched strong baselines in their original, non-terminal-aware form."""

from __future__ import annotations

import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)

from .common import (
    BUDGETS,
    CHECKPOINT_DIR,
    RAW_DIR,
    atomic_json,
    atomic_torch_save,
    driver_for,
    parameter_count,
    seed_everything,
    stable_seed,
    synchronize,
)
from .models import OriginalAmortizedDeepBSDE, RawTanhPINN, architecture_payload


def _cpu_state(model: nn.Module) -> dict[str, torch.Tensor]:
    return {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}


def _gradient_norm(parameters: Any) -> float:
    square = 0.0
    for parameter in parameters:
        if parameter.grad is not None:
            square += float(parameter.grad.detach().double().square().sum().item())
    return math.sqrt(square)


def _elapsed(started: float, device: torch.device) -> float:
    synchronize(device)
    return time.perf_counter() - started


def _set_timed_lr(
    optimizer: torch.optim.Optimizer, base_lr: float, seconds: float
) -> float:
    factor = 1.0 if seconds < 480.0 else (0.2 if seconds < 800.0 else 0.04)
    value = base_lr * factor
    for group in optimizer.param_groups:
        group["lr"] = value
    return value


def _clipped_square(delta: torch.Tensor, clip: float = 50.0) -> torch.Tensor:
    absolute = delta.abs()
    return torch.where(absolute < clip, delta.square(), 2 * clip * absolute - clip**2).mean()


@torch.no_grad()
def deep_bsde_validation(
    model: OriginalAmortizedDeepBSDE,
    driver: ConcaveSemilinearDriver,
    x: torch.Tensor,
    normals: torch.Tensor,
) -> tuple[float, int]:
    model.eval()
    terminal_y, terminal_x, nonfinite = model.rollout(driver, x, normals)
    return float(_clipped_square(terminal_y - driver.terminal(terminal_x)).item()), nonfinite


def train_deep_bsde_original(
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    budgets: tuple[float, ...] = BUDGETS,
    time_steps: int = 16,
    batch_size: int = 128,
    learning_rate: float = 5e-4,
    validation_size: int = 512,
    validation_every: int = 25,
) -> dict[str, Any]:
    """Official nonshared Deep-BSDE core, amortized only at raw Y0."""
    driver = driver_for(benchmark)
    full_seed = stable_seed("high_budget", "deep_bsde_original", benchmark, seed)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = OriginalAmortizedDeepBSDE(driver.d, time_steps=time_steps).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, eps=1e-8)
    valid_x = torch.randn((validation_size, driver.d), device=device, generator=valid_generator)
    valid_normals = torch.randn(
        (validation_size, time_steps, driver.d), device=device, generator=valid_generator
    )
    maximum = max(budgets)
    next_budget = step = nonfinite_count = 0
    history: list[dict[str, Any]] = []
    checkpoints: list[dict[str, Any]] = []
    last_validation = float("inf")
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    while next_budget < len(budgets):
        step += 1
        model.train()
        x0 = torch.randn((batch_size, driver.d), device=device, generator=generator)
        normals = torch.randn(
            (batch_size, time_steps, driver.d), device=device, generator=generator
        )
        terminal_y, terminal_x, rollout_nonfinite = model.rollout(driver, x0, normals)
        loss = _clipped_square(terminal_y - driver.terminal(terminal_x))
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        norm = _gradient_norm(model.parameters())
        finite = bool(torch.isfinite(loss).item() and math.isfinite(norm))
        if finite:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 50.0)
            optimizer.step()
        else:
            optimizer.zero_grad(set_to_none=True)
            nonfinite_count += 1
        nonfinite_count += rollout_nonfinite
        elapsed = _elapsed(started, device)
        lr = _set_timed_lr(optimizer, learning_rate, elapsed)
        validate = step == 1 or step % validation_every == 0
        if validate:
            last_validation, validation_nonfinite = deep_bsde_validation(
                model, driver, valid_x, valid_normals
            )
            nonfinite_count += validation_nonfinite
            elapsed = _elapsed(started, device)
            history.append(
                {
                    "step": step, "elapsed_seconds": elapsed,
                    "training_terminal_clipped_mse": float(loss.detach().item()),
                    "validation_terminal_clipped_mse": last_validation,
                    "gradient_norm": norm, "learning_rate": lr,
                    "nonfinite_count_cumulative": nonfinite_count,
                }
            )
        while next_budget < len(budgets) and elapsed >= budgets[next_budget] and validate:
            target = budgets[next_budget]
            path = CHECKPOINT_DIR / f"deep_bsde_original_{benchmark}_s{seed}_b{int(target)}.pt"
            atomic_torch_save(
                path,
                {
                    "state_dict": _cpu_state(model), "architecture": architecture_payload(model),
                    "benchmark": benchmark, "seed": seed, "step": step,
                },
            )
            actual = _elapsed(started, device)
            checkpoints.append(
                {
                    "method": "deep_bsde_original", "benchmark": benchmark, "seed": seed,
                    "target_budget_seconds": target, "training_seconds": actual,
                    "step": step, "validation_objective": last_validation,
                    "checkpoint": str(path), "learning_rate": lr,
                    "time_steps": time_steps, "batch_size": batch_size,
                    "parameter_count": parameter_count(model),
                    "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
                    if device.type == "cuda" else None,
                    "nonfinite_count": nonfinite_count,
                    "selection_criterion": "predeclared elapsed-time checkpoint",
                    "final_test_used_for_selection": False,
                }
            )
            next_budget += 1
            elapsed = actual
        if elapsed > maximum * 1.35 and next_budget < len(budgets):
            raise RuntimeError("Deep-BSDE checkpoint timing missed the budget ladder")
    result = {
        "method": "deep_bsde_original", "benchmark": benchmark, "seed": seed,
        "configuration": {
            "adaptation": "raw amortized Y0(x0); official-style nonshared BN/ReLU Z_n",
            "terminal_aware_features": False, "time_steps": time_steps,
            "batch_size": batch_size, "learning_rate": learning_rate,
            "z_width": model.width, "validation_size": validation_size,
            "official_source_commit": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
        },
        "history": history, "checkpoints": checkpoints,
    }
    atomic_json(RAW_DIR / f"deep_bsde_original_{benchmark}_s{seed}.json", result)
    return result


def _pinn_base(
    model: RawTanhPINN,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    value = model.value(driver, t, x)
    u_t = torch.autograd.grad(value.sum(), t, create_graph=True, retain_graph=True)[0]
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True, retain_graph=True)[0]
    return u_t + driver.driver(value, driver.sigma * gradient), gradient, value


def _sampled_laplacian(
    gradient: torch.Tensor,
    x: torch.Tensor,
    coordinates: torch.Tensor,
    d: int,
    *,
    create_graph: bool,
) -> torch.Tensor:
    result = torch.zeros(len(x), dtype=x.dtype, device=x.device)
    coordinate_list = coordinates.tolist()
    for index, coordinate in enumerate(coordinate_list):
        second = torch.autograd.grad(
            gradient[:, coordinate].sum(), x, create_graph=create_graph,
            retain_graph=create_graph or index + 1 < len(coordinate_list),
        )[0][:, coordinate]
        result = result + second
    return (d / len(coordinate_list)) * result


def _terminal_penalty(
    model: RawTanhPINN,
    driver: ConcaveSemilinearDriver,
    terminal_x: torch.Tensor,
    *,
    create_graph: bool,
) -> tuple[torch.Tensor, dict[str, float]]:
    local_x = terminal_x.detach().requires_grad_(True)
    t = torch.full((len(local_x),), driver.T, dtype=local_x.dtype, device=local_x.device)
    value = model.value(driver, t, local_x)
    gradient = torch.autograd.grad(
        value.sum(), local_x, create_graph=create_graph, retain_graph=create_graph
    )[0]
    value_mse = torch.mean((value - driver.terminal(local_x)) ** 2)
    z_mse = torch.mean(
        (driver.sigma * gradient - driver.sigma * driver.terminal_gradient_torch(local_x)) ** 2
    )
    penalty = 20.0 * value_mse + z_mse
    return penalty, {
        "terminal_value_rmse": float(torch.sqrt(value_mse.detach()).item()),
        "terminal_gradient_component_rmse": float(torch.sqrt(z_mse.detach()).item()),
    }


def sdgd_objective(
    model: RawTanhPINN,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    terminal_x: torch.Tensor,
    *,
    coordinate_budget: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, dict[str, float]]:
    local_t = t.detach().requires_grad_(True)
    local_x = x.detach().requires_grad_(True)
    base, gradient, _ = _pinn_base(model, driver, local_t, local_x)
    budget = min(coordinate_budget, driver.d)
    first = torch.randperm(driver.d, device=x.device, generator=generator)[:budget]
    second = torch.randperm(driver.d, device=x.device, generator=generator)[:budget]
    residual_first = base + _sampled_laplacian(
        gradient, local_x, first, driver.d, create_graph=True
    )
    residual_second = base + _sampled_laplacian(
        gradient, local_x, second, driver.d, create_graph=True
    )
    pde_objective = torch.mean(
        residual_first.detach() * residual_second
        + residual_second.detach() * residual_first
    )
    terminal, terminal_metrics = _terminal_penalty(
        model, driver, terminal_x, create_graph=True
    )
    objective = pde_objective + terminal
    return objective, {
        "pde_unbiased_product": float(
            torch.mean(residual_first.detach() * residual_second.detach()).item()
        ),
        "sampled_residual_rms": float(
            torch.sqrt(0.5 * torch.mean(residual_first.detach() ** 2 + residual_second.detach() ** 2)).item()
        ),
        "terminal_penalty": float(terminal.detach().item()),
        **terminal_metrics,
    }


def pinn_validation_proxy(
    model: RawTanhPINN,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    terminal_x: torch.Tensor,
    coordinate_pairs: list[tuple[torch.Tensor, torch.Tensor]],
) -> tuple[float, dict[str, float]]:
    model.eval()
    products: list[float] = []
    rms_values: list[float] = []
    for first, second in coordinate_pairs:
        residuals: list[torch.Tensor] = []
        for coordinates in (first, second):
            local_t = t.detach().clone().requires_grad_(True)
            local_x = x.detach().clone().requires_grad_(True)
            base, gradient, _ = _pinn_base(model, driver, local_t, local_x)
            residuals.append(
                base + _sampled_laplacian(
                    gradient, local_x, coordinates, driver.d, create_graph=False
                )
            )
        products.append(float(torch.mean(residuals[0].detach() * residuals[1].detach()).item()))
        rms_values.append(float(torch.sqrt(0.5 * torch.mean(residuals[0].detach() ** 2 + residuals[1].detach() ** 2)).item()))
    terminal, terminal_metrics = _terminal_penalty(
        model, driver, terminal_x, create_graph=False
    )
    residual_rms = float(np.mean(rms_values))
    objective = residual_rms + float(torch.sqrt(torch.clamp(terminal.detach(), min=0.0)).item())
    return objective, {
        "sampled_residual_rms": residual_rms,
        "unbiased_residual_product": float(np.mean(products)),
        "terminal_penalty": float(terminal.detach().item()),
        **terminal_metrics,
    }


def train_sdgd_pinn_original(
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    budgets: tuple[float, ...] = BUDGETS,
    width: int = 192,
    depth: int = 3,
    batch_size: int = 128,
    coordinate_budget: int = 16,
    learning_rate: float = 5e-4,
    validation_points: int = 64,
    validation_replicates: int = 4,
    validation_every: int = 25,
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    full_seed = stable_seed("high_budget", "sdgd_pinn_original", benchmark, seed)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = RawTanhPINN(driver.d, width=width, depth=depth).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, eps=1e-8)
    valid_t = torch.rand(validation_points, device=device, generator=valid_generator) * driver.T
    valid_x = torch.randn(
        (validation_points, driver.d), device=device, generator=valid_generator
    ) * torch.sqrt(1.0 + driver.sigma**2 * valid_t[:, None])
    valid_terminal_x = torch.randn(
        (validation_points, driver.d), device=device, generator=valid_generator
    ) * math.sqrt(1.0 + driver.sigma**2 * driver.T)
    budget = min(coordinate_budget, driver.d)
    coordinate_pairs = [
        (
            torch.randperm(driver.d, device=device, generator=valid_generator)[:budget],
            torch.randperm(driver.d, device=device, generator=valid_generator)[:budget],
        )
        for _ in range(validation_replicates)
    ]
    next_budget = step = nonfinite_count = 0
    maximum = max(budgets)
    history: list[dict[str, Any]] = []
    checkpoints: list[dict[str, Any]] = []
    last_validation = float("inf")
    last_metrics: dict[str, float] = {}
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    while next_budget < len(budgets):
        step += 1
        model.train()
        t = torch.rand(batch_size, device=device, generator=generator) * driver.T
        x = torch.randn((batch_size, driver.d), device=device, generator=generator)
        x = x * torch.sqrt(1.0 + driver.sigma**2 * t[:, None])
        terminal_x = torch.randn((batch_size, driver.d), device=device, generator=generator)
        terminal_x = terminal_x * math.sqrt(1.0 + driver.sigma**2 * driver.T)
        objective, diagnostics = sdgd_objective(
            model, driver, t, x, terminal_x,
            coordinate_budget=coordinate_budget, generator=generator,
        )
        optimizer.zero_grad(set_to_none=True)
        objective.backward()
        norm = _gradient_norm(model.parameters())
        finite = bool(torch.isfinite(objective).item() and math.isfinite(norm))
        if finite:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            optimizer.step()
        else:
            optimizer.zero_grad(set_to_none=True)
            nonfinite_count += 1
        elapsed = _elapsed(started, device)
        lr = _set_timed_lr(optimizer, learning_rate, elapsed)
        validate = step == 1 or step % validation_every == 0
        if validate:
            last_validation, last_metrics = pinn_validation_proxy(
                model, driver, valid_t, valid_x, valid_terminal_x, coordinate_pairs
            )
            elapsed = _elapsed(started, device)
            history.append(
                {
                    "step": step, "elapsed_seconds": elapsed,
                    "training_objective": float(objective.detach().item()),
                    **{f"training_{key}": value for key, value in diagnostics.items()},
                    **{f"validation_{key}": value for key, value in last_metrics.items()},
                    "validation_objective": last_validation,
                    "gradient_norm": norm, "learning_rate": lr,
                    "nonfinite_count_cumulative": nonfinite_count,
                }
            )
        while next_budget < len(budgets) and elapsed >= budgets[next_budget] and validate:
            target = budgets[next_budget]
            path = CHECKPOINT_DIR / f"sdgd_pinn_original_{benchmark}_s{seed}_b{int(target)}.pt"
            atomic_torch_save(
                path,
                {
                    "state_dict": _cpu_state(model), "architecture": architecture_payload(model),
                    "benchmark": benchmark, "seed": seed, "step": step,
                },
            )
            actual = _elapsed(started, device)
            checkpoints.append(
                {
                    "method": "sdgd_pinn_original", "benchmark": benchmark, "seed": seed,
                    "target_budget_seconds": target, "training_seconds": actual,
                    "step": step, "validation_objective": last_validation,
                    **{f"validation_{key}": value for key, value in last_metrics.items()},
                    "checkpoint": str(path), "learning_rate": lr,
                    "coordinate_budget": coordinate_budget, "batch_size": batch_size,
                    "parameter_count": parameter_count(model),
                    "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
                    if device.type == "cuda" else None,
                    "nonfinite_count": nonfinite_count,
                    "selection_criterion": "predeclared elapsed-time checkpoint",
                    "final_test_used_for_selection": False,
                }
            )
            next_budget += 1
            elapsed = actual
        if elapsed > maximum * 1.35 and next_budget < len(budgets):
            raise RuntimeError("SDGD checkpoint timing missed the budget ladder")
    result = {
        "method": "sdgd_pinn_original", "benchmark": benchmark, "seed": seed,
        "configuration": {
            "network": "ordinary tanh MLP (t,x)->u", "terminal_aware_features": False,
            "terminal_training": "official soft penalty 20*MSE(u_T,g)+MSE(z_T,z_g)",
            "width": width, "depth": depth, "batch_size": batch_size,
            "coordinate_budget": coordinate_budget, "learning_rate": learning_rate,
            "official_source_commit": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
        },
        "history": history, "checkpoints": checkpoints,
    }
    atomic_json(RAW_DIR / f"sdgd_pinn_original_{benchmark}_s{seed}.json", result)
    return result


def load_deep_bsde_original(
    path: Path, device: torch.device
) -> tuple[OriginalAmortizedDeepBSDE, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = dict(payload["architecture"])
    architecture.pop("kind")
    model = OriginalAmortizedDeepBSDE(**architecture).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


def load_sdgd_pinn_original(
    path: Path, device: torch.device
) -> tuple[RawTanhPINN, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = dict(payload["architecture"])
    architecture.pop("kind")
    model = RawTanhPINN(**architecture).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


def predict_deep_bsde_original(
    model: OriginalAmortizedDeepBSDE,
    driver: ConcaveSemilinearDriver,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    state = torch.as_tensor(x, dtype=torch.float32, device=device)
    synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        value = model.initial_value(driver, state)
    elapsed = _elapsed(started, device)
    return value.cpu().numpy().astype(np.float64), elapsed


def predict_sdgd_pinn_original(
    model: RawTanhPINN,
    driver: ConcaveSemilinearDriver,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    state = torch.as_tensor(x, dtype=torch.float32, device=device)
    t = torch.zeros(len(state), dtype=state.dtype, device=device)
    synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        value = model.value(driver, t, state)
    elapsed = _elapsed(started, device)
    return value.cpu().numpy().astype(np.float64), elapsed


__all__ = [
    "load_deep_bsde_original", "load_sdgd_pinn_original",
    "pinn_validation_proxy", "predict_deep_bsde_original",
    "predict_sdgd_pinn_original", "sdgd_objective", "train_deep_bsde_original",
    "train_sdgd_pinn_original",
]
