"""Timed amortized Deep-BSDE and SDGD-PINN baselines for generic f(u,z)."""

from __future__ import annotations

import copy
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
from .models import AmortizedDeepBSDE, HardTerminalValue, architecture_payload


def _cpu_state(model: nn.Module) -> dict[str, torch.Tensor]:
    return {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}


def _gradient_norm(parameters: Any) -> float:
    square = 0.0
    for parameter in parameters:
        if parameter.grad is not None:
            square += float(parameter.grad.detach().double().square().sum().item())
    return math.sqrt(square)


@torch.no_grad()
def deep_bsde_validation(
    model: AmortizedDeepBSDE,
    driver: ConcaveSemilinearDriver,
    x: torch.Tensor,
    normals: torch.Tensor,
) -> tuple[float, int]:
    model.eval()
    terminal_y, terminal_x, nonfinite = model.rollout(driver, x, normals)
    target = driver.terminal(terminal_x)
    loss = torch.nn.functional.huber_loss(terminal_y, target, delta=1.0)
    return float(loss.item()), nonfinite


def _elapsed(started: float, device: torch.device) -> float:
    synchronize(device)
    return time.perf_counter() - started


def _set_timed_lr(optimizer: torch.optim.Optimizer, base_lr: float, seconds: float) -> float:
    factor = 1.0 if seconds < 480.0 else (0.2 if seconds < 800.0 else 0.04)
    value = base_lr * factor
    for group in optimizer.param_groups:
        group["lr"] = value
    return value


def train_deep_bsde_timed(
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    budgets: tuple[float, ...] = BUDGETS,
    width: int = 192,
    depth: int = 3,
    time_steps: int = 16,
    batch_size: int = 128,
    learning_rate: float = 5e-4,
    validation_size: int = 512,
    validation_every: int = 50,
    calibration_only_seconds: float | None = None,
    trial_tag: str = "",
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    full_seed = stable_seed("strong_baseline", "deep_bsde", benchmark, seed)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = AmortizedDeepBSDE(
        driver.d, width=width, depth=depth, activation="silu"
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, eps=1e-8)
    valid_x = torch.randn(
        (validation_size, driver.d), device=device, generator=validation_generator
    )
    valid_normals = torch.randn(
        (validation_size, time_steps, driver.d),
        device=device,
        generator=validation_generator,
    )
    effective_budgets = (
        (float(calibration_only_seconds),)
        if calibration_only_seconds is not None
        else budgets
    )
    maximum = max(effective_budgets)
    next_budget = 0
    step = 0
    history: list[dict[str, Any]] = []
    checkpoints: list[dict[str, Any]] = []
    nonfinite_count = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    last_validation = float("inf")
    while next_budget < len(effective_budgets):
        step += 1
        model.train()
        x0 = torch.randn((batch_size, driver.d), device=device, generator=generator)
        normals = torch.randn(
            (batch_size, time_steps, driver.d), device=device, generator=generator
        )
        terminal_y, terminal_x, rollout_nonfinite = model.rollout(driver, x0, normals)
        target = driver.terminal(terminal_x)
        loss = torch.nn.functional.huber_loss(terminal_y, target, delta=1.0)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        gradient_norm = _gradient_norm(model.parameters())
        finite = bool(torch.isfinite(loss).item() and math.isfinite(gradient_norm))
        if finite:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            optimizer.step()
        else:
            nonfinite_count += 1
            optimizer.zero_grad(set_to_none=True)
        nonfinite_count += rollout_nonfinite
        elapsed = _elapsed(started, device)
        lr = _set_timed_lr(optimizer, learning_rate, elapsed)
        should_validate = step == 1 or step % validation_every == 0
        if should_validate:
            last_validation, validation_nonfinite = deep_bsde_validation(
                model, driver, valid_x, valid_normals
            )
            nonfinite_count += validation_nonfinite
            elapsed = _elapsed(started, device)
            history.append(
                {
                    "step": step,
                    "elapsed_seconds": elapsed,
                    "training_terminal_huber": float(loss.detach().item()),
                    "validation_terminal_huber": last_validation,
                    "gradient_norm": gradient_norm,
                    "learning_rate": lr,
                    "nonfinite_count_cumulative": nonfinite_count,
                }
            )
        while (
            next_budget < len(effective_budgets)
            and elapsed >= effective_budgets[next_budget]
            and should_validate
        ):
            target_budget = effective_budgets[next_budget]
            state = _cpu_state(model)
            checkpoint_path = CHECKPOINT_DIR / (
                f"deep_bsde_{benchmark}_s{seed}{trial_tag}_b{int(target_budget)}.pt"
            )
            atomic_torch_save(
                checkpoint_path,
                {
                    "state_dict": state,
                    "architecture": architecture_payload(model),
                    "benchmark": benchmark,
                    "seed": seed,
                    "step": step,
                },
            )
            elapsed_after_save = _elapsed(started, device)
            record = {
                "method": "deep_bsde",
                "benchmark": benchmark,
                "seed": seed,
                "target_budget_seconds": target_budget,
                "training_seconds": elapsed_after_save,
                "step": step,
                "validation_objective": last_validation,
                "checkpoint": str(checkpoint_path),
                "learning_rate": lr,
                "time_steps": time_steps,
                "batch_size": batch_size,
                "parameter_count": parameter_count(model),
                "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
                if device.type == "cuda"
                else None,
                "nonfinite_count": nonfinite_count,
                "selection_criterion": "predeclared elapsed-time checkpoint; validation terminal Huber reported",
                "final_test_used_for_selection": False,
            }
            checkpoints.append(record)
            next_budget += 1
            elapsed = elapsed_after_save
        if elapsed > maximum * 1.35 and next_budget < len(effective_budgets):
            raise RuntimeError("Deep BSDE failed to emit a checkpoint near a time budget")
    result = {
        "method": "deep_bsde",
        "benchmark": benchmark,
        "seed": seed,
        "configuration": {
            "width": width,
            "depth": depth,
            "time_steps": time_steps,
            "batch_size": batch_size,
            "learning_rate": learning_rate,
            "validation_size": validation_size,
            "validation_every": validation_every,
            "adaptation": "amortized Y0(x0), shared terminal-aware Z(t,x)",
            "bsde_sign": "dY=-f(Y,Z)dt+Z dW",
        },
        "history": history,
        "checkpoints": checkpoints,
    }
    suffix = f"cal{int(calibration_only_seconds)}{trial_tag}" if calibration_only_seconds else "main"
    atomic_json(RAW_DIR / f"deep_bsde_{benchmark}_s{seed}_{suffix}.json", result)
    return result


def _pinn_base_gradient(
    model: HardTerminalValue,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    value = model.value(driver, t, x)
    u_t = torch.autograd.grad(value.sum(), t, create_graph=True, retain_graph=True)[0]
    gradient = torch.autograd.grad(
        value.sum(), x, create_graph=True, retain_graph=True
    )[0]
    z = driver.sigma * gradient
    base = u_t + driver.driver(value, z)
    return base, gradient, value


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
            gradient[:, coordinate].sum(),
            x,
            create_graph=create_graph,
            retain_graph=create_graph or index + 1 < len(coordinate_list),
        )[0][:, coordinate]
        result = result + second
    return (d / len(coordinate_list)) * result


def sdgd_objective(
    model: HardTerminalValue,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    coordinate_budget: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, dict[str, float]]:
    local_t = t.detach().requires_grad_(True)
    local_x = x.detach().requires_grad_(True)
    base, gradient, _ = _pinn_base_gradient(model, driver, local_t, local_x)
    budget = min(coordinate_budget, driver.d)
    first = torch.randperm(driver.d, device=x.device, generator=generator)[:budget]
    second = torch.randperm(driver.d, device=x.device, generator=generator)[:budget]
    residual_first = base + _sampled_laplacian(
        gradient, local_x, first, driver.d, create_graph=True
    )
    residual_second = base + _sampled_laplacian(
        gradient, local_x, second, driver.d, create_graph=True
    )
    objective = torch.mean(
        residual_first.detach() * residual_second
        + residual_second.detach() * residual_first
    )
    proxy = torch.mean(residual_first.detach() * residual_second.detach())
    return objective, {
        "residual_product": float(proxy.item()),
        "residual_first_rms": float(torch.sqrt(torch.mean(residual_first.detach() ** 2)).item()),
        "residual_second_rms": float(torch.sqrt(torch.mean(residual_second.detach() ** 2)).item()),
    }


def pinn_validation_proxy(
    model: HardTerminalValue,
    driver: ConcaveSemilinearDriver,
    t: torch.Tensor,
    x: torch.Tensor,
    coordinate_pairs: list[tuple[torch.Tensor, torch.Tensor]],
) -> tuple[float, float]:
    model.eval()
    estimates: list[float] = []
    residual_rms: list[float] = []
    for first, second in coordinate_pairs:
        local_t = t.detach().clone().requires_grad_(True)
        local_x = x.detach().clone().requires_grad_(True)
        base, gradient, _ = _pinn_base_gradient(model, driver, local_t, local_x)
        residual_first = base + _sampled_laplacian(
            gradient, local_x, first, driver.d, create_graph=False
        )
        # Rebuild the graph because create_graph=False releases it.
        local_t = t.detach().clone().requires_grad_(True)
        local_x = x.detach().clone().requires_grad_(True)
        base, gradient, _ = _pinn_base_gradient(model, driver, local_t, local_x)
        residual_second = base + _sampled_laplacian(
            gradient, local_x, second, driver.d, create_graph=False
        )
        estimates.append(float(torch.mean(residual_first.detach() * residual_second.detach()).item()))
        residual_rms.append(
            float(
                torch.sqrt(
                    0.5
                    * (
                        torch.mean(residual_first.detach() ** 2)
                        + torch.mean(residual_second.detach() ** 2)
                    )
                ).item()
            )
        )
    return float(np.mean(estimates)), float(np.mean(residual_rms))


def train_sdgd_pinn_timed(
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
    calibration_only_seconds: float | None = None,
    trial_tag: str = "",
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    full_seed = stable_seed("strong_baseline", "sdgd_pinn", benchmark, seed)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = HardTerminalValue(
        driver.d, width=width, depth=depth, activation="silu"
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, eps=1e-8)
    valid_t = torch.rand(
        validation_points, device=device, generator=validation_generator
    ) * (driver.T - 1e-5)
    valid_x = torch.randn(
        (validation_points, driver.d), device=device, generator=validation_generator
    )
    coordinate_pairs = []
    budget = min(coordinate_budget, driver.d)
    for _ in range(validation_replicates):
        coordinate_pairs.append(
            (
                torch.randperm(driver.d, device=device, generator=validation_generator)[:budget],
                torch.randperm(driver.d, device=device, generator=validation_generator)[:budget],
            )
        )
    effective_budgets = (
        (float(calibration_only_seconds),)
        if calibration_only_seconds is not None
        else budgets
    )
    maximum = max(effective_budgets)
    next_budget = step = nonfinite_count = 0
    history: list[dict[str, Any]] = []
    checkpoints: list[dict[str, Any]] = []
    last_validation = last_validation_rms = float("inf")
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    while next_budget < len(effective_budgets):
        step += 1
        model.train()
        t = torch.rand(batch_size, device=device, generator=generator) * (driver.T - 1e-5)
        x = torch.randn((batch_size, driver.d), device=device, generator=generator)
        objective, diagnostics = sdgd_objective(
            model,
            driver,
            t,
            x,
            coordinate_budget=coordinate_budget,
            generator=generator,
        )
        optimizer.zero_grad(set_to_none=True)
        objective.backward()
        gradient_norm = _gradient_norm(model.parameters())
        finite = bool(torch.isfinite(objective).item() and math.isfinite(gradient_norm))
        if finite:
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            optimizer.step()
        else:
            nonfinite_count += 1
            optimizer.zero_grad(set_to_none=True)
        elapsed = _elapsed(started, device)
        lr = _set_timed_lr(optimizer, learning_rate, elapsed)
        should_validate = step == 1 or step % validation_every == 0
        if should_validate:
            last_validation, last_validation_rms = pinn_validation_proxy(
                model, driver, valid_t, valid_x, coordinate_pairs
            )
            elapsed = _elapsed(started, device)
            history.append(
                {
                    "step": step,
                    "elapsed_seconds": elapsed,
                    "training_objective": float(objective.detach().item()),
                    "training_residual_product": diagnostics["residual_product"],
                    "validation_residual_product": last_validation,
                    "validation_sampled_residual_rms": last_validation_rms,
                    "gradient_norm": gradient_norm,
                    "learning_rate": lr,
                    "nonfinite_count_cumulative": nonfinite_count,
                }
            )
        while (
            next_budget < len(effective_budgets)
            and elapsed >= effective_budgets[next_budget]
            and should_validate
        ):
            target_budget = effective_budgets[next_budget]
            checkpoint_path = CHECKPOINT_DIR / (
                f"sdgd_pinn_{benchmark}_s{seed}{trial_tag}_b{int(target_budget)}.pt"
            )
            atomic_torch_save(
                checkpoint_path,
                {
                    "state_dict": _cpu_state(model),
                    "architecture": architecture_payload(model),
                    "benchmark": benchmark,
                    "seed": seed,
                    "step": step,
                },
            )
            elapsed_after_save = _elapsed(started, device)
            checkpoints.append(
                {
                    "method": "sdgd_pinn",
                    "benchmark": benchmark,
                    "seed": seed,
                    "target_budget_seconds": target_budget,
                    "training_seconds": elapsed_after_save,
                    "step": step,
                    # Use the nonnegative frozen sampled-residual RMS for
                    # checkpoint diagnostics.  Retain the unbiased two-draw
                    # product separately; it can be negative at finite sample
                    # size and is therefore a poor ranking scalar.
                    "validation_objective": last_validation_rms,
                    "validation_unbiased_residual_product": last_validation,
                    "validation_sampled_residual_rms": last_validation_rms,
                    "checkpoint": str(checkpoint_path),
                    "learning_rate": lr,
                    "coordinate_budget": coordinate_budget,
                    "batch_size": batch_size,
                    "parameter_count": parameter_count(model),
                    "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
                    if device.type == "cuda"
                    else None,
                    "nonfinite_count": nonfinite_count,
                    "selection_criterion": "predeclared elapsed-time checkpoint; frozen reference-free sampled-residual RMS reported",
                    "final_test_used_for_selection": False,
                }
            )
            next_budget += 1
            elapsed = elapsed_after_save
        if elapsed > maximum * 1.35 and next_budget < len(effective_budgets):
            raise RuntimeError("SDGD-PINN failed to emit a checkpoint near a time budget")
    result = {
        "method": "sdgd_pinn",
        "benchmark": benchmark,
        "seed": seed,
        "configuration": {
            "width": width,
            "depth": depth,
            "batch_size": batch_size,
            "coordinate_budget": coordinate_budget,
            "learning_rate": learning_rate,
            "validation_points": validation_points,
            "validation_replicates": validation_replicates,
            "validation_every": validation_every,
            "terminal_condition": "hard enforced",
            "objective": "two-independent-coordinate-set unbiased squared-residual gradient",
        },
        "history": history,
        "checkpoints": checkpoints,
    }
    suffix = f"cal{int(calibration_only_seconds)}{trial_tag}" if calibration_only_seconds else "main"
    atomic_json(RAW_DIR / f"sdgd_pinn_{benchmark}_s{seed}_{suffix}.json", result)
    return result


def load_deep_bsde(path: Path, device: torch.device) -> tuple[AmortizedDeepBSDE, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = AmortizedDeepBSDE(**architecture).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


def load_sdgd_pinn(path: Path, device: torch.device) -> tuple[HardTerminalValue, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = HardTerminalValue(**architecture).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


def predict_deep_bsde(
    model: AmortizedDeepBSDE,
    driver: ConcaveSemilinearDriver,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    state = torch.as_tensor(np.asarray(x), dtype=torch.float32, device=device)
    synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        value = model.initial_value(driver, state)
    elapsed = _elapsed(started, device)
    return value.cpu().numpy().astype(np.float64), elapsed


def predict_sdgd_pinn(
    model: HardTerminalValue,
    driver: ConcaveSemilinearDriver,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, float]:
    state = torch.as_tensor(np.asarray(x), dtype=torch.float32, device=device)
    t = torch.zeros(len(state), dtype=state.dtype, device=device)
    synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        value = model.value(driver, t, state)
    elapsed = _elapsed(started, device)
    return value.cpu().numpy().astype(np.float64), elapsed


__all__ = [
    "deep_bsde_validation", "load_deep_bsde", "load_sdgd_pinn",
    "pinn_validation_proxy", "predict_deep_bsde", "predict_sdgd_pinn",
    "sdgd_objective", "train_deep_bsde_timed", "train_sdgd_pinn_timed",
]
