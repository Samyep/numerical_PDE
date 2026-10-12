"""PINN baseline with published Stochastic Dimension Gradient Descent.

Training uses two independent coordinate subsets in each step.  Their product
gives an unbiased stochastic gradient of the squared full PDE residual,
avoiding the extra variance term from naively squaring one Laplacian estimate.
The terminal condition is hard-enforced and its residual is still measured.
"""

from __future__ import annotations

import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .common import load_config, seed_everything, stable_seed
from .equation import B1Equation
from .neural_common import DenseNet, clone_state_dict, count_parameters, tx_input


def _terminal_torch(equation: B1Equation, x: torch.Tensor) -> torch.Tensor:
    A = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
    return -torch.logsumexp(x @ A.T - math.log(equation.K), dim=-1)


class HardTerminalPINN(nn.Module):
    def __init__(self, d: int, *, width: int, depth: int, activation: str) -> None:
        super().__init__()
        self.d = d
        self.correction = DenseNet(
            d + 1, 1, width=width, depth=depth, activation=activation
        )

    def value(self, equation: B1Equation, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        correction = self.correction(tx_input(t, x, equation.T)).reshape(-1)
        return _terminal_torch(equation, x) + (equation.T - t) * correction


def _base_and_gradient(
    model: HardTerminalPINN,
    equation: B1Equation,
    t: torch.Tensor,
    x: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    value = model.value(equation, t, x)
    u_t = torch.autograd.grad(value.sum(), t, create_graph=True, retain_graph=True)[0]
    gradient = torch.autograd.grad(
        value.sum(), x, create_graph=True, retain_graph=True
    )[0]
    base = u_t - torch.sum(gradient.square(), dim=-1)
    return base, gradient, value


def _sampled_laplacian(
    gradient: torch.Tensor,
    x: torch.Tensor,
    coordinates: torch.Tensor,
    d: int,
) -> torch.Tensor:
    diagonal = torch.zeros(len(x), dtype=x.dtype, device=x.device)
    for coordinate in coordinates.tolist():
        second = torch.autograd.grad(
            gradient[:, coordinate].sum(),
            x,
            create_graph=True,
            retain_graph=True,
        )[0][:, coordinate]
        diagonal = diagonal + second
    return (d / len(coordinates)) * diagonal


def sdgd_loss(
    model: HardTerminalPINN,
    equation: B1Equation,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    coordinate_budget: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, dict[str, float]]:
    t = t.detach().requires_grad_(True)
    x = x.detach().requires_grad_(True)
    base, gradient, _ = _base_and_gradient(model, equation, t, x)
    budget = min(coordinate_budget, equation.d)
    first = torch.randperm(equation.d, device=x.device, generator=generator)[:budget]
    second = torch.randperm(equation.d, device=x.device, generator=generator)[:budget]
    residual_first = base + _sampled_laplacian(gradient, x, first, equation.d)
    residual_second = base + _sampled_laplacian(gradient, x, second, equation.d)
    # Independent subsets make E[2 R_1 stopgrad(R_2)] the gradient of E[R]^2.
    physics = torch.mean(2.0 * residual_first.detach() * residual_second)
    terminal_x = x.detach()
    terminal_t = torch.full_like(t, equation.T)
    terminal_residual = model.value(equation, terminal_t, terminal_x) - _terminal_torch(
        equation, terminal_x
    )
    terminal_loss = torch.mean(terminal_residual.square())
    return physics + terminal_loss, {
        "physics_proxy": float(torch.mean(residual_first.detach() * residual_second.detach()).item()),
        "terminal_loss": float(terminal_loss.detach().item()),
    }


def _validation_proxy(
    model: HardTerminalPINN,
    equation: B1Equation,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    coordinate_budget: int,
    generator: torch.Generator,
) -> float:
    model.eval()
    estimates: list[float] = []
    # Averaging four independent products reduces selection noise without an
    # O(d) full Laplacian validation pass.
    for _ in range(4):
        loss, diagnostics = sdgd_loss(
            model,
            equation,
            t,
            x,
            coordinate_budget=coordinate_budget,
            generator=generator,
        )
        estimates.append(diagnostics["physics_proxy"] + diagnostics["terminal_loss"])
    return float(np.mean(estimates))


def sdgd_checkpoint_path(root: Path, equation: B1Equation, seed: int, stage: str) -> Path:
    return root / f"sdgd_pinn_{stage}_d{equation.d}_k{equation.k}_s{seed}.pt"


def train_sdgd_pinn(
    equation: B1Equation,
    *,
    seed: int,
    stage: str,
    device: torch.device,
    checkpoint_path: Path,
    force: bool = False,
) -> tuple[HardTerminalPINN, dict[str, Any]]:
    config = load_config()["sdgd_pinn"]
    if checkpoint_path.is_file() and not force:
        payload = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model = HardTerminalPINN(
            equation.d,
            width=config["width"],
            depth=config["depth"],
            activation=config["activation"],
        ).to(device)
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return model, payload["metadata"]

    full_seed = stable_seed("B1", "sdgd_pinn", equation.d, equation.k, seed, stage)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = HardTerminalPINN(
        equation.d,
        width=config["width"],
        depth=config["depth"],
        activation=config["activation"],
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
    steps = int(config[f"{stage}_steps"])
    valid_t = torch.rand(
        config["validation_points"], device=device, generator=validation_generator
    ) * (equation.T - 1e-6)
    valid_x = torch.randn(
        (config["validation_points"], equation.d),
        device=device,
        generator=validation_generator,
    )
    best_loss = float("inf")
    best_step = 0
    best_state: dict[str, torch.Tensor] | None = None
    patience = 0
    history: list[dict[str, float | int]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    executed = 0
    for step in range(1, steps + 1):
        executed = step
        model.train()
        t = torch.rand(config["batch_size"], device=device, generator=generator) * (
            equation.T - 1e-6
        )
        x = torch.randn(
            (config["batch_size"], equation.d), device=device, generator=generator
        )
        loss, diagnostics = sdgd_loss(
            model,
            equation,
            t,
            x,
            coordinate_budget=config["coordinate_budget"],
            generator=generator,
        )
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), config["gradient_clip"])
        optimizer.step()
        if step % config["validation_every"] == 0 or step == steps:
            validation = _validation_proxy(
                model,
                equation,
                valid_t,
                valid_x,
                coordinate_budget=config["coordinate_budget"],
                generator=validation_generator,
            )
            history.append(
                {
                    "step": step,
                    "training_objective": float(loss.detach().item()),
                    "training_physics_proxy": diagnostics["physics_proxy"],
                    "validation_residual_product": validation,
                }
            )
            # The unbiased product can be slightly negative; select by its
            # absolute magnitude as a fixed, reference-free validation proxy.
            selection_loss = abs(validation)
            if selection_loss < best_loss * (1.0 - 1e-5):
                best_loss = selection_loss
                best_step = step
                best_state = clone_state_dict(model)
                patience = 0
            else:
                patience += 1
            if patience >= config["early_stopping_patience"]:
                break
    if best_state is None:
        raise RuntimeError("SDGD-PINN produced no checkpoint")
    model.load_state_dict(best_state)
    model.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    elapsed = time.perf_counter() - started
    metadata: dict[str, Any] = {
        "method": "pinn_sdgd",
        "d": equation.d,
        "k": equation.k,
        "seed": seed,
        "full_seed": full_seed,
        "stage": stage,
        "architecture": f"hard-terminal PINN {config['depth']}x{config['width']} {config['activation']}",
        "parameter_count": count_parameters(model),
        "optimizer_steps": executed,
        "best_step": best_step,
        "coordinate_budget": min(config["coordinate_budget"], equation.d),
        "coordinate_sampling": "uniform without replacement; two independent subsets",
        "laplacian_scale": equation.d / min(config["coordinate_budget"], equation.d),
        "training_objective": "independent-product unbiased gradient of squared full residual",
        "terminal_condition": "hard-enforced; residual also measured",
        "best_validation_residual_product_abs": best_loss,
        "training_seconds": elapsed,
        "peak_gpu_bytes": peak,
        "history": history,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
        "source_commit": load_config()["external_sources"]["sdgd_pinn"]["commit"],
    }
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": best_state, "metadata": metadata}, checkpoint_path)
    return model, metadata


def predict_sdgd_pinn(
    model: HardTerminalPINN,
    equation: B1Equation,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, float]:
    started = time.perf_counter()
    state = torch.as_tensor(x, dtype=torch.float32, device=device).requires_grad_(True)
    t = torch.zeros(len(state), dtype=torch.float32, device=device)
    value = model.value(equation, t, state)
    gradient = torch.autograd.grad(value.sum(), state, create_graph=False)[0]
    z = math.sqrt(2.0) * gradient
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    return value.detach().cpu().numpy(), z.detach().cpu().numpy(), time.perf_counter() - started
