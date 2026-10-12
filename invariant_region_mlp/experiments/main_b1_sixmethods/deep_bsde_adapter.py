"""Global/amortized Deep BSDE baseline for the common initial-state task."""

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


class GlobalDeepBSDE(nn.Module):
    """Amortized Y_0(x_0) plus a shared Z(t,x) network."""

    def __init__(
        self,
        d: int,
        *,
        width: int,
        depth: int,
        activation: str,
        z_bound: float,
    ) -> None:
        super().__init__()
        self.d = d
        self.u0 = DenseNet(d, 1, width=width, depth=depth, activation=activation)
        self.z = DenseNet(
            d + 1,
            d,
            width=width,
            depth=depth,
            activation=activation,
            output_bound=z_bound,
        )

    def initial_value(self, x: torch.Tensor) -> torch.Tensor:
        return self.u0(x).reshape(-1)

    def z_value(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        return self.z(tx_input(t, x, T))

    def rollout(
        self,
        equation: B1Equation,
        x0: torch.Tensor,
        normals: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        steps = normals.shape[1]
        dt = equation.T / steps
        sqrt_dt = math.sqrt(dt)
        x = x0
        y = self.initial_value(x0)
        for step in range(steps):
            t = torch.full(
                (len(x0),), step * dt, dtype=x0.dtype, device=x0.device
            )
            z = self.z_value(t, x, equation.T)
            dw = sqrt_dt * normals[:, step, :]
            # dY = -f(Z)dt + Z dW and f=-|Z|^2/2.
            y = y + 0.5 * torch.sum(z.square(), dim=-1) * dt + torch.sum(z * dw, dim=-1)
            x = x + equation.sigma * dw
        return y, x


def deep_bsde_checkpoint_path(
    root: Path, equation: B1Equation, seed: int, stage: str
) -> Path:
    return root / f"deep_bsde_{stage}_d{equation.d}_k{equation.k}_s{seed}.pt"


def _validation_loss(
    model: GlobalDeepBSDE,
    equation: B1Equation,
    x0: torch.Tensor,
    normals: torch.Tensor,
) -> float:
    model.eval()
    with torch.no_grad():
        terminal_y, terminal_x = model.rollout(equation, x0, normals)
        difference = terminal_y - _terminal_torch(equation, terminal_x)
        return float(torch.mean(torch.nn.functional.huber_loss(
            terminal_y,
            _terminal_torch(equation, terminal_x),
            reduction="none",
            delta=10.0,
        )).item())


def train_deep_bsde(
    equation: B1Equation,
    *,
    seed: int,
    stage: str,
    device: torch.device,
    checkpoint_path: Path,
    force: bool = False,
) -> tuple[GlobalDeepBSDE, dict[str, Any]]:
    config = load_config()["deep_bsde"]
    if checkpoint_path.is_file() and not force:
        payload = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model = GlobalDeepBSDE(
            equation.d,
            width=config["width"],
            depth=config["depth"],
            activation=config["activation"],
            z_bound=config["z_output_bound"],
        ).to(device)
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return model, payload["metadata"]

    full_seed = stable_seed("B1", "deep_bsde", equation.d, equation.k, seed, stage)
    seed_everything(full_seed)
    generator = torch.Generator(device=device).manual_seed(full_seed)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    model = GlobalDeepBSDE(
        equation.d,
        width=config["width"],
        depth=config["depth"],
        activation=config["activation"],
        z_bound=config["z_output_bound"],
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
    steps = int(config[f"{stage}_steps"])
    valid_x = torch.randn(
        (config["validation_size"], equation.d), device=device, generator=valid_generator
    )
    valid_normals = torch.randn(
        (config["validation_size"], config["time_steps"], equation.d),
        device=device,
        generator=valid_generator,
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
        x0 = torch.randn(
            (config["batch_size"], equation.d), device=device, generator=generator
        )
        normals = torch.randn(
            (config["batch_size"], config["time_steps"], equation.d),
            device=device,
            generator=generator,
        )
        terminal_y, terminal_x = model.rollout(equation, x0, normals)
        target = _terminal_torch(equation, terminal_x)
        loss = torch.nn.functional.huber_loss(terminal_y, target, delta=10.0)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), config["gradient_clip"])
        optimizer.step()
        if step % config["validation_every"] == 0 or step == steps:
            valid_loss = _validation_loss(model, equation, valid_x, valid_normals)
            history.append(
                {"step": step, "training_terminal_loss": loss.item(), "validation_terminal_loss": valid_loss}
            )
            if valid_loss < best_loss * (1.0 - 1e-5):
                best_loss = valid_loss
                best_step = step
                best_state = clone_state_dict(model)
                patience = 0
            else:
                patience += 1
            if patience >= config["early_stopping_patience"]:
                break
    if best_state is None:
        raise RuntimeError("Deep BSDE produced no checkpoint")
    model.load_state_dict(best_state)
    model.eval()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    elapsed = time.perf_counter() - started
    metadata: dict[str, Any] = {
        "method": "deep_bsde",
        "protocol": "global/amortized over x0~N(0,I_d)",
        "d": equation.d,
        "k": equation.k,
        "seed": seed,
        "full_seed": full_seed,
        "stage": stage,
        "architecture": f"u0(x) and shared z(t,x), each {config['depth']}x{config['width']} {config['activation']}",
        "parameter_count": count_parameters(model),
        "z_protocol": "learned independently; not constrained to sqrt(2)*grad(u0)",
        "time_steps": config["time_steps"],
        "optimizer_steps": executed,
        "best_step": best_step,
        "best_validation_terminal_loss": best_loss,
        "training_seconds": elapsed,
        "peak_gpu_bytes": peak,
        "history": history,
        "bsde_sign": "dY=+0.5*|Z|^2 dt + Z dW",
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
        "source_commit": load_config()["external_sources"]["deep_bsde"]["commit"],
    }
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": best_state, "metadata": metadata}, checkpoint_path)
    return model, metadata


def predict_deep_bsde(
    model: GlobalDeepBSDE,
    equation: B1Equation,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, float]:
    started = time.perf_counter()
    state = torch.as_tensor(x, dtype=torch.float32, device=device)
    t = torch.zeros(len(state), dtype=torch.float32, device=device)
    with torch.no_grad():
        value = model.initial_value(state)
        z = model.z_value(t, state, equation.T)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    return value.cpu().numpy(), z.cpu().numpy(), time.perf_counter() - started
