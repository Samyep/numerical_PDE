"""Deep Picard Iteration adapted to the B1 Brownian quadratic HJB.

This preserves the published method's Picard rounds, Feynman--Kac value
labels, BEL gradient labels with control variates, joint value/gradient
regression, and warm starts.  Exact B1 references never enter training or
checkpoint selection.
"""

from __future__ import annotations

import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from .common import load_config, seed_everything, stable_seed
from .equation import B1Equation
from .neural_common import DenseNet, clone_state_dict, count_parameters, tx_input, value_and_z


class DPINet(DenseNet):
    def __init__(self, d: int, width: int, depth: int, activation: str) -> None:
        super().__init__(d + 1, 1, width=width, depth=depth, activation=activation)
        self.d = d

    def value(self, t: torch.Tensor, x: torch.Tensor, T: float = 0.25) -> torch.Tensor:
        return self(tx_input(t, x, T)).reshape(-1)


def _terminal_torch(equation: B1Equation, x: torch.Tensor) -> torch.Tensor:
    A = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
    return -torch.logsumexp(x @ A.T - math.log(equation.K), dim=-1)


def _previous_state(
    previous: DPINet | None,
    t: torch.Tensor,
    x: torch.Tensor,
    equation: B1Equation,
) -> tuple[torch.Tensor, torch.Tensor]:
    if previous is None:
        return (
            torch.zeros(len(t), dtype=x.dtype, device=x.device),
            torch.zeros_like(x),
        )
    with torch.enable_grad():
        local_x = x.detach().requires_grad_(True)
        value, z = value_and_z(
            previous, t.detach(), local_x, T=equation.T, create_graph=False
        )
    return value.detach(), z.detach()


def generate_picard_labels(
    equation: B1Equation,
    previous: DPINet | None,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    mc_samples: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Published one-layer value/z formula with finite-variance control variates."""

    batch, d = x.shape
    tau = torch.clamp(equation.T - t, min=1e-7)
    terminal_normal = torch.randn(
        (batch, mc_samples, d), dtype=x.dtype, device=x.device, generator=generator
    )
    terminal_x = x[:, None, :] + torch.sqrt(2.0 * tau)[:, None, None] * terminal_normal
    terminal_value = _terminal_torch(equation, terminal_x.reshape(-1, d)).reshape(
        batch, mc_samples
    )
    gx = _terminal_torch(equation, x)
    value_label = terminal_value.mean(dim=1)
    z_label = torch.mean(
        (terminal_value - gx[:, None])[:, :, None]
        * terminal_normal
        / torch.sqrt(tau)[:, None, None],
        dim=1,
    )

    if previous is not None:
        uniform = torch.rand(
            (batch, mc_samples), dtype=x.dtype, device=x.device, generator=generator
        )
        delta = tau[:, None] * torch.clamp(uniform, min=1e-6)
        integral_normal = torch.randn(
            (batch, mc_samples, d), dtype=x.dtype, device=x.device, generator=generator
        )
        sample_t = t[:, None] + delta
        sample_x = x[:, None, :] + torch.sqrt(2.0 * delta)[:, :, None] * integral_normal
        _, sample_z = _previous_state(
            previous, sample_t.reshape(-1), sample_x.reshape(-1, d), equation
        )
        sample_driver = -0.5 * torch.sum(sample_z.square(), dim=-1).reshape(
            batch, mc_samples
        )
        _, parent_z = _previous_state(previous, t, x, equation)
        parent_driver = -0.5 * torch.sum(parent_z.square(), dim=-1)
        value_label = value_label + tau * sample_driver.mean(dim=1)
        z_label = z_label + tau[:, None] * torch.mean(
            (sample_driver - parent_driver[:, None])[:, :, None]
            * integral_normal
            / torch.sqrt(delta)[:, :, None],
            dim=1,
        )
    return value_label.detach(), z_label.detach()


def dpi_checkpoint_path(root: Path, equation: B1Equation, seed: int, stage: str) -> Path:
    return root / f"dpi_{stage}_d{equation.d}_k{equation.k}_s{seed}.pt"


def train_dpi(
    equation: B1Equation,
    *,
    seed: int,
    stage: str,
    device: torch.device,
    checkpoint_path: Path,
    force: bool = False,
) -> tuple[DPINet, dict[str, Any]]:
    config = load_config()["dpi"]
    rounds = int(config[f"{stage}_picard_rounds"])
    steps_per_round = int(config[f"{stage}_steps_per_round"])
    if checkpoint_path.is_file() and not force:
        payload = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model = DPINet(
            equation.d, config["width"], config["depth"], config["activation"]
        ).to(device)
        model.load_state_dict(payload["state_dict"])
        model.eval()
        return model, payload["metadata"]

    full_seed = stable_seed("B1", "dpi", equation.d, equation.k, seed, stage)
    seed_everything(full_seed)
    data_generator = torch.Generator(device=device).manual_seed(full_seed)
    validation_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    previous: DPINet | None = None
    round_history: list[dict[str, Any]] = []
    total_steps = 0
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for picard_round in range(1, rounds + 1):
        model = DPINet(
            equation.d, config["width"], config["depth"], config["activation"]
        ).to(device)
        if previous is not None:
            model.load_state_dict(previous.state_dict())
        optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
        valid_t = torch.rand(
            config["validation_points"], device=device, generator=validation_generator
        ) * (equation.T - 1e-6)
        valid_x = torch.randn(
            (config["validation_points"], equation.d),
            device=device,
            generator=validation_generator,
        )
        valid_u, valid_z = generate_picard_labels(
            equation,
            previous,
            valid_t,
            valid_x,
            mc_samples=config["validation_mc_samples"],
            generator=validation_generator,
        )
        best_loss = float("inf")
        best_state: dict[str, torch.Tensor] | None = None
        best_step = 0
        patience = 0
        executed = 0
        for step in range(1, steps_per_round + 1):
            executed = step
            t = torch.rand(config["batch_size"], device=device, generator=data_generator) * (
                equation.T - 1e-6
            )
            x = torch.randn(
                (config["batch_size"], equation.d), device=device, generator=data_generator
            )
            target_u, target_z = generate_picard_labels(
                equation,
                previous,
                t,
                x,
                mc_samples=config["mc_samples"],
                generator=data_generator,
            )
            train_x = x.detach().requires_grad_(True)
            predicted_u, predicted_z = value_and_z(
                model, t, train_x, T=equation.T, create_graph=True
            )
            value_loss = torch.mean((predicted_u - target_u) ** 2)
            gradient_loss = torch.mean((predicted_z - target_z) ** 2)
            loss = value_loss + config["gradient_loss_weight"] * gradient_loss
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            total_steps += 1
            if step % config["validation_every"] == 0 or step == steps_per_round:
                valid_input = valid_x.detach().requires_grad_(True)
                predicted_valid_u, predicted_valid_z = value_and_z(
                    model, valid_t, valid_input, T=equation.T, create_graph=False
                )
                valid_loss = (
                    torch.mean((predicted_valid_u.detach() - valid_u) ** 2)
                    + config["gradient_loss_weight"]
                    * torch.mean((predicted_valid_z.detach() - valid_z) ** 2)
                ).item()
                if valid_loss < best_loss * (1.0 - 1e-5):
                    best_loss = valid_loss
                    best_state = clone_state_dict(model)
                    best_step = step
                    patience = 0
                else:
                    patience += 1
                if patience >= config["early_stopping_patience"]:
                    break
        if best_state is None:
            raise RuntimeError("DPI round produced no checkpoint")
        model.load_state_dict(best_state)
        model.eval()
        previous = model
        round_history.append(
            {
                "picard_round": picard_round,
                "steps_executed": executed,
                "best_step": best_step,
                "best_validation_label_loss": best_loss,
            }
        )
    if previous is None:
        raise RuntimeError("DPI ran zero Picard rounds")
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        peak = int(torch.cuda.max_memory_allocated(device))
    else:
        peak = None
    elapsed = time.perf_counter() - started
    metadata: dict[str, Any] = {
        "method": "dpi",
        "d": equation.d,
        "k": equation.k,
        "seed": seed,
        "full_seed": full_seed,
        "stage": stage,
        "architecture": f"4-layer-style value MLP {config['depth']}x{config['width']} {config['activation']}",
        "parameter_count": count_parameters(previous),
        "picard_rounds": rounds,
        "optimizer_steps": total_steps,
        "mc_samples_per_label": config["mc_samples"],
        "gradient_loss_weight": config["gradient_loss_weight"],
        "training_seconds": elapsed,
        "peak_gpu_bytes": peak,
        "round_history": round_history,
        "analytic_reference_used_for_training": False,
        "analytic_reference_used_for_checkpoint_selection": False,
        "source_commit": load_config()["external_sources"]["dpi"]["commit"],
    }
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"state_dict": clone_state_dict(previous), "metadata": metadata}, checkpoint_path)
    return previous, metadata


def predict_dpi(
    model: DPINet,
    equation: B1Equation,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, float]:
    started = time.perf_counter()
    state = torch.tensor(x, dtype=torch.float32, device=device).requires_grad_(True)
    t = torch.zeros(len(state), dtype=torch.float32, device=device)
    value, z = value_and_z(model, t, state, T=equation.T, create_graph=False)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    return value.detach().cpu().numpy(), z.detach().cpu().numpy(), time.perf_counter() - started
