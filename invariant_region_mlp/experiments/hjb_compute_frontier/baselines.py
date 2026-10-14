"""Frozen strong external baselines for the two HJB benchmarks.

The module intentionally exposes one configuration per baseline.  Exact
nonterminal references are absent from every training and checkpoint-selection
function.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .common import atomic_torch_save, seed_everything
from .equations import Equation, torch_terminal, torch_terminal_z
from .models import DenseNet, clone_state_dict, count_parameters, tx_input


def _peak(device: torch.device) -> int | None:
    if device.type != "cuda":
        return None
    torch.cuda.synchronize(device)
    return int(torch.cuda.max_memory_allocated(device))


class HardTerminalPINN(nn.Module):
    def __init__(self, d: int, *, width: int = 512, depth: int = 4) -> None:
        super().__init__()
        self.d = d
        self.width = width
        self.depth = depth
        self.correction = DenseNet(
            d + 1, 1, width=width, depth=depth, activation="tanh"
        )

    def value(self, equation: Equation, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return torch_terminal(equation, x) + (equation.T - t) * self.correction(
            tx_input(t, x, equation.T)
        ).reshape(-1)


def _pinn_base_gradient(
    model: HardTerminalPINN,
    equation: Equation,
    t: torch.Tensor,
    x: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    value = model.value(equation, t, x)
    u_t = torch.autograd.grad(value.sum(), t, create_graph=True, retain_graph=True)[0]
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True, retain_graph=True)[0]
    base = u_t - torch.sum(gradient.square(), dim=-1)
    return base, gradient, value


def _sampled_laplacian(
    gradient: torch.Tensor,
    x: torch.Tensor,
    coordinates: torch.Tensor,
    d: int,
    *,
    create_graph: bool,
) -> torch.Tensor:
    output = torch.zeros(len(x), dtype=x.dtype, device=x.device)
    coordinate_list = coordinates.tolist()
    for index, coordinate in enumerate(coordinate_list):
        diagonal = torch.autograd.grad(
            gradient[:, coordinate].sum(),
            x,
            create_graph=create_graph,
            retain_graph=create_graph or index + 1 < len(coordinate_list),
        )[0][:, coordinate]
        output = output + diagonal
    return (d / len(coordinate_list)) * output


def sdgd_objective(
    model: HardTerminalPINN,
    equation: Equation,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    coordinate_budget: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, dict[str, float]]:
    local_t = t.detach().requires_grad_(True)
    local_x = x.detach().requires_grad_(True)
    base, gradient, _ = _pinn_base_gradient(model, equation, local_t, local_x)
    budget = min(coordinate_budget, equation.d)
    first = torch.randperm(equation.d, device=x.device, generator=generator)[:budget]
    second = torch.randperm(equation.d, device=x.device, generator=generator)[:budget]
    residual_first = base + _sampled_laplacian(
        gradient, local_x, first, equation.d, create_graph=True
    )
    residual_second = base + _sampled_laplacian(
        gradient, local_x, second, equation.d, create_graph=True
    )
    # Two independent coordinate subsets give an unbiased gradient of the
    # squared full residual.  Symmetrization uses both stochastic graphs.
    objective = torch.mean(
        residual_first.detach() * residual_second
        + residual_second.detach() * residual_first
    )
    return objective, {
        "residual_product": float(
            torch.mean(residual_first.detach() * residual_second.detach()).item()
        )
    }


def full_pinn_residual_mse(
    model: HardTerminalPINN,
    equation: Equation,
    t: torch.Tensor,
    x: torch.Tensor,
) -> float:
    """Stable selection loss using the exact 400-coordinate Laplacian."""

    model.eval()
    local_t = t.detach().clone().requires_grad_(True)
    local_x = x.detach().clone().requires_grad_(True)
    base, gradient, _ = _pinn_base_gradient(model, equation, local_t, local_x)
    coordinates = torch.arange(equation.d, device=x.device)
    laplacian = _sampled_laplacian(
        gradient, local_x, coordinates, equation.d, create_graph=False
    )
    residual = base.detach() + laplacian.detach()
    return float(torch.mean(residual.square()).item())


def train_pinn_sdgd(
    equation: Equation,
    *,
    seed: int,
    device: torch.device,
    output_dir: Path,
    max_steps: int = 10_000,
    checkpoints: tuple[int, ...] = (2_000, 5_000, 10_000),
    width: int = 512,
    batch_size: int = 128,
    coordinate_budget: int = 32,
) -> dict[str, Any]:
    seed_everything(seed)
    model = HardTerminalPINN(equation.d, width=width, depth=4).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[5_000, 8_000], gamma=0.2
    )
    generator = torch.Generator(device=device).manual_seed(seed)
    valid_generator = torch.Generator(device=device).manual_seed(seed + 1)
    valid_t = torch.rand(8, device=device, generator=valid_generator) * (equation.T - 1e-5)
    valid_x = torch.randn((8, equation.d), device=device, generator=valid_generator)
    audit_t = torch.rand((4, 8), device=device, generator=valid_generator) * (
        equation.T - 1e-5
    )
    audit_x = torch.randn(
        (4, 8, equation.d), device=device, generator=valid_generator
    )
    history: list[dict[str, Any]] = []
    saved: list[dict[str, Any]] = []
    cumulative_started = time.perf_counter()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    for step in range(1, max_steps + 1):
        model.train()
        t = torch.rand(batch_size, device=device, generator=generator) * (equation.T - 1e-5)
        x = torch.randn((batch_size, equation.d), device=device, generator=generator)
        objective, diagnostics = sdgd_objective(
            model,
            equation,
            t,
            x,
            coordinate_budget=coordinate_budget,
            generator=generator,
        )
        optimizer.zero_grad(set_to_none=True)
        objective.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step % 250 == 0:
            history.append(
                {
                    "step": step,
                    "training_objective": float(objective.detach().item()),
                    "training_residual_product": diagnostics["residual_product"],
                    "learning_rate": float(optimizer.param_groups[0]["lr"]),
                    "elapsed_seconds": time.perf_counter() - cumulative_started,
                }
            )
        if step in checkpoints:
            validation = full_pinn_residual_mse(model, equation, valid_t, valid_x)
            validation_audit = [
                full_pinn_residual_mse(model, equation, audit_t[index], audit_x[index])
                for index in range(len(audit_t))
            ]
            elapsed = time.perf_counter() - cumulative_started
            metadata = {
                "method": "pinn_sdgd",
                "seed": seed,
                "step": step,
                "training_seconds": elapsed,
                "validation_full_residual_mse": validation,
                "validation_points": len(valid_t),
                "validation_laplacian_coordinates": equation.d,
                "validation_audit_residual_mse_banks": validation_audit,
                "validation_audit_residual_mse_mean": float(np.mean(validation_audit)),
                "validation_audit_residual_mse_std": float(
                    np.std(validation_audit, ddof=1)
                ),
                "validation_audit_used_for_selection": False,
                "coordinate_budget": coordinate_budget,
                "parameter_count": count_parameters(model),
                "peak_gpu_bytes": _peak(device),
                "analytic_reference_used_for_training_or_selection": False,
            }
            path = output_dir / f"pinn_sdgd_step{step}.pt"
            atomic_torch_save(
                path,
                {
                    "state_dict": clone_state_dict(model),
                    "architecture": {"d": equation.d, "width": width, "depth": 4},
                    "metadata": metadata,
                },
            )
            saved.append({"path": str(path), **metadata})
    return {
        "method": "pinn_sdgd",
        "seed": seed,
        "history": history,
        "checkpoints": saved,
        "selected_step": min(saved, key=lambda item: item["validation_full_residual_mse"])["step"],
        "selection_metric": "minimum frozen full-Laplacian squared PDE residual",
    }


def load_pinn(path: Path, device: torch.device) -> tuple[HardTerminalPINN, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = HardTerminalPINN(
        int(architecture["d"]),
        width=int(architecture["width"]),
        depth=int(architecture["depth"]),
    ).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload["metadata"]


class AmortizedDeepBSDE(nn.Module):
    def __init__(self, d: int, *, width: int = 256, depth: int = 3) -> None:
        super().__init__()
        self.d = d
        self.width = width
        self.depth = depth
        self.u0 = DenseNet(d, 1, width=width, depth=depth, activation="silu")
        self.z = DenseNet(d + 1, d, width=width, depth=depth, activation="silu", final_gain=0.1)

    def initial_value(self, x: torch.Tensor) -> torch.Tensor:
        return self.u0(x).reshape(-1)

    def z_value(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        return self.z(tx_input(t, x, T))

    def rollout(
        self,
        equation: Equation,
        x0: torch.Tensor,
        normals: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        steps = normals.shape[1]
        dt = equation.T / steps
        root_dt = math.sqrt(dt)
        state = x0
        value = self.initial_value(x0)
        for step in range(steps):
            t = torch.full((len(x0),), step * dt, dtype=x0.dtype, device=x0.device)
            z = self.z_value(t, state, equation.T)
            dw = root_dt * normals[:, step, :]
            value = value + 0.5 * torch.sum(z.square(), dim=-1) * dt + torch.sum(z * dw, dim=-1)
            state = state + equation.sigma * dw
        return value, state


def _bsde_validation(
    model: AmortizedDeepBSDE,
    equation: Equation,
    x: torch.Tensor,
    normal: torch.Tensor,
) -> float:
    model.eval()
    with torch.no_grad():
        terminal_y, terminal_x = model.rollout(equation, x, normal)
        error = terminal_y - torch_terminal(equation, terminal_x)
        return float(torch.mean(error.square()).item())


def train_deep_bsde(
    equation: Equation,
    *,
    seed: int,
    device: torch.device,
    output_dir: Path,
    max_steps: int = 10_000,
    checkpoints: tuple[int, ...] = (2_000, 5_000, 10_000),
    width: int = 256,
    time_steps: int = 32,
    batch_size: int = 128,
) -> dict[str, Any]:
    seed_everything(seed)
    model = AmortizedDeepBSDE(equation.d, width=width, depth=3).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[5_000, 8_000], gamma=0.2
    )
    generator = torch.Generator(device=device).manual_seed(seed)
    validation_generator = torch.Generator(device=device).manual_seed(seed + 1)
    valid_x = torch.randn((256, equation.d), device=device, generator=validation_generator)
    valid_normal = torch.randn(
        (256, time_steps, equation.d), device=device, generator=validation_generator
    )
    history: list[dict[str, Any]] = []
    saved: list[dict[str, Any]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for step in range(1, max_steps + 1):
        model.train()
        x0 = torch.randn((batch_size, equation.d), device=device, generator=generator)
        normal = torch.randn(
            (batch_size, time_steps, equation.d), device=device, generator=generator
        )
        terminal_y, terminal_x = model.rollout(equation, x0, normal)
        target = torch_terminal(equation, terminal_x)
        loss = torch.nn.functional.huber_loss(terminal_y, target, delta=10.0)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        if step % 250 == 0:
            history.append(
                {
                    "step": step,
                    "training_terminal_huber": float(loss.detach().item()),
                    "learning_rate": float(optimizer.param_groups[0]["lr"]),
                    "elapsed_seconds": time.perf_counter() - started,
                }
            )
        if step in checkpoints:
            validation = _bsde_validation(model, equation, valid_x, valid_normal)
            metadata = {
                "method": "deep_bsde",
                "adaptation": "amortized u0(x0) with shared z(t,x)",
                "seed": seed,
                "step": step,
                "time_steps": time_steps,
                "training_seconds": time.perf_counter() - started,
                "validation_terminal_mse": validation,
                "parameter_count": count_parameters(model),
                "peak_gpu_bytes": _peak(device),
                "analytic_reference_used_for_training_or_selection": False,
            }
            path = output_dir / f"deep_bsde_step{step}.pt"
            atomic_torch_save(
                path,
                {
                    "state_dict": clone_state_dict(model),
                    "architecture": {"d": equation.d, "width": width, "depth": 3},
                    "metadata": metadata,
                },
            )
            saved.append({"path": str(path), **metadata})
    return {
        "method": "deep_bsde",
        "seed": seed,
        "history": history,
        "checkpoints": saved,
        "selected_step": min(saved, key=lambda item: item["validation_terminal_mse"])["step"],
        "selection_metric": "minimum frozen terminal mismatch MSE",
    }


def load_deep_bsde(path: Path, device: torch.device) -> tuple[AmortizedDeepBSDE, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = AmortizedDeepBSDE(
        int(architecture["d"]),
        width=int(architecture["width"]),
        depth=int(architecture["depth"]),
    ).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload["metadata"]


class TerminalAwareDPI(nn.Module):
    """PISGRADNET-style joint value/z output with exact terminal blending."""

    def __init__(self, d: int, *, width: int = 512, depth: int = 4) -> None:
        super().__init__()
        self.d = d
        self.width = width
        self.depth = depth
        self.network = DenseNet(
            d + 1, d + 1, width=width, depth=depth, activation="elu"
        )

    def value_z(
        self, equation: Equation, t: torch.Tensor, x: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        correction = self.network(tx_input(t, x, equation.T))
        tau = (equation.T - t).reshape(-1, 1)
        value = torch_terminal(equation, x) + tau[:, 0] * correction[:, 0]
        z = torch_terminal_z(equation, x) + tau * correction[:, 1:]
        return value, z


def _dpi_previous(
    previous: TerminalAwareDPI | None,
    equation: Equation,
    t: torch.Tensor,
    x: torch.Tensor,
) -> torch.Tensor:
    if previous is None:
        return torch.zeros_like(x)
    with torch.no_grad():
        return previous.value_z(equation, t, x)[1]


def generate_dpi_labels(
    equation: Equation,
    previous: TerminalAwareDPI | None,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    paths: int,
    generator: torch.Generator,
) -> tuple[torch.Tensor, torch.Tensor]:
    batch, d = x.shape
    tau = torch.clamp(equation.T - t, min=1e-6)
    normal = torch.randn(
        (batch, paths, d), dtype=x.dtype, device=x.device, generator=generator
    )
    terminal_x = x[:, None, :] + torch.sqrt(2.0 * tau)[:, None, None] * normal
    terminal = torch_terminal(equation, terminal_x.reshape(-1, d)).reshape(batch, paths)
    parent_terminal = torch_terminal(equation, x)
    value = terminal.mean(dim=1)
    z = torch.mean(
        (terminal - parent_terminal[:, None])[:, :, None]
        * normal
        / torch.sqrt(tau)[:, None, None],
        dim=1,
    )
    if previous is not None:
        uniform = torch.clamp(
            torch.rand((batch, paths), device=x.device, generator=generator), min=1e-5
        )
        delta = tau[:, None] * uniform
        integral_normal = torch.randn(
            (batch, paths, d), dtype=x.dtype, device=x.device, generator=generator
        )
        sample_t = t[:, None] + delta
        sample_x = x[:, None, :] + torch.sqrt(2.0 * delta)[:, :, None] * integral_normal
        sample_z = _dpi_previous(
            previous,
            equation,
            sample_t.reshape(-1),
            sample_x.reshape(-1, d),
        ).reshape(batch, paths, d)
        sample_driver = -0.5 * torch.sum(sample_z.square(), dim=-1)
        parent_z = _dpi_previous(previous, equation, t, x)
        parent_driver = -0.5 * torch.sum(parent_z.square(), dim=-1)
        value = value + tau * sample_driver.mean(dim=1)
        z = z + tau[:, None] * torch.mean(
            (sample_driver - parent_driver[:, None])[:, :, None]
            * integral_normal
            / torch.sqrt(delta)[:, :, None],
            dim=1,
        )
    return value.detach(), z.detach()


def _dpi_dataset(
    equation: Equation,
    previous: TerminalAwareDPI | None,
    *,
    points: int,
    paths: int,
    generator: torch.Generator,
    device: torch.device,
    parent_batch: int = 8,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, float]:
    t = torch.rand(points, device=device, generator=generator) * (equation.T - 1e-5)
    x = torch.randn((points, equation.d), device=device, generator=generator)
    value = torch.empty(points, device=device)
    z = torch.empty((points, equation.d), device=device)
    started = time.perf_counter()
    for start in range(0, points, parent_batch):
        stop = min(start + parent_batch, points)
        value[start:stop], z[start:stop] = generate_dpi_labels(
            equation,
            previous,
            t[start:stop],
            x[start:stop],
            paths=paths,
            generator=generator,
        )
    return t, x, value, z, time.perf_counter() - started


def train_dpi(
    equation: Equation,
    *,
    seed: int,
    device: torch.device,
    output_dir: Path,
    rounds: int = 20,
    training_points: int = 4096,
    paths: int = 512,
    epochs: int = 16,
    batch_size: int = 512,
    width: int = 512,
    checkpoint_rounds: tuple[int, ...] = (1, 2, 5, 10, 20),
) -> dict[str, Any]:
    seed_everything(seed)
    generator = torch.Generator(device=device).manual_seed(seed)
    valid_generator = torch.Generator(device=device).manual_seed(seed + 1)
    previous: TerminalAwareDPI | None = None
    history: list[dict[str, Any]] = []
    saved: list[dict[str, Any]] = []
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    label_seconds_total = 0.0
    for picard_round in range(1, rounds + 1):
        train_t, train_x, train_u, train_z, label_seconds = _dpi_dataset(
            equation,
            previous,
            points=training_points,
            paths=paths,
            generator=generator,
            device=device,
        )
        valid_t, valid_x, valid_u, valid_z, valid_label_seconds = _dpi_dataset(
            equation,
            previous,
            points=256,
            paths=paths,
            generator=valid_generator,
            device=device,
        )
        label_seconds_total += label_seconds + valid_label_seconds
        model = TerminalAwareDPI(equation.d, width=width, depth=4).to(device)
        if previous is not None:
            model.load_state_dict(previous.state_dict())
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        round_started = time.perf_counter()
        last_loss = float("nan")
        for _epoch in range(epochs):
            order = torch.randperm(training_points, device=device, generator=generator)
            for start in range(0, training_points, batch_size):
                index = order[start : start + batch_size]
                predicted_u, predicted_z = model.value_z(
                    equation, train_t[index], train_x[index]
                )
                value_loss = torch.mean((predicted_u - train_u[index]).square())
                gradient_loss = torch.mean((predicted_z - train_z[index]).square())
                loss = value_loss + 0.1 * gradient_loss
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
                optimizer.step()
                last_loss = float(loss.detach().item())
        model.eval()
        with torch.no_grad():
            predicted_u, predicted_z = model.value_z(equation, valid_t, valid_x)
            validation = float(
                (
                    torch.mean((predicted_u - valid_u).square())
                    + 0.1 * torch.mean((predicted_z - valid_z).square())
                ).item()
            )
        round_record = {
            "round": picard_round,
            "training_loss_last": last_loss,
            "validation_target_loss": validation,
            "label_seconds": label_seconds + valid_label_seconds,
            "supervised_seconds": time.perf_counter() - round_started,
            "cumulative_seconds": time.perf_counter() - started,
        }
        history.append(round_record)
        previous = model
        if picard_round in checkpoint_rounds:
            metadata = {
                "method": "dpi",
                "seed": seed,
                "round": picard_round,
                "training_seconds": round_record["cumulative_seconds"],
                "validation_target_loss": validation,
                "paths_per_label": paths,
                "training_points_per_round": training_points,
                "epochs_per_round": epochs,
                "parameter_count": count_parameters(model),
                "peak_gpu_bytes": _peak(device),
                "analytic_reference_used_for_training_or_selection": False,
            }
            path = output_dir / f"dpi_round{picard_round}.pt"
            atomic_torch_save(
                path,
                {
                    "state_dict": clone_state_dict(model),
                    "architecture": {"d": equation.d, "width": width, "depth": 4},
                    "metadata": metadata,
                },
            )
            saved.append({"path": str(path), **metadata})
    return {
        "method": "dpi",
        "seed": seed,
        "history": history,
        "checkpoints": saved,
        # Round-local target losses are not comparable across rounds.  The
        # preset ceiling is therefore the eligible final checkpoint.
        "selected_round": rounds,
        "selection_metric": "preset Picard-round ceiling; within-round targets only",
        "label_seconds_total": label_seconds_total,
    }


def load_dpi(path: Path, device: torch.device) -> tuple[TerminalAwareDPI, dict[str, Any]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    architecture = payload["architecture"]
    model = TerminalAwareDPI(
        int(architecture["d"]),
        width=int(architecture["width"]),
        depth=int(architecture["depth"]),
    ).to(device)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload["metadata"]


def predict_baseline(
    method: str,
    model: nn.Module,
    equation: Equation,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray | None, float]:
    state = torch.as_tensor(
        np.array(x, dtype=np.float32, copy=True, order="C"),
        dtype=torch.float32,
        device=device,
    )
    t = torch.zeros(len(state), dtype=torch.float32, device=device)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    started = time.perf_counter()
    if method == "pinn_sdgd":
        local = state.detach().requires_grad_(True)
        value = model.value(equation, t, local)
        gradient = torch.autograd.grad(value.sum(), local)[0]
        z = equation.sigma * gradient
    elif method == "deep_bsde":
        with torch.no_grad():
            value = model.initial_value(state)
            z = model.z_value(t, state, equation.T)
    elif method == "dpi":
        with torch.no_grad():
            value, z = model.value_z(equation, t, state)
    else:
        raise ValueError(method)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - started
    return value.detach().cpu().numpy(), z.detach().cpu().numpy(), elapsed
