"""Preregistered PINN and float64 frozen-network differential evaluator."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np
import torch
from torch import nn

from lqg import LQGEquation, seed_sequence


CHECKPOINTS = (500, 1000, 2500)


def _torch_seed(sequence: np.random.SeedSequence) -> int:
    return int(sequence.generate_state(1, dtype=np.uint64)[0] % np.uint64(2**63 - 1))


class PINN(nn.Module):
    """Five-hidden-layer, width-50 tanh network with Glorot-normal weights."""

    def __init__(self, d: int, width: int = 50, depth: int = 5) -> None:
        super().__init__()
        sizes = [d + 1] + [width] * depth + [1]
        self.layers = nn.ModuleList(
            [nn.Linear(sizes[index], sizes[index + 1]) for index in range(len(sizes) - 1)]
        )
        for layer in self.layers:
            nn.init.xavier_normal_(layer.weight)
            nn.init.zeros_(layer.bias)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        value = inputs
        for layer in self.layers[:-1]:
            value = torch.tanh(layer(value))
        return self.layers[-1](value)


def _unit_ball_torch(generator: torch.Generator, count: int, d: int, device: torch.device) -> torch.Tensor:
    direction = torch.randn((count, d), generator=generator, device=device)
    direction = direction / direction.norm(dim=1, keepdim=True).clamp_min(1e-30)
    radius = torch.rand((count, 1), generator=generator, device=device).pow(1.0 / d)
    return direction * radius


def _terminal_torch(eq: LQGEquation, x: torch.Tensor, c1: torch.Tensor, c2: torch.Tensor) -> torch.Tensor:
    quadratic = torch.sum(
        c1 * (x[:, :-1] - x[:, 1:]).square() + c2 * x[:, 1:].square(),
        dim=1,
    )
    return torch.log((1.0 + quadratic) / 2.0)


def _residual(
    model: PINN,
    inputs: torch.Tensor,
    probes: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return PDE residual and its components using batched Hutchinson VJPs."""

    inputs = inputs.requires_grad_(True)
    u = model(inputs)[:, 0]
    first = torch.autograd.grad(u.sum(), inputs, create_graph=True)[0]
    grad_x = first[:, :-1]
    u_t = first[:, -1]
    # Each row of the pointwise network is independent.  A batched VJP with
    # shape (probe,batch,d) therefore returns every per-point Hessian-vector
    # product without a Python loop over probes.
    hvp = torch.autograd.grad(
        grad_x,
        inputs,
        grad_outputs=probes,
        create_graph=True,
        is_grads_batched=True,
    )[0][..., :-1]
    laplacian = torch.mean(torch.sum(hvp * probes, dim=-1), dim=0)
    residual = u_t + laplacian - torch.sum(grad_x.square(), dim=1)
    return residual, u, grad_x, laplacian


def train_surrogate(
    eq: LQGEquation,
    network_seed: int,
    output_dir: Path,
    *,
    device: str = "cuda",
    iterations: int = 2500,
    interior_count: int = 100,
    terminal_count: int = 1000,
    progress_every: int = 50,
) -> dict[str, Any]:
    """Train one preregistered PINN, writing immutable resume checkpoints."""

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = output_dir / "training.json"
    final_path = output_dir / f"checkpoint_{iterations}.pt"
    if final_path.exists() and metadata_path.exists():
        return json.loads(metadata_path.read_text(encoding="utf-8"))

    torch_device = torch.device(device if device == "cpu" or torch.cuda.is_available() else "cpu")
    torch.manual_seed(_torch_seed(seed_sequence(eq.d, network_seed, 10)))
    if torch_device.type == "cuda":
        torch.cuda.manual_seed_all(_torch_seed(seed_sequence(eq.d, network_seed, 10)))
    model = PINN(eq.d).to(torch_device, dtype=torch.float32)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, betas=(0.9, 0.99))
    sample_generator = torch.Generator(device=torch_device)
    sample_generator.manual_seed(_torch_seed(seed_sequence(eq.d, network_seed, 11)))
    probe_generator = torch.Generator(device=torch_device)
    probe_generator.manual_seed(_torch_seed(seed_sequence(eq.d, network_seed, 12)))
    c1 = torch.as_tensor(eq.c1, dtype=torch.float32, device=torch_device)
    c2 = torch.as_tensor(eq.c2, dtype=torch.float32, device=torch_device)
    probe_count = math.ceil(eq.d / 4)
    history: list[dict[str, float | int]] = []
    started = time.perf_counter()
    forward_points = 0
    backward_steps = 0
    for step in range(1, iterations + 1):
        x_interior = _unit_ball_torch(sample_generator, interior_count, eq.d, torch_device)
        t_interior = torch.rand(
            (interior_count, 1), generator=sample_generator, device=torch_device
        ) * eq.T
        inputs = torch.cat([x_interior, t_interior], dim=1)
        probes = torch.randint(
            0,
            2,
            (probe_count, interior_count, eq.d),
            generator=probe_generator,
            device=torch_device,
            dtype=torch.int8,
        ).to(torch.float32)
        probes.mul_(2.0).sub_(1.0)
        residual, _, _, _ = _residual(model, inputs, probes)
        x_terminal = _unit_ball_torch(sample_generator, terminal_count, eq.d, torch_device)
        terminal_inputs = torch.cat(
            [x_terminal, torch.full((terminal_count, 1), eq.T, device=torch_device)],
            dim=1,
        )
        terminal_prediction = model(terminal_inputs)[:, 0]
        terminal_truth = _terminal_torch(eq, x_terminal, c1, c2)
        pde_loss = torch.mean(residual.square())
        terminal_loss = torch.mean((terminal_prediction - terminal_truth).square())
        loss = pde_loss + terminal_loss
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        forward_points += interior_count + terminal_count
        backward_steps += 1
        if step == 1 or step % progress_every == 0 or step in CHECKPOINTS:
            row = {
                "step": step,
                "loss": float(loss.detach().cpu()),
                "pde_loss": float(pde_loss.detach().cpu()),
                "terminal_loss": float(terminal_loss.detach().cpu()),
                "elapsed_seconds": time.perf_counter() - started,
            }
            history.append(row)
            print(
                f"TRAIN d={eq.d} seed={network_seed} step={step}/{iterations} "
                f"loss={row['loss']:.6g} pde={row['pde_loss']:.6g} "
                f"terminal={row['terminal_loss']:.6g} elapsed={row['elapsed_seconds']:.1f}s",
                flush=True,
            )
        if step in CHECKPOINTS or step == iterations:
            temporary = output_dir / f"checkpoint_{step}.tmp.pt"
            torch.save(
                {
                    "state_dict": model.state_dict(),
                    "d": eq.d,
                    "network_seed": network_seed,
                    "step": step,
                    "architecture": [eq.d + 1] + [50] * 5 + [1],
                    "dtype_training": "float32",
                },
                temporary,
            )
            temporary.replace(output_dir / f"checkpoint_{step}.pt")
    elapsed = time.perf_counter() - started
    payload: dict[str, Any] = {
        "dimension": eq.d,
        "network_seed": network_seed,
        "iterations": iterations,
        "interior_per_iteration": interior_count,
        "terminal_per_iteration": terminal_count,
        "hutchinson_probes": probe_count,
        "optimizer": {"name": "Adam", "lr": 1e-3, "betas": [0.9, 0.99]},
        "network": {"hidden_layers": 5, "width": 50, "activation": "tanh", "init": "Glorot normal"},
        "device": str(torch_device),
        "torch_version": torch.__version__,
        "wall_clock_seconds": elapsed,
        "network_forward_points": forward_points,
        "network_backward_steps": backward_steps,
        "history": history,
        "seed_scheme": "SeedSequence([20261101,1,1,d,network_seed,component,...])",
        "components": {"network_init": 10, "collocation": 11, "hutchinson": 12},
    }
    temporary_json = metadata_path.with_suffix(".tmp.json")
    temporary_json.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    temporary_json.replace(metadata_path)
    return payload


@dataclass
class NetworkCounters:
    forward_calls: int = 0
    forward_points: int = 0
    derivative_calls: int = 0
    derivative_points: int = 0
    laplacian_probe_points: int = 0

    def as_dict(self) -> dict[str, int]:
        return {
            "forward_calls": self.forward_calls,
            "forward_points": self.forward_points,
            "backward_calls": self.derivative_calls,
            "backward_points": self.derivative_points,
            "laplacian_probe_points": self.laplacian_probe_points,
        }


@dataclass
class FrozenMLP:
    """Float64 manual differentiator for a trained tanh network.

    Manual directional propagation avoids building autograd graphs inside the
    recursive Monte Carlo solver while still evaluating the cast float64
    checkpoint exactly.  Spatial Laplacians use either the coordinate basis or
    fixed component-isolated Rademacher probes.
    """

    weights: list[torch.Tensor]
    biases: list[torch.Tensor]
    d: int
    device: torch.device
    probes: torch.Tensor
    exact_laplacian: bool = False
    max_batch: int = 4096
    counters: NetworkCounters = field(default_factory=NetworkCounters)

    @classmethod
    def load(
        cls,
        checkpoint: Path,
        *,
        d: int,
        network_seed: int,
        device: str = "cuda",
        exact_laplacian: bool | None = None,
        max_batch: int = 4096,
    ) -> "FrozenMLP":
        target = torch.device(device if device == "cpu" or torch.cuda.is_available() else "cpu")
        payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
        state = payload["state_dict"]
        indices = sorted(
            int(key.split(".")[1]) for key in state if key.startswith("layers.") and key.endswith(".weight")
        )
        weights = [state[f"layers.{index}.weight"].to(device=target, dtype=torch.float64) for index in indices]
        biases = [state[f"layers.{index}.bias"].to(device=target, dtype=torch.float64) for index in indices]
        if exact_laplacian is None:
            exact_laplacian = d <= 50
        probe_count = d if exact_laplacian else math.ceil(d / 4)
        if exact_laplacian:
            probes = torch.eye(d, dtype=torch.float64, device=target)
        else:
            generator = torch.Generator(device=target)
            generator.manual_seed(_torch_seed(seed_sequence(d, network_seed, 30, int(payload["step"]))))
            probes = torch.randint(
                0, 2, (probe_count, d), generator=generator, device=target, dtype=torch.int8
            ).to(torch.float64)
            probes.mul_(2.0).sub_(1.0)
        return cls(
            weights=weights,
            biases=biases,
            d=d,
            device=target,
            probes=probes,
            exact_laplacian=bool(exact_laplacian),
            max_batch=max_batch,
        )

    def _forward_cache(self, inputs: torch.Tensor) -> tuple[torch.Tensor, list[torch.Tensor], list[torch.Tensor]]:
        hidden = inputs
        derivatives: list[torch.Tensor] = []
        second_derivatives: list[torch.Tensor] = []
        for weight, bias in zip(self.weights[:-1], self.biases[:-1]):
            hidden = torch.tanh(hidden @ weight.T + bias)
            first = 1.0 - hidden.square()
            derivatives.append(first)
            second_derivatives.append(-2.0 * hidden * first)
        value = hidden @ self.weights[-1].T + self.biases[-1]
        return value[:, 0], derivatives, second_derivatives

    def _one_chunk(
        self,
        t: np.ndarray,
        x: np.ndarray,
        *,
        need_laplacian: bool,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
        xt = torch.as_tensor(
            np.concatenate([x, t[:, None]], axis=1), dtype=torch.float64, device=self.device
        )
        with torch.no_grad():
            hidden = xt
            hiddens: list[torch.Tensor] = []
            first_phi: list[torch.Tensor] = []
            second_phi: list[torch.Tensor] = []
            for weight, bias in zip(self.weights[:-1], self.biases[:-1]):
                hidden = torch.tanh(hidden @ weight.T + bias)
                hiddens.append(hidden)
                first = 1.0 - hidden.square()
                first_phi.append(first)
                second_phi.append(-2.0 * hidden * first)
            value = (hidden @ self.weights[-1].T + self.biases[-1])[:, 0]
            adjoint = self.weights[-1][0].expand(len(xt), -1)
            for layer_index in range(len(self.weights) - 2, -1, -1):
                delta = adjoint * first_phi[layer_index]
                adjoint = delta @ self.weights[layer_index]
            gradient = adjoint
            laplacian: torch.Tensor | None = None
            if need_laplacian:
                p = len(self.probes)
                direction = torch.zeros(
                    (len(xt), p, self.d + 1), dtype=torch.float64, device=self.device
                )
                direction[:, :, : self.d] = self.probes[None, :, :]
                second = torch.zeros_like(direction)
                for layer_index, weight in enumerate(self.weights[:-1]):
                    direction_affine = torch.matmul(direction, weight.T)
                    second_affine = torch.matmul(second, weight.T)
                    second = (
                        first_phi[layer_index][:, None, :] * second_affine
                        + second_phi[layer_index][:, None, :] * direction_affine.square()
                    )
                    direction = first_phi[layer_index][:, None, :] * direction_affine
                second_out = torch.matmul(second, self.weights[-1].T)[..., 0]
                # Coordinate-basis directions partition the trace and must be
                # summed.  Rademacher directions are independent trace
                # estimates and must be averaged.
                laplacian = (
                    torch.sum(second_out, dim=1)
                    if self.exact_laplacian
                    else torch.mean(second_out, dim=1)
                )
                self.counters.laplacian_probe_points += len(xt) * p
            self.counters.forward_calls += 1
            self.counters.forward_points += len(xt)
            self.counters.derivative_calls += 1
            self.counters.derivative_points += len(xt)
            return (
                value.cpu().numpy(),
                (math.sqrt(2.0) * gradient[:, : self.d]).cpu().numpy(),
                gradient[:, self.d].cpu().numpy(),
                None if laplacian is None else laplacian.cpu().numpy(),
            )

    def evaluate(
        self,
        t: np.ndarray,
        x: np.ndarray,
        *,
        need_laplacian: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.d)
        pieces: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]] = []
        # The directional tensor is O(batch * probes * width), so high-d
        # Laplacians use smaller chunks than value/gradient-only calls.
        chunk = min(self.max_batch, 512 if need_laplacian else self.max_batch)
        for begin in range(0, len(t), chunk):
            end = min(begin + chunk, len(t))
            pieces.append(self._one_chunk(t[begin:end], x[begin:end], need_laplacian=need_laplacian))
        value = np.concatenate([piece[0] for piece in pieces])
        z = np.concatenate([piece[1] for piece in pieces])
        u_t = np.concatenate([piece[2] for piece in pieces])
        lap = None
        if need_laplacian:
            lap = np.concatenate([piece[3] for piece in pieces if piece[3] is not None])
        return value, z, u_t, lap

    def value_z(self, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        value, z, _, _ = self.evaluate(t, x, need_laplacian=False)
        return value, z

    def residual(self, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        value, z, u_t, lap = self.evaluate(t, x, need_laplacian=True)
        if lap is None:
            raise RuntimeError("laplacian was not evaluated")
        residual = u_t + lap - 0.5 * np.sum(z * z, axis=1)
        return residual, value, z

    def terminal_defect(self, eq: LQGEquation, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        t = np.full(len(x), eq.T, dtype=np.float64)
        value, z = self.value_z(t, x)
        return eq.terminal(x) - value, eq.sigma_grad_terminal(x) - z


def checkpoint_paths(root: Path, d: int, network_seed: int) -> Iterable[tuple[int, Path]]:
    directory = root / "networks" / f"d{d}" / f"seed{network_seed}"
    for step in CHECKPOINTS:
        yield step, directory / f"checkpoint_{step}.pt"
