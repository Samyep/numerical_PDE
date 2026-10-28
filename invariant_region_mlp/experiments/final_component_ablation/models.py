"""Matched no-hard-terminal-anchor Actors for the controlled ablation."""

from __future__ import annotations

import copy
from typing import Literal

import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.benchmarks import (
    OriginalBenchmark,
)


Activation = Literal["silu", "elu", "tanh"]


def _activation(name: Activation) -> nn.Module:
    if name == "silu":
        return nn.SiLU()
    if name == "elu":
        return nn.ELU()
    if name == "tanh":
        return nn.Tanh()
    raise KeyError(name)


def _network(
    input_dim: int,
    output_dim: int,
    width: int,
    depth: int,
    activation: Activation = "silu",
) -> nn.Sequential:
    layers: list[nn.Module] = []
    previous = input_dim
    for _ in range(depth):
        linear = nn.Linear(previous, width)
        nn.init.xavier_uniform_(linear.weight)
        nn.init.zeros_(linear.bias)
        layers.extend((linear, _activation(activation)))
        previous = width
    final = nn.Linear(previous, output_dim)
    nn.init.zeros_(final.weight)
    nn.init.zeros_(final.bias)
    layers.append(final)
    return nn.Sequential(*layers)


class UnanchoredOriginalActor(nn.Module):
    """Original-benchmark Actor retaining terminal offsets but not `(T-t)`."""

    actor_kind = "terminal_offset_unanchored"

    def __init__(
        self,
        benchmark: OriginalBenchmark,
        *,
        width: int,
        depth: int,
        activation: Activation,
    ) -> None:
        super().__init__()
        self.benchmark = benchmark
        self.width = width
        self.depth = depth
        self.activation_name = activation
        output_dim = benchmark.d if benchmark.actor_kind == "gradient" else 1
        self.network = _network(
            benchmark.d + 1, output_dim, width, depth, activation
        )

    def _raw(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        inputs = torch.cat(
            (t.reshape(-1, 1) / self.benchmark.T, x / self.benchmark.input_scale),
            dim=-1,
        )
        return self.network(inputs)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor
    ) -> tuple[torch.Tensor | None, torch.Tensor | None]:
        raw = self._raw(t, x)
        if self.benchmark.actor_kind == "gradient":
            terminal_z = self.benchmark.sigma * self.benchmark.terminal_gradient(x)
            return None, terminal_z + raw
        bounds = self.benchmark.value_bounds
        terminal = self.benchmark.terminal(x)
        if bounds is None:
            return terminal + raw[:, 0], None
        lower, upper = bounds
        unit = (terminal - lower) / (upper - lower)
        epsilon = max(torch.finfo(x.dtype).eps, 1e-7)
        unit = unit.clamp(epsilon, 1.0 - epsilon)
        logit = torch.log(unit) - torch.log1p(-unit)
        value = lower + (upper - lower) * torch.sigmoid(logit + raw[:, 0])
        return value, None

    def clone(self) -> "UnanchoredOriginalActor":
        return copy.deepcopy(self)

    def specification(self) -> dict[str, object]:
        return {
            "kind": self.actor_kind,
            "width": self.width,
            "depth": self.depth,
            "activation": self.activation_name,
            "input_dimension": self.benchmark.d + 1,
            "output_dimension": self.benchmark.d
            if self.benchmark.actor_kind == "gradient"
            else 1,
            "parameters": sum(parameter.numel() for parameter in self.parameters()),
            "terminal_offset_retained": True,
            "hard_terminal_identity": False,
            "oracle_projection": False,
        }


class UnanchoredJointActor(nn.Module):
    """Joint Actor with identical residual net and terminal offset, no tau."""

    actor_kind = "terminal_offset_unanchored"

    def __init__(
        self, driver: ConcaveSemilinearDriver, *, width: int, depth: int
    ) -> None:
        super().__init__()
        self.driver = driver
        self.width = width
        self.depth = depth
        output_dim = 1 + (driver.d if driver.has_gradient else 0)
        self.network = _network(driver.d + 1, output_dim, width, depth)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        output = self.network(inputs)
        terminal = self.driver.terminal(x)
        span = self.driver.r_high - self.driver.r_low
        unit_terminal = (terminal - self.driver.r_low) / span
        eps = torch.finfo(x.dtype).eps
        unit_terminal = torch.clamp(unit_terminal, min=eps, max=1.0 - eps)
        terminal_logit = torch.log(unit_terminal) - torch.log1p(-unit_terminal)
        r = self.driver.r_low + span * torch.sigmoid(
            terminal_logit + output[:, 0]
        )
        if self.driver.has_gradient:
            terminal_q = self.driver.sigma * self.driver.terminal_gradient_torch(x)
            q = terminal_q + output[:, 1:]
        else:
            q = None
        return r, q


__all__ = ["UnanchoredJointActor", "UnanchoredOriginalActor"]
