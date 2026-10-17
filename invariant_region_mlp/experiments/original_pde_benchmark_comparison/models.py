"""Full-ambient Terminal-Aware Actors used only by our method."""

from __future__ import annotations

import copy
from typing import Literal

import torch
from torch import nn

from .benchmarks import OriginalBenchmark


Activation = Literal["silu", "elu", "tanh"]


def _activation(name: Activation) -> nn.Module:
    if name == "silu":
        return nn.SiLU()
    if name == "elu":
        return nn.ELU()
    if name == "tanh":
        return nn.Tanh()
    raise KeyError(name)


class TerminalAwareActor(nn.Module):
    """Terminal value/gradient plus a vanishing learned residual.

    The network always receives all 100 ambient coordinates.  No norm,
    projection, GMM responsibility, or exact nonterminal feature is exposed.
    """

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
        layers: list[nn.Module] = []
        previous = benchmark.d + 1
        for _ in range(depth):
            linear = nn.Linear(previous, width)
            nn.init.xavier_uniform_(linear.weight)
            nn.init.zeros_(linear.bias)
            layers.extend((linear, _activation(activation)))
            previous = width
        final = nn.Linear(previous, output_dim)
        # A newly constructed model is the deterministic terminal Actor a0.
        nn.init.zeros_(final.weight)
        nn.init.zeros_(final.bias)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def _raw(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        inputs = torch.cat(
            (t.reshape(-1, 1) / self.benchmark.T, x / self.benchmark.input_scale),
            dim=-1,
        )
        return self.network(inputs)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor
    ) -> tuple[torch.Tensor | None, torch.Tensor | None]:
        tau = torch.clamp(self.benchmark.T - t.reshape(-1), min=0.0)
        raw = self._raw(t, x)
        if self.benchmark.actor_kind == "gradient":
            terminal_z = self.benchmark.sigma * self.benchmark.terminal_gradient(x)
            return None, terminal_z + tau[:, None] * raw
        bounds = self.benchmark.value_bounds
        if bounds is None:
            return self.benchmark.terminal(x) + tau * raw[:, 0], None
        lower, upper = bounds
        terminal = self.benchmark.terminal(x)
        unit = (terminal - lower) / (upper - lower)
        epsilon = max(torch.finfo(x.dtype).eps, 1e-7)
        unit = unit.clamp(epsilon, 1.0 - epsilon)
        logit = torch.log(unit) - torch.log1p(-unit)
        r = lower + (upper - lower) * torch.sigmoid(logit + tau * raw[:, 0])
        return r, None

    def clone(self) -> "TerminalAwareActor":
        return copy.deepcopy(self)

    def specification(self) -> dict[str, object]:
        return {
            "kind": "terminal_aware_full_ambient_tangent_actor",
            "width": self.width,
            "depth": self.depth,
            "activation": self.activation_name,
            "input_dimension": self.benchmark.d + 1,
            "output_dimension": self.benchmark.d
            if self.benchmark.actor_kind == "gradient"
            else 1,
            "parameters": sum(parameter.numel() for parameter in self.parameters()),
            "oracle_projection": False,
            "analytic_nonterminal_features": False,
        }


def load_actor(
    path: str, benchmark: OriginalBenchmark, device: torch.device
) -> tuple[TerminalAwareActor, dict[str, object]]:
    payload = torch.load(path, map_location=device, weights_only=False)
    actor = TerminalAwareActor(
        benchmark,
        width=int(payload["width"]),
        depth=int(payload["depth"]),
        activation=str(payload["activation"]),
    ).to(device)
    actor.load_state_dict(payload["state_dict"])
    actor.eval()
    return actor, payload


__all__ = ["Activation", "TerminalAwareActor", "load_actor"]
