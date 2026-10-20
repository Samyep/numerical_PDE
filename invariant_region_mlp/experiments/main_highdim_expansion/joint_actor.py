"""Preregistered unbounded value-gradient terminal-aware actor."""

from __future__ import annotations

import copy

import torch
from torch import nn

from .benchmarks import DenseFourierDriver


class JointTerminalActor(nn.Module):
    """r=g+(T-t)R and q=sqrt(2) grad(g)+(T-t)Q exactly."""

    actor_kind = "terminal_aware_joint_unbounded"

    def __init__(self, driver: DenseFourierDriver, *, width: int, depth: int) -> None:
        super().__init__()
        self.driver = driver
        self.width, self.depth = int(width), int(depth)
        layers: list[nn.Module] = []
        previous = driver.d + 1
        for _ in range(self.depth):
            linear = nn.Linear(previous, self.width)
            nn.init.xavier_uniform_(linear.weight)
            nn.init.zeros_(linear.bias)
            layers.extend((linear, nn.SiLU()))
            previous = self.width
        final = nn.Linear(previous, 1 + driver.d)
        nn.init.zeros_(final.weight)
        nn.init.zeros_(final.bias)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float | None = None
    ) -> tuple[torch.Tensor, torch.Tensor]:
        terminal_time = self.driver.T if T is None else float(T)
        if abs(terminal_time - self.driver.T) > 1e-15:
            raise ValueError("actor terminal time differs from its driver")
        time = t.reshape(-1)
        raw = self.network(torch.cat((time[:, None] / terminal_time, x), dim=-1))
        tau = torch.clamp(terminal_time - time, min=0.0)
        r = self.driver.terminal(x) + tau * raw[:, 0]
        q = (
            self.driver.sigma * self.driver.terminal_gradient_torch(x)
            + tau[:, None] * raw[:, 1:]
        )
        return r, q

    def clone(self) -> "JointTerminalActor":
        return copy.deepcopy(self)

    def specification(self) -> dict[str, object]:
        return {
            "kind": self.actor_kind,
            "width": self.width,
            "depth": self.depth,
            "activation": "silu",
            "input_dimension": self.driver.d + 1,
            "output_dimension": self.driver.d + 1,
            "parameters": sum(parameter.numel() for parameter in self.parameters()),
            "value_anchor": "g(x)+(T-t)R_theta(t,x)",
            "gradient_anchor": "sqrt(2)grad(g)(x)+(T-t)Q_theta(t,x)",
            "value_projection": False,
            "oracle_projection": False,
        }


__all__ = ["JointTerminalActor"]
