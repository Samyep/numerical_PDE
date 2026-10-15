"""Standard and terminal-aware full-ambient Tangent Actors."""

from __future__ import annotations

from typing import Literal

import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)


ActorKind = Literal["standard", "terminal_aware"]


def _residual_network(input_dim: int, output_dim: int, width: int, depth: int) -> nn.Sequential:
    layers: list[nn.Module] = []
    previous = input_dim
    for _ in range(depth):
        layer = nn.Linear(previous, width)
        nn.init.xavier_uniform_(layer.weight)
        nn.init.zeros_(layer.bias)
        layers.extend((layer, nn.SiLU()))
        previous = width
    final = nn.Linear(previous, output_dim)
    nn.init.zeros_(final.weight)
    nn.init.zeros_(final.bias)
    layers.append(final)
    return nn.Sequential(*layers)


class StandardActor(nn.Module):
    actor_kind: ActorKind = "standard"

    def __init__(self, driver: ConcaveSemilinearDriver, *, width: int, depth: int) -> None:
        super().__init__()
        self.driver = driver
        self.width, self.depth = width, depth
        output_dim = 1 + (driver.d if driver.has_gradient else 0)
        self.network = _residual_network(driver.d + 1, output_dim, width, depth)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        output = self.network(inputs)
        r = self.driver.constrain_r(output[:, 0])
        q = output[:, 1:] if self.driver.has_gradient else None
        return r, q


class TerminalAwareActor(nn.Module):
    actor_kind: ActorKind = "terminal_aware"

    def __init__(self, driver: ConcaveSemilinearDriver, *, width: int, depth: int) -> None:
        super().__init__()
        self.driver = driver
        self.width, self.depth = width, depth
        output_dim = 1 + (driver.d if driver.has_gradient else 0)
        self.network = _residual_network(driver.d + 1, output_dim, width, depth)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        output = self.network(inputs)
        tau = torch.clamp(T - t.reshape(-1), min=0.0)
        terminal = self.driver.terminal(x)
        span = self.driver.r_high - self.driver.r_low
        unit_terminal = (terminal - self.driver.r_low) / span
        # The benchmark terminal values lie strictly inside the tangent domain.
        # Clamp only guards roundoff in generic use and is inactive here.
        eps = torch.finfo(x.dtype).eps
        unit_terminal = torch.clamp(unit_terminal, min=eps, max=1.0 - eps)
        terminal_logit = torch.log(unit_terminal) - torch.log1p(-unit_terminal)
        r = self.driver.r_low + span * torch.sigmoid(
            terminal_logit + tau * output[:, 0]
        )
        if self.driver.has_gradient:
            terminal_q = self.driver.sigma * self.driver.terminal_gradient_torch(x)
            q = terminal_q + tau[:, None] * output[:, 1:]
        else:
            q = None
        return r, q


def make_actor(
    kind: ActorKind,
    driver: ConcaveSemilinearDriver,
    *,
    width: int,
    depth: int,
) -> StandardActor | TerminalAwareActor:
    if kind == "standard":
        return StandardActor(driver, width=width, depth=depth)
    if kind == "terminal_aware":
        return TerminalAwareActor(driver, width=width, depth=depth)
    raise KeyError(kind)


def count_parameters(module: nn.Module) -> int:
    return sum(parameter.numel() for parameter in module.parameters())

