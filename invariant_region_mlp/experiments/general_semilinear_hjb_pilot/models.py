"""Generic actor networks with separate state and action dimensions."""

from __future__ import annotations

import torch
from torch import nn

from .problems import ControlledHJB


class ActorNet(nn.Module):
    def __init__(self, problem: ControlledHJB, *, width: int, depth: int) -> None:
        super().__init__()
        self.d, self.m = problem.d, problem.m
        self.width, self.depth = width, depth
        self.problem = problem
        layers: list[nn.Module] = []
        previous = problem.d + 1
        for _ in range(depth):
            layer = nn.Linear(previous, width)
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)
            layers.extend((layer, nn.SiLU()))
            previous = width
        final = nn.Linear(previous, problem.m)
        nn.init.zeros_(final.weight)
        nn.init.zeros_(final.bias)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def raw(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        return self.network(inputs)

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        return self.problem.feasible_action(self.raw(t, x, T))


class ZeroPolicy(nn.Module):
    def __init__(self, problem: ControlledHJB) -> None:
        super().__init__()
        self.problem = problem
        self.d, self.m = problem.d, problem.m

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        del t, T
        return torch.zeros((*x.shape[:-1], self.m), dtype=x.dtype, device=x.device)


def count_parameters(module: nn.Module) -> int:
    return sum(parameter.numel() for parameter in module.parameters())

