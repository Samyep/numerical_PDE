"""Small network and differentiation utilities shared by neural methods."""

from __future__ import annotations

import copy
import math
from typing import Iterable

import numpy as np
import torch
from torch import nn


class DenseNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        *,
        width: int,
        depth: int,
        activation: str = "tanh",
        output_bound: float | None = None,
    ) -> None:
        super().__init__()
        activation_class = {"tanh": nn.Tanh, "silu": nn.SiLU, "relu": nn.ReLU}[activation]
        layers: list[nn.Module] = []
        previous = input_dim
        for _ in range(depth):
            linear = nn.Linear(previous, width)
            nn.init.xavier_uniform_(linear.weight)
            nn.init.zeros_(linear.bias)
            layers.extend((linear, activation_class()))
            previous = width
        final = nn.Linear(previous, output_dim)
        nn.init.xavier_uniform_(final.weight, gain=0.25)
        nn.init.zeros_(final.bias)
        layers.append(final)
        if output_bound is not None:
            layers.append(nn.Hardtanh(-output_bound, output_bound))
        self.model = nn.Sequential(*layers)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        return self.model(value)


def tx_input(t: torch.Tensor, x: torch.Tensor, T: float = 0.25) -> torch.Tensor:
    return torch.cat((t.reshape(-1, 1) / T, x), dim=-1)


def value_and_z(
    model: nn.Module,
    t: torch.Tensor,
    x: torch.Tensor,
    *,
    T: float = 0.25,
    create_graph: bool,
) -> tuple[torch.Tensor, torch.Tensor]:
    x_input = x if x.requires_grad else x.detach().requires_grad_(True)
    value = model(tx_input(t, x_input, T)).reshape(-1)
    gradient = torch.autograd.grad(
        value.sum(), x_input, create_graph=create_graph, retain_graph=create_graph
    )[0]
    return value, math.sqrt(2.0) * gradient


def count_parameters(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}


def batches(size: int, batch_size: int, generator: torch.Generator) -> Iterable[torch.Tensor]:
    order = torch.randperm(size, generator=generator)
    for start in range(0, size, batch_size):
        yield order[start : start + batch_size]
