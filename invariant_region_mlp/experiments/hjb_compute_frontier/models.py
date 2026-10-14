"""Networks shared by the four methods."""

from __future__ import annotations

import math

import torch
from torch import nn


def make_activation(name: str) -> nn.Module:
    return {
        "tanh": nn.Tanh,
        "silu": nn.SiLU,
        "elu": nn.ELU,
        "relu": nn.ReLU,
    }[name]()


class DenseNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        *,
        width: int,
        depth: int,
        activation: str,
        final_gain: float = 0.25,
    ) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        previous = input_dim
        for _ in range(depth):
            layer = nn.Linear(previous, width)
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)
            layers.extend((layer, make_activation(activation)))
            previous = width
        final = nn.Linear(previous, output_dim)
        nn.init.xavier_uniform_(final.weight, gain=final_gain)
        nn.init.zeros_(final.bias)
        layers.append(final)
        self.model = nn.Sequential(*layers)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        return self.model(value)


def tx_input(t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
    return torch.cat((t.reshape(-1, 1) / T, x), dim=-1)


class ActorNet(nn.Module):
    def __init__(self, d: int, *, width: int, depth: int) -> None:
        super().__init__()
        self.d = d
        self.width = width
        self.depth = depth
        self.network = DenseNet(
            d + 1, d, width=width, depth=depth, activation="silu"
        )

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        return self.network(tx_input(t, x, T))


class ScalarSpaceTimeNet(nn.Module):
    def __init__(self, d: int, *, width: int, depth: int, activation: str) -> None:
        super().__init__()
        self.d = d
        self.network = DenseNet(
            d + 1, 1, width=width, depth=depth, activation=activation
        )

    def value(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        return self.network(tx_input(t, x, T)).reshape(-1)


def value_and_z(
    model: ScalarSpaceTimeNet,
    t: torch.Tensor,
    x: torch.Tensor,
    T: float,
    *,
    create_graph: bool,
) -> tuple[torch.Tensor, torch.Tensor]:
    state = x if x.requires_grad else x.detach().requires_grad_(True)
    value = model.value(t, state, T)
    gradient = torch.autograd.grad(
        value.sum(), state, create_graph=create_graph, retain_graph=create_graph
    )[0]
    return value, math.sqrt(2.0) * gradient


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}


def count_parameters(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters())
