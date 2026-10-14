"""Tangent Actors with a bounded Value head and optional full Gradient head."""

from __future__ import annotations

import torch
from torch import nn

from .drivers import ConcaveSemilinearDriver


class TangentActor(nn.Module):
    def __init__(self, driver: ConcaveSemilinearDriver, *, width: int, depth: int) -> None:
        super().__init__()
        self.driver = driver
        self.width, self.depth = width, depth
        output_dim = 1 + (driver.d if driver.has_gradient else 0)
        layers: list[nn.Module] = []
        previous = driver.d + 1
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
        self.network = nn.Sequential(*layers)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        inputs = torch.cat((t.reshape(-1, 1) / T, x), dim=-1)
        output = self.network(inputs)
        r = self.driver.constrain_r(output[:, 0])
        q = output[:, 1:] if self.driver.has_gradient else None
        return r, q


class ConstantTangentActor(nn.Module):
    def __init__(self, driver: ConcaveSemilinearDriver) -> None:
        super().__init__()
        self.driver = driver
        self.initial_r = 0.5 * (driver.r_low + driver.r_high)

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        del t, T
        r = torch.full(x.shape[:-1], self.initial_r, dtype=x.dtype, device=x.device)
        q = torch.zeros_like(x) if self.driver.has_gradient else None
        return r, q


def count_parameters(module: nn.Module) -> int:
    return sum(parameter.numel() for parameter in module.parameters())
