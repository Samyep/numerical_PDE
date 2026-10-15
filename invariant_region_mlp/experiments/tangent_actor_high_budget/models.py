"""Original-form competitor models without Tangent-Actor terminal features."""

from __future__ import annotations

import math

import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)


class RawTanhPINN(nn.Module):
    """The ordinary tanh PINN used by SDGD: ``(t,x) -> u`` directly."""

    def __init__(self, d: int, *, width: int = 192, depth: int = 3) -> None:
        super().__init__()
        self.d, self.width, self.depth = d, width, depth
        layers: list[nn.Module] = []
        previous = d + 1
        for _ in range(depth):
            linear = nn.Linear(previous, width)
            nn.init.xavier_uniform_(linear.weight)
            nn.init.zeros_(linear.bias)
            layers.extend((linear, nn.Tanh()))
            previous = width
        final = nn.Linear(previous, 1)
        nn.init.xavier_uniform_(final.weight, gain=0.1)
        nn.init.constant_(final.bias, 0.5)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def value(
        self, driver: ConcaveSemilinearDriver, t: torch.Tensor, x: torch.Tensor
    ) -> torch.Tensor:
        inputs = torch.cat((t.reshape(-1, 1) / driver.T, x), dim=1)
        return self.network(inputs).reshape(-1)


class OfficialStyleZNet(nn.Module):
    """Deep-BSDE's nonshared BN/ReLU Z subnet, adapted to PyTorch."""

    def __init__(self, d: int, width: int) -> None:
        super().__init__()
        self.input_bn = nn.BatchNorm1d(d, eps=1e-6, momentum=0.01)
        self.linear1 = nn.Linear(d, width, bias=False)
        self.bn1 = nn.BatchNorm1d(width, eps=1e-6, momentum=0.01)
        self.linear2 = nn.Linear(width, width, bias=False)
        self.bn2 = nn.BatchNorm1d(width, eps=1e-6, momentum=0.01)
        self.output = nn.Linear(width, d)
        self.output_bn = nn.BatchNorm1d(d, eps=1e-6, momentum=0.01)
        for layer in (self.linear1, self.linear2, self.output):
            nn.init.xavier_uniform_(layer.weight)
            if layer.bias is not None:
                nn.init.zeros_(layer.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = torch.relu(self.bn1(self.linear1(self.input_bn(x))))
        h = torch.relu(self.bn2(self.linear2(h)))
        return self.output_bn(self.output(h))


class RawInitialValue(nn.Module):
    """Amortization-only adaptation of scalar Deep-BSDE Y0."""

    def __init__(self, d: int, width: int) -> None:
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(d, width), nn.Tanh(), nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 1),
        )
        for layer in self.network:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_uniform_(layer.weight)
                nn.init.zeros_(layer.bias)
        nn.init.constant_(self.network[-1].bias, 0.5)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).reshape(-1)


class OriginalAmortizedDeepBSDE(nn.Module):
    """Raw Y0 and nonshared Z_n networks, with no terminal-aware ansatz."""

    def __init__(self, d: int, *, time_steps: int = 16, width: int | None = None) -> None:
        super().__init__()
        self.d = d
        self.time_steps = time_steps
        self.width = int(width if width is not None else min(d + 10, 256))
        self.initial_net = RawInitialValue(d, self.width)
        self.z_nets = nn.ModuleList(
            OfficialStyleZNet(d, self.width) for _ in range(time_steps)
        )

    def initial_value(
        self, driver: ConcaveSemilinearDriver, x: torch.Tensor
    ) -> torch.Tensor:
        del driver
        return self.initial_net(x)

    def rollout(
        self,
        driver: ConcaveSemilinearDriver,
        x0: torch.Tensor,
        normals: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, int]:
        if normals.shape[1] != self.time_steps:
            raise ValueError("time-step mismatch")
        dt = driver.T / self.time_steps
        root_dt = math.sqrt(dt)
        x = x0
        y = self.initial_value(driver, x0)
        nonfinite = 0
        for step, z_net in enumerate(self.z_nets):
            # The official implementation divides each subnet output by d.
            z = z_net(x) / driver.d
            dw = root_dt * normals[:, step, :]
            y = y - driver.driver(y, z) * dt + torch.sum(z * dw, dim=-1)
            x = x + driver.sigma * dw
            nonfinite += int((~torch.isfinite(y)).sum().item())
        return y, x, nonfinite


def architecture_payload(model: nn.Module) -> dict[str, int | str | None]:
    if isinstance(model, RawTanhPINN):
        return {"kind": "raw_tanh_pinn", "d": model.d, "width": model.width, "depth": model.depth}
    if isinstance(model, OriginalAmortizedDeepBSDE):
        return {
            "kind": "original_amortized_deep_bsde", "d": model.d,
            "time_steps": model.time_steps, "width": model.width,
        }
    raise TypeError(type(model))


__all__ = [
    "OriginalAmortizedDeepBSDE", "RawTanhPINN", "architecture_payload",
]
