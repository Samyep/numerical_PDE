"""Terminal-aware neural models used by external baseline adaptations."""

from __future__ import annotations

import math

import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)


def _activation(name: str) -> nn.Module:
    if name == "silu":
        return nn.SiLU()
    if name == "tanh":
        return nn.Tanh()
    if name == "elu":
        return nn.ELU()
    raise KeyError(name)


class DenseNet(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        *,
        width: int,
        depth: int,
        activation: str = "silu",
        zero_final: bool = False,
    ) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        previous = input_dim
        for _ in range(depth):
            layer = nn.Linear(previous, width)
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)
            layers.extend((layer, _activation(activation)))
            previous = width
        final = nn.Linear(previous, output_dim)
        if zero_final:
            nn.init.zeros_(final.weight)
            nn.init.zeros_(final.bias)
        else:
            nn.init.xavier_uniform_(final.weight, gain=0.1)
            nn.init.zeros_(final.bias)
        layers.append(final)
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.network(inputs)


def tx_input(t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
    return torch.cat((t.reshape(-1, 1) / T, x), dim=-1)


class HardTerminalValue(nn.Module):
    """PINN ansatz that enforces u(T,x)=g(x) exactly."""

    def __init__(
        self,
        d: int,
        *,
        width: int = 192,
        depth: int = 3,
        activation: str = "silu",
    ) -> None:
        super().__init__()
        self.d, self.width, self.depth, self.activation = d, width, depth, activation
        self.correction = DenseNet(
            d + 1, 1, width=width, depth=depth, activation=activation
        )

    def value(
        self,
        driver: ConcaveSemilinearDriver,
        t: torch.Tensor,
        x: torch.Tensor,
    ) -> torch.Tensor:
        correction = self.correction(tx_input(t, x, driver.T)).reshape(-1)
        return driver.terminal(x) + (driver.T - t.reshape(-1)) * correction


class AmortizedDeepBSDE(nn.Module):
    """Multiple-initial-state Deep BSDE with shared terminal-aware Z network.

    The original algorithm learns one scalar Y0.  Here Y0(x0) is a network so
    the comparison predicts the same 960 held-out states in one amortized fit.
    """

    def __init__(
        self,
        d: int,
        *,
        width: int = 192,
        depth: int = 3,
        activation: str = "silu",
    ) -> None:
        super().__init__()
        self.d, self.width, self.depth, self.activation = d, width, depth, activation
        self.initial_correction = DenseNet(
            d, 1, width=width, depth=depth, activation=activation, zero_final=True
        )
        self.z_correction = DenseNet(
            d + 1,
            d,
            width=width,
            depth=depth,
            activation=activation,
            zero_final=True,
        )

    def initial_value(
        self, driver: ConcaveSemilinearDriver, x: torch.Tensor
    ) -> torch.Tensor:
        return driver.terminal(x) + driver.T * self.initial_correction(x).reshape(-1)

    def z_value(
        self,
        driver: ConcaveSemilinearDriver,
        t: torch.Tensor,
        x: torch.Tensor,
    ) -> torch.Tensor:
        terminal_z = driver.sigma * driver.terminal_gradient_torch(x)
        tau = (driver.T - t.reshape(-1)).reshape(-1, 1)
        return terminal_z + tau * self.z_correction(tx_input(t, x, driver.T))

    def rollout(
        self,
        driver: ConcaveSemilinearDriver,
        x0: torch.Tensor,
        normals: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, int]:
        steps = normals.shape[1]
        dt = driver.T / steps
        root_dt = math.sqrt(dt)
        x = x0
        y = self.initial_value(driver, x0)
        nonfinite = 0
        for step in range(steps):
            t = torch.full(
                (len(x0),), step * dt, dtype=x0.dtype, device=x0.device
            )
            z = self.z_value(driver, t, x)
            dw = root_dt * normals[:, step, :]
            # PDE: u_t + Delta u + f(u,z)=0, z=sqrt(2) grad u.
            # Therefore dY=-f(Y,Z)dt+Z dW along dX=sqrt(2)dW.
            y = y - driver.driver(y, z) * dt + torch.sum(z * dw, dim=-1)
            x = x + driver.sigma * dw
            nonfinite += int((~torch.isfinite(y)).sum().item())
        return y, x, nonfinite


def architecture_payload(model: nn.Module) -> dict[str, int | str]:
    return {
        "d": int(getattr(model, "d")),
        "width": int(getattr(model, "width")),
        "depth": int(getattr(model, "depth")),
        "activation": str(getattr(model, "activation")),
    }


__all__ = [
    "AmortizedDeepBSDE", "DenseNet", "HardTerminalValue", "architecture_payload",
    "tx_input",
]
