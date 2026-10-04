"""Same-data 1D Euler adaptation of the published RoeNet architecture.

The two learned maps and the wave-splitting formula follow the authors'
``nontrivial_3c/train_net.py`` at commit
``ef877957c1c0ddb16eac17006d75bdf7bd786d45``.  This module is intentionally
labelled an adaptation: the grid, training distribution, periodic boundary,
normalization, and validation-based stopping rule are matched to the HCFL
benchmark rather than to the paper's single Sod trajectory experiment.
"""

from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as functional


class ResidualPointBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.first = nn.Conv1d(in_channels, out_channels, 1)
        self.second = nn.Conv1d(out_channels, out_channels, 1)
        self.shortcut = (
            nn.Identity()
            if in_channels == out_channels
            else nn.Conv1d(in_channels, out_channels, 1)
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        residual = self.second(functional.relu(self.first(values)))
        return functional.relu(residual + self.shortcut(values))


def point_network(in_channels: int, out_channels: int) -> nn.Sequential:
    return nn.Sequential(
        ResidualPointBlock(in_channels, 16),
        ResidualPointBlock(16, 16),
        ResidualPointBlock(16, 32),
        ResidualPointBlock(32, 64),
        ResidualPointBlock(64, 64),
        ResidualPointBlock(64, 64),
        nn.Conv1d(64, out_channels, 1),
    )


class RoeNetEuler1d(nn.Module):
    """RoeNet wave splitting with a learned 64-dimensional hidden basis."""

    def __init__(
        self,
        state_mean: torch.Tensor,
        state_std: torch.Tensor,
        cells: int = 64,
        hidden_waves: int = 64,
        internal_steps: int = 1,
        regularization: float = 1.0e-5,
    ) -> None:
        super().__init__()
        self.channels = 3
        self.cells = cells
        self.hidden_waves = hidden_waves
        self.internal_steps = internal_steps
        self.regularization = regularization
        self.lambda_net = point_network(2 * self.channels, hidden_waves)
        self.left_net = point_network(
            2 * self.channels, hidden_waves * self.channels
        )
        self.register_buffer(
            "state_mean", state_mean.reshape(1, 1, self.channels).float()
        )
        self.register_buffer(
            "state_std", state_std.reshape(1, 1, self.channels).float()
        )

    def _decomposition(
        self, first: torch.Tensor, second: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        normalized_first = (first - self.state_mean) / self.state_std
        normalized_second = (second - self.state_mean) / self.state_std
        features = torch.cat((normalized_first, normalized_second), dim=-1)
        channels_first = features.transpose(1, 2)
        eigenvalues = self.lambda_net(channels_first).transpose(1, 2) / 10.0
        left = self.left_net(channels_first).transpose(1, 2)
        left = left.reshape(
            left.shape[0], left.shape[1], self.hidden_waves, self.channels
        )

        # The official code forms (L^T L)^(-1)L^T.  A small Tikhonov term
        # keeps the same least-squares inverse well-defined during early
        # training, where a random learned L can be nearly rank deficient.
        transpose = left.transpose(-2, -1)
        gram = transpose @ left
        identity = torch.eye(
            self.channels, dtype=left.dtype, device=left.device
        ).reshape(1, 1, self.channels, self.channels)
        right = torch.linalg.solve(
            gram + self.regularization * identity, transpose
        )
        return eigenvalues, (left, right)

    @staticmethod
    def _apply_split(
        eigenvalues: torch.Tensor,
        left: torch.Tensor,
        right: torch.Tensor,
        jump: torch.Tensor,
        sign: int,
    ) -> torch.Tensor:
        characteristic = torch.einsum("bnhc,bnc->bnh", left, jump)
        if sign < 0:
            speeds = eigenvalues - eigenvalues.abs()
        else:
            speeds = eigenvalues + eigenvalues.abs()
        return torch.einsum("bnch,bnh->bnc", right, speeds * characteristic)

    def spatial_derivative(self, state: torch.Tensor) -> torch.Tensor:
        left_state = torch.roll(state, 1, dims=-2)
        right_state = torch.roll(state, -1, dims=-2)
        lambda_right, (left_right, right_right) = self._decomposition(
            state, right_state
        )
        lambda_left, (left_left, right_left) = self._decomposition(
            left_state, state
        )
        right_going = self._apply_split(
            lambda_right,
            left_right,
            right_right,
            right_state - state,
            -1,
        )
        left_going = self._apply_split(
            lambda_left,
            left_left,
            right_left,
            state - left_state,
            +1,
        )
        dx = 1.0 / self.cells
        return -(right_going + left_going) / (2.0 * dx)

    def forward(self, state: torch.Tensor, dt: float = 4.0e-4) -> torch.Tensor:
        step = dt / self.internal_steps
        output = state
        for _ in range(self.internal_steps):
            output = output + step * self.spatial_derivative(output)
        return output


