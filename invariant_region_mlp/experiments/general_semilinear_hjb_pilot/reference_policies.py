"""Post-training analytic/numerical reference policies for evaluation only."""

from __future__ import annotations

import numpy as np
import torch
from torch import nn

from .problems import DirectionalProblem, LQGProblem
from .reference import RiccatiReference, ScalarReference


class LQGReferencePolicy(nn.Module):
    def __init__(self, problem: LQGProblem, reference: RiccatiReference) -> None:
        super().__init__()
        self.problem = problem
        grid = np.linspace(0.0, problem.T, 4097)
        P, _ = reference.coefficients(grid)
        self.register_buffer("P_grid", torch.as_tensor(P, dtype=torch.float32))

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        position = torch.clamp(t / T * (len(self.P_grid) - 1), 0.0, len(self.P_grid) - 1.0)
        left = torch.floor(position).long().clamp(max=len(self.P_grid) - 2)
        weight = position - left
        P = self.P_grid[left] * (1.0 - weight) + self.P_grid[left + 1] * weight
        return -(P / self.problem.R)[:, None] * x


class ScalarReferencePolicy(nn.Module):
    def __init__(
        self, problem: DirectionalProblem, reference: ScalarReference
    ) -> None:
        super().__init__()
        self.problem = problem
        self.t0 = float(reference.times[0])
        self.dt = float(reference.times[1] - reference.times[0])
        self.x0 = float(reference.space[0])
        self.dx = float(reference.space[1] - reference.space[0])
        self.register_buffer("actions", torch.as_tensor(reference.actions, dtype=torch.float32))
        self.register_buffer("beta", torch.as_tensor(problem.beta, dtype=torch.float32))

    def predict(self, t: torch.Tensor, x: torch.Tensor, T: float) -> torch.Tensor:
        del T
        s = x @ self.beta.to(dtype=x.dtype)
        ft = torch.clamp((t - self.t0) / self.dt, 0.0, self.actions.shape[0] - 1.0)
        fs = torch.clamp((s - self.x0) / self.dx, 0.0, self.actions.shape[1] - 1.0)
        it = torch.floor(ft).long().clamp(max=self.actions.shape[0] - 2)
        ix = torch.floor(fs).long().clamp(max=self.actions.shape[1] - 2)
        wt, wx = ft - it, fs - ix
        table = self.actions
        a0 = table[it, ix] * (1.0 - wx) + table[it, ix + 1] * wx
        a1 = table[it + 1, ix] * (1.0 - wx) + table[it + 1, ix + 1] * wx
        return (a0 * (1.0 - wt) + a1 * wt).unsqueeze(-1)

