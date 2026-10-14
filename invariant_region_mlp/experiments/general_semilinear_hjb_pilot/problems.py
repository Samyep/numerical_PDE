"""Controlled HJB definitions used by the general semilinear pilot."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np
import torch


def quartic_root_numpy(p: np.ndarray, lam: float = 1.0, iterations: int = 60) -> np.ndarray:
    """Solve a + lam*a^3 = -p using monotone vectorized bisection."""
    p = np.asarray(p, dtype=np.float64)
    radius = np.abs(p) + 1.0
    lo, hi = -radius, radius
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        residual = mid + lam * mid**3 + p
        hi = np.where(residual > 0.0, mid, hi)
        lo = np.where(residual <= 0.0, mid, lo)
    return 0.5 * (lo + hi)


def quartic_root_torch(p: torch.Tensor, lam: float = 1.0, iterations: int = 36) -> torch.Tensor:
    """Torch counterpart of :func:`quartic_root_numpy` for detached labels."""
    radius = p.abs() + 1.0
    lo, hi = -radius, radius
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        residual = mid + lam * mid.pow(3) + p
        hi = torch.where(residual > 0.0, mid, hi)
        lo = torch.where(residual <= 0.0, mid, lo)
    return 0.5 * (lo + hi)


class ControlledHJB(ABC):
    """One interface for state dimension d and independent action dimension m."""

    name: str
    d: int
    m: int
    T: float
    sigma: float = math.sqrt(2.0)

    @abstractmethod
    def terminal(self, x: torch.Tensor) -> torch.Tensor: ...

    @abstractmethod
    def terminal_numpy(self, x: np.ndarray) -> np.ndarray: ...

    @abstractmethod
    def drift(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor: ...

    def diffusion(self, t: torch.Tensor, x: torch.Tensor) -> float:
        del t, x
        return self.sigma

    @abstractmethod
    def running_cost(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor: ...

    @abstractmethod
    def greedy_action(
        self, t: torch.Tensor, x: torch.Tensor, value_gradient: torch.Tensor
    ) -> torch.Tensor: ...

    @abstractmethod
    def feasible_action(self, raw_actor_output: torch.Tensor) -> torch.Tensor: ...

    @abstractmethod
    def hamiltonian_numpy(self, action: np.ndarray, projected_gradient: np.ndarray) -> np.ndarray: ...

    def metadata(self) -> dict[str, Any]:
        payload = dict(self.__dict__)
        if "beta" in payload:
            beta = np.asarray(payload.pop("beta"))
            payload["beta_norm"] = float(np.linalg.norm(beta))
            payload["beta_seed"] = getattr(self, "beta_seed", None)
        payload.update(name=self.name, d=self.d, m=self.m, T=self.T, sigma=self.sigma)
        return payload


@dataclass
class LQGProblem(ControlledHJB):
    d: int = 100
    T: float = 0.5
    gamma: float = 0.2
    Q: float = 0.5
    R: float = 1.0
    G: float = 1.0
    name: str = "lqg"

    def __post_init__(self) -> None:
        self.m = self.d

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return 0.5 * self.G * x.square().sum(dim=-1)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        return 0.5 * self.G * np.square(x).sum(axis=-1)

    def drift(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        del t
        return -self.gamma * x + a

    def running_cost(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        del t
        return 0.5 * self.Q * x.square().sum(dim=-1) + 0.5 * self.R * a.square().sum(dim=-1)

    def greedy_action(
        self, t: torch.Tensor, x: torch.Tensor, value_gradient: torch.Tensor
    ) -> torch.Tensor:
        del t, x
        return -value_gradient / self.R

    def feasible_action(self, raw_actor_output: torch.Tensor) -> torch.Tensor:
        return raw_actor_output

    def hamiltonian_numpy(self, action: np.ndarray, projected_gradient: np.ndarray) -> np.ndarray:
        return np.sum(action * projected_gradient, axis=-1) + 0.5 * self.R * np.sum(action**2, axis=-1)


class DirectionalProblem(ControlledHJB):
    """Base for scalar controls acting along a fixed dense physical direction."""

    def __init__(self, *, d: int, T: float, beta_seed: int, name: str) -> None:
        self.d, self.m, self.T, self.name = d, 1, T, name
        self.beta_seed = beta_seed
        rng = np.random.default_rng(beta_seed)
        beta = rng.standard_normal(d)
        self.beta = beta / np.linalg.norm(beta)

    def _beta_torch(self, x: torch.Tensor) -> torch.Tensor:
        return torch.as_tensor(self.beta, dtype=x.dtype, device=x.device)

    def projection(self, x: torch.Tensor) -> torch.Tensor:
        return x @ self._beta_torch(x)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        s = self.projection(x)
        return 1.5 * torch.logaddexp(s, -s) - 1.5 * math.log(2.0) + 0.2 * torch.sin(2.0 * s)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        s = np.asarray(x) @ self.beta
        return 1.5 * np.logaddexp(s, -s) - 1.5 * math.log(2.0) + 0.2 * np.sin(2.0 * s)

    def terminal_1d(self, s: np.ndarray) -> np.ndarray:
        return 1.5 * np.logaddexp(s, -s) - 1.5 * math.log(2.0) + 0.2 * np.sin(2.0 * s)

    def drift(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        del t
        return a * self._beta_torch(x)


class QuarticControlProblem(DirectionalProblem):
    def __init__(self, *, d: int = 100, T: float = 0.5, lam: float = 1.0, beta_seed: int = 9137) -> None:
        super().__init__(d=d, T=T, beta_seed=beta_seed, name="quartic")
        self.lam = lam

    def running_cost(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        del t, x
        scalar = a[..., 0]
        return 0.5 * scalar.square() + 0.25 * self.lam * scalar.pow(4)

    def greedy_action(
        self, t: torch.Tensor, x: torch.Tensor, value_gradient: torch.Tensor
    ) -> torch.Tensor:
        del t, x
        p = value_gradient @ self._beta_torch(value_gradient)
        return quartic_root_torch(p, self.lam).unsqueeze(-1)

    def feasible_action(self, raw_actor_output: torch.Tensor) -> torch.Tensor:
        return raw_actor_output

    def greedy_numpy(self, p: np.ndarray) -> np.ndarray:
        return quartic_root_numpy(p, self.lam)

    def hamiltonian_numpy(self, action: np.ndarray, projected_gradient: np.ndarray) -> np.ndarray:
        a = np.asarray(action)[..., 0] if np.asarray(action).ndim > np.asarray(projected_gradient).ndim else np.asarray(action)
        return a * projected_gradient + 0.5 * a**2 + 0.25 * self.lam * a**4


class BoundedControlProblem(DirectionalProblem):
    def __init__(self, *, d: int = 100, T: float = 0.5, A: float = 0.35, beta_seed: int = 9137) -> None:
        super().__init__(d=d, T=T, beta_seed=beta_seed, name="bounded")
        self.A = A

    def running_cost(self, t: torch.Tensor, x: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        del t, x
        return 0.5 * a[..., 0].square()

    def greedy_action(
        self, t: torch.Tensor, x: torch.Tensor, value_gradient: torch.Tensor
    ) -> torch.Tensor:
        del t, x
        p = value_gradient @ self._beta_torch(value_gradient)
        return torch.clamp(-p, -self.A, self.A).unsqueeze(-1)

    def feasible_action(self, raw_actor_output: torch.Tensor) -> torch.Tensor:
        return self.A * torch.tanh(raw_actor_output)

    def greedy_numpy(self, p: np.ndarray) -> np.ndarray:
        return np.clip(-np.asarray(p), -self.A, self.A)

    def hamiltonian_numpy(self, action: np.ndarray, projected_gradient: np.ndarray) -> np.ndarray:
        a = np.asarray(action)[..., 0] if np.asarray(action).ndim > np.asarray(projected_gradient).ndim else np.asarray(action)
        return a * projected_gradient + 0.5 * a**2


def make_problem(name: str, *, d: int = 100) -> ControlledHJB:
    if name == "lqg":
        return LQGProblem(d=d)
    if name == "quartic":
        return QuarticControlProblem(d=d)
    if name == "bounded":
        return BoundedControlProblem(d=d)
    raise KeyError(name)

