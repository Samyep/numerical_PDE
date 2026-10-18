"""Dynamic-dimension definitions for the original and nonsmooth benchmarks."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Literal

import numpy as np
import torch


ActorKind = Literal["gradient", "value"]


class Benchmark:
    pde: str
    name: str
    d: int
    T: float
    sigma: float
    actor_kind: ActorKind
    input_scale: float = 1.0
    value_bounds: tuple[float, float] | None = None

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        x_var = x.detach().clone().requires_grad_(True)
        with torch.enable_grad():
            value = self.terminal(x_var)
            return torch.autograd.grad(value.sum(), x_var)[0]

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return torch.zeros_like(x)

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        raise NotImplementedError

    def tangent_coefficients(
        self, r: torch.Tensor | None, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        raise NotImplementedError

    def sample_states(
        self, count: int, *, seed: int, purpose: str = "function"
    ) -> tuple[np.ndarray, np.ndarray]:
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        x = rng.standard_normal((count, self.d))
        return t.astype(np.float64), x.astype(np.float64)


class HJBLogBenchmark(Benchmark):
    pde = "hjb_log"
    T = 1.0
    sigma = math.sqrt(2.0)
    actor_kind: ActorKind = "gradient"

    def __init__(self, d: int) -> None:
        self.d = int(d)
        self.name = f"hjb_log_d{d}"

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return torch.log((1.0 + x.square().sum(dim=-1)) / 2.0)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        return 2.0 * x / (1.0 + x.square().sum(dim=-1, keepdim=True))

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del r
        assert q is not None
        return -0.5 * q.square().sum(dim=-1)

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("HJB-Log needs a Gradient Actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        return c, -q, 0.5 * q.square().sum(dim=-1)

    def sample_states(self, count: int, *, seed: int, purpose: str = "function"):
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        x = rng.standard_normal((count, self.d))
        if purpose == "native":
            x *= np.sqrt(2.0 * t[:, None])
        elif purpose == "mixed":
            mask = np.arange(count) % 2 == 1
            x[mask] *= np.sqrt(2.0 * t[mask, None])
        elif purpose not in {"function", "official"}:
            raise KeyError(purpose)
        return t.astype(np.float64), x.astype(np.float64)


class AllenCahnBenchmark(Benchmark):
    pde = "allen_cahn"
    T = 0.3
    sigma = math.sqrt(2.0)
    actor_kind: ActorKind = "value"
    input_scale = 4.0
    value_bounds = (0.0, 1.0)

    def __init__(self, d: int) -> None:
        self.d = int(d)
        self.name = f"allen_cahn_d{d}"

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return 0.5 / (1.0 + 0.2 * x.square().sum(dim=-1))

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        denominator = 1.0 + 0.2 * x.square().sum(dim=-1, keepdim=True)
        return -0.2 * x / denominator.square()

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del q
        return r - r.pow(3)

    def tangent_coefficients(self, r, q):
        del q
        if r is None:
            raise ValueError("Allen-Cahn needs a Value Actor")
        c = 1.0 - 3.0 * r.square()
        b = torch.zeros((*r.shape, self.d), dtype=r.dtype, device=r.device)
        return c, b, 2.0 * r.pow(3)

    def sample_states(self, count: int, *, seed: int, purpose: str = "function"):
        if purpose not in {"function", "official", "native", "mixed"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        x = rng.standard_normal((count, self.d)) * np.sqrt(2.0 * t[:, None])
        return t.astype(np.float64), x.astype(np.float64)


class OUGMMBenchmark(Benchmark):
    pde = "ou_gmm"
    T = 1.0
    sigma = 1.0
    actor_kind: ActorKind = "gradient"
    input_scale = math.sqrt(5.0)
    initial_variance = 2.0
    evaluation_initial_variance = 4.0

    def __init__(self, d: int, mean: np.ndarray, pi: np.ndarray) -> None:
        self.d = int(d)
        self.name = f"ou_gmm_d{d}"
        self.mean = torch.as_tensor(np.asarray(mean), dtype=torch.float64)
        self.pi = torch.as_tensor(np.asarray(pi), dtype=torch.float64)
        if self.mean.shape != (5, self.d) or self.pi.shape != (5,):
            raise ValueError((self.mean.shape, self.pi.shape))
        if not torch.isclose(self.pi.sum(), torch.tensor(1.0, dtype=torch.float64)):
            raise ValueError("mixture weights must sum to one")

    @classmethod
    def from_npz(cls, path: Path) -> "OUGMMBenchmark":
        with np.load(path) as data:
            return cls(int(data["dimension"]), data["mean"], data["pi"])

    def _parameters(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.mean.to(x), self.pi.to(x)

    def negative_log_gmm(self, x: torch.Tensor, tau: torch.Tensor | float) -> torch.Tensor:
        mean, pi = self._parameters(x)
        local = torch.as_tensor(tau, dtype=x.dtype, device=x.device)
        while local.ndim < x.ndim:
            local = local.unsqueeze(-1)
        decay = torch.exp(-local)
        variance = self.initial_variance * decay.square() + 0.5 * (1.0 - decay.square())
        shifted = mean.unsqueeze(0) * decay.unsqueeze(-2)
        difference = x.unsqueeze(-2) - shifted
        logits = (
            torch.log(pi).unsqueeze(0)
            - 0.5 * self.d * torch.log(2.0 * torch.pi * variance)
            - 0.5 * difference.square().sum(dim=-1) / variance
        )
        return -torch.logsumexp(logits, dim=-1)

    def exact_value(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return self.negative_log_gmm(x, self.T - t)

    def exact_gradient(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        x_var = x.detach().clone().requires_grad_(True)
        with torch.enable_grad():
            value = self.exact_value(t, x_var)
            return torch.autograd.grad(value.sum(), x_var)[0]

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.negative_log_gmm(x, 0.0)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        t = torch.full(x.shape[:-1], self.T, dtype=x.dtype, device=x.device)
        return self.exact_gradient(t, x)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return x

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del r
        assert q is not None
        return -0.5 * q.square().sum(dim=-1) - self.d

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("OU/GMM needs a Gradient Actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        return c, -q, 0.5 * q.square().sum(dim=-1) - self.d

    def sample_states(self, count: int, *, seed: int, purpose: str = "official"):
        if purpose not in {"function", "official", "native", "mixed"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        x = rng.standard_normal((count, self.d)) * np.sqrt(4.0 + t[:, None])
        return t.astype(np.float64), x.astype(np.float64)


class MLPFailureBenchmark(Benchmark):
    pde = "mlp_failure"
    T = 1.0
    sigma = 1.0
    actor_kind: ActorKind = "gradient"

    def __init__(self, d: int, *, zero_tolerance: float = 1e-12) -> None:
        if d < 2:
            raise ValueError("inactive-gradient norm requires d>=2")
        self.d = int(d)
        self.name = f"mlp_failure_d{d}"
        self.zero_tolerance = float(zero_tolerance)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return x[..., 0].abs()

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        result = torch.zeros_like(x)
        result[..., 0] = torch.sign(x[..., 0])
        return result

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del r
        assert q is not None
        return torch.linalg.vector_norm(q[..., 1:], dim=-1)

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("MLP failure benchmark needs a Gradient Actor")
        inactive = q[..., 1:]
        norm = torch.linalg.vector_norm(inactive, dim=-1)
        direction = inactive / norm.clamp_min(self.zero_tolerance)[..., None]
        direction = torch.where((norm > self.zero_tolerance)[..., None], direction, 0.0)
        b = torch.zeros_like(q)
        b[..., 1:] = direction
        c = torch.zeros_like(norm)
        h = norm - (b * q).sum(dim=-1)
        return c, b, h

    def exact_value_numpy(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        times = np.asarray(t, dtype=np.float64).reshape(-1)
        mu = np.asarray(x, dtype=np.float64).reshape(len(times), self.d)[:, 0]
        variance = np.maximum(self.T - times, 0.0)
        sigma = np.sqrt(variance)
        result = np.abs(mu)
        active = sigma > 1e-14
        a = mu[active] / sigma[active]
        result[active] = (
            sigma[active] * math.sqrt(2.0 / math.pi) * np.exp(-0.5 * np.square(a))
            + mu[active] * np.vectorize(math.erf)(a / math.sqrt(2.0))
        )
        return result

    def exact_z_numpy(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        times = np.asarray(t, dtype=np.float64).reshape(-1)
        states = np.asarray(x, dtype=np.float64).reshape(len(times), self.d)
        sigma = np.sqrt(np.maximum(self.T - times, 0.0))
        result = np.zeros_like(states)
        active = sigma > 1e-14
        result[active, 0] = np.vectorize(math.erf)(
            states[active, 0] / (math.sqrt(2.0) * sigma[active])
        )
        result[~active, 0] = np.sign(states[~active, 0])
        return result


__all__ = [
    "AllenCahnBenchmark", "Benchmark", "HJBLogBenchmark", "MLPFailureBenchmark",
    "OUGMMBenchmark",
]
