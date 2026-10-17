"""Exact definitions of the three original-paper semilinear benchmarks.

The tangent solver uses ``z = sigma * grad(u)``.  For the DPI benchmark we
place the known linear term ``-F dot grad(u)`` in the forward generator and
leave the concave quadratic term in the nonlinear driver.  This is algebraic
bookkeeping only; :meth:`residual` checks the complete published PDE.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Literal

import numpy as np
import torch

from .common import PACKAGE_DIR


ActorKind = Literal["gradient", "value"]


class OriginalBenchmark:
    name: str
    d: int = 100
    T: float
    sigma: float
    actor_kind: ActorKind
    value_bounds: tuple[float, float] | None = None
    input_scale: float = 1.0

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return torch.zeros_like(x)

    def tangent_coefficients(
        self, r: torch.Tensor | None, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Return ``c,b,h`` for ``f ~= c*y + b dot z + h``."""
        raise NotImplementedError

    def sample_states(
        self, count: int, *, seed: int, purpose: str
    ) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError

    def metadata(self) -> dict[str, object]:
        return {
            "name": self.name,
            "dimension": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "z_convention": "z=sigma*grad(u)",
            "actor_kind": self.actor_kind,
            "value_bounds": self.value_bounds,
            "input_scale": self.input_scale,
        }


class HJBLogBenchmark(OriginalBenchmark):
    name = "hjb_log_d100"
    T = 1.0
    sigma = math.sqrt(2.0)
    actor_kind: ActorKind = "gradient"
    input_scale = 1.0

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return torch.log((1.0 + x.square().sum(dim=-1)) / 2.0)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        return 2.0 * x / (1.0 + x.square().sum(dim=-1, keepdim=True))

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("HJB-Log requires a Gradient Actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        b = -q
        h = 0.5 * q.square().sum(dim=-1)
        return c, b, h

    def sample_states(self, count: int, *, seed: int, purpose: str):
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        if purpose in {"function", "sdgd", "official"}:
            x = rng.standard_normal((count, self.d))
        elif purpose == "native":
            x = rng.standard_normal((count, self.d)) * np.sqrt(2.0 * t[:, None])
        elif purpose == "mixed":
            mask = np.arange(count) % 2 == 0
            x = rng.standard_normal((count, self.d))
            x[~mask] *= np.sqrt(2.0 * t[~mask, None])
        else:
            raise KeyError(purpose)
        return t.astype(np.float64), x.astype(np.float64)

    def metadata(self):
        return {
            **super().metadata(),
            "equation": "u_t + Delta u - ||grad u||^2 = 0",
            "driver": "f(z)=-||z||^2/2",
            "terminal": "log((1+||x||^2)/2)",
            "native_initial_state": "x0=0",
            "native_horizon": 1.0,
            "reference": "independent Hopf-Cole one-dimensional quadrature",
        }


class AllenCahnBenchmark(OriginalBenchmark):
    name = "allen_cahn_d100"
    T = 0.3
    sigma = math.sqrt(2.0)
    actor_kind: ActorKind = "value"
    value_bounds = (0.0, 1.0)
    input_scale = 4.0

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return 0.5 / (1.0 + 0.2 * x.square().sum(dim=-1))

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        denominator = 1.0 + 0.2 * x.square().sum(dim=-1, keepdim=True)
        return -0.2 * x / denominator.square()

    def tangent_coefficients(self, r, q):
        del q
        if r is None:
            raise ValueError("Allen-Cahn requires a Value Actor")
        c = 1.0 - 3.0 * r.square()
        b = torch.zeros((*r.shape, self.d), dtype=r.dtype, device=r.device)
        h = 2.0 * r.pow(3)
        return c, b, h

    def sample_states(self, count: int, *, seed: int, purpose: str):
        if purpose not in {"native", "mixed", "function", "official"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        # The original single-point task starts at zero.  These are exactly the
        # marginal states visited by its uncontrolled Brownian forward process.
        x = rng.standard_normal((count, self.d)) * np.sqrt(2.0 * t[:, None])
        return t.astype(np.float64), x.astype(np.float64)

    def metadata(self):
        return {
            **super().metadata(),
            "equation": "u_t + Delta u + u - u^3 = 0",
            "driver": "f(u)=u-u^3",
            "terminal": "1/(2+0.4||x||^2)",
            "native_initial_state": "x0=0",
            "native_horizon": 0.3,
            "tangent_domain_justification": "comparison principle preserves 0<=u<=1",
            "reference": "independent radial squared-radius PDE convergence study",
        }


class DPIOUBenchmark(OriginalBenchmark):
    name = "dpi_ou_gmm_d100"
    T = 1.0
    sigma = 1.0
    actor_kind: ActorKind = "gradient"
    input_scale = math.sqrt(5.0)
    theta = 1.0
    mu = 0.0
    alpha = 1.0
    components = 5
    initial_variance = 2.0
    evaluation_initial_variance = 4.0

    def __init__(self, input_dir: Path | None = None) -> None:
        directory = input_dir or PACKAGE_DIR / "official_inputs"
        mean_path = directory / "mean_100d_ms=1.0_vs=2.0_5.pt"
        pi_path = directory / "pi_100d_ms=1.0_vs=2.0_5.pt"
        if not mean_path.exists() or not pi_path.exists():
            raise FileNotFoundError(
                "official DPI GMM tensors are absent; run prepare_official_inputs.py"
            )
        self.mean = torch.load(mean_path, map_location="cpu", weights_only=False).double()
        self.pi = torch.load(pi_path, map_location="cpu", weights_only=False).double()
        if self.mean.shape != (self.components, self.d) or self.pi.shape != (self.components,):
            raise ValueError("unexpected official DPI GMM shape")
        if not torch.isclose(self.pi.sum(), torch.tensor(1.0, dtype=torch.float64)):
            raise ValueError("DPI mixture weights do not sum to one")

    def _parameters(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.mean.to(dtype=x.dtype, device=x.device), self.pi.to(
            dtype=x.dtype, device=x.device
        )

    def negative_log_gmm(self, x: torch.Tensor, tau: torch.Tensor | float) -> torch.Tensor:
        means, weights = self._parameters(x)
        local_tau = torch.as_tensor(tau, dtype=x.dtype, device=x.device)
        while local_tau.ndim < x.ndim:
            local_tau = local_tau.unsqueeze(-1)
        decay = torch.exp(-self.theta * local_tau)
        variance = self.initial_variance * decay.square() + self.alpha / (
            2.0 * self.theta
        ) * (1.0 - decay.square())
        component_means = means.unsqueeze(0) * decay.unsqueeze(-2)
        difference = x.unsqueeze(-2) - component_means
        log_component = (
            torch.log(weights).unsqueeze(0)
            - 0.5 * self.d * torch.log(2.0 * torch.pi * variance)
            - 0.5 * difference.square().sum(dim=-1) / variance
        )
        return -torch.logsumexp(log_component, dim=-1)

    def exact_value(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return self.negative_log_gmm(x, self.T - t)

    def exact_gradient(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        means, weights = self._parameters(x)
        tau = self.T - t
        while tau.ndim < x.ndim:
            tau = tau.unsqueeze(-1)
        decay = torch.exp(-self.theta * tau)
        variance = self.initial_variance * decay.square() + self.alpha / (
            2.0 * self.theta
        ) * (1.0 - decay.square())
        component_means = means.unsqueeze(0) * decay.unsqueeze(-2)
        difference = x.unsqueeze(-2) - component_means
        logits = (
            torch.log(weights).unsqueeze(0)
            - 0.5 * self.d * torch.log(2.0 * torch.pi * variance)
            - 0.5 * difference.square().sum(dim=-1) / variance
        )
        responsibilities = torch.softmax(logits, dim=-1)
        return (
            responsibilities.unsqueeze(-1) * difference / variance.unsqueeze(-1)
        ).sum(dim=-2)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.negative_log_gmm(x, 0.0)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        t = torch.full(x.shape[:-1], self.T, dtype=x.dtype, device=x.device)
        return self.exact_gradient(t, x)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        # The complete PDE contains -F dot grad(u), F=theta(mu-x).
        return self.theta * (x - self.mu)

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("DPI OU/GMM requires a Gradient Actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        b = -q
        h = 0.5 * q.square().sum(dim=-1) - self.d * self.theta
        return c, b, h

    def sample_states(self, count: int, *, seed: int, purpose: str):
        if purpose not in {"official", "function", "mixed", "native"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), count)
        # Exact distribution in the authors' Equation.sample_x(t): X0~N(0,4I)
        # followed by the Brownian sampling generator with alpha=1.
        x = rng.standard_normal((count, self.d)) * np.sqrt(
            self.evaluation_initial_variance + self.alpha * t[:, None]
        )
        return t.astype(np.float64), x.astype(np.float64)

    def metadata(self):
        return {
            **super().metadata(),
            "equation": (
                "u_t + alpha/2 Delta u - theta(mu-x).grad u "
                "- alpha/2||grad u||^2 - d theta = 0"
            ),
            "driver_z_form": "f(z,x)=-F.grad(u)-||z||^2/2-d theta",
            "base_drift_used_by_tangent_solver": "-F=theta(x-mu)",
            "terminal": "negative log density of official five-component GMM",
            "exact_solution": "negative log density of GMM evolved by forward OU process",
            "evaluation_distribution": "t~U[0,1], x|t~N(0,(4+t)I)",
            "official_missing_variance_recovery": "five covariance matrices fixed to 2I per source constructor",
        }


BENCHMARK_NAMES = ("hjb_log_d100", "allen_cahn_d100", "dpi_ou_gmm_d100")


def make_benchmark(name: str) -> OriginalBenchmark:
    if name == "hjb_log_d100":
        return HJBLogBenchmark()
    if name == "allen_cahn_d100":
        return AllenCahnBenchmark()
    if name == "dpi_ou_gmm_d100":
        return DPIOUBenchmark()
    raise KeyError(name)


__all__ = [
    "ActorKind",
    "AllenCahnBenchmark",
    "BENCHMARK_NAMES",
    "DPIOUBenchmark",
    "HJBLogBenchmark",
    "OriginalBenchmark",
    "make_benchmark",
]
