"""Dense-Fourier joint PDE and recreated published coupled-LQG benchmark."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from typing import Any

import numpy as np
import torch
from scipy.special import roots_laguerre

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)


SQRT2 = math.sqrt(2.0)


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(str(value.shape).encode("ascii"))
        digest.update(value.tobytes())
    return digest.hexdigest()


class DenseFourierDriver(ConcaveSemilinearDriver):
    """Joint f(u,z)=0.5*tanh(u)*|z|^2 with an analytic sinh transform."""

    family = "dense_fourier_joint"
    has_gradient = True
    sigma = SQRT2
    reconstruction_sign = 1.0
    # These are metadata-only guards: the preregistered actor below is
    # intentionally unbounded and does not project or sigmoid-transform r.
    r_low, r_high = -20.0, 20.0

    def __init__(
        self,
        omega: np.ndarray,
        signs: np.ndarray,
        phases: np.ndarray,
        *,
        T: float = 0.5,
        parameter_seed: tuple[int, ...] = (),
    ) -> None:
        omega = np.asarray(omega, dtype=np.float64)
        signs = np.asarray(signs, dtype=np.float64)
        phases = np.asarray(phases, dtype=np.float64)
        if omega.ndim != 2 or len(signs) != len(omega) or len(phases) != len(omega):
            raise ValueError("incompatible Dense Fourier parameters")
        self.k, self.d = omega.shape
        if self.k != 2 * self.d:
            raise ValueError("the preregistered benchmark requires k=2d")
        norms = np.linalg.norm(omega, axis=1)
        if not np.allclose(norms, 1.0, rtol=0.0, atol=2e-14):
            raise ValueError("omega rows must be unit vectors")
        if not np.all(np.isin(signs, (-1.0, 1.0))):
            raise ValueError("a_j must be signs")
        self.T = float(T)
        self.name = f"dense_fourier_joint_d{self.d}_k{self.k}"
        self.omega = omega.copy()
        self.signs = signs.copy()
        self.phases = phases.copy()
        self.omega_norm_squared = np.sum(self.omega * self.omega, axis=1)
        self.scale = 2.0 / math.sqrt(2.0 * self.d)
        self.parameter_seed = tuple(int(value) for value in parameter_seed)
        self.basis_seed = -1
        # Some inherited diagnostics expect B.  It is the complete dense
        # direction matrix, not a low-rank oracle feature exposed to the actor.
        self.B = self.omega.T
        self._tensor_cache: dict[tuple[str, str, torch.dtype], torch.Tensor] = {}

    def _tensor(self, name: str, value: np.ndarray, like: torch.Tensor) -> torch.Tensor:
        key = (name, str(like.device), like.dtype)
        result = self._tensor_cache.get(key)
        if result is None:
            result = torch.as_tensor(value, dtype=like.dtype, device=like.device)
            self._tensor_cache[key] = result
        return result

    def _h_torch(
        self, x: torch.Tensor, attenuation: torch.Tensor | float = 1.0
    ) -> tuple[torch.Tensor, torch.Tensor]:
        omega = self._tensor("omega", self.omega, x)
        signs = self._tensor("signs", self.signs, x)
        phases = self._tensor("phases", self.phases, x)
        angle = x @ omega.T + phases
        factor = signs * attenuation
        h = 0.5 + self.scale * torch.sum(factor * torch.sin(angle), dim=-1)
        grad_h = self.scale * ((factor * torch.cos(angle)) @ omega)
        return h, grad_h

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """Forward diffusion drift for the semilinear PDE (zero here)."""
        del t
        return torch.zeros_like(x)

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        h, _ = self._h_torch(x)
        return torch.asinh(h)

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        h, grad_h = self._h_torch(x)
        return grad_h / torch.sqrt(1.0 + h.square()).unsqueeze(-1)

    def terminal_numpy(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        angle = state @ self.omega.T + self.phases
        h = 0.5 + self.scale * np.sum(self.signs * np.sin(angle), axis=-1)
        return np.arcsinh(h)

    def terminal_gradient_numpy(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        angle = state @ self.omega.T + self.phases
        h = 0.5 + self.scale * np.sum(self.signs * np.sin(angle), axis=-1)
        grad_h = self.scale * ((self.signs * np.cos(angle)) @ self.omega)
        return grad_h / np.sqrt(1.0 + h * h)[..., None]

    def terminal_from_coordinates_numpy(self, s: np.ndarray) -> np.ndarray:
        angle = np.asarray(s, dtype=np.float64) + self.phases
        return np.arcsinh(
            0.5 + self.scale * np.sum(self.signs * np.sin(angle), axis=-1)
        )

    def terminal_coordinate_gradient_numpy(self, s: np.ndarray) -> np.ndarray:
        angle = np.asarray(s, dtype=np.float64) + self.phases
        h = 0.5 + self.scale * np.sum(self.signs * np.sin(angle), axis=-1)
        return (
            self.scale
            * self.signs
            * np.cos(angle)
            / np.sqrt(1.0 + h * h)[..., None]
        )

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading)
        tau = np.maximum(self.T - time, 0.0)
        attenuation = np.exp(-tau[..., None] * self.omega_norm_squared)
        angle = state @ self.omega.T + self.phases
        weighted = attenuation * self.signs
        w = 0.5 + self.scale * np.sum(weighted * np.sin(angle), axis=-1)
        grad_w = self.scale * ((weighted * np.cos(angle)) @ self.omega)
        u = np.arcsinh(w)
        z = self.sigma * grad_w / np.sqrt(1.0 + w * w)[..., None]
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 0]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 1:]

    def exact_u_torch(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        tau = torch.clamp(self.T - t.reshape(-1), min=0.0)
        norms = self._tensor("omega_n2", self.omega_norm_squared, x)
        attenuation = torch.exp(-tau[:, None] * norms)
        w, _ = self._h_torch(x, attenuation)
        return torch.asinh(w)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("Dense Fourier is jointly value-gradient nonlinear")
        return 0.5 * torch.tanh(y) * torch.sum(z.square(), dim=-1)

    def driver_y(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("Dense Fourier is jointly value-gradient nonlinear")
        return 0.5 * (1.0 - torch.tanh(y).square()) * torch.sum(
            z.square(), dim=-1
        )

    def driver_z(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        if z is None:
            raise ValueError("Dense Fourier is jointly value-gradient nonlinear")
        return torch.tanh(y).unsqueeze(-1) * z

    def driver_reduced_numpy(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        return 0.5 * np.tanh(y) * np.sum(np.square(z), axis=-1)

    def metadata(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "family": self.family,
            "d": self.d,
            "k": self.k,
            "T": self.T,
            "sigma": self.sigma,
            "has_gradient": True,
            "actor_input": "full ambient x",
            "actor_parameterization": "r=g+tau*R; q=sqrt(2)*grad(g)+tau*Q",
            "analytically_integrable_high_dimensional_stress_test": True,
            "parameter_seed": list(self.parameter_seed),
            "parameter_sha256": array_sha256(self.omega, self.signs, self.phases),
            "omega_sha256": array_sha256(self.omega),
            "signs_sha256": array_sha256(self.signs),
            "phases_sha256": array_sha256(self.phases),
            "omega_norm_max_error": float(
                np.max(np.abs(self.omega_norm_squared - 1.0))
            ),
            "oracle_projection_used": False,
        }


@lru_cache(maxsize=None)
def _laguerre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = roots_laguerre(int(order))
    return nodes.astype(np.float64), weights.astype(np.float64)


@dataclass
class CoupledLQGEquation:
    """Paper-compatible LQG with a frozen, globally defined coefficient draw."""

    c1: np.ndarray
    c2: np.ndarray
    T: float = 0.5
    coefficient_seed: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        self.c1 = np.asarray(self.c1, dtype=np.float64).copy()
        self.c2 = np.asarray(self.c2, dtype=np.float64).copy()
        if self.c1.shape != self.c2.shape or self.c1.ndim != 1:
            raise ValueError("c1/c2 must be matching vectors")
        if np.any((self.c1 < 0.5) | (self.c1 > 1.5)) or np.any(
            (self.c2 < 0.5) | (self.c2 > 1.5)
        ):
            raise ValueError("published coefficients lie in [0.5,1.5]")
        self.d = len(self.c1) + 1
        self.sigma = SQRT2
        self.family = "published_coupled_lqg_recreated"
        self.name = f"published_coupled_lqg_d{self.d}"
        A = np.zeros((self.d, self.d), dtype=np.float64)
        index = np.arange(self.d - 1)
        A[index, index] += self.c1
        A[index + 1, index + 1] += self.c1 + self.c2
        A[index, index + 1] -= self.c1
        A[index + 1, index] -= self.c1
        eigenvalues, eigenvectors = np.linalg.eigh(A)
        if eigenvalues[0] <= 0.0:
            raise ValueError("constructed LQG matrix is not SPD")
        self.A = A
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors
        self.coefficient_seed = tuple(int(value) for value in self.coefficient_seed)
        self._tensor_cache: dict[tuple[str, str, torch.dtype], torch.Tensor] = {}

    def _tensor(self, name: str, value: np.ndarray, like: torch.Tensor) -> torch.Tensor:
        key = (name, str(like.device), like.dtype)
        result = self._tensor_cache.get(key)
        if result is None:
            result = torch.as_tensor(value, dtype=like.dtype, device=like.device)
            self._tensor_cache[key] = result
        return result

    def quadratic(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        return np.einsum("...i,ij,...j->...", state, self.A, state)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return np.log((1.0 + self.quadratic(x)) / 2.0)

    def terminal_gradient(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        q = self.quadratic(state)
        return 2.0 * (state @ self.A) / (1.0 + q)[..., None]

    def terminal_torch(self, x: torch.Tensor) -> torch.Tensor:
        c1 = self._tensor("c1", self.c1, x)
        c2 = self._tensor("c2", self.c2, x)
        q = torch.sum(
            c1 * (x[..., :-1] - x[..., 1:]).square()
            + c2 * x[..., 1:].square(),
            dim=-1,
        )
        return torch.log((1.0 + q) / 2.0)

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        A = self._tensor("A", self.A, x)
        ax = x @ A
        q = torch.sum(ax * x, dim=-1)
        return 2.0 * ax / (1.0 + q).unsqueeze(-1)

    def exact_state(
        self, t: np.ndarray | float, x: np.ndarray, *, order: int = 192
    ) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        flat = state.reshape(-1, self.d)
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading).reshape(-1)
        tau = np.maximum(self.T - time, 0.0)
        y = flat @ self.eigenvectors
        nodes, weights = _laguerre_rule(order)
        q_zero = np.sum(self.eigenvalues[None, :] * np.square(y), axis=1)
        # In high dimension the integral is concentrated at s=O(1/d), while
        # vanilla Laguerre nodes live at O(1).  Match the exact logarithmic
        # slope at s=0 and integrate y=rho*s against exp(-y).  This is an
        # algebraic change of variables, not an approximation or fitted knob.
        rho = 1.0 + 2.0 * tau * np.sum(self.eigenvalues) + q_zero
        scaled_nodes = nodes[None, :] / rho[:, None]
        denominator = (
            1.0
            + 4.0
            * tau[:, None, None]
            * scaled_nodes[:, :, None]
            * self.eigenvalues[None, None, :]
        )
        log_factor = (
            nodes[None, :]
            - scaled_nodes
            - 0.5 * np.sum(np.log(denominator), axis=-1)
        )
        quadratic = np.sum(
            self.eigenvalues[None, None, :]
            * np.square(y[:, None, :])
            / denominator,
            axis=-1,
        )
        with np.errstate(divide="ignore"):
            log_weights = np.log(weights)
        weighted = np.exp(
            log_weights[None, :]
            + log_factor
            - scaled_nodes * quadratic
            - np.log(rho)[:, None]
        )
        w = 2.0 * np.sum(weighted, axis=1)
        eig_integral = np.sum(
            weighted[:, :, None]
            * scaled_nodes[:, :, None]
            * self.eigenvalues[None, None, :]
            * y[:, None, :]
            / denominator,
            axis=1,
        )
        grad_u = (4.0 / w[:, None]) * (eig_integral @ self.eigenvectors.T)
        u = -np.log(w)
        z = self.sigma * grad_u
        terminal_rows = tau <= 4.0 * np.finfo(np.float64).eps
        if np.any(terminal_rows):
            u[terminal_rows] = self.terminal(flat[terminal_rows])
            z[terminal_rows] = self.sigma * self.terminal_gradient(
                flat[terminal_rows]
            )
        return np.concatenate((u[:, None], z), axis=1).reshape(
            *leading, self.d + 1
        )

    def exact_u(self, t: np.ndarray | float, x: np.ndarray, *, order: int = 192) -> np.ndarray:
        return self.exact_state(t, x, order=order)[..., 0]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray, *, order: int = 192) -> np.ndarray:
        return self.exact_state(t, x, order=order)[..., 1:]

    def exact_u_torch(
        self, t: torch.Tensor, x: torch.Tensor, *, order: int = 192
    ) -> torch.Tensor:
        eigenvalues = self._tensor("eigenvalues", self.eigenvalues, x)
        eigenvectors = self._tensor("eigenvectors", self.eigenvectors, x)
        node_np, weight_np = _laguerre_rule(order)
        nodes = self._tensor(f"nodes_{order}", node_np, x)
        weights = self._tensor(f"weights_{order}", weight_np, x)
        tau = torch.clamp(self.T - t.reshape(-1), min=0.0)
        y = x @ eigenvectors
        q_zero = torch.sum(eigenvalues[None, :] * y.square(), dim=1)
        rho = 1.0 + 2.0 * tau * torch.sum(eigenvalues) + q_zero
        scaled_nodes = nodes[None, :] / rho[:, None]
        denominator = (
            1.0
            + 4.0
            * tau[:, None, None]
            * scaled_nodes[:, :, None]
            * eigenvalues[None, None, :]
        )
        log_factor = (
            nodes[None, :]
            - scaled_nodes
            - 0.5 * torch.sum(torch.log(denominator), dim=-1)
        )
        quadratic = torch.sum(
            eigenvalues[None, None, :] * y[:, None, :].square() / denominator,
            dim=-1,
        )
        weighted = torch.exp(
            torch.log(weights)[None, :]
            + log_factor
            - scaled_nodes * quadratic
            - torch.log(rho)[:, None]
        )
        w = 2.0 * torch.sum(weighted, dim=1)
        return -torch.log(w)

    def metadata(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "family": self.family,
            "label": "RECREATED_PUBLISHED_BENCHMARK",
            "d": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "coefficient_seed": list(self.coefficient_seed),
            "coefficient_sha256": array_sha256(self.c1, self.c2),
            "matrix_sha256": array_sha256(self.A),
            "eigenvalue_min": float(self.eigenvalues[0]),
            "eigenvalue_max": float(self.eigenvalues[-1]),
            "condition_number": float(self.eigenvalues[-1] / self.eigenvalues[0]),
            "paper_time_interval": [0.0, self.T],
            "state_distribution": "uniform volume measure on the unit ball",
            "official_repository_commit_audited": "c2412b942b1720084accf9ac9651c6085e073919",
            "official_source_conflict": (
                "ICLR Section 3.3 uses T=0.5; repository LQG.md/class defaults to T=1. "
                "The paper interval is used, and coefficients are one frozen global draw."
            ),
        }


class LQGActorBenchmark:
    """Gradient-only adapter for the mature terminal-aware Tangent-Actor."""

    actor_kind = "gradient"
    input_scale = 1.0
    value_bounds = None
    sigma = SQRT2

    def __init__(self, equation: CoupledLQGEquation) -> None:
        self.equation = equation
        self.d, self.T = equation.d, equation.T
        self.name = equation.name

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.equation.terminal_torch(x)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        return self.equation.terminal_gradient_torch(x)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return torch.zeros_like(x)

    def driver(self, r: torch.Tensor | None, q: torch.Tensor | None) -> torch.Tensor:
        del r
        if q is None:
            raise ValueError("LQG requires a gradient actor")
        return -0.5 * torch.sum(q.square(), dim=-1)

    def tangent_coefficients(
        self, r: torch.Tensor | None, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        del r
        if q is None:
            raise ValueError("LQG requires a gradient actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        return c, -q, 0.5 * torch.sum(q.square(), dim=-1)


def sample_unit_ball(rng: np.random.Generator, count: int, d: int) -> np.ndarray:
    direction = rng.standard_normal((int(count), int(d)))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    radius = rng.random(int(count)) ** (1.0 / int(d))
    return (direction * radius[:, None]).astype(np.float64)


__all__ = [
    "CoupledLQGEquation",
    "DenseFourierDriver",
    "LQGActorBenchmark",
    "SQRT2",
    "array_sha256",
    "sample_unit_ball",
]
