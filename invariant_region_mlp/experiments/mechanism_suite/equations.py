"""Exact equations, references, certificates, and fixed test distributions.

The MLP state is always ``(u, z)`` with ``z = sigma * grad(u)`` and every
array is evaluated in float64.  The equation classes intentionally expose the
same small interface consumed by the unchanged active-VB ``FullHistoryMLP``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.special import logsumexp, ndtri
from scipy.stats import qmc


HERE = Path(__file__).resolve().parent
VB_DIR = HERE.parent / "active_vb_high_budget"
if str(VB_DIR) not in sys.path:
    sys.path.insert(0, str(VB_DIR))

from vb_equation import make_test_points as make_vb_test_points  # noqa: E402


BASE_SEED = 20261006


def stable_seed(*parts: object) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def stable_sigmoid(value: np.ndarray | float) -> np.ndarray:
    value = np.asarray(value, dtype=np.float64)
    return np.exp(-np.logaddexp(0.0, -value))


def _as_points(t: np.ndarray | float, x: np.ndarray, d: int) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(t, dtype=np.float64)
    state = np.asarray(x, dtype=np.float64)
    if state.shape[-1] != d:
        raise ValueError(f"last state dimension must be {d}, got {state.shape}")
    return time, state


def make_uniform_points(
    equation: Any,
    *,
    n_points: int,
    seed: int,
    validation_fraction: float = 0.2,
    spatial_low: float = -1.0,
    spatial_high: float = 1.0,
) -> dict[str, np.ndarray]:
    """Fixed interior points with a deterministic held-out split."""

    rng = np.random.default_rng(
        np.random.SeedSequence([seed, equation.d, n_points, stable_seed(equation.name)])
    )
    t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), size=n_points).astype(np.float64)
    x = rng.uniform(spatial_low, spatial_high, size=(n_points, equation.d)).astype(np.float64)
    is_validation = np.zeros(n_points, dtype=bool)
    count = int(round(validation_fraction * n_points))
    if count:
        is_validation[rng.choice(n_points, size=count, replace=False)] = True
    return {
        "t": t,
        "x": x,
        "is_validation": is_validation,
        "is_boundary": np.zeros(n_points, dtype=bool),
    }


def make_points(equation: Any, *, n_points: int = 1200, seed: int = BASE_SEED) -> dict[str, np.ndarray]:
    """Dispatch to the pre-registered test geometry for an equation."""

    if equation.family in {"vba", "burgers_fisher", "published_vb"}:
        boundary = int(round(n_points / 6.0))
        interior = n_points - boundary
        points = make_vb_test_points(
            equation.d,
            n_interior=interior,
            n_boundary=boundary,
            seed=seed,
            validation_fraction=0.2,
        )
        # The inherited point generator uses T=0.5, which is exactly the
        # horizon of P2, P3, and N4.
        return {key: np.asarray(value) for key, value in points.items()}
    return make_uniform_points(equation, n_points=n_points, seed=seed)


class ExactEquation:
    """Small protocol shared by all analytic equations."""

    name: str
    family: str
    d: int
    T: float
    sigma: float
    mu: float | np.ndarray

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = np.asarray(self.exact_u(t, x), dtype=np.float64)
        z = np.asarray(self.exact_z(t, x), dtype=np.float64)
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_generator(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        state = self.exact_state(t, x)
        return self.generator(state[..., 0], state[..., 1:])

    @property
    def has_u_certificate(self) -> bool:
        return False

    @property
    def u_interval(self) -> tuple[float, float] | None:
        return None

    def certificate_labels(self) -> dict[str, str]:
        raise NotImplementedError

    def certificate_parameters(self) -> dict[str, Any]:
        raise NotImplementedError

    def provenance(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "family": self.family,
            "dimension": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "mu": np.asarray(self.mu, dtype=np.float64).tolist(),
            "dtype": "float64",
            "state_convention": "z=sigma*grad(u)",
            "certificate_labels": self.certificate_labels(),
            "certificate_parameters": self.certificate_parameters(),
        }


@dataclass(frozen=True)
class RidgeLSEHJB(ExactEquation):
    """P1: ridge log-sum-exp HJB with an analytic Hopf--Cole solution."""

    d: int
    kappa: float = 1.0
    direction_seed: int = BASE_SEED
    T: float = 0.25
    sigma: float = math.sqrt(2.0)

    name: str = "P1_ridge_lse_hjb"
    family: str = "ridge_lse"

    def __post_init__(self) -> None:
        if self.d <= 0 or self.kappa <= 0.0:
            raise ValueError("dimension and kappa must be positive")
        rng = np.random.default_rng(
            np.random.SeedSequence([self.direction_seed, self.d, 1])
        )
        direction = rng.standard_normal(self.d, dtype=np.float64)
        direction /= np.linalg.norm(direction)
        object.__setattr__(self, "w", direction)
        object.__setattr__(self, "lambdas", self.kappa * np.array([0.5, 1.5, 3.0], dtype=np.float64))
        object.__setattr__(self, "mu", 0.0)

    def _logits(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        s = np.einsum("...d,d->...", state, self.w)
        h = self.T - time
        return s[..., None] * self.lambdas + h[..., None] * self.lambdas**2

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return -(logsumexp(self._logits(t, x), axis=-1) - math.log(len(self.lambdas)))

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        normalization = logsumexp(logits, axis=-1, keepdims=True)
        u = -(normalization[..., 0] - math.log(len(self.lambdas)))
        weights = np.exp(logits - normalization)
        coefficient = np.sum(weights * self.lambdas, axis=-1)
        z = -self.sigma * coefficient[..., None] * self.w
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_grad(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
        coefficient = np.sum(weights * self.lambdas, axis=-1)
        return -coefficient[..., None] * self.w

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.exact_grad(t, x)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return self.exact_u(self.T, x)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
        first = np.sum(weights * self.lambdas, axis=-1)
        second = np.sum(weights * self.lambdas**2, axis=-1)
        u_t = second
        laplacian = -(second - first**2)
        z = -self.sigma * first[..., None] * self.w
        return u_t + laplacian + self.generator(np.zeros_like(first), z)

    def fzero_reference(self, t: np.ndarray, x: np.ndarray, quadrature_order: int = 80) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        s = state @ self.w
        h = np.maximum(self.T - time, 0.0)
        shifted = s[:, None] + 2.0 * np.sqrt(h)[:, None] * nodes[None, :]
        logits = shifted[..., None] * self.lambdas
        terminal = -(logsumexp(logits, axis=-1) - math.log(len(self.lambdas)))
        return np.sum(terminal * weights[None, :], axis=1) / math.sqrt(math.pi)

    @property
    def segment_coefficients(self) -> tuple[float, float]:
        return float(self.lambdas.min()), float(self.lambdas.max())

    @property
    def segment_endpoints(self) -> tuple[np.ndarray, np.ndarray]:
        low_c, high_c = self.segment_coefficients
        return -self.sigma * low_c * self.w, -self.sigma * high_c * self.w

    @property
    def box_low(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.minimum(first, second)

    @property
    def box_high(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.maximum(first, second)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma * float(self.lambdas.max())

    @property
    def z_center(self) -> np.ndarray:
        low_c, high_c = self.segment_coefficients
        return -self.sigma * (low_c + high_c) * 0.5 * self.w

    @property
    def dose_scale(self) -> float:
        low_c, high_c = self.segment_coefficients
        return self.sigma * (high_c - low_c) * 0.5

    @property
    def picard_proxy(self) -> float:
        return self.z_ball_radius * self.T

    def certificate_labels(self) -> dict[str, str]:
        return {
            "segment": "PDE-derived: Hopf-Cole tilted-expectation convex hull of grad(g)",
            "box": "PDE-derived: coordinatewise hull of the certified segment",
            "ball": "PDE-derived: Euclidean hull bound from the certified segment",
            "span": "PDE-derived: translation invariance orthogonal to w",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "w_seed": self.direction_seed,
            "w_sha256": hashlib.sha256(self.w.tobytes()).hexdigest(),
            "lambdas": self.lambdas.tolist(),
            "segment_coefficients": list(self.segment_coefficients),
            "box_low": self.box_low.tolist(),
            "box_high": self.box_high.tolist(),
            "ball_radius": self.z_ball_radius,
        }


@dataclass(frozen=True)
class VBa(ExactEquation):
    """P2 and N4: logistic product-driver equation."""

    d: int
    a: float = 4.0
    sigma: float = 0.5
    T: float = 0.5
    name: str = "P2_vba"
    family: str = "vba"

    def __post_init__(self) -> None:
        if self.d <= 0 or self.sigma <= 0.0:
            raise ValueError("dimension and sigma must be positive")
        object.__setattr__(self, "mu", -(self.a / self.d + 0.5 * self.sigma**2))

    @property
    def has_u_certificate(self) -> bool:
        return True

    @property
    def u_interval(self) -> tuple[float, float]:
        return (0.0, 1.0)

    @property
    def z_upper(self) -> float:
        return self.sigma / 4.0

    @property
    def box_low(self) -> np.ndarray:
        return np.zeros(self.d, dtype=np.float64)

    @property
    def box_high(self) -> np.ndarray:
        return np.full(self.d, self.z_upper, dtype=np.float64)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma * math.sqrt(self.d) / 4.0

    @property
    def z_center(self) -> np.ndarray:
        return np.full(self.d, self.z_upper / 2.0, dtype=np.float64)

    @property
    def dose_scale(self) -> float:
        # Pre-registered scale in the protocol (the certified coordinate cap).
        return self.z_upper

    @property
    def picard_proxy(self) -> float:
        # max |df/du| and ||df/dz|| are bounded jointly by this conservative
        # scalar proxy on [0,1] x [0,sigma/4]^d.
        return self.T * max(
            self.sigma * self.d * self.z_upper,
            self.sigma * math.sqrt(self.d),
        )

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        return stable_sigmoid(self.a * time + np.sum(state, axis=-1))

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        u = stable_sigmoid(self.a * time + np.sum(state, axis=-1))
        coordinate = self.sigma * u * (1.0 - u)
        z = np.broadcast_to(coordinate[..., None], state.shape)
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        coordinate = self.sigma * u * (1.0 - u)
        return np.broadcast_to(coordinate[..., None], np.asarray(x).shape).copy()

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return self.exact_u(self.T, x)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        return self.sigma * np.asarray(u, dtype=np.float64) * np.sum(
            np.asarray(z, dtype=np.float64), axis=-1
        )

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        q = u * (1.0 - u)
        u_t = self.a * q
        sum_grad = self.d * q
        laplacian = self.d * q * (1.0 - 2.0 * u)
        z = np.broadcast_to((self.sigma * q)[..., None], np.asarray(x).shape)
        return (
            u_t
            + self.mu * sum_grad
            + 0.5 * self.sigma**2 * laplacian
            + self.generator(u, z)
        )

    def fzero_reference(self, t: np.ndarray, x: np.ndarray, quadrature_order: int = 80) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        h = np.maximum(self.T - time, 0.0)
        mean_argument = self.a * self.T + np.sum(state, axis=-1) + self.d * self.mu * h
        std_argument = self.sigma * np.sqrt(self.d * h)
        arguments = mean_argument[:, None] + math.sqrt(2.0) * std_argument[:, None] * nodes[None, :]
        return np.sum(stable_sigmoid(arguments) * weights[None, :], axis=1) / math.sqrt(math.pi)

    def certificate_labels(self) -> dict[str, str]:
        return {
            "u_[0,1]": "PDE-derived: comparison with stationary constants 0 and 1",
            "z_nonnegative": "PDE-derived: coordinatewise monotonicity and comparison",
            "z_upper_sigma_over_4": "solution-informed: logistic profile u(1-u)<=1/4",
            "ball": "solution-informed: the coordinate cap implies |z|<=sigma*sqrt(d)/4",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "u_interval": [0.0, 1.0],
            "z_lower": 0.0,
            "z_upper": self.z_upper,
            "z_ball_radius": self.z_ball_radius,
            "a": self.a,
        }


@dataclass(frozen=True)
class BurgersFisher(VBa):
    """P3: product driver plus Fisher reaction."""

    rho: float = 1.0
    a: float = 4.0
    sigma: float = 0.5
    name: str = "P3_burgers_fisher"
    family: str = "burgers_fisher"

    def __post_init__(self) -> None:
        if self.rho <= 0.0:
            raise ValueError("rho must be positive")
        object.__setattr__(
            self,
            "mu",
            -((self.a + self.rho) / self.d + 0.5 * self.sigma**2),
        )

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        value = np.asarray(u, dtype=np.float64)
        return (
            self.sigma * value * np.sum(np.asarray(z, dtype=np.float64), axis=-1)
            + self.rho * value * (1.0 - value)
        )

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        q = u * (1.0 - u)
        u_t = self.a * q
        sum_grad = self.d * q
        laplacian = self.d * q * (1.0 - 2.0 * u)
        z = np.broadcast_to((self.sigma * q)[..., None], np.asarray(x).shape)
        return (
            u_t
            + self.mu * sum_grad
            + 0.5 * self.sigma**2 * laplacian
            + self.generator(u, z)
        )

    @property
    def picard_proxy(self) -> float:
        reaction_bound = self.rho
        return self.T * max(
            self.sigma * self.d * self.z_upper + reaction_bound,
            self.sigma * math.sqrt(self.d),
        )

    def certificate_parameters(self) -> dict[str, Any]:
        return {**super().certificate_parameters(), "rho": self.rho}


@dataclass(frozen=True)
class MultiDirectionLSEHJB(ExactEquation):
    """N3 reconstruction: equal-norm multi-direction Hopf--Cole HJB.

    The attachment referred to ``reference_code/lse_mlp.py`` but did not
    include it.  This is the mathematically canonical ``make_equation(d,k,s)``
    construction consistent with the specified exact Hopf--Cole family.
    """

    d: int = 100
    n_directions: int = 5
    strength: float = 0.5
    direction_seed: int = BASE_SEED
    T: float = 0.5
    sigma: float = math.sqrt(2.0)
    name: str = "N3_multidirection_lse"
    family: str = "multidirection_lse"

    def __post_init__(self) -> None:
        rng = np.random.default_rng(
            np.random.SeedSequence(
                [self.direction_seed, self.d, self.n_directions, int(1000 * self.strength)]
            )
        )
        directions = rng.standard_normal((self.n_directions, self.d), dtype=np.float64)
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        vectors = self.strength * directions
        object.__setattr__(self, "vectors", vectors)
        object.__setattr__(self, "mu", 0.0)

    def _logits(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        base = np.einsum("...d,kd->...k", state, self.vectors)
        norms2 = np.sum(self.vectors**2, axis=1)
        return base + (self.T - time)[..., None] * norms2

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return -(logsumexp(self._logits(t, x), axis=-1) - math.log(self.n_directions))

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        normalization = logsumexp(logits, axis=-1, keepdims=True)
        u = -(normalization[..., 0] - math.log(self.n_directions))
        weights = np.exp(logits - normalization)
        z = -self.sigma * np.einsum("...k,kd->...d", weights, self.vectors)
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_grad(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
        return -np.einsum("...k,kd->...d", weights, self.vectors)

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.exact_grad(t, x)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return self.exact_u(self.T, x)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
        norms2 = np.sum(self.vectors**2, axis=1)
        mean_vector = np.einsum("...k,kd->...d", weights, self.vectors)
        mean_norm2 = np.sum(weights * norms2, axis=-1)
        u_t = mean_norm2
        laplacian = -(mean_norm2 - np.sum(mean_vector**2, axis=-1))
        z = -self.sigma * mean_vector
        return u_t + laplacian + self.generator(np.zeros_like(u_t), z)

    def fzero_reference(self, t: np.ndarray, x: np.ndarray, qmc_power: int = 12) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        sampler = qmc.Sobol(
            d=self.n_directions,
            scramble=True,
            seed=stable_seed(self.name, self.d, self.n_directions, self.strength),
        )
        uniform = sampler.random_base2(qmc_power)
        normal = ndtri(np.clip(uniform, 1e-12, 1.0 - 1e-12))
        covariance_root = np.linalg.cholesky(
            self.vectors @ self.vectors.T + 1e-14 * np.eye(self.n_directions)
        )
        correlated = normal @ covariance_root.T
        base_logits = state @ self.vectors.T
        h = np.maximum(self.T - time, 0.0)
        result = np.empty(len(time), dtype=np.float64)
        for begin in range(0, len(time), 32):
            end = min(begin + 32, len(time))
            logits = (
                base_logits[begin:end, None, :]
                + self.sigma * np.sqrt(h[begin:end])[:, None, None] * correlated[None, :, :]
            )
            terminal = -(logsumexp(logits, axis=-1) - math.log(self.n_directions))
            result[begin:end] = np.mean(terminal, axis=1)
        return result

    @property
    def box_low(self) -> np.ndarray:
        return np.min(-self.sigma * self.vectors, axis=0)

    @property
    def box_high(self) -> np.ndarray:
        return np.max(-self.sigma * self.vectors, axis=0)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma * self.strength

    @property
    def z_center(self) -> np.ndarray:
        return 0.5 * (self.box_low + self.box_high)

    @property
    def dose_scale(self) -> float:
        return self.z_ball_radius

    @property
    def picard_proxy(self) -> float:
        return self.z_ball_radius * self.T

    def certificate_labels(self) -> dict[str, str]:
        return {
            "convex_hull": "PDE-derived: Hopf-Cole tilted-expectation hull",
            "box": "PDE-derived: coordinatewise hull of exact gradient vertices",
            "ball": "PDE-derived: equal-norm vertex hull",
            "span": "PDE-derived: translation invariance outside the direction span",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "direction_seed": self.direction_seed,
            "directions_sha256": hashlib.sha256(self.vectors.tobytes()).hexdigest(),
            "n_directions": self.n_directions,
            "strength": self.strength,
            "box_low": self.box_low.tolist(),
            "box_high": self.box_high.tolist(),
            "ball_radius": self.z_ball_radius,
        }


def published_vb(d: int) -> VBa:
    return VBa(
        d=d,
        a=1.0,
        sigma=math.sqrt(2.0),
        name="N4_published_vb",
        family="published_vb",
    )


def equation_registry() -> dict[str, list[ExactEquation]]:
    """All pre-registered dimension/parameter instances except P4."""

    equations: dict[str, list[ExactEquation]] = {
        "P1": [RidgeLSEHJB(d=d) for d in (20, 50, 100, 200)],
        "P2_a4": [VBa(d=d, a=4.0, name="P2_vba_a4") for d in (20, 50, 100)],
        "P2_a8": [VBa(d=d, a=8.0, name="P2_vba_a8") for d in (20, 50, 100)],
        "P3_rho1": [BurgersFisher(d=d, rho=1.0, name="P3_burgers_fisher_rho1") for d in (20, 50)],
        "P3_rho2": [BurgersFisher(d=d, rho=2.0, name="P3_burgers_fisher_rho2") for d in (20, 50)],
        "N3": [MultiDirectionLSEHJB()],
        "N4": [published_vb(d) for d in (20, 50, 100)],
    }
    return equations
