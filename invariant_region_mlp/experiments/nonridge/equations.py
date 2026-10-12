"""Spectral-separable full-rank HJB equations and frozen point sets."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from typing import Any

import numpy as np
from scipy.fft import dct, idct
from scipy.special import logsumexp

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    POINT_COUNT,
)
from invariant_region_mlp.experiments.mechanism_suite.equations import (
    ExactEquation,
    stable_seed,
)


SQRT2 = math.sqrt(2.0)
PROFILE_SLOPES = np.array([0.5, 1.5, 3.0], dtype=np.float64)
FAMILIES: dict[str, float] = {"S0": 0.0, "S05": 0.5, "S1": 1.0}


@lru_cache(maxsize=None)
def hermite_rule(order: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.hermite.hermgauss(int(order))
    nodes = np.asarray(nodes, dtype=np.float64)
    normalized = np.asarray(weights / math.sqrt(math.pi), dtype=np.float64)
    log_weights = np.log(normalized)
    nodes.setflags(write=False)
    normalized.setflags(write=False)
    log_weights.setflags(write=False)
    return nodes, normalized, log_weights


def profile(y: np.ndarray) -> np.ndarray:
    values = np.asarray(y, dtype=np.float64)
    return math.log(3.0) - logsumexp(values[..., None] * PROFILE_SLOPES, axis=-1)


def profile_prime(y: np.ndarray) -> np.ndarray:
    values = np.asarray(y, dtype=np.float64)
    logits = values[..., None] * PROFILE_SLOPES
    weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
    return -np.sum(weights * PROFILE_SLOPES, axis=-1)


def profile_second(y: np.ndarray) -> np.ndarray:
    values = np.asarray(y, dtype=np.float64)
    logits = values[..., None] * PROFILE_SLOPES
    weights = np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))
    first = np.sum(weights * PROFILE_SLOPES, axis=-1)
    second = np.sum(weights * PROFILE_SLOPES.square(), axis=-1)
    return -(second - first.square())


@dataclass(frozen=True)
class SpectralSeparableHJB(ExactEquation):
    """Hopf--Cole HJB with curvature spread across a fixed DCT spectrum."""

    d: int
    alpha: float
    pde_id: str
    T: float = 0.25
    sigma: float = SQRT2
    quadrature_order: int = 96
    family: str = "spectral_separable"

    def __post_init__(self) -> None:
        if self.d <= 0 or self.alpha < 0.0:
            raise ValueError("dimension must be positive and alpha nonnegative")
        if self.pde_id not in FAMILIES or FAMILIES[self.pde_id] != self.alpha:
            raise ValueError(f"inconsistent family parameters: {self.pde_id}, {self.alpha}")
        indices = np.arange(1, self.d + 1, dtype=np.float64)
        spectrum = indices ** (-self.alpha)
        spectrum /= np.linalg.norm(spectrum)
        rng = np.random.default_rng(
            np.random.SeedSequence([BASE_SEED, self.d, 31])
        )
        signs = 2 * rng.integers(0, 2, self.d) - 1
        object.__setattr__(self, "c", spectrum.astype(np.float64))
        object.__setattr__(self, "signs", signs.astype(np.float64))
        object.__setattr__(self, "mu", 0.0)
        object.__setattr__(self, "name", f"{self.pde_id}_spectral_separable_hjb")

    def rotate(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        if state.shape[-1] != self.d:
            raise ValueError(f"expected last dimension {self.d}, got {state.shape}")
        return dct(state * self.signs, type=2, norm="ortho", axis=-1)

    def inverse_rotate(self, y: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(y, dtype=np.float64)
        if coordinates.shape[-1] != self.d:
            raise ValueError(f"expected last dimension {self.d}, got {coordinates.shape}")
        return self.signs * idct(coordinates, type=2, norm="ortho", axis=-1)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        y = self.rotate(x)
        return np.sum(self.c * profile(y), axis=-1)

    def terminal_gradient_z(self, x: np.ndarray) -> np.ndarray:
        y = self.rotate(x)
        gradient_y = self.c * profile_prime(y)
        return self.sigma * self.inverse_rotate(gradient_y)

    def terminal_and_z(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        y = self.rotate(x)
        value = np.sum(self.c * profile(y), axis=-1)
        z = self.sigma * self.inverse_rotate(self.c * profile_prime(y))
        return value, z

    @staticmethod
    def one_dimensional_state(
        c: np.ndarray,
        tau: np.ndarray,
        y: np.ndarray,
        *,
        order: int = 96,
    ) -> tuple[np.ndarray, np.ndarray]:
        coefficient, horizon, location = np.broadcast_arrays(
            np.asarray(c, dtype=np.float64),
            np.asarray(tau, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
        )
        nodes, _weights, log_weights = hermite_rule(order)
        shifted = (
            location[..., None]
            + 2.0 * np.sqrt(np.maximum(horizon, 0.0))[..., None] * nodes
        )
        log_integrand = (
            -coefficient[..., None] * profile(shifted) + log_weights
        )
        log_denominator = logsumexp(log_integrand, axis=-1)
        tilted = np.exp(log_integrand - log_denominator[..., None])
        value = -log_denominator
        derivative = np.sum(
            tilted
            * coefficient[..., None]
            * profile_prime(shifted),
            axis=-1,
        )
        return value, derivative

    def _exact_flat(
        self,
        t: np.ndarray,
        x: np.ndarray,
        *,
        order: int,
        batch_size: int = 64,
    ) -> tuple[np.ndarray, np.ndarray]:
        values = np.empty(len(x), dtype=np.float64)
        gradients_y = np.empty_like(x)
        for begin in range(0, len(x), batch_size):
            end = min(begin + batch_size, len(x))
            y = self.rotate(x[begin:end])
            tau = np.maximum(self.T - t[begin:end], 0.0)
            v, derivative = self.one_dimensional_state(
                self.c[None, :],
                tau[:, None],
                y,
                order=order,
            )
            values[begin:end] = np.sum(v, axis=-1)
            gradients_y[begin:end] = derivative
        return values, self.sigma * self.inverse_rotate(gradients_y)

    def exact_state_order(
        self,
        t: np.ndarray | float,
        x: np.ndarray,
        *,
        order: int,
    ) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        if state.shape[-1] != self.d:
            raise ValueError(f"expected last dimension {self.d}, got {state.shape}")
        leading = state.shape[:-1]
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading)
        flat_x = state.reshape(-1, self.d)
        flat_t = time.reshape(-1)
        value, z = self._exact_flat(flat_t, flat_x, order=order)
        return np.concatenate((value[:, None], z), axis=-1).reshape(*leading, self.d + 1)

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state_order(t, x, order=self.quadrature_order)

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 0]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), state.shape[:-1])
        if np.all(np.abs(time - self.T) <= 4.0 * np.finfo(np.float64).eps):
            return self.terminal_gradient_z(state)
        return self.exact_state(time, state)[..., 1:]

    def heat_iterate(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading).reshape(-1)
        flat = state.reshape(-1, self.d)
        nodes, weights, _log_weights = hermite_rule(self.quadrature_order)
        output = np.empty(len(flat), dtype=np.float64)
        for begin in range(0, len(flat), 64):
            end = min(begin + 64, len(flat))
            y = self.rotate(flat[begin:end])
            tau = np.maximum(self.T - time[begin:end], 0.0)
            shifted = y[..., None] + 2.0 * np.sqrt(tau)[:, None, None] * nodes
            expectation = np.sum(profile(shifted) * weights, axis=-1)
            output[begin:end] = np.sum(self.c * expectation, axis=-1)
        return output.reshape(leading)

    def terminal_hessian_eigenvalues(self, x: np.ndarray) -> np.ndarray:
        return self.c * profile_second(self.rotate(x))

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    @staticmethod
    def double_generator(first: np.ndarray, second: np.ndarray) -> np.ndarray:
        return -0.5 * np.sum(
            np.asarray(first, dtype=np.float64)
            * np.asarray(second, dtype=np.float64),
            axis=-1,
        )

    @property
    def box_low(self) -> np.ndarray:
        return np.full(self.d, -np.inf, dtype=np.float64)

    @property
    def box_high(self) -> np.ndarray:
        return np.full(self.d, np.inf, dtype=np.float64)

    @property
    def z_ball_radius(self) -> float:
        return float("inf")

    @property
    def z_center(self) -> np.ndarray:
        return np.zeros(self.d, dtype=np.float64)

    @property
    def dose_scale(self) -> float:
        return 1.0

    @property
    def picard_proxy(self) -> float:
        return float("nan")

    def certificate_labels(self) -> dict[str, str]:
        return {"none": "diagnostic-only unbounded certificate attributes"}

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "pde_id": self.pde_id,
            "alpha": self.alpha,
            "quadrature_order": self.quadrature_order,
            "spectrum_sha256": hashlib.sha256(self.c.tobytes()).hexdigest(),
            "signs_sha256": hashlib.sha256(self.signs.tobytes()).hexdigest(),
        }


@lru_cache(maxsize=None)
def make_equation(pde_id: str, d: int) -> SpectralSeparableHJB:
    if pde_id not in FAMILIES:
        raise ValueError(f"unknown non-ridge PDE {pde_id!r}")
    return SpectralSeparableHJB(d=int(d), alpha=FAMILIES[pde_id], pde_id=pde_id)


@lru_cache(maxsize=None)
def fixed_points(pde_id: str, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    equation = make_equation(pde_id, d)
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, d, POINT_COUNT, stable_seed(pde_id)])
    )
    t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), POINT_COUNT).astype(np.float64)
    x = rng.uniform(-1.0, 1.0, (POINT_COUNT, d)).astype(np.float64)
    is_validation = np.zeros(POINT_COUNT, dtype=bool)
    selected = rng.choice(POINT_COUNT, size=POINT_COUNT // 5, replace=False)
    is_validation[selected] = True
    t.setflags(write=False)
    x.setflags(write=False)
    is_validation.setflags(write=False)
    return t, x, is_validation


__all__ = [
    "FAMILIES",
    "PROFILE_SLOPES",
    "SQRT2",
    "SpectralSeparableHJB",
    "fixed_points",
    "hermite_rule",
    "make_equation",
    "profile",
    "profile_prime",
    "profile_second",
]
