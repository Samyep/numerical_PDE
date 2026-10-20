"""Verified B6 and Dense-Fourier benchmark adapters."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from functools import lru_cache
from typing import Any

import numpy as np
import torch
from scipy.special import logsumexp

from invariant_region_mlp.experiments.double_estimator.equations import POINT_COUNT
from invariant_region_mlp.experiments.mechanism_suite.equations import (
    ExactEquation,
    stable_seed as point_stable_seed,
)
from invariant_region_mlp.experiments.main_highdim_pilot.equations import (
    DenseFourierDriver,
    array_sha256,
)

from .protocol import BASE_SEED, CALIBRATION_PATH


SQRT2 = math.sqrt(2.0)


@dataclass(frozen=True)
class B6Equation(ExactEquation):
    """Full-rank nonseparable log-sum-exp Hopf--Cole benchmark."""

    d: int
    beta: float
    gamma: float
    T: float = 0.25
    sigma: float = SQRT2
    family: str = "b6_fullrank_lse"
    name: str = "B6_fullrank_nonseparable_hjb"

    def __post_init__(self) -> None:
        if self.d <= 0 or self.beta <= 0.0 or self.gamma < 0.0:
            raise ValueError((self.d, self.beta, self.gamma))
        count = 2 * int(self.d)
        rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, self.d, 66]))
        directions = rng.standard_normal((count, self.d))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        multipliers = rng.uniform(0.5, 1.5, count)
        common = rng.standard_normal(self.d)
        common /= np.linalg.norm(common)
        vectors = (
            self.gamma * common[None, :]
            + self.beta * multipliers[:, None] * directions
        ).astype(np.float64)
        squared_norms = np.sum(vectors * vectors, axis=1)
        for value in (directions, multipliers, common, vectors, squared_norms):
            value.setflags(write=False)
        object.__setattr__(self, "K", count)
        object.__setattr__(self, "e", directions)
        object.__setattr__(self, "U", multipliers)
        object.__setattr__(self, "v", common)
        object.__setattr__(self, "b", vectors)
        object.__setattr__(self, "n2", squared_norms)
        object.__setattr__(self, "mu", 0.0)

    def _flat_logits(
        self, t: np.ndarray | float, x: np.ndarray
    ) -> tuple[np.ndarray, tuple[int, ...]]:
        state = np.asarray(x, dtype=np.float64)
        if state.shape[-1] != self.d:
            raise ValueError((state.shape, self.d))
        leading = state.shape[:-1]
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading).reshape(-1)
        flat = state.reshape(-1, self.d)
        logits = flat @ self.b.T + np.maximum(self.T - time, 0.0)[:, None] * self.n2
        return logits, leading

    def terminal(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        logits = state.reshape(-1, self.d) @ self.b.T
        return (-(logsumexp(logits, axis=1) - math.log(self.K))).reshape(leading)

    def terminal_gradient_z(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        logits = state.reshape(-1, self.d) @ self.b.T
        weights = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
        return (-self.sigma * (weights @ self.b)).reshape(*leading, self.d)

    def terminal_and_z(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return self.terminal(x), self.terminal_gradient_z(x)

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits, leading = self._flat_logits(t, x)
        normalizer = logsumexp(logits, axis=1, keepdims=True)
        weights = np.exp(logits - normalizer)
        value = -(normalizer[:, 0] - math.log(self.K))
        z = -self.sigma * (weights @ self.b)
        return np.concatenate((value[:, None], z), axis=1).reshape(*leading, self.d + 1)

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 0]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 1:]

    def probabilities(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits, leading = self._flat_logits(t, x)
        result = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
        return result.reshape(*leading, self.K)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        return -0.5 * np.sum(np.square(z), axis=-1)

    @staticmethod
    def double_generator(first: np.ndarray, second: np.ndarray) -> np.ndarray:
        return -0.5 * np.sum(np.asarray(first) * np.asarray(second), axis=-1)

    @property
    def box_low(self) -> np.ndarray:
        return np.full(self.d, -np.inf)

    @property
    def box_high(self) -> np.ndarray:
        return np.full(self.d, np.inf)

    @property
    def z_ball_radius(self) -> float:
        return float("inf")

    @property
    def z_center(self) -> np.ndarray:
        return np.zeros(self.d)

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
            "pde_id": "B6",
            "beta": self.beta,
            "gamma": self.gamma,
            "K": self.K,
            "construction_seed": [BASE_SEED, self.d, 66],
            "directions_sha256": hashlib.sha256(self.e.tobytes()).hexdigest(),
            "multipliers_sha256": hashlib.sha256(self.U.tobytes()).hexdigest(),
            "common_vector_sha256": hashlib.sha256(self.v.tobytes()).hexdigest(),
            "vectors_sha256": hashlib.sha256(self.b.tobytes()).hexdigest(),
        }


@lru_cache(maxsize=None)
def make_b6_equation(d: int, beta: float, gamma: float) -> B6Equation:
    return B6Equation(int(d), float(beta), float(gamma))


@lru_cache(maxsize=None)
def fixed_b6_points(d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, int(d), POINT_COUNT, point_stable_seed("B6")])
    )
    t = rng.uniform(0.0, np.nextafter(0.25, 0.0), POINT_COUNT).astype(np.float64)
    x = rng.uniform(-1.0, 1.0, (POINT_COUNT, int(d))).astype(np.float64)
    is_validation = np.zeros(POINT_COUNT, dtype=bool)
    selected = rng.choice(POINT_COUNT, size=POINT_COUNT // 5, replace=False)
    is_validation[selected] = True
    for value in (t, x, is_validation):
        value.setflags(write=False)
    return t, x, is_validation


class HJBActorBenchmark:
    T = 0.25
    sigma = SQRT2
    actor_kind = "gradient"
    input_scale = 1.0
    value_bounds = None

    def __init__(self, equation: B6Equation):
        self.pde = "b6_fullrank"
        self.name = f"B6_d{equation.d}"
        self.d = equation.d
        self.equation = equation
        self._tensor_cache: dict[tuple[str, str, torch.dtype], torch.Tensor] = {}

    def _vectors(self, x: torch.Tensor) -> torch.Tensor:
        key = ("vectors", str(x.device), x.dtype)
        value = self._tensor_cache.get(key)
        if value is None:
            value = torch.as_tensor(self.equation.b, dtype=x.dtype, device=x.device)
            self._tensor_cache[key] = value
        return value

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        logits = x @ self._vectors(x).T
        return -(torch.logsumexp(logits, dim=-1) - math.log(self.equation.K))

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        vectors = self._vectors(x)
        return -(torch.softmax(x @ vectors.T, dim=-1) @ vectors)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return torch.zeros_like(x)

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del r
        if q is None:
            raise ValueError("B6 needs a Gradient Actor")
        return -0.5 * q.square().sum(-1)

    def tangent_coefficients(self, r, q):
        del r
        if q is None:
            raise ValueError("B6 needs a Gradient Actor")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        return c, -q, 0.5 * q.square().sum(-1)

    def sample_states(self, count: int, *, seed: int, purpose: str = "function"):
        if purpose not in {"function", "official"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), int(count))
        x = rng.uniform(-1.0, 1.0, (int(count), self.d))
        return t.astype(np.float64), x.astype(np.float64)


def dense_parameters(d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce the pilot's exact draw order and seed namespace."""
    d = int(d)
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, d, 41]))
    omega = rng.standard_normal((2 * d, d)).astype(np.float64)
    omega /= np.linalg.norm(omega, axis=1, keepdims=True)
    signs = (2 * rng.integers(0, 2, 2 * d) - 1).astype(np.float64)
    phases = rng.uniform(0.0, 2.0 * np.pi, 2 * d).astype(np.float64)
    return omega, signs, phases


@lru_cache(maxsize=None)
def make_dense_fourier(d: int) -> DenseFourierDriver:
    omega, signs, phases = dense_parameters(int(d))
    return DenseFourierDriver(
        omega,
        signs,
        phases,
        T=0.5,
        parameter_seed=(BASE_SEED, int(d), 41),
    )


def _b6_parameters(d: int) -> tuple[float, float]:
    if not CALIBRATION_PATH.is_file():
        raise RuntimeError(
            f"B6 parameters are not frozen: missing {CALIBRATION_PATH}. "
            "Run calibration before any training."
        )
    payload = json.loads(CALIBRATION_PATH.read_text(encoding="utf-8"))
    row = payload["selections"][str(int(d))]
    return float(row["selected_beta"]), float(row["selected_gamma"])


@lru_cache(maxsize=None)
def make_b6(d: int) -> HJBActorBenchmark:
    beta, gamma = _b6_parameters(int(d))
    return HJBActorBenchmark(make_b6_equation(int(d), beta=beta, gamma=gamma))


def make_benchmark(pde: str, d: int):
    if pde == "b6_fullrank":
        return make_b6(int(d))
    if pde == "dense_fourier":
        return make_dense_fourier(int(d))
    raise KeyError(pde)


def exact_state(benchmark: Any, t: np.ndarray, x: np.ndarray) -> np.ndarray:
    if isinstance(benchmark, HJBActorBenchmark):
        return benchmark.equation.exact_state(t, x)
    if isinstance(benchmark, DenseFourierDriver):
        return benchmark.exact_state(t, x)
    raise TypeError(type(benchmark))


def exact_u(benchmark: Any, t: np.ndarray, x: np.ndarray) -> np.ndarray:
    return exact_state(benchmark, t, x)[..., 0]


def exact_z(benchmark: Any, t: np.ndarray, x: np.ndarray) -> np.ndarray:
    return exact_state(benchmark, t, x)[..., 1:]


def terminal_gradient_torch(benchmark: Any, x: torch.Tensor) -> torch.Tensor:
    """Return grad g through either validated benchmark interface."""
    if isinstance(benchmark, HJBActorBenchmark):
        return benchmark.terminal_gradient(x)
    if isinstance(benchmark, DenseFourierDriver):
        return benchmark.terminal_gradient_torch(x)
    raise TypeError(type(benchmark))


def sample_cube_states(
    pde: str, d: int, count: int, *, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    benchmark = make_benchmark(pde, d)
    rng = np.random.default_rng(int(seed))
    t = rng.uniform(0.0, np.nextafter(benchmark.T, 0.0), int(count))
    x = rng.uniform(-1.0, 1.0, (int(count), int(d)))
    return t.astype(np.float64), x.astype(np.float64)


def metadata(benchmark: Any) -> dict[str, Any]:
    if isinstance(benchmark, HJBActorBenchmark):
        equation: B6Equation = benchmark.equation
        result = equation.certificate_parameters()
        result.update(
            {
                "name": benchmark.name,
                "family": equation.family,
                "d": equation.d,
                "T": equation.T,
                "sigma": equation.sigma,
                "beta": equation.beta,
                "gamma": equation.gamma,
                "actor_kind": "gradient",
                "actor_input": "full ambient x",
                "oracle_projection_used": False,
                "parameter_sha256": array_sha256(
                    equation.e, equation.U, equation.v, equation.b, equation.n2
                ),
            }
        )
        return result
    if isinstance(benchmark, DenseFourierDriver):
        return benchmark.metadata()
    raise TypeError(type(benchmark))


__all__ = [
    "B6Equation",
    "DenseFourierDriver",
    "array_sha256",
    "dense_parameters",
    "exact_state",
    "exact_u",
    "exact_z",
    "fixed_b6_points",
    "make_b6",
    "make_benchmark",
    "make_dense_fourier",
    "metadata",
    "sample_cube_states",
]
