"""S05 and B6 equations plus Tangent-Actor benchmark adapters."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from typing import Any

import numpy as np
import torch
from scipy.special import logsumexp

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    POINT_COUNT,
)
from invariant_region_mlp.experiments.mechanism_suite.equations import (
    ExactEquation,
    stable_seed as point_stable_seed,
)
from invariant_region_mlp.experiments.nonridge.equations import (
    PROFILE_SLOPES,
    SpectralSeparableHJB,
    fixed_points as nonridge_fixed_points,
    make_equation as make_nonridge_equation,
)
from invariant_region_mlp.experiments.nonridge.solver import (
    torch_dct_ortho,
    torch_idct_ortho,
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


@dataclass(frozen=True)
class B6Equation(ExactEquation):
    """Full-rank non-separable log-sum-exp Hopf--Cole benchmark."""

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
        rng = np.random.default_rng(
            np.random.SeedSequence([BASE_SEED, int(self.d), 66])
        )
        directions = rng.standard_normal((count, self.d))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        multipliers = rng.uniform(0.5, 1.5, count)
        common = rng.standard_normal(self.d)
        common /= np.linalg.norm(common)
        vectors = (
            float(self.gamma) * common[None, :]
            + float(self.beta) * multipliers[:, None] * directions
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
            raise ValueError(f"expected last dimension {self.d}, got {state.shape}")
        leading = state.shape[:-1]
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), leading).reshape(-1)
        flat = state.reshape(-1, self.d)
        tau = np.maximum(self.T - time, 0.0)
        logits = flat @ self.b.T + tau[:, None] * self.n2
        return logits, leading

    def terminal(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        logits = state.reshape(-1, self.d) @ self.b.T
        value = -(logsumexp(logits, axis=1) - math.log(self.K))
        return value.reshape(leading)

    def terminal_gradient_z(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        leading = state.shape[:-1]
        logits = state.reshape(-1, self.d) @ self.b.T
        weights = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
        z = -self.sigma * (weights @ self.b)
        return z.reshape(*leading, self.d)

    def terminal_and_z(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return self.terminal(x), self.terminal_gradient_z(x)

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits, leading = self._flat_logits(t, x)
        normalizer = logsumexp(logits, axis=1, keepdims=True)
        weights = np.exp(logits - normalizer)
        value = -(normalizer[:, 0] - math.log(self.K))
        z = -self.sigma * (weights @ self.b)
        return np.concatenate((value[:, None], z), axis=1).reshape(
            *leading, self.d + 1
        )

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 0]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.exact_state(t, x)[..., 1:]

    def probabilities(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits, leading = self._flat_logits(t, x)
        weights = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
        return weights.reshape(*leading, self.K)

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
def fixed_points(pde_id: str, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if pde_id == "S05":
        return nonridge_fixed_points("S05", int(d))
    if pde_id != "B6":
        raise KeyError(pde_id)
    rng = np.random.default_rng(
        np.random.SeedSequence(
            [BASE_SEED, int(d), POINT_COUNT, point_stable_seed("B6")]
        )
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
    """Duck-typed adapter for the frozen Terminal-Aware actor machinery."""

    T = 0.25
    sigma = SQRT2
    actor_kind = "gradient"
    input_scale = 1.0
    value_bounds = None

    def __init__(self, pde_id: str, equation: SpectralSeparableHJB | B6Equation):
        self.pde = pde_id
        self.name = f"{pde_id}_d{equation.d}"
        self.d = equation.d
        self.equation = equation
        self._tensor_cache: dict[tuple[str, str, torch.dtype], torch.Tensor] = {}

    def _tensor(
        self, name: str, array: np.ndarray, like: torch.Tensor
    ) -> torch.Tensor:
        key = (name, str(like.device), like.dtype)
        value = self._tensor_cache.get(key)
        if value is None:
            value = torch.tensor(array, dtype=like.dtype, device=like.device)
            self._tensor_cache[key] = value
        return value

    def _s05_parameters(
        self, x: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        equation = self.equation
        if not isinstance(equation, SpectralSeparableHJB):
            raise TypeError(type(equation))
        signs = self._tensor("signs", equation.signs, x)
        coefficients = self._tensor("coefficients", equation.c, x)
        slopes = self._tensor("slopes", PROFILE_SLOPES, x)
        return signs, coefficients, slopes

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        if isinstance(self.equation, B6Equation):
            vectors = self._tensor("vectors", self.equation.b, x)
            logits = torch.matmul(x, vectors.T)
            return -(torch.logsumexp(logits, dim=-1) - math.log(self.equation.K))
        signs, coefficients, slopes = self._s05_parameters(x)
        y = torch_dct_ortho(x * signs)
        values = math.log(3.0) - torch.logsumexp(
            y.unsqueeze(-1) * slopes, dim=-1
        )
        return torch.sum(coefficients * values, dim=-1)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        if isinstance(self.equation, B6Equation):
            vectors = self._tensor("vectors", self.equation.b, x)
            weights = torch.softmax(torch.matmul(x, vectors.T), dim=-1)
            return -torch.matmul(weights, vectors)
        signs, coefficients, slopes = self._s05_parameters(x)
        y = torch_dct_ortho(x * signs)
        logits = y.unsqueeze(-1) * slopes
        weights = torch.softmax(logits, dim=-1)
        derivative = -torch.sum(weights * slopes, dim=-1)
        return signs * torch_idct_ortho(coefficients * derivative)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return torch.zeros_like(x)

    def driver(self, r: torch.Tensor, q: torch.Tensor | None) -> torch.Tensor:
        del r
        if q is None:
            raise ValueError("gradient actor required")
        return -0.5 * torch.sum(q.square(), dim=-1)

    def tangent_coefficients(
        self, r: torch.Tensor | None, q: torch.Tensor | None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        del r
        if q is None:
            raise ValueError("gradient actor required")
        c = torch.zeros(q.shape[:-1], dtype=q.dtype, device=q.device)
        return c, -q, 0.5 * torch.sum(q.square(), dim=-1)

    def sample_states(
        self, count: int, *, seed: int, purpose: str = "function"
    ) -> tuple[np.ndarray, np.ndarray]:
        if purpose not in {"function", "official"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(seed)
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), int(count))
        x = rng.uniform(-1.0, 1.0, (int(count), self.d))
        return t.astype(np.float64), x.astype(np.float64)


def make_s05_equation(d: int) -> SpectralSeparableHJB:
    return make_nonridge_equation("S05", int(d))


def make_actor_benchmark(
    pde_id: str,
    d: int,
    *,
    beta: float | None = None,
    gamma: float | None = None,
) -> HJBActorBenchmark:
    if pde_id == "S05":
        equation: SpectralSeparableHJB | B6Equation = make_s05_equation(d)
    elif pde_id == "B6":
        if beta is None or gamma is None:
            raise ValueError("B6 requires frozen beta and gamma")
        equation = make_b6_equation(d, beta, gamma)
    else:
        raise KeyError(pde_id)
    return HJBActorBenchmark(pde_id, equation)


__all__ = [
    "B6Equation",
    "HJBActorBenchmark",
    "SQRT2",
    "array_sha256",
    "fixed_points",
    "make_actor_benchmark",
    "make_b6_equation",
    "make_s05_equation",
]
