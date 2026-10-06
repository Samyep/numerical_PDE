"""Source-faithful Allen--Cahn equation and shared random-tree utilities.

The numerical instance is Becker et al. (2020), CICP 28(5), Section 3.1:

    u_t + Delta u + u - u^3 = 0,       T = 1,
    u(T,x) = 1 / (2 + (2/5) ||x||^2).

The forward diffusion is ``x + sqrt(2) (W_s-W_t)``.  Beck et al.'s
published numerical code truncates the scalar Picard state to [-4, 4].
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Iterable, Sequence

import numpy as np


HORIZON = 1.0
BECK_NUMERICAL_RADIUS = 4.0
CERTIFIED_INTERVAL = (0.0, 1.0)
TERMINAL_QUADRATIC_COEFFICIENT = 2.0 / 5.0

# Table 1 of Becker et al. (2020).  These are numerical references, not exact
# values: each is an average of five V_{8,8,4}(0,0) evaluations.
PUBLISHED_MLP_REFERENCES = {10: 0.29555, 100: 0.03373, 1000: 0.00340}
PUBLISHED_DS_REFERENCES = {10: 0.29614, 100: 0.03376, 1000: 0.00339}


def reaction(value: float | np.ndarray) -> float | np.ndarray:
    """Allen--Cahn reaction f(u)=u-u^3, evaluated in float64."""

    u = np.asarray(value, dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        result = u - u * u * u
    if result.ndim == 0:
        return float(result)
    return result


def reaction_derivative(value: float | np.ndarray) -> float | np.ndarray:
    u = np.asarray(value, dtype=np.float64)
    result = 1.0 - 3.0 * u * u
    if result.ndim == 0:
        return float(result)
    return result


def sharp_clipped_driver_lipschitz(radius: float) -> float:
    """Sharp global Lipschitz constant of u -> f(P_r(u))."""

    if radius < 0.0:
        raise ValueError("radius must be nonnegative")
    return max(1.0, abs(1.0 - 3.0 * radius * radius))


def theorem_example_radius(sample_parameter: int) -> float:
    """The explicit rho_M=log(1+log(M)) example used in Beck Theorem 1.1."""

    if sample_parameter < 1:
        raise ValueError("sample_parameter must be positive")
    return math.log1p(math.log(float(sample_parameter)))


@dataclass(frozen=True)
class AllenCahnEquation:
    dimension: int
    horizon: float = HORIZON

    def __post_init__(self) -> None:
        if self.dimension < 1:
            raise ValueError("dimension must be positive")
        if self.horizon <= 0.0:
            raise ValueError("horizon must be positive")

    def terminal(self, state: np.ndarray) -> float | np.ndarray:
        x = np.asarray(state, dtype=np.float64)
        squared_norm = np.sum(x * x, axis=-1, dtype=np.float64)
        value = 1.0 / (2.0 + TERMINAL_QUADRATIC_COEFFICIENT * squared_norm)
        if np.ndim(value) == 0:
            return float(value)
        return value

    def transition(
        self, state: np.ndarray, elapsed: float | np.ndarray, normal: np.ndarray
    ) -> np.ndarray:
        x = np.asarray(state, dtype=np.float64)
        z = np.asarray(normal, dtype=np.float64)
        dt = np.asarray(elapsed, dtype=np.float64)
        scale = np.sqrt(2.0 * dt)
        if scale.ndim > 0:
            scale = scale[..., None]
        return x + scale * z


@dataclass(frozen=True)
class IntervalStateProjector:
    """Euclidean projector onto [low,high] x R^d."""

    low: float
    high: float

    def __post_init__(self) -> None:
        if not self.low <= self.high:
            raise ValueError("low must not exceed high")

    def project_value(self, value: float) -> float:
        return float(min(self.high, max(self.low, float(value))))

    def project_state(
        self, value: float, gradient_state: np.ndarray
    ) -> tuple[float, np.ndarray]:
        z = np.asarray(gradient_state, dtype=np.float64)
        return self.project_value(value), z.copy()

    def overshoot(self, value: float) -> float:
        u = float(value)
        return max(self.low - u, 0.0, u - self.high)


@dataclass
class RunningMoments:
    count: int = 0
    mean: float = 0.0
    m2: float = 0.0
    nonfinite: int = 0

    def add(self, value: float) -> None:
        x = float(value)
        if not math.isfinite(x):
            self.nonfinite += 1
            return
        self.count += 1
        delta = x - self.mean
        self.mean += delta / self.count
        self.m2 += delta * (x - self.mean)

    @property
    def variance(self) -> float:
        return self.m2 / (self.count - 1) if self.count > 1 else 0.0


@dataclass
class WorkDiagnostics:
    """Actual work and online mechanism diagnostics for one root estimate."""

    capture_trace: bool = False
    recursive_states: int = 0
    terminal_g_evals: int = 0
    f_evals: int = 0
    normal_scalar_draws: int = 0
    uniform_draws: int = 0
    state_transitions: int = 0
    checked_states: int = 0
    violating_states: int = 0
    projection_activations: int = 0
    overshoot_square_sum: float = 0.0
    nonfinite_states: int = 0
    pre_projection_min: float = math.inf
    pre_projection_max: float = -math.inf
    generator_before_abs: RunningMoments = field(default_factory=RunningMoments)
    generator_after_abs: RunningMoments = field(default_factory=RunningMoments)
    correction_moments: RunningMoments = field(default_factory=RunningMoments)
    correction_trace: list[float] = field(default_factory=list)

    def observe_driver(
        self,
        value_before: float,
        value_after: float,
        interval: IntervalStateProjector,
        generator_before: float,
        generator_after: float,
    ) -> None:
        before = float(value_before)
        after = float(value_after)
        self.checked_states += 1
        overshoot = interval.overshoot(before)
        if overshoot > 0.0:
            self.violating_states += 1
        self.overshoot_square_sum += overshoot * overshoot
        if after != before:
            self.projection_activations += 1
        if not math.isfinite(before) or not math.isfinite(after):
            self.nonfinite_states += 1
        else:
            self.pre_projection_min = min(self.pre_projection_min, before)
            self.pre_projection_max = max(self.pre_projection_max, before)
        self.generator_before_abs.add(abs(float(generator_before)))
        self.generator_after_abs.add(abs(float(generator_after)))
        self.f_evals += 1

    def observe_correction(self, correction: float) -> None:
        value = float(correction)
        self.correction_moments.add(value)
        if self.capture_trace:
            self.correction_trace.append(value)

    def as_dict(self) -> dict[str, float | int]:
        checked = max(self.checked_states, 1)
        return {
            "recursive_states": self.recursive_states,
            "terminal_g_evals": self.terminal_g_evals,
            "f_evals": self.f_evals,
            "normal_scalar_draws": self.normal_scalar_draws,
            "uniform_draws": self.uniform_draws,
            "total_stochastic_samples": self.normal_scalar_draws + self.uniform_draws,
            "state_transitions": self.state_transitions,
            "checked_states": self.checked_states,
            "violating_states": self.violating_states,
            "pre_truncation_violation_rate": self.violating_states / checked,
            "projection_activations": self.projection_activations,
            "truncation_activation_rate": self.projection_activations / checked,
            "mean_squared_overshoot": self.overshoot_square_sum / checked,
            "generator_abs_before_mean": self.generator_before_abs.mean,
            "generator_abs_after_mean": self.generator_after_abs.mean,
            "generator_before_nonfinite": self.generator_before_abs.nonfinite,
            "generator_after_nonfinite": self.generator_after_abs.nonfinite,
            "nonlinear_correction_variance": self.correction_moments.variance,
            "nonfinite_corrections": self.correction_moments.nonfinite,
            "nonfinite_states": self.nonfinite_states,
            "pre_projection_min": (
                self.pre_projection_min
                if math.isfinite(self.pre_projection_min)
                else 0.0
            ),
            "pre_projection_max": (
                self.pre_projection_max
                if math.isfinite(self.pre_projection_max)
                else 0.0
            ),
        }


@dataclass(frozen=True)
class KeyedRandomTree:
    """Order-independent random tree shared pathwise by all methods.

    Each logical node/stream obtains a NumPy ``Generator`` from a SeedSequence
    composed of the experiment seed, dimension, path, and stream identifier.
    The implementations can therefore be independent without relying on mutable
    RNG call order.
    """

    seed: int
    dimension: int

    def _generator(self, path: Sequence[int], stream: int) -> np.random.Generator:
        if any(int(item) < 0 for item in path):
            raise ValueError("random-tree path components must be nonnegative")
        entropy = [
            int(self.seed) & 0xFFFFFFFF,
            (int(self.seed) >> 32) & 0xFFFFFFFF,
            int(self.dimension),
            int(stream),
            *(int(item) for item in path),
        ]
        return np.random.default_rng(np.random.SeedSequence(entropy))

    def normal(
        self, path: Sequence[int], stream: int, shape: tuple[int, ...]
    ) -> np.ndarray:
        return self._generator(path, stream).standard_normal(shape, dtype=np.float64)

    def uniform(
        self, path: Sequence[int], stream: int, count: int
    ) -> np.ndarray:
        return self._generator(path, stream).random(count, dtype=np.float64)

    @property
    def fingerprint(self) -> str:
        payload = json.dumps(
            {"schema": 1, "seed": int(self.seed), "dimension": self.dimension},
            sort_keys=True,
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


def child_path(
    parent: Sequence[int], level: int, sample: int, branch: int
) -> tuple[int, ...]:
    """Stable nonnegative key for a recursive plus/minus child."""

    if branch not in (0, 1):
        raise ValueError("branch must be 0 (fine) or 1 (coarse)")
    return (*parent, level + 1, sample + 1, branch + 1)


def deterministic_test_point(dimension: int, kind: str) -> np.ndarray:
    if kind == "zero":
        return np.zeros(dimension, dtype=np.float64)
    if kind == "ramp":
        values = np.arange(1, dimension + 1, dtype=np.float64)
        return 0.1 * values / np.linalg.norm(values)
    raise ValueError(f"unknown test-point kind {kind!r}")


def aggregate_fingerprint(items: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for item in items:
        digest.update(item.encode("ascii"))
    return digest.hexdigest()
