"""Direct Beck-style truncated full-history recursive MLP implementation."""

from __future__ import annotations

from dataclasses import dataclass
import math
import time
from typing import Sequence

import numpy as np

from .allen_cahn_equation import (
    AllenCahnEquation,
    IntervalStateProjector,
    KeyedRandomTree,
    WorkDiagnostics,
    child_path,
)


@dataclass(frozen=True)
class MLPRun:
    value: float
    wall_clock_seconds: float
    work: dict[str, float | int]
    correction_trace: np.ndarray
    draw_fingerprint: str


def _direct_cubic(value: float) -> float:
    """Local copy of f(u)=u-u^3 for an independent Beck code path."""

    u = np.float64(value)
    with np.errstate(over="ignore", invalid="ignore"):
        return float(u - u * u * u)


class BeckTruncatedMLP:
    """Equation (3) of Becker et al. with direct scalar truncation.

    This implementation deliberately computes ``f(P_r(u))`` locally instead
    of calling the generic interval projector used by ``IntervalIRMLP``.
    """

    method_name = "beck_truncated"

    def __init__(
        self,
        equation: AllenCahnEquation,
        sample_size: int,
        radius: float,
        random_tree: KeyedRandomTree,
        *,
        capture_trace: bool = False,
    ) -> None:
        if sample_size < 1:
            raise ValueError("sample_size must be positive")
        if radius <= 0.0:
            raise ValueError("radius must be positive")
        self.equation = equation
        self.sample_size = int(sample_size)
        self.radius = float(radius)
        self.random_tree = random_tree
        self.interval = IntervalStateProjector(-self.radius, self.radius)
        self.capture_trace = capture_trace
        self.stats = WorkDiagnostics(capture_trace=capture_trace)

    def _driver(self, value: float) -> float:
        before = float(value)
        after = min(self.radius, max(-self.radius, before))
        generator_before = _direct_cubic(before)
        generator_after = _direct_cubic(after)
        self.stats.observe_driver(
            before, after, self.interval, generator_before, generator_after
        )
        return generator_after

    def _solve(
        self,
        current_time: float,
        state: np.ndarray,
        depth: int,
        path: Sequence[int],
    ) -> float:
        self.stats.recursive_states += 1
        if depth <= 0:
            return 0.0

        remaining = self.equation.horizon - float(current_time)
        terminal_count = self.sample_size**depth
        terminal_normals = self.random_tree.normal(
            path, 1, (terminal_count, self.equation.dimension)
        )
        self.stats.normal_scalar_draws += terminal_normals.size
        self.stats.state_transitions += terminal_count
        terminal_states = self.equation.transition(
            state, remaining, terminal_normals
        )
        terminal_values = self.equation.terminal(terminal_states)
        self.stats.terminal_g_evals += terminal_count
        result = float(np.mean(terminal_values, dtype=np.float64))

        for level in range(depth):
            count = self.sample_size ** (depth - level)
            uniforms = self.random_tree.uniform(path, 2 + 2 * level, count)
            normals = self.random_tree.normal(
                path,
                3 + 2 * level,
                (count, self.equation.dimension),
            )
            self.stats.uniform_draws += count
            self.stats.normal_scalar_draws += normals.size
            self.stats.state_transitions += count
            child_times = current_time + remaining * uniforms
            child_states = self.equation.transition(
                state, child_times - current_time, normals
            )
            correction_sum = np.float64(0.0)
            for sample in range(count):
                fine = self._solve(
                    float(child_times[sample]),
                    child_states[sample],
                    level,
                    child_path(path, level, sample, 0),
                )
                correction = self._driver(fine)
                if level > 0:
                    coarse = self._solve(
                        float(child_times[sample]),
                        child_states[sample],
                        level - 1,
                        child_path(path, level, sample, 1),
                    )
                    correction -= self._driver(coarse)
                self.stats.observe_correction(correction)
                correction_sum += np.float64(correction)
            result += remaining * float(correction_sum / np.float64(count))
        return float(result)

    def run(self, depth: int, state: np.ndarray, current_time: float = 0.0) -> MLPRun:
        if depth < 0:
            raise ValueError("depth must be nonnegative")
        x = np.asarray(state, dtype=np.float64)
        if x.shape != (self.equation.dimension,):
            raise ValueError(
                f"state must have shape {(self.equation.dimension,)}, got {x.shape}"
            )
        self.stats = WorkDiagnostics(capture_trace=self.capture_trace)
        started = time.perf_counter()
        value = self._solve(float(current_time), x, int(depth), ())
        elapsed = time.perf_counter() - started
        return MLPRun(
            value=value,
            wall_clock_seconds=elapsed,
            work=self.stats.as_dict(),
            correction_trace=np.asarray(self.stats.correction_trace, dtype=np.float64),
            draw_fingerprint=self.random_tree.fingerprint,
        )


class RawFullHistoryMLP(BeckTruncatedMLP):
    """Untruncated full-history MLP on the same exact random tree."""

    method_name = "raw"

    def _driver(self, value: float) -> float:
        before = float(value)
        generator = _direct_cubic(before)
        self.stats.observe_driver(
            before, before, self.interval, generator, generator
        )
        return generator
