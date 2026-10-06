"""Generic Samplewise interval IR implementation for Allen--Cahn MLP."""

from __future__ import annotations

from typing import Sequence
import time

import numpy as np

from .allen_cahn_equation import (
    AllenCahnEquation,
    IntervalStateProjector,
    KeyedRandomTree,
    WorkDiagnostics,
    child_path,
    reaction,
)
from .beck_truncated_mlp import MLPRun


class IntervalIRMLP:
    """Samplewise IR using Pi_[a,b]xR^d before every nonlinear reuse.

    The recursion is intentionally implemented independently of
    ``BeckTruncatedMLP``.  Only the immutable keyed random tree and equation
    primitives are shared.  Because the Allen--Cahn driver is z-independent,
    the implementation only needs the scalar component of the generic product
    projector during the recursion; ``project_state`` separately verifies that
    every gradient component is unchanged.
    """

    method_name = "interval_ir"

    def __init__(
        self,
        equation: AllenCahnEquation,
        sample_size: int,
        interval: tuple[float, float],
        random_tree: KeyedRandomTree,
        *,
        capture_trace: bool = False,
        method_name: str | None = None,
    ) -> None:
        if sample_size < 1:
            raise ValueError("sample_size must be positive")
        self.equation = equation
        self.sample_size = int(sample_size)
        self.projector = IntervalStateProjector(*interval)
        self.random_tree = random_tree
        self.capture_trace = capture_trace
        self.stats = WorkDiagnostics(capture_trace=capture_trace)
        if method_name is not None:
            self.method_name = method_name

    def _projected_driver(self, value: float) -> float:
        before = float(value)
        projected = self.projector.project_value(before)
        generator_before = float(reaction(before))
        generator_after = float(reaction(projected))
        self.stats.observe_driver(
            before,
            projected,
            self.projector,
            generator_before,
            generator_after,
        )
        return generator_after

    def _evaluate(
        self,
        time_now: float,
        position: np.ndarray,
        picard_depth: int,
        node_key: Sequence[int],
    ) -> float:
        self.stats.recursive_states += 1
        if picard_depth == 0:
            return 0.0

        duration = self.equation.horizon - float(time_now)
        leaf_count = self.sample_size**picard_depth
        leaf_noise = self.random_tree.normal(
            node_key, 1, (leaf_count, self.equation.dimension)
        )
        self.stats.normal_scalar_draws += leaf_noise.size
        self.stats.state_transitions += leaf_count
        leaves = self.equation.transition(position, duration, leaf_noise)
        leaf_payoffs = self.equation.terminal(leaves)
        self.stats.terminal_g_evals += leaf_count
        estimate = float(np.sum(leaf_payoffs, dtype=np.float64) / leaf_count)

        for picard_level in range(picard_depth):
            batch_size = self.sample_size ** (picard_depth - picard_level)
            random_times = self.random_tree.uniform(
                node_key, 2 + 2 * picard_level, batch_size
            )
            brownian_noise = self.random_tree.normal(
                node_key,
                3 + 2 * picard_level,
                (batch_size, self.equation.dimension),
            )
            self.stats.uniform_draws += batch_size
            self.stats.normal_scalar_draws += brownian_noise.size
            self.stats.state_transitions += batch_size
            evaluation_times = time_now + duration * random_times
            evaluation_points = self.equation.transition(
                position, evaluation_times - time_now, brownian_noise
            )
            level_sum = np.float64(0.0)
            for index in range(batch_size):
                fine_value = self._evaluate(
                    float(evaluation_times[index]),
                    evaluation_points[index],
                    picard_level,
                    child_path(node_key, picard_level, index, 0),
                )
                increment = self._projected_driver(fine_value)
                if picard_level:
                    coarse_value = self._evaluate(
                        float(evaluation_times[index]),
                        evaluation_points[index],
                        picard_level - 1,
                        child_path(node_key, picard_level, index, 1),
                    )
                    increment -= self._projected_driver(coarse_value)
                self.stats.observe_correction(increment)
                level_sum += np.float64(increment)
            estimate += duration * float(level_sum / np.float64(batch_size))
        return float(estimate)

    def run(
        self, depth: int, state: np.ndarray, current_time: float = 0.0
    ) -> MLPRun:
        if depth < 0:
            raise ValueError("depth must be nonnegative")
        x = np.asarray(state, dtype=np.float64)
        if x.shape != (self.equation.dimension,):
            raise ValueError(
                f"state must have shape {(self.equation.dimension,)}, got {x.shape}"
            )
        self.stats = WorkDiagnostics(capture_trace=self.capture_trace)
        started = time.perf_counter()
        value = self._evaluate(float(current_time), x, int(depth), ())
        elapsed = time.perf_counter() - started
        return MLPRun(
            value=value,
            wall_clock_seconds=elapsed,
            work=self.stats.as_dict(),
            correction_trace=np.asarray(self.stats.correction_trace, dtype=np.float64),
            draw_fingerprint=self.random_tree.fingerprint,
        )
