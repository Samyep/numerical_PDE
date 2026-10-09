from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import pytest

from invariant_region_mlp.experiments.double_estimator.equations import make_equation
from invariant_region_mlp.experiments.double_estimator.solver import (
    DoubleEstimatorMLP,
    METHODS as FROZEN_METHODS,
)

from .protocol import EXPECTED_NEW_ROWS, EXPECTED_ROWS, all_tasks, reusable_source
from .solver import CorrectionBaselineMLP, METHODS


def _solver(cls, equation, method, seed_parts, M=4, trace=False):
    return cls(
        equation,
        M,
        method,
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        trace_draws=trace,
    )


@pytest.mark.parametrize("method_name", ["raw", "path", "double", "double_path"])
def test_nonbaseline_subclass_matches_frozen_solver(method_name: str) -> None:
    equation = make_equation("P1", 12)
    t = np.array([0.03, 0.17], dtype=np.float64)
    x = np.linspace(-0.8, 0.9, 24, dtype=np.float64).reshape(2, 12)
    seed_parts = [20261201, 12, 3, 4, 7, 0]
    frozen = _solver(
        DoubleEstimatorMLP, equation, FROZEN_METHODS[method_name], seed_parts, trace=True
    )
    baseline = _solver(
        CorrectionBaselineMLP, equation, METHODS[method_name], seed_parts, trace=True
    )
    expected = frozen.solve(3, t, x)
    observed = baseline.solve(3, t, x)
    np.testing.assert_array_equal(observed, expected)
    assert baseline.draw_fingerprint == frozen.draw_fingerprint
    assert baseline.auxiliary_draw_fingerprint == frozen.auxiliary_draw_fingerprint


def test_constant_baseline_preserves_mean_within_five_standard_errors() -> None:
    rng = np.random.default_rng(20261201)
    parents, siblings, d = 4000, 32, 3
    remaining = rng.uniform(0.02, 0.25, (parents, 1))
    uniform = rng.random((parents, siblings))
    r = np.maximum(uniform, np.finfo(np.float64).tiny) ** 2
    dt = remaining * r
    normal = rng.standard_normal((parents, siblings, d))
    driver = rng.normal(1.7, 0.8, (parents, siblings))
    importance = 2.0 * remaining * np.sqrt(r)
    ebl = normal / np.sqrt(dt)[:, :, None]
    constant = 3.25
    raw = np.mean(importance[:, :, None] * driver[:, :, None] * ebl, axis=1)
    centred = np.mean(
        importance[:, :, None] * (driver - constant)[:, :, None] * ebl,
        axis=1,
    )
    paired = centred - raw
    mean = np.mean(paired, axis=0)
    standard_error = np.std(paired, axis=0, ddof=1) / math.sqrt(parents)
    assert np.all(np.abs(mean) <= 5.0 * standard_error)


@dataclass(frozen=True)
class LinearTerminalHJB:
    d: int = 100
    T: float = 0.5
    sigma: float = math.sqrt(2.0)
    mu: float = 0.0
    family: str = "linear_test"
    name: str = "linear_terminal_test"

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, dtype=np.float64)[..., 0]

    def exact_u(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, dtype=np.float64)[..., 0] - (self.T - np.asarray(t))

    def exact_z(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        del t
        z = np.zeros_like(np.asarray(x, dtype=np.float64))
        z[..., 0] = self.sigma
        return z

    def exact_state(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        return np.concatenate(
            (self.exact_u(t, x)[..., None], self.exact_z(t, x)), axis=-1
        )

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        return -0.5 * np.sum(np.asarray(z, dtype=np.float64) ** 2, axis=-1)

    @property
    def box_low(self) -> np.ndarray:
        return np.zeros(self.d, dtype=np.float64)

    @property
    def box_high(self) -> np.ndarray:
        return np.full(self.d, self.sigma, dtype=np.float64)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma


def test_linear_terminal_level1_increment_is_exactly_zero_at_every_node() -> None:
    equation = LinearTerminalHJB()
    seed_parts = [20261201, equation.d, 3, 4, 0, 0]
    solver = _solver(
        CorrectionBaselineMLP,
        equation,
        METHODS["path_base"],
        seed_parts,
        M=4,
    )
    t = np.array([0.0, 0.11, 0.23], dtype=np.float64)
    x = np.zeros((len(t), equation.d), dtype=np.float64)
    solver.solve(3, t, x, collect_root=True)
    assert solver.max_abs_level1_gradient_integrand == 0.0
    assert solver.root_level1_gradient_increment is not None
    np.testing.assert_array_equal(
        solver.root_level1_gradient_increment,
        np.zeros_like(solver.root_level1_gradient_increment),
    )


@pytest.mark.parametrize(
    ("baseline_method", "counterpart"),
    [
        ("path_base", "path"),
        ("path_base_oracle", "path"),
        ("double_path_base", "double_path"),
    ],
)
def test_baseline_consumes_no_draws(baseline_method: str, counterpart: str) -> None:
    equation = make_equation("MR", 20)
    t = np.array([0.02], dtype=np.float64)
    x = np.linspace(-1.0, 1.0, 20, dtype=np.float64)[None, :]
    seed_parts = [20261201, 20, 3, 4, 2, 0]
    baseline = _solver(
        CorrectionBaselineMLP, equation, METHODS[baseline_method], seed_parts, trace=True
    )
    plain = _solver(
        CorrectionBaselineMLP, equation, METHODS[counterpart], seed_parts, trace=True
    )
    baseline.solve(3, t, x)
    plain.solve(3, t, x)
    assert baseline.draw_fingerprint == plain.draw_fingerprint
    assert baseline.auxiliary_draw_fingerprint == plain.auxiliary_draw_fingerprint
    assert baseline.stats.standard_normal_variates == plain.stats.standard_normal_variates
    assert baseline.stats.time_uniform_variates == plain.stats.time_uniform_variates


def test_registered_task_counts_and_reuse() -> None:
    tasks = all_tasks()
    assert len(tasks) == EXPECTED_ROWS == 880
    new = [task for task in tasks if reusable_source(task) is None]
    assert len(new) == EXPECTED_NEW_ROWS == 600

