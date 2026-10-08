from __future__ import annotations

import numpy as np

from invariant_region_mlp.experiments.candidate_screen.candidates import LQGameHJB
from invariant_region_mlp.experiments.mechanism_suite.equations import RidgeLSEHJB

from .equations import MultiRidgeHJB
from .protocol import all_tasks, task_identity
from .solver import DoubleEstimatorMLP, METHODS


def _solver(equation, method: str, *, trace: bool = False) -> DoubleEstimatorMLP:
    seed = [20261201, equation.d, 2, 2, 0, 0]
    return DoubleEstimatorMLP(
        equation,
        2,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed)),
        aux_seed_sequence=np.random.SeedSequence(seed + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed + [0xD05E])),
        trace_draws=trace,
    )


def test_quadratic_cross_estimator_removes_jensen_bias() -> None:
    rng = np.random.default_rng(123)
    truth_z = np.array([0.5, -0.25, 0.75])
    first = truth_z + rng.standard_normal((400_000, 3))
    second = truth_z + rng.standard_normal((400_000, 3))
    truth = -0.5 * np.dot(truth_z, truth_z)
    raw_bias = float(np.mean(-0.5 * np.sum(first * first, axis=1)) - truth)
    double_bias = float(np.mean(-0.5 * np.sum(first * second, axis=1)) - truth)
    assert raw_bias < -1.45
    assert abs(double_bias) < 0.01


def test_game_cross_formula() -> None:
    equation = LQGameHJB(d=8, a=3.0, b=1.0, frac_A=0.25)
    solver = _solver(equation, "double")
    first = np.zeros((1, 1, 9))
    second = np.zeros((1, 1, 9))
    first[..., 1:] = np.arange(1.0, 9.0)
    second[..., 1:] = np.arange(2.0, 10.0)
    expected = (
        -1.5 * np.sum(first[..., 1:3] * second[..., 1:3], axis=-1)
        + 0.5 * np.sum(first[..., 3:] * second[..., 3:], axis=-1)
    )
    np.testing.assert_allclose(solver._double_value(first, second), expected)


def test_primary_draw_stream_is_paired_with_raw() -> None:
    equation = RidgeLSEHJB(d=3)
    t = np.array([0.03, 0.12])
    x = np.zeros((2, 3))
    raw = _solver(equation, "raw", trace=True)
    double = _solver(equation, "double", trace=True)
    # n=3 exercises auxiliary calls nested inside the primary recursive tree.
    raw.solve(3, t, x)
    double.solve(3, t, x)
    assert raw.draw_fingerprint == double.draw_fingerprint
    assert double.auxiliary_draw_fingerprint is not None
    assert double.auxiliary_recursive_invocations > 0
    assert double.recursive_invocations > raw.recursive_invocations


def test_pathwise_terminal_preserves_value_draws() -> None:
    equation = RidgeLSEHJB(d=4)
    t = np.array([0.02, 0.08])
    x = np.zeros((2, 4))
    raw = _solver(equation, "raw", trace=True)
    path = _solver(equation, "path", trace=True)
    raw_state = raw._terminal_estimate(2, t, x)
    path_state = path._terminal_estimate(2, t, x)
    np.testing.assert_allclose(raw_state[:, 0], path_state[:, 0])
    assert raw.draw_fingerprint == path.draw_fingerprint
    assert np.all(np.isfinite(path_state))


def test_multiridge_subbox_projection_is_feasible() -> None:
    equation = MultiRidgeHJB(d=12, n_grad_samples=100)
    solver = _solver(equation, "sub_box")
    rng = np.random.default_rng(7)
    state = rng.standard_normal((2, 3, equation.d + 1)) * 20.0
    exact = equation.exact_state(np.zeros((2, 3)), np.zeros((2, 3, equation.d)))
    corrected = solver._transform(state, exact)
    coordinates = corrected[..., 1:] @ equation.Qh
    assert np.all(coordinates >= equation.sub_lo - 1e-12)
    assert np.all(coordinates <= equation.sub_hi + 1e-12)
    orthogonal = corrected[..., 1:] - coordinates @ equation.Qh.T
    assert np.max(np.abs(orthogonal)) < 1e-10


def test_task_matrix_is_unique_and_complete() -> None:
    tasks = all_tasks()
    identities = [task_identity(task) for task in tasks]
    assert len(tasks) == 6080
    assert len(identities) == len(set(identities))
    assert {task["pde_id"] for task in tasks} == {
        "P1",
        "MR",
        "C2-convex",
        "C2-cancel",
        "C2-flip",
        "P4",
    }
