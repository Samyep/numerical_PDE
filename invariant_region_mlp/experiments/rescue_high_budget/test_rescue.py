"""Unit and paired-stream tests for the high-budget rescue implementations."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from funding_rescue import (  # noqa: E402
    FundingEquation,
    FundingFullHistoryMLP,
    RAW as FUNDING_RAW,
    SAMPLEWISE as FUNDING_SAMPLEWISE,
)
from hjb_rescue import (  # noqa: E402
    CERTIFIED_RADIUS,
    HJBFullHistoryMLP,
    RAW as HJB_RAW,
    RosenbrockHJB,
    SAMPLEWISE as HJB_SAMPLEWISE,
)


class FundingTests(unittest.TestCase):
    def test_published_protocol_parameters(self) -> None:
        equation = FundingEquation()
        self.assertEqual(equation.d, 100)
        self.assertEqual(equation.T, 0.5)
        self.assertEqual(equation.sigma, 0.2)
        self.assertEqual(equation.mu, 0.06)

    def test_projection_is_feasible_and_identity_on_feasible_state(self) -> None:
        equation = FundingEquation()
        t = np.asarray([0.0, 0.25])
        x = np.full((2, equation.d), 100.0)
        delta = np.zeros_like(x)
        delta[:, 0] = 0.5 * equation.delta_radius(t)
        state = np.zeros((2, equation.d + 1))
        state[:, 0] = [20.0, 21.0]
        state[:, 1:] = equation.sigma * x * delta
        solver = FundingFullHistoryMLP(
            M=2, method=FUNDING_SAMPLEWISE, seed=7
        )
        corrected = solver._correct(state, t, x)
        np.testing.assert_array_equal(corrected, state)

        bad = state.copy()
        bad[:, 1] *= 4.0
        corrected = solver._correct(bad, t, x)
        norm = np.linalg.norm(
            corrected[:, 1:] / (equation.sigma * x), axis=1
        )
        self.assertTrue(np.all(norm <= equation.delta_radius(t) * (1.0 + 1e-14)))

    def test_corrected_terminal_ebl_and_paired_tree(self) -> None:
        equation = FundingEquation()
        t = np.zeros(3)
        x = np.full((3, equation.d), 100.0)
        raw = FundingFullHistoryMLP(M=2, method=FUNDING_RAW, seed=11)
        projected = FundingFullHistoryMLP(
            M=2, method=FUNDING_SAMPLEWISE, seed=11
        )
        raw_state = raw.solve(2, t, x, collect_root=True)
        projected_state = projected.solve(2, t, x, collect_root=True)
        self.assertEqual(raw.draws.fingerprint, projected.draws.fingerprint)
        np.testing.assert_array_equal(raw.root_terminal, projected.root_terminal)
        self.assertTrue(np.all(np.isfinite(raw_state)))
        self.assertTrue(np.all(np.isfinite(projected_state)))

        seed = 23
        solver = FundingFullHistoryMLP(M=2, method=FUNDING_RAW, seed=seed)
        estimate = solver._terminal_estimate(1, t[:1], x[:1])
        rng = np.random.default_rng(seed)
        normal = rng.standard_normal((1, 2, equation.d))
        remaining = equation.T - t[:1]
        brownian = np.sqrt(remaining)[:, None, None] * normal
        xt = x[:1, None, :] * np.exp(
            (equation.mu - 0.5 * equation.sigma**2)
            * remaining[:, None, None]
            + equation.sigma * brownian
        )
        difference = equation.terminal(xt) - equation.terminal(x[:1])[:, None]
        expected_z = np.mean(
            difference[:, :, None]
            * normal
            / np.sqrt(remaining)[:, None, None],
            axis=1,
        )
        np.testing.assert_allclose(estimate[:, 1:], expected_z, rtol=0.0, atol=0.0)

    def test_base_cases(self) -> None:
        equation = FundingEquation()
        t = np.zeros(1)
        x = np.full((1, equation.d), 100.0)
        solver = FundingFullHistoryMLP(M=2, method=FUNDING_RAW, seed=3)
        np.testing.assert_array_equal(solver.solve(0, t, x), np.zeros((1, 101)))
        one = solver.solve(1, t, x)
        self.assertEqual(one.shape, (1, 101))
        self.assertEqual(solver.work.f_evals, 0)


class HJBTests(unittest.TestCase):
    def test_existing_matrix_construction(self) -> None:
        equation = RosenbrockHJB(100)
        self.assertAlmostEqual(equation.trace, 304.8856089115143, places=12)
        self.assertAlmostEqual(equation.lambda_max, 6.383491595137865, places=12)

    def test_exact_state_is_certified_and_projection_identity(self) -> None:
        equation = RosenbrockHJB(100)
        rng = np.random.default_rng(17)
        x = rng.normal(size=(12, 100))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        t = rng.random(12)
        u, z = equation.exact_state(t, x, n_quad=48)
        self.assertTrue(np.all(np.linalg.norm(z, axis=1) <= CERTIFIED_RADIUS))
        state = np.column_stack([u, z])
        solver = HJBFullHistoryMLP(
            equation, M=2, method=HJB_SAMPLEWISE, seed=5
        )
        corrected = solver._correct(state)
        np.testing.assert_allclose(corrected, state, rtol=0.0, atol=0.0)

    def test_hjb_projection_and_paired_tree(self) -> None:
        equation = RosenbrockHJB(20)
        state = np.zeros((3, 21))
        state[:, 1:] = 10.0
        solver = HJBFullHistoryMLP(
            equation, M=2, method=HJB_SAMPLEWISE, seed=9
        )
        corrected = solver._correct(state)
        self.assertTrue(
            np.all(np.linalg.norm(corrected[:, 1:], axis=1) <= CERTIFIED_RADIUS)
        )

        rng = np.random.default_rng(8)
        x = rng.normal(size=(3, 20)) * 0.1
        t = rng.random(3) * 0.5
        raw = HJBFullHistoryMLP(equation, M=2, method=HJB_RAW, seed=13)
        projected = HJBFullHistoryMLP(
            equation, M=2, method=HJB_SAMPLEWISE, seed=13
        )
        raw_state = raw.solve(2, t, x, collect_root=True)
        projected_state = projected.solve(2, t, x, collect_root=True)
        self.assertEqual(raw.draws.fingerprint, projected.draws.fingerprint)
        np.testing.assert_array_equal(raw.root_terminal, projected.root_terminal)
        self.assertTrue(np.all(np.isfinite(raw_state)))
        self.assertTrue(np.all(np.isfinite(projected_state)))

    def test_hjb_corrected_ebl_and_base_cases(self) -> None:
        equation = RosenbrockHJB(10)
        t = np.asarray([0.2])
        x = np.zeros((1, 10))
        solver = HJBFullHistoryMLP(equation, M=2, method=HJB_RAW, seed=29)
        np.testing.assert_array_equal(solver.solve(0, t, x), np.zeros((1, 11)))

        solver = HJBFullHistoryMLP(equation, M=2, method=HJB_RAW, seed=29)
        estimate = solver._terminal_estimate(1, t, x)
        rng = np.random.default_rng(29)
        normal = rng.standard_normal((1, 2, 10))
        remaining = equation.T - t
        xt = x[:, None, :] + math.sqrt(2.0) * np.sqrt(remaining)[
            :, None, None
        ] * normal
        gt = equation.terminal(xt)
        expected_z = np.mean(
            (gt - equation.terminal(x)[:, None])[:, :, None]
            * normal
            / np.sqrt(remaining)[:, None, None],
            axis=1,
        )
        np.testing.assert_allclose(estimate[:, 1:], expected_z, rtol=0.0, atol=0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
