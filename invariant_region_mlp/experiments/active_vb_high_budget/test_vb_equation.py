"""Scientific and implementation tests for the active VB experiment."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from vb_equation import PUBLISHED_SIGMA, ViscousBurgersEquation  # noqa: E402
from vb_mlp_methods import (  # noqa: E402
    BATCH_BOX,
    RAW,
    SAMPLE_BALL,
    SAMPLE_BOX,
    Z_ONLY_BOX,
    FullHistoryMLP,
    terminal_ebl_weight,
)


class EquationTests(unittest.TestCase):
    def test_exact_analytic_pde_residual(self) -> None:
        rng = np.random.default_rng(7)
        for d in (2, 20, 80):
            equation = ViscousBurgersEquation(d)
            x = rng.uniform(-0.5, 0.5, size=(1000, d))
            t = rng.uniform(0.0, equation.T, size=1000)
            residual = equation.analytic_residual(t, x)
            self.assertLess(float(np.max(np.abs(residual))), 2e-13)

    def test_exact_numerical_finite_difference_residual(self) -> None:
        equation = ViscousBurgersEquation(5)
        rng = np.random.default_rng(71)
        x = rng.uniform(-0.35, 0.35, size=(12, equation.d))
        t = rng.uniform(0.1, 0.4, size=12)
        h = 2.0e-4
        u0 = equation.exact_u(t, x)
        ut = (equation.exact_u(t + h, x) - equation.exact_u(t - h, x)) / (2.0 * h)
        gradient = np.empty_like(x)
        laplacian = np.zeros(len(t))
        for coordinate in range(equation.d):
            xp = x.copy()
            xm = x.copy()
            xp[:, coordinate] += h
            xm[:, coordinate] -= h
            up = equation.exact_u(t, xp)
            um = equation.exact_u(t, xm)
            gradient[:, coordinate] = (up - um) / (2.0 * h)
            laplacian += (up - 2.0 * u0 + um) / h**2
        z = equation.sigma * gradient
        residual = (
            ut
            + equation.mu * np.sum(gradient, axis=1)
            + 0.5 * equation.sigma**2 * laplacian
            + equation.generator(u0, z)
        )
        self.assertLess(float(np.max(np.abs(residual))), 2e-7)

    def test_exact_state_satisfies_box_and_ball_certificates(self) -> None:
        rng = np.random.default_rng(8)
        for d in (20, 40, 60, 80):
            equation = ViscousBurgersEquation(d)
            x = rng.uniform(-3.0, 3.0, size=(500, d))
            t = rng.uniform(0.0, equation.T, size=500)
            state = equation.exact_state(t, x)
            self.assertTrue(np.all(state[:, 0] >= 0.0))
            self.assertTrue(np.all(state[:, 0] <= 1.0))
            self.assertTrue(np.all(state[:, 1:] >= 0.0))
            self.assertTrue(np.all(state[:, 1:] <= equation.z_upper + 2e-16))
            self.assertTrue(
                np.all(np.linalg.norm(state[:, 1:], axis=1) <= equation.z_ball_radius + 2e-15)
            )

    def test_published_parameters(self) -> None:
        equation = ViscousBurgersEquation(20)
        self.assertEqual(equation.sigma, PUBLISHED_SIGMA)
        self.assertAlmostEqual(equation.mu, -(1.0 / 20.0 + 1.0))
        self.assertAlmostEqual(equation.z_upper, math.sqrt(2.0) / 4.0)


class ProjectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.equation = ViscousBurgersEquation(4)
        self.bad = np.array(
            [
                [[-0.2, -1.0, 0.1, 2.0, 0.2], [1.4, 0.2, -0.3, 0.1, 5.0]],
                [[0.5, 1.0, 1.0, 1.0, 1.0], [0.1, -2.0, -1.0, -0.5, -0.1]],
            ],
            dtype=np.float64,
        )

    def solver(self, method) -> FullHistoryMLP:
        return FullHistoryMLP(self.equation, 2, method, np.random.default_rng(0))

    def test_box_projections_are_feasible(self) -> None:
        for method in (SAMPLE_BOX, Z_ONLY_BOX, BATCH_BOX):
            corrected = self.solver(method)._correct(self.bad)
            self.assertTrue(np.all(corrected[..., 1:] >= 0.0))
            self.assertTrue(np.all(corrected[..., 1:] <= self.equation.z_upper + 1e-15))
            if method != Z_ONLY_BOX:
                self.assertTrue(np.all(corrected[..., 0] >= 0.0))
                self.assertTrue(np.all(corrected[..., 0] <= 1.0))

    def test_ball_projection_is_feasible(self) -> None:
        corrected = self.solver(SAMPLE_BALL)._correct(self.bad)
        self.assertTrue(np.all(corrected[..., 0] >= 0.0))
        self.assertTrue(np.all(corrected[..., 0] <= 1.0))
        self.assertTrue(
            np.all(
                np.linalg.norm(corrected[..., 1:], axis=-1)
                <= self.equation.z_ball_radius + 1e-15
            )
        )

    def test_all_certified_projections_are_identity_on_exact_state(self) -> None:
        rng = np.random.default_rng(9)
        x = rng.normal(size=(6, 3, self.equation.d))
        t = rng.uniform(0.0, self.equation.T, size=(6, 3))
        exact = self.equation.exact_state(t, x)
        for method in (SAMPLE_BOX, Z_ONLY_BOX, SAMPLE_BALL, BATCH_BOX):
            corrected = self.solver(method)._correct(exact)
            np.testing.assert_allclose(corrected, exact, rtol=0.0, atol=2e-16)

    def test_batch_box_uses_one_upper_contraction_per_sibling_group(self) -> None:
        corrected = self.solver(BATCH_BOX)._correct(self.bad)
        positive = np.maximum(self.bad[..., 1:], 0.0)
        # Wherever both entries are positive, a common alpha preserves ratios.
        ratio0 = corrected[0, 0, 2] / positive[0, 0, 1]
        ratio1 = corrected[0, 1, 1] / positive[0, 1, 0]
        self.assertAlmostEqual(float(ratio0), float(ratio1))


class MLPTests(unittest.TestCase):
    def test_n0_and_n1_base_cases(self) -> None:
        equation = ViscousBurgersEquation(3)
        t = np.array([0.1, 0.2])
        x = np.arange(6, dtype=np.float64).reshape(2, 3) / 10.0
        zero_solver = FullHistoryMLP(equation, 3, RAW, np.random.default_rng(10))
        np.testing.assert_array_equal(zero_solver.solve(0, t, x), np.zeros((2, 4)))
        solver = FullHistoryMLP(equation, 3, RAW, np.random.default_rng(11))
        level_one = solver.solve(1, t, x)
        direct = FullHistoryMLP(equation, 3, RAW, np.random.default_rng(11))
        terminal = direct._terminal_estimate(1, t, x)
        np.testing.assert_array_equal(level_one, terminal)
        self.assertEqual(solver.stats.f_evals, 0)

    def test_corrected_terminal_ebl_scaling(self) -> None:
        normal = np.array([[[1.0, -2.0]], [[3.0, 4.0]]])
        remaining = np.array([[0.25], [0.04]])
        weight = terminal_ebl_weight(normal, remaining)
        np.testing.assert_allclose(weight[0], normal[0] / 0.5)
        np.testing.assert_allclose(weight[1], normal[1] / 0.2)
        # Explicitly distinguish the corrected normalization from the public
        # implementation's standard_normal/(T-t).
        self.assertNotAlmostEqual(float(weight[0, 0, 0]), 1.0 / 0.25)

    def test_paired_methods_consume_identical_random_trees(self) -> None:
        equation = ViscousBurgersEquation(3)
        t = np.array([0.05, 0.25])
        x = np.array([[0.1, -0.2, 0.3], [0.2, 0.1, -0.1]])
        fingerprints = []
        for method in (RAW, SAMPLE_BOX, BATCH_BOX):
            solver = FullHistoryMLP(
                equation,
                2,
                method,
                np.random.default_rng(12),
                trace_draws=True,
            )
            solver.solve(3, t, x)
            fingerprints.append(solver.draw_fingerprint)
        self.assertEqual(len(set(fingerprints)), 1)

    def test_more_terminal_samples_reduce_smoke_mse(self) -> None:
        equation = ViscousBurgersEquation(3)
        t = np.array([0.15])
        x = np.array([[0.2, -0.1, 0.05]])
        # A level-one estimate contains only the terminal expectation, hence
        # its reference is the f=0 linear-PDE solution rather than nonlinear u*.
        truth = equation.fzero_reference(t, x, quadrature_order=120)[0]
        mse: dict[int, float] = {}
        for M in (2, 32):
            errors = []
            for rep in range(96):
                solver = FullHistoryMLP(
                    equation,
                    M,
                    RAW,
                    np.random.default_rng(np.random.SeedSequence([13, rep])),
                )
                estimate = solver.solve(1, t, x)[0, 0]
                errors.append((estimate - truth) ** 2)
            mse[M] = float(np.mean(errors))
        self.assertLess(mse[32], 0.35 * mse[2])


if __name__ == "__main__":
    unittest.main(verbosity=2)
