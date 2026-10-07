"""Scientific and pairing tests for the mechanism benchmark suite."""

from __future__ import annotations

import math
import unittest

import numpy as np

from .equations import BurgersFisher, MultiDirectionLSEHJB, RidgeLSEHJB, VBa
from .gates import finite_difference_residual, numerical_reference_finite_difference_residual
from .mechanism_mlp import (
    BALL,
    BATCH_BOX,
    BOX,
    RAW,
    SEGMENT,
    SIGN_ONLY,
    SPAN_ONLY,
    MechanismMLP,
    run_single_repetition,
)
from .norm_hjb import NormDriverHJB, load_norm_reference
from .run_suite import P4_REFERENCE, e5_tasks, e6_tasks, reuse_source_task


class AnalyticEquationTests(unittest.TestCase):
    def test_analytic_residuals(self) -> None:
        rng = np.random.default_rng(101)
        equations = (
            RidgeLSEHJB(7),
            VBa(7, a=4.0),
            VBa(7, a=8.0),
            BurgersFisher(7, rho=1.0),
            BurgersFisher(7, rho=2.0),
            MultiDirectionLSEHJB(d=12, n_directions=5, strength=0.5),
        )
        for equation in equations:
            x = rng.uniform(-0.5, 0.5, size=(128, equation.d))
            t = rng.uniform(0.0, equation.T, size=128)
            maximum = float(np.max(np.abs(equation.analytic_residual(t, x))))
            self.assertLess(maximum, 2e-13, equation.name)

    def test_fourth_order_finite_difference_residuals(self) -> None:
        for equation in (RidgeLSEHJB(5), VBa(5, a=4.0), BurgersFisher(5, rho=2.0)):
            result = finite_difference_residual(equation, n_points=8)
            self.assertLess(result["max_abs_residual"], 1e-6, equation.name)

    def test_terminal_consistency_and_exact_gradient(self) -> None:
        rng = np.random.default_rng(102)
        for equation in (RidgeLSEHJB(8), VBa(8, a=4.0), BurgersFisher(8, rho=1.0)):
            x = rng.uniform(-1.0, 1.0, size=(32, equation.d))
            np.testing.assert_allclose(
                equation.exact_u(equation.T, x), equation.terminal(x), rtol=0.0, atol=1e-14
            )
            t = rng.uniform(0.05, equation.T - 0.05, size=32)
            coordinate = 2
            step = 2e-5
            plus = x.copy()
            minus = x.copy()
            plus[:, coordinate] += step
            minus[:, coordinate] -= step
            numerical = (equation.exact_u(t, plus) - equation.exact_u(t, minus)) / (2 * step)
            exact = equation.exact_z(t, x)[:, coordinate] / equation.sigma
            np.testing.assert_allclose(numerical, exact, rtol=2e-8, atol=2e-9)


@unittest.skipUnless(P4_REFERENCE.exists(), "P4 reference has not been built")
class NumericalReferenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.equation = NormDriverHJB(20, str(P4_REFERENCE.resolve()), T=0.5)

    def test_refinement_and_reference_residual(self) -> None:
        metadata = load_norm_reference(str(P4_REFERENCE.resolve())).metadata
        self.assertTrue(metadata["refinement_passed"])
        self.assertLess(metadata["richardson_refinement_max_difference"], 1e-6)
        rng = np.random.default_rng(103)
        x = rng.uniform(-1.0, 1.0, size=(64, self.equation.d))
        t = rng.uniform(0.005, self.equation.T - 0.005, size=64)
        self.assertLess(float(np.max(np.abs(self.equation.reference_residual(t, x)))), 1e-6)

    def test_independent_finite_difference_residual(self) -> None:
        result = numerical_reference_finite_difference_residual(self.equation)
        self.assertLess(result["max_abs_residual"], 1e-6)

    def test_terminal_gradient_is_analytic(self) -> None:
        rng = np.random.default_rng(109)
        x = rng.uniform(-1.0, 1.0, size=(64, self.equation.d))
        s = x @ self.equation.w
        expected = np.tanh(self.equation.beta * s)[:, None] * self.equation.w
        actual = self.equation.exact_grad(self.equation.T, x)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2e-15)


class ProjectionAndPairingTests(unittest.TestCase):
    @staticmethod
    def solver(equation: object, method: object, seed: int = 1) -> MechanismMLP:
        return MechanismMLP(
            equation,
            2,
            method,
            np.random.default_rng(seed),
            dose_rng=np.random.default_rng(seed + 1000),
        )

    def test_box_is_identity_on_exact_states(self) -> None:
        rng = np.random.default_rng(104)
        for equation in (RidgeLSEHJB(8), VBa(8, a=4.0), BurgersFisher(8, rho=2.0)):
            x = rng.uniform(-1.0, 1.0, size=(5, 3, equation.d))
            t = rng.uniform(0.0, equation.T, size=(5, 3))
            exact = equation.exact_state(t, x)
            corrected = self.solver(equation, BOX)._transform(exact, exact)
            np.testing.assert_allclose(corrected, exact, rtol=0.0, atol=3e-15)

    def test_ridge_geometry_projections_are_feasible(self) -> None:
        equation = RidgeLSEHJB(8)
        rng = np.random.default_rng(107)
        state = rng.normal(size=(4, 5, equation.d + 1)) * 8.0
        segment = self.solver(equation, SEGMENT)._transform(state, state)
        coefficient = -np.einsum("...d,d->...", segment[..., 1:], equation.w) / equation.sigma
        low, high = equation.segment_coefficients
        self.assertTrue(np.all(coefficient >= low - 1e-14))
        self.assertTrue(np.all(coefficient <= high + 1e-14))
        parallel = -equation.sigma * coefficient[..., None] * equation.w
        np.testing.assert_allclose(segment[..., 1:], parallel, atol=2e-14, rtol=0.0)

        ball = self.solver(equation, BALL)._transform(state, state)
        self.assertTrue(
            np.all(np.linalg.norm(ball[..., 1:], axis=-1) <= equation.z_ball_radius + 2e-14)
        )
        span = self.solver(equation, SPAN_ONLY)._transform(state, state)
        parallel = np.einsum("...d,d->...", span[..., 1:], equation.w)[..., None] * equation.w
        np.testing.assert_allclose(span[..., 1:], parallel, atol=2e-14, rtol=0.0)

    def test_vb_sign_and_batch_rules(self) -> None:
        equation = VBa(6, a=4.0)
        rng = np.random.default_rng(108)
        state = rng.normal(size=(3, 4, equation.d + 1))
        signed = self.solver(equation, SIGN_ONLY)._transform(state, state)
        self.assertTrue(np.all((0.0 <= signed[..., 0]) & (signed[..., 0] <= 1.0)))
        self.assertTrue(np.all(signed[..., 1:] >= 0.0))

        batched = self.solver(equation, BATCH_BOX)._transform(state, state)
        self.assertTrue(np.all((0.0 <= batched[..., 0]) & (batched[..., 0] <= 1.0)))
        self.assertTrue(np.all(batched[..., 1:] >= 0.0))
        self.assertTrue(np.all(batched[..., 1:] <= equation.z_upper + 2e-15))
        positive = np.maximum(state[..., 1:], 0.0)
        maxima = np.max(positive, axis=(1, 2))
        alpha = np.minimum(1.0, equation.z_upper / np.maximum(maxima, 1e-300))
        np.testing.assert_allclose(
            batched[..., 1:], positive * alpha[:, None, None], rtol=0.0, atol=2e-15
        )

    def test_root_is_not_clipped_at_level_one(self) -> None:
        equation = VBa(4, a=4.0)
        t = np.array([0.1, 0.2], dtype=np.float64)
        x = np.array([[0.1, -0.2, 0.3, 0.0], [0.2, 0.1, -0.1, 0.4]])
        raw = self.solver(equation, RAW, 105).solve(1, t, x)
        box = self.solver(equation, BOX, 105).solve(1, t, x)
        np.testing.assert_array_equal(raw, box)

    def test_paired_methods_have_identical_tree_fingerprints(self) -> None:
        equation = VBa(4, a=4.0)
        rng = np.random.default_rng(106)
        t = rng.uniform(0.0, equation.T, size=12)
        x = rng.uniform(-0.5, 0.5, size=(12, equation.d))
        validation = np.zeros(12, dtype=bool)
        validation[:3] = True
        fingerprints = []
        for method in (RAW, BOX):
            result = run_single_repetition(
                pde_id="test",
                equation=equation,
                method=method,
                n=3,
                M=2,
                repetition=0,
                t=t,
                x=x,
                is_validation=validation,
                base_seed=20261006,
                chunk_size=4,
                trace_draws=True,
            )
            fingerprints.append(result["metadata"]["draw_fingerprint"])
        self.assertIsNotNone(fingerprints[0])
        self.assertEqual(fingerprints[0], fingerprints[1])

    def test_constant_predictor_skill_is_one(self) -> None:
        truth = np.linspace(-2.0, 3.0, 101)
        prediction = np.full_like(truth, np.mean(truth))
        skill = np.sqrt(np.mean((prediction - truth) ** 2)) / np.std(truth, ddof=0)
        self.assertEqual(float(skill), 1.0)


class RunnerReuseTests(unittest.TestCase):
    def test_only_identical_chunked_tasks_are_reused(self) -> None:
        e5_d100_rep0 = next(
            task
            for task in e5_tasks()
            if task["spec"]["d"] == 100
            and task["method"]["name"] == "centre"
            and task["repetition"] == 0
        )
        e5_d20_rep0 = next(
            task
            for task in e5_tasks()
            if task["spec"]["d"] == 20
            and task["method"]["name"] == "centre"
            and task["repetition"] == 0
        )
        self.assertEqual(reuse_source_task(e5_d100_rep0)["stage"], "e0")
        self.assertIsNone(reuse_source_task(e5_d20_rep0))

        e6_d20_rep0 = next(
            task
            for task in e6_tasks()
            if task["spec"]["d"] == 20
            and task["method"]["name"] == "centre"
            and task["repetition"] == 0
        )
        e6_d100_rep0 = next(
            task
            for task in e6_tasks()
            if task["spec"]["d"] == 100
            and task["method"]["name"] == "centre"
            and task["repetition"] == 0
        )
        e6_d100_rep3 = next(
            task
            for task in e6_tasks()
            if task["spec"]["d"] == 100
            and task["method"]["name"] == "centre"
            and task["repetition"] == 3
        )
        self.assertEqual(reuse_source_task(e6_d20_rep0)["stage"], "e3")
        self.assertEqual(reuse_source_task(e6_d100_rep0)["stage"], "e0")
        self.assertIsNone(reuse_source_task(e6_d100_rep3))


if __name__ == "__main__":
    unittest.main(verbosity=2)
