"""Unit and pathwise tests for the Allen--Cahn recovery experiment."""

from __future__ import annotations

import math
import unittest

import numpy as np

from .allen_cahn_equation import (
    AllenCahnEquation,
    IntervalStateProjector,
    KeyedRandomTree,
    reaction,
    sharp_clipped_driver_lipschitz,
    theorem_example_radius,
)
from .beck_truncated_mlp import BeckTruncatedMLP, RawFullHistoryMLP
from .ir_interval_mlp import IntervalIRMLP


class AllenCahnEquationTests(unittest.TestCase):
    def test_source_reaction_convention(self) -> None:
        values = np.array([-2.0, -0.5, 0.0, 0.5, 2.0], dtype=np.float64)
        np.testing.assert_array_equal(reaction(values), values - values**3)

    def test_source_terminal_and_diffusion_convention(self) -> None:
        equation = AllenCahnEquation(2)
        self.assertEqual(equation.horizon, 1.0)
        self.assertEqual(equation.terminal(np.zeros(2)), 0.5)
        point = equation.transition(
            np.zeros(2), 0.5, np.ones(2, dtype=np.float64)
        )
        # sigma=sqrt(2), so sqrt(2*0.5)=1.
        np.testing.assert_array_equal(point, np.ones(2))
        self.assertEqual(equation.terminal(point), 1.0 / 2.8)

    def test_scalar_projector(self) -> None:
        projector = IntervalStateProjector(-4.0, 4.0)
        self.assertEqual(projector.project_value(-9.0), -4.0)
        self.assertEqual(projector.project_value(2.0), 2.0)
        self.assertEqual(projector.project_value(9.0), 4.0)

    def test_product_projection_changes_only_u(self) -> None:
        projector = IntervalStateProjector(-1.0, 1.0)
        z = np.array([1.5, -2.0, 3.25], dtype=np.float64)
        projected_u, projected_z = projector.project_state(5.0, z)
        self.assertEqual(projected_u, 1.0)
        np.testing.assert_array_equal(projected_z, z)
        self.assertIsNot(projected_z, z)

    def test_exact_feasible_state_is_identity(self) -> None:
        projector = IntervalStateProjector(0.0, 1.0)
        z = np.array([-1.0, 2.0], dtype=np.float64)
        for value in (0.0, 0.25, 1.0):
            projected_u, projected_z = projector.project_state(value, z)
            self.assertEqual(projected_u, value)
            np.testing.assert_array_equal(projected_z, z)

    def test_theorem_radius_schedule(self) -> None:
        self.assertAlmostEqual(theorem_example_radius(2), math.log1p(math.log(2.0)))
        self.assertLess(theorem_example_radius(10), math.log(math.log(10.0)) + 1.0)


class ContainmentTests(unittest.TestCase):
    def test_direct_beck_driver_equals_generic_ir_driver(self) -> None:
        equation = AllenCahnEquation(1)
        tree = KeyedRandomTree(11, 1)
        beck = BeckTruncatedMLP(equation, 2, 4.0, tree)
        generic = IntervalIRMLP(equation, 2, (-4.0, 4.0), tree)
        for value in (-20.0, -4.0, -0.2, 0.0, 0.8, 4.0, 20.0):
            self.assertEqual(beck._driver(value), generic._projected_driver(value))

    def test_matched_keyed_random_trees(self) -> None:
        first = KeyedRandomTree(123456789, 10)
        second = KeyedRandomTree(123456789, 10)
        np.testing.assert_array_equal(
            first.normal((1, 2, 3), 7, (9, 10)),
            second.normal((1, 2, 3), 7, (9, 10)),
        )
        np.testing.assert_array_equal(
            first.uniform((4, 5), 8, 20), second.uniform((4, 5), 8, 20)
        )
        self.assertEqual(first.fingerprint, second.fingerprint)

    def test_pathwise_beck_ir_equality(self) -> None:
        for dimension in (1, 10):
            equation = AllenCahnEquation(dimension)
            point = np.linspace(0.0, 0.1, dimension, dtype=np.float64)
            for seed in (1, 99):
                tree = KeyedRandomTree(seed, dimension)
                beck = BeckTruncatedMLP(
                    equation, 3, 4.0, tree, capture_trace=True
                ).run(3, point)
                generic = IntervalIRMLP(
                    equation,
                    3,
                    (-4.0, 4.0),
                    tree,
                    capture_trace=True,
                ).run(3, point)
                self.assertEqual(beck.value, generic.value)
                np.testing.assert_array_equal(
                    beck.correction_trace, generic.correction_trace
                )
                self.assertEqual(beck.draw_fingerprint, generic.draw_fingerprint)

    def test_base_recursion(self) -> None:
        equation = AllenCahnEquation(3)
        tree = KeyedRandomTree(7, 3)
        point = np.zeros(3)
        for solver in (
            BeckTruncatedMLP(equation, 2, 4.0, tree),
            IntervalIRMLP(equation, 2, (-4.0, 4.0), tree),
            RawFullHistoryMLP(equation, 2, 4.0, tree),
        ):
            result = solver.run(0, point)
            self.assertEqual(result.value, 0.0)
            self.assertEqual(result.work["recursive_states"], 1)
            self.assertEqual(result.work["terminal_g_evals"], 0)
            self.assertEqual(result.work["f_evals"], 0)

    def test_raw_has_no_hidden_clipping(self) -> None:
        equation = AllenCahnEquation(1)
        tree = KeyedRandomTree(5, 1)
        raw = RawFullHistoryMLP(equation, 2, 4.0, tree)
        beck = BeckTruncatedMLP(equation, 2, 4.0, tree)
        self.assertEqual(raw._driver(5.0), -120.0)
        self.assertEqual(beck._driver(5.0), -60.0)

    def test_float64_outputs_are_finite(self) -> None:
        equation = AllenCahnEquation(10)
        point = np.zeros(10, dtype=np.float64)
        tree = KeyedRandomTree(2020, 10)
        for solver in (
            RawFullHistoryMLP(equation, 2, 4.0, tree),
            BeckTruncatedMLP(equation, 2, 4.0, tree),
            IntervalIRMLP(equation, 2, (-4.0, 4.0), tree),
        ):
            result = solver.run(3, point)
            self.assertTrue(math.isfinite(result.value))
            self.assertIsInstance(result.value, float)

    def test_clipped_driver_lipschitz_bound(self) -> None:
        grid = np.linspace(-8.0, 8.0, 20001)
        for radius in (1.0, 4.0):
            projector = IntervalStateProjector(-radius, radius)
            transformed = reaction(np.clip(grid, projector.low, projector.high))
            slopes = np.abs(np.diff(transformed) / np.diff(grid))
            self.assertLessEqual(
                float(np.max(slopes)), sharp_clipped_driver_lipschitz(radius) + 1e-10
            )


if __name__ == "__main__":
    unittest.main()
