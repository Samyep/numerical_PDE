"""Focused tests for the frozen mechanism-suite round-2 implementation."""

from __future__ import annotations

import math
from unittest import mock
import unittest

import numpy as np

from .equations import make_points
from .mechanism_mlp import BOX, RAW
from .round2_certificate import (
    ROUND2_BASE_SEED,
    p4_derivative_bounds,
    p4_tightness_interval,
)
from .round2_mlp import (
    _project_p4_tight,
    run_single_repetition_r2,
)
from .run_round2 import (
    ILLEGAL_R2_FACTORS,
    R1_DEEP_CONFIGS,
    R1_SHALLOW_CONFIGS,
    RESULTS_ROOT,
    primary_method,
    r1_tasks,
    r2_methods,
    r2_tasks,
    r3_tasks,
    r4_tasks,
    r5_tasks,
    r6_tasks,
)
from .run_suite import P4_REFERENCE, make_equation


class Round2CertificateTests(unittest.TestCase):
    def test_gauss_hermite_spot_values_match_preregistration(self) -> None:
        expected = {
            0.1: 0.12859138271504306,
            0.25: 0.23637707301493457,
            0.5: 0.3523427414907869,
        }
        for tau, magnitude in expected.items():
            lower, upper = p4_derivative_bounds(tau, 0.0)
            self.assertAlmostEqual(float(lower), -magnitude, places=13)
            self.assertAlmostEqual(float(upper), magnitude, places=13)

    def test_tightness_endpoints_are_registered_certificates(self) -> None:
        tau = np.array([0.1, 0.25, 0.5])
        s = np.array([-0.2, 0.0, 0.3])
        loose_low, loose_high = p4_tightness_interval(tau, s, 0.0)
        tight_low, tight_high = p4_tightness_interval(tau, s, 1.0)
        direct_low, direct_high = p4_derivative_bounds(tau, s)
        np.testing.assert_array_equal(loose_low, -np.ones(3))
        np.testing.assert_array_equal(loose_high, np.ones(3))
        np.testing.assert_allclose(tight_low, direct_low, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(tight_high, direct_high, rtol=0.0, atol=0.0)

    def test_projection_uses_interval_without_inward_margin(self) -> None:
        equation = make_equation({"pde_id": "P4", "d": 20, "kind": "p4"})
        state = np.zeros((1, 3, 21), dtype=np.float64)
        state[0, :, 1:] = (
            equation.sigma
            * np.array([-2.0, 0.1, 3.0])[:, None]
            * equation.w
        )
        time = np.full((1, 3), 0.25)
        x = np.zeros((1, 3, 20), dtype=np.float64)
        with mock.patch(
            "invariant_region_mlp.experiments.mechanism_suite.round2_mlp.p4_derivative_bounds_cached",
            return_value=(
                np.full((1, 3), -0.2),
                np.full((1, 3), 0.3),
            ),
        ):
            projected = _project_p4_tight(
                equation, state, time, x, theta=1.0
            )
        coefficient = projected[0, :, 1:] @ equation.w / equation.sigma
        np.testing.assert_allclose(coefficient, [-0.2, 0.1, 0.3], atol=2e-15)


class Round2DesignTests(unittest.TestCase):
    def test_task_counts_and_new_configurations_are_fixed(self) -> None:
        self.assertEqual(len(r1_tasks()), 4720)
        self.assertEqual(len(r2_tasks()), 6000)
        self.assertEqual(len(r3_tasks()), 200)
        self.assertEqual(len(r4_tasks()), 160)
        self.assertEqual(len(r5_tasks()), 1920)
        self.assertEqual(len(r6_tasks()), 250)
        self.assertEqual(
            set(R1_DEEP_CONFIGS + R1_SHALLOW_CONFIGS),
            {
                (4, 3),
                (4, 4),
                (4, 6),
                (5, 2),
                (3, 6),
                (3, 10),
                (3, 16),
                (2, 32),
            },
        )

    def test_illegal_point_nine_is_excluded(self) -> None:
        self.assertEqual(ILLEGAL_R2_FACTORS, (0.25, 0.5, 0.75))
        p4 = {"pde_id": "P4", "d": 20, "kind": "p4"}
        methods = r2_methods(p4)
        illegal = [item.factor for item in methods if "illegal" in item.transform]
        self.assertEqual(illegal, [0.25, 0.5, 0.75])

    def test_primary_method_is_tight_segment_only_for_p4(self) -> None:
        self.assertEqual(
            primary_method({"kind": "p4"}).name, "tight_segment"
        )
        self.assertEqual(primary_method({"kind": "p2"}).name, "box")

    def test_round2_output_root_is_separate(self) -> None:
        self.assertEqual(RESULTS_ROOT.name, "mechanism_suite_r2")
        self.assertNotEqual(RESULTS_ROOT.name, "mechanism_suite")

    def test_fresh_points_but_unchanged_direction(self) -> None:
        spec = {"pde_id": "P1", "d": 20, "kind": "p1"}
        first = make_equation(spec)
        second = make_equation(spec)
        np.testing.assert_array_equal(first.w, second.w)
        old_points = make_points(first, n_points=40, seed=20261006)
        new_points = make_points(first, n_points=40, seed=ROUND2_BASE_SEED)
        self.assertFalse(np.array_equal(old_points["t"], new_points["t"]))
        self.assertFalse(np.array_equal(old_points["x"], new_points["x"]))
        self.assertFalse(
            np.array_equal(old_points["is_validation"], new_points["is_validation"])
        )

    def test_tree_draws_remain_paired_across_methods(self) -> None:
        equation = make_equation({"pde_id": "P1", "d": 20, "kind": "p1"})
        points = make_points(equation, n_points=8, seed=ROUND2_BASE_SEED)
        fingerprints = []
        for method in (RAW, BOX):
            result = run_single_repetition_r2(
                pde_id="P1",
                equation=equation,
                method=method,
                n=1,
                M=2,
                repetition=0,
                t=points["t"],
                x=points["x"],
                is_validation=points["is_validation"],
                base_seed=ROUND2_BASE_SEED,
                chunk_size=4,
                study="unit_test",
                trace_draws=True,
            )
            fingerprints.append(result["metadata"]["draw_fingerprint"])
        self.assertIsNotNone(fingerprints[0])
        self.assertEqual(fingerprints[0], fingerprints[1])


if __name__ == "__main__":
    unittest.main()
