from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import BOX, RAW

from .equations import (
    BASE_SEED,
    LQGEquation,
    ProductionCubicHJB,
    make_equation,
    make_points,
)
from .mlp import _generic_segment, _lqg_geometry, run_single_repetition
from .protocol import S1_CONFIGS, regular_methods, s0_tasks, s1_tasks
from .reference_c3 import build_and_audit_reference


class Round3Tests(unittest.TestCase):
    def test_registered_matrix_and_fresh_seed(self) -> None:
        self.assertEqual(BASE_SEED, 20261207)
        self.assertEqual(len(S1_CONFIGS), 19)
        self.assertEqual(len(s0_tasks()), 300)
        identities = {
            (
                task["pde_id"], task["d"], task["n"], task["M"],
                task["method"]["name"], task["repetition"],
            )
            for task in s1_tasks()
        }
        self.assertEqual(len(identities), len(s1_tasks()))

    def test_c2_prediction_signs_and_exact_residual(self) -> None:
        for d in (20, 100, 400):
            convex = make_equation("C2-convex", d)
            cancel = make_equation("C2-cancel", d)
            flip = make_equation("C2-flip", d)
            self.assertLess(convex.predicted_orthogonal_bias_factor(), 0.0)
            self.assertLessEqual(abs(cancel.predicted_orthogonal_bias_factor()), 1.5)
            self.assertGreater(flip.predicted_orthogonal_bias_factor(), 0.0)
            points = make_points(convex, pde_id="C2-convex", n_points=12)
            self.assertLess(
                np.max(np.abs(convex.analytic_residual(points["t"], points["x"]))),
                1e-12,
            )

    def test_lqg_formula_and_dynamic_geometries(self) -> None:
        equation = LQGEquation(d=10)
        points = make_points(equation, pde_id="LQG", n_points=64)
        self.assertLess(
            np.max(np.abs(equation.analytic_residual(points["t"], points["x"]))),
            1e-12,
        )
        self.assertLess(
            np.max(
                np.abs(
                    equation.exact_u(equation.T, points["x"])
                    - equation.terminal(points["x"])
                )
            ),
            1e-12,
        )
        exact = equation.exact_state(points["t"], points["x"])
        noisy = exact.copy()
        noisy[:, 1:] += 10.0
        for transform in ("segment", "box", "ball", "span_only", "centre"):
            corrected = _lqg_geometry(equation, noisy, points["x"], transform)
            self.assertTrue(np.all(np.isfinite(corrected)))
        segment = _lqg_geometry(equation, noisy, points["x"], "segment")
        x = points["x"]
        coefficient = np.sum(segment[:, 1:] * x, axis=1) / np.sum(x * x, axis=1)
        low, high = equation.p_interval
        self.assertTrue(np.all(coefficient >= 2.0 * equation.sigma * low - 1e-12))
        self.assertTrue(np.all(coefficient <= 2.0 * equation.sigma * high + 1e-12))

    def test_c3_reference_and_segment(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "c3.npz"
            metadata = build_and_audit_reference(path)
            self.assertTrue(metadata["passed"])
            self.assertLessEqual(metadata["finest_richardson_difference_L8"], 1e-6)
            self.assertLessEqual(metadata["domain_L6_L8_max_difference"], 1e-7)
            self.assertLessEqual(
                metadata["independent_fd_residual"]["max_abs_residual"], 1e-5
            )
            equation = ProductionCubicHJB(d=20, reference_path=str(path))
            points = make_points(equation, pde_id="C3", n_points=32)
            self.assertLessEqual(
                float(np.max(np.abs(
                    equation.exact_u(equation.T, points["x"])
                    - equation.terminal(points["x"])
                ))),
                1e-14,
            )
            state = equation.exact_state(points["t"], points["x"])
            perturbed = state.copy()
            perturbed[:, 1:] += 50.0
            projected = _generic_segment(equation, perturbed)
            coefficient = projected[:, 1:] @ equation.w / equation.sigma
            self.assertLessEqual(np.max(np.abs(coefficient)), equation.A + 1e-12)
            orthogonal = projected[:, 1:] - (
                projected[:, 1:] @ equation.w
            )[:, None] * equation.w
            self.assertLess(np.max(np.linalg.norm(orthogonal, axis=1)), 1e-12)

    def test_tree_draws_are_paired_across_methods(self) -> None:
        equation = make_equation("P1", 20)
        points = make_points(equation, pde_id="P1", n_points=16)
        common = dict(
            pde_id="P1", equation=equation, n=2, M=2, repetition=0,
            t=points["t"], x=points["x"], is_validation=points["is_validation"],
            base_seed=BASE_SEED, chunk_size=8, study="test", code_commit=None,
            trace_draws=True,
        )
        raw = run_single_repetition(method=RAW, **common)
        box = run_single_repetition(method=BOX, **common)
        self.assertEqual(
            raw["metadata"]["draw_fingerprint"],
            box["metadata"]["draw_fingerprint"],
        )
        self.assertEqual(
            raw["metadata"]["work"]["f_evals"],
            box["metadata"]["work"]["f_evals"],
        )
        self.assertEqual(raw["metadata"]["dtype"], "float64")

    def test_regular_method_names_are_unique(self) -> None:
        for pde in ("P1", "P4", "C1", "C2-convex", "C3", "LQG"):
            names = [method.name for method in regular_methods(pde)]
            self.assertEqual(len(names), len(set(names)))


if __name__ == "__main__":
    unittest.main()
