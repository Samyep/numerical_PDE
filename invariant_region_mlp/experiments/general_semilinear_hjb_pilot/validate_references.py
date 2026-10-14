"""Run and record independent numerical-reference refinement checks."""

from __future__ import annotations

import numpy as np

from .common import RESULTS_DIR, atomic_json, ensure_dirs
from .problems import BoundedControlProblem, LQGProblem, QuarticControlProblem
from .reference import RiccatiReference, validate_scalar_reference


def main() -> None:
    ensure_dirs()
    problem = LQGProblem()
    reference = RiccatiReference(problem)
    grid = np.linspace(0.0, problem.T, 1001)
    P, c = reference.coefficients(grid)
    derivative = np.gradient(P, grid, edge_order=2)
    residual = derivative - 2.0 * problem.gamma * P - P**2 / problem.R + problem.Q
    payload = {
        "lqg": {
            "solver": "scipy.solve_ivp RK45 with dense output",
            "rtol": 1e-11,
            "atol": 1e-13,
            "terminal_P_error": float(abs(P[-1] - problem.G)),
            "terminal_c_error": float(abs(c[-1])),
            "riccati_residual_grid_max_abs": float(np.max(np.abs(residual[2:-2]))),
            "P_min": float(P.min()),
            "P_max": float(P.max()),
        }
    }
    for scalar_problem in (QuarticControlProblem(), BoundedControlProblem()):
        _, metrics = validate_scalar_reference(scalar_problem)
        payload[scalar_problem.name] = metrics
    payload["training_leakage"] = {
        "reference_imported_by_training_core": False,
        "reference_used_for_checkpoint_selection": False,
        "reference_used_for_policy_iteration_decision": False,
    }
    atomic_json(RESULTS_DIR / "reference_validation.json", payload)
    print(RESULTS_DIR / "reference_validation.json")


if __name__ == "__main__":
    main()

