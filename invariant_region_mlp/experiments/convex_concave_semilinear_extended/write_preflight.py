"""Run CPU correctness checks and save their reproducibility record."""

from __future__ import annotations

import subprocess
import time

from .common import RESULTS_DIR, atomic_json, ensure_dirs


def main() -> None:
    ensure_dirs()
    started = time.perf_counter()
    command = [
        "/sw/user/python/miniforge3-pytorch-2.10.0/bin/python",
        "-m",
        "pytest",
        "invariant_region_mlp/experiments/convex_concave_semilinear_extended/test_correctness.py",
        "-q",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    atomic_json(
        RESULTS_DIR / "numerical_correctness.json",
        {
            "passed": completed.returncode == 0,
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
            "wall_seconds": time.perf_counter() - started,
            "coverage": [
                "generic analytic tangent coefficients versus driver autodiff",
                "joint tangent support for all four concave working drivers",
                "nonseparable mixed derivative and Hessian semidefiniteness on u in [0,1]",
                "dense orthonormal bases and k-dimensional terminal-gradient rank",
                "controlled terminal signal scaling",
                "full ambient Actor input with no basis projection",
                "generic weighted evaluator terminal recovery",
                "independent 2D periodic reference refinement, residual, terminal, and domain",
                "training core has no numerical-reference or basis dependency",
            ],
            "external_baselines_run": False,
            "rank_projection_used": False,
            "analytic_nonterminal_labels_used": False,
        },
    )
    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
