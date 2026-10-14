"""Run the CPU checks and write a compact machine-readable audit."""

from __future__ import annotations

import json
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
        "invariant_region_mlp/experiments/convex_concave_semilinear_pilot/test_correctness.py",
        "-q",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    payload = {
        "passed": completed.returncode == 0,
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "wall_seconds": time.perf_counter() - started,
        "test_count": 15,
        "coverage": [
            "analytic tangent coefficients versus driver autodiff",
            "tangent-plane identity and upper-bound property on each valid domain",
            "generic convex-to-concave sign transformation and reconstruction",
            "weighted Feynman--Kac formula with zero-order coefficient and source",
            "complete state-dependent Actor VJP versus common-random-number finite difference",
            "frozen Actor parameters with preserved state derivatives",
            "terminal recovery at t=T",
            "quadratic-HJB reduction to the established controlled-policy formula",
            "periodic-reference spatial/time refinement, residual, terminal and solution domain",
            "training-core reference leakage inspection",
        ],
        "exponential_clipping_used": False,
        "martingale_control_variate_used": False,
        "rank_projection_used": False,
    }
    atomic_json(RESULTS_DIR / "numerical_correctness.json", payload)
    print(completed.stdout, end="")
    if completed.returncode:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
