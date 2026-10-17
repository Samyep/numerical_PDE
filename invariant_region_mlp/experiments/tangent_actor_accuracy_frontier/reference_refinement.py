"""Independent refinement of the two-dimensional primary numerical reference."""

from __future__ import annotations

import argparse
import json
import time

import numpy as np

from .common import (
    OLD_PREDICTION_DIR,
    OLD_REFERENCE_DIR,
    RAW_DIR,
    REFERENCE_DIR,
    array_sha256,
    atomic_json,
    ensure_dirs,
)
from .drivers import SineQuadraticDriver
from .reference2d import load_reference, save_reference, solve_strang_reference


def _comparison(left, right, t: np.ndarray, coordinates: np.ndarray) -> dict[str, float]:
    value = left.value(t, coordinates) - right.value(t, coordinates)
    gradient = left.gradient(t, coordinates) - right.gradient(t, coordinates)
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(value)))),
        "value_max_abs": float(np.max(np.abs(value))),
        "gradient_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient).sum(axis=1)))
        ),
        "gradient_max_norm": float(np.max(np.linalg.norm(gradient, axis=1))),
    }


def run(*, n: int, dt: float, query_points: int) -> dict:
    ensure_dirs()
    driver = SineQuadraticDriver(d=100, lam=0.2, eta=0.0)
    old_path = OLD_REFERENCE_DIR / f"{driver.name}.npz"
    old = load_reference(old_path)
    output_path = REFERENCE_DIR / f"{driver.name}_refined.npz"
    previous = load_reference(output_path) if output_path.exists() else None
    started = time.perf_counter()
    refined = solve_strang_reference(driver, n=n, dt=dt, store_dt=0.002)
    wall = time.perf_counter() - started
    rng = np.random.default_rng(2026101707)
    t = rng.uniform(0.0, driver.T, query_points)
    coordinates = rng.uniform(0.0, 2.0 * np.pi, (query_points, 2))
    random_comparison = _comparison(old, refined, t, coordinates)

    bank_path = OLD_PREDICTION_DIR / (
        "post_freeze_primary_d100_k2_lam0p2_eta0_main_mcv_s0.npz"
    )
    with np.load(bank_path) as bank:
        test_t = bank["test_t"].astype(np.float64)
        test_x = bank["test_x"].astype(np.float64)
    frozen_comparison = _comparison(old, refined, test_t, test_x @ driver.B)
    previous_random = None if previous is None else _comparison(
        previous, refined, t, coordinates
    )
    previous_frozen = None if previous is None else _comparison(
        previous, refined, test_t, test_x @ driver.B
    )
    save_reference(output_path, refined)
    # Once a third level exists, use the complete previous-to-current
    # discrepancy without assuming a Richardson order. Otherwise retain the
    # old-to-current discrepancy as a conservative two-level proxy.
    uncertainty_source = (
        (previous_random, previous_frozen)
        if previous is not None else (random_comparison, frozen_comparison)
    )
    uncertainty = max(item["value_rmse"] for item in uncertainty_source)
    result = {
        "driver": driver.metadata(),
        "old_reference": str(old_path),
        "old_n": old.n,
        "old_dt": old.dt,
        "refined_reference": str(output_path),
        "refined_n": refined.n,
        "refined_dt": refined.dt,
        "random_periodic_query": random_comparison,
        "frozen_test_bank": frozen_comparison,
        "previous_refined_n": None if previous is None else previous.n,
        "previous_refined_dt": None if previous is None else previous.dt,
        "previous_refined_random_query": previous_random,
        "previous_refined_frozen_test_bank": previous_frozen,
        "reference_uncertainty_value_rmse": uncertainty,
        "uncertainty_rule": (
            "max previous-refined-vs-current RMSE; no Richardson-order assumption"
            if previous is not None
            else "max old-vs-refined RMSE; conservative two-level proxy"
        ),
        "query_points": query_points,
        "frozen_test_bank_sha256": array_sha256(test_t, test_x),
        "solution_grid_min": float(refined.values.min()),
        "solution_grid_max": float(refined.values.max()),
        "inside_invariant_interval": bool(
            refined.values.min() >= -2e-8 and refined.values.max() <= 1.0 + 2e-8
        ),
        "wall_seconds": wall,
    }
    atomic_json(REFERENCE_DIR / "primary_refinement_audit.json", result)
    atomic_json(RAW_DIR / "reference_refinement_primary.json", result)
    print(json.dumps(result, indent=2), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=256)
    parser.add_argument("--dt", type=float, default=0.00003125)
    parser.add_argument("--query-points", type=int, default=32768)
    args = parser.parse_args()
    run(n=args.n, dt=args.dt, query_points=args.query_points)


if __name__ == "__main__":
    main()

