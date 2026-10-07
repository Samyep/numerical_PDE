"""Round-2 P4 certificate and its pre-registered containment audit."""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline

from .equations import make_points
from .norm_hjb import (
    NormDriverHJB,
    _bilinear_uniform,
    _monotone_half_line,
    load_norm_reference,
)


ROUND2_BASE_SEED = 20261107
GAUSS_HERMITE_ORDER = 80
HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
BOUND_CACHE_PATH = (
    PACKAGE_ROOT
    / "results"
    / "mechanism_suite_r2"
    / "reference_cache"
    / "p4_tanh_heat_gh80.npz"
)
BOUND_CACHE_PATH_TEXT = str(BOUND_CACHE_PATH.resolve())
CONTAINMENT_POINT_COUNT = 100_000
NEAR_TERMINAL_POINT_COUNT = 10_000
NEAR_TERMINAL_TAU_MAX = 0.1
NEAR_ZERO_S_HALF_WIDTH = 0.25

# Four successive monotone grids are fixed before the audit is run.  Their
# spacings on [0,12] are 0.005, 0.0025, 0.00125, and 0.000625.
CONTAINMENT_LEVELS = (
    (2401, 512),
    (4801, 1024),
    (9601, 2048),
    (19201, 4096),
)


@lru_cache(maxsize=4)
def _hermite_rule(order: int = GAUSS_HERMITE_ORDER) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    return (
        np.asarray(nodes, dtype=np.float64),
        np.asarray(weights / math.sqrt(math.pi), dtype=np.float64),
    )


def p4_derivative_bounds(
    tau: np.ndarray | float,
    s: np.ndarray | float,
    *,
    beta: float = 2.0,
    lambda_f: float = 1.0,
    order: int = GAUSS_HERMITE_ORDER,
    block_size: int = 65_536,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the constant-drift comparison bounds with Gauss--Hermite."""

    tau_array, s_array = np.broadcast_arrays(
        np.asarray(tau, dtype=np.float64), np.asarray(s, dtype=np.float64)
    )
    if np.any(tau_array < 0.0):
        raise ValueError("tau must be nonnegative")
    nodes, weights = _hermite_rule(order)
    flat_tau = tau_array.reshape(-1)
    flat_s = s_array.reshape(-1)
    lower = np.empty_like(flat_s)
    upper = np.empty_like(flat_s)
    for begin in range(0, len(flat_s), block_size):
        end = min(begin + block_size, len(flat_s))
        local_tau = flat_tau[begin:end]
        local_s = flat_s[begin:end]
        gaussian = 2.0 * np.sqrt(local_tau)[:, None] * nodes[None, :]
        drift = lambda_f * local_tau
        lower[begin:end] = np.sum(
            np.tanh(beta * (local_s[:, None] - drift[:, None] + gaussian))
            * weights[None, :],
            axis=1,
        )
        upper[begin:end] = np.sum(
            np.tanh(beta * (local_s[:, None] + drift[:, None] + gaussian))
            * weights[None, :],
            axis=1,
        )
    return lower.reshape(tau_array.shape), upper.reshape(tau_array.shape)


def _heat_tanh_direct(
    tau: np.ndarray | float,
    location: np.ndarray | float,
    *,
    beta: float = 2.0,
    order: int = GAUSS_HERMITE_ORDER,
    block_size: int = 65_536,
) -> np.ndarray:
    """Evaluate E[tanh(beta*(location+sqrt(2*tau)*xi))] directly."""

    tau_array, location_array = np.broadcast_arrays(
        np.asarray(tau, dtype=np.float64),
        np.asarray(location, dtype=np.float64),
    )
    nodes, weights = _hermite_rule(order)
    flat_tau = tau_array.reshape(-1)
    flat_location = location_array.reshape(-1)
    output = np.empty_like(flat_location)
    for begin in range(0, len(output), block_size):
        end = min(begin + block_size, len(output))
        gaussian = (
            2.0 * np.sqrt(flat_tau[begin:end])[:, None] * nodes[None, :]
        )
        output[begin:end] = np.sum(
            np.tanh(
                beta * (flat_location[begin:end, None] + gaussian)
            )
            * weights[None, :],
            axis=1,
        )
    return output.reshape(tau_array.shape)


def build_p4_bound_cache(path: Path = BOUND_CACHE_PATH) -> dict[str, Any]:
    """Tabulate the 80-node quadrature once for fast recursive projection."""

    path.parent.mkdir(parents=True, exist_ok=True)
    tau = np.linspace(0.0, 0.5, 257, dtype=np.float64)
    location = np.linspace(0.0, 12.0, 6001, dtype=np.float64)
    values = np.empty((len(tau), len(location)), dtype=np.float64)
    for index, local_tau in enumerate(tau):
        values[index] = _heat_tanh_direct(local_tau, location)
    values[:, 0] = 0.0

    spline = RectBivariateSpline(tau, location, values, kx=3, ky=3, s=0.0)
    rng = np.random.default_rng(
        np.random.SeedSequence([ROUND2_BASE_SEED, 4, 0xCA5E])
    )
    audit_tau = np.concatenate(
        (
            rng.uniform(0.0, 0.5, size=20_000),
            rng.uniform(0.0, 0.02, size=5_000),
        )
    )
    audit_location = np.concatenate(
        (
            rng.uniform(-8.0, 8.0, size=20_000),
            rng.uniform(-0.05, 0.05, size=5_000),
        )
    )
    direct = _heat_tanh_direct(audit_tau, audit_location)
    cached = np.sign(audit_location) * spline.ev(
        audit_tau, np.abs(audit_location)
    )
    maximum_error = float(np.max(np.abs(cached - direct)))
    metadata = {
        "schema_version": 1,
        "definition": "E[tanh(beta*(y+sqrt(2*tau)*xi))]",
        "beta": 2.0,
        "gauss_hermite_nodes": GAUSS_HERMITE_ORDER,
        "tau_points": len(tau),
        "tau_interval": [float(tau[0]), float(tau[-1])],
        "nonnegative_location_points": len(location),
        "location_interval": [float(location[0]), float(location[-1])],
        "interpolation": "odd cubic RectBivariateSpline",
        "audit_points": len(audit_tau),
        "audit_maximum_absolute_error": maximum_error,
        "audit_threshold": 2e-6,
        "audit_passed": maximum_error <= 2e-6,
    }
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        tau=tau,
        location=location,
        values=values,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    temporary.replace(path)
    load_p4_bound_cache.cache_clear()
    return metadata


@lru_cache(maxsize=4)
def load_p4_bound_cache(
    path_text: str = BOUND_CACHE_PATH_TEXT,
) -> tuple[RectBivariateSpline, dict[str, Any], float]:
    path = Path(path_text)
    with np.load(path, allow_pickle=False) as data:
        tau = np.asarray(data["tau"], dtype=np.float64)
        location = np.asarray(data["location"], dtype=np.float64)
        values = np.asarray(data["values"], dtype=np.float64)
        metadata = json.loads(str(data["metadata_json"]))
    if not metadata.get("audit_passed", False):
        raise RuntimeError("P4 bound cache failed its interpolation audit")
    return (
        RectBivariateSpline(tau, location, values, kx=3, ky=3, s=0.0),
        metadata,
        float(location[-1]),
    )


def ensure_p4_bound_cache(path: Path = BOUND_CACHE_PATH) -> dict[str, Any]:
    if not path.exists():
        return build_p4_bound_cache(path)
    _, metadata, _ = load_p4_bound_cache(str(path.resolve()))
    return metadata


def p4_derivative_bounds_cached(
    tau: np.ndarray | float,
    s: np.ndarray | float,
    *,
    lambda_f: float = 1.0,
    cache_path: Path = BOUND_CACHE_PATH,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate bounds from the audited table generated by 80-node GH."""

    tau_array, s_array = np.broadcast_arrays(
        np.asarray(tau, dtype=np.float64), np.asarray(s, dtype=np.float64)
    )
    path_text = (
        BOUND_CACHE_PATH_TEXT
        if cache_path == BOUND_CACHE_PATH
        else str(cache_path.resolve())
    )
    spline, _, domain = load_p4_bound_cache(path_text)

    def evaluate(location: np.ndarray) -> np.ndarray:
        result = np.empty_like(location, dtype=np.float64)
        inside = np.abs(location) <= domain
        if np.any(inside):
            result[inside] = np.sign(location[inside]) * spline.ev(
                tau_array[inside], np.abs(location[inside])
            )
        if np.any(~inside):
            result[~inside] = _heat_tanh_direct(
                tau_array[~inside], location[~inside]
            )
        return result

    drift = lambda_f * tau_array
    return evaluate(s_array - drift), evaluate(s_array + drift)


def p4_tightness_interval(
    tau: np.ndarray | float,
    s: np.ndarray | float,
    theta: float,
    *,
    beta: float = 2.0,
    lambda_f: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate between the round-1 segment and the tight certificate."""

    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must lie in [0,1]")
    v_minus, v_plus = p4_derivative_bounds(
        tau, s, beta=beta, lambda_f=lambda_f
    )
    lower = (1.0 - theta) * -1.0 + theta * v_minus
    upper = (1.0 - theta) * 1.0 + theta * v_plus
    return lower, upper


def _maximum_violation(
    value: np.ndarray, lower: np.ndarray, upper: np.ndarray
) -> tuple[float, int, float]:
    violation = np.maximum(
        np.maximum(lower - value, value - upper), 0.0
    )
    return (
        float(np.max(violation)),
        int(np.count_nonzero(violation > 0.0)),
        float(np.quantile(violation, 0.999)),
    )


def run_p4_containment_audit(reference_path: Path) -> dict[str, Any]:
    """Run the pre-registered sharp-bound audit on four successive grids."""

    equation = NormDriverHJB(
        d=20,
        reference_path=str(reference_path.resolve()),
        beta=2.0,
        lambda_f=1.0,
        T=0.5,
    )
    points = make_points(
        equation,
        n_points=CONTAINMENT_POINT_COUNT,
        seed=ROUND2_BASE_SEED,
    )
    test_tau = equation.T - np.asarray(points["t"], dtype=np.float64)
    test_s = np.asarray(points["x"], dtype=np.float64) @ equation.w

    targeted_rng = np.random.default_rng(
        np.random.SeedSequence([ROUND2_BASE_SEED, 4, 0xA503])
    )
    targeted_tau = targeted_rng.uniform(
        0.0, NEAR_TERMINAL_TAU_MAX, size=NEAR_TERMINAL_POINT_COUNT
    ).astype(np.float64)
    targeted_s = targeted_rng.uniform(
        -NEAR_ZERO_S_HALF_WIDTH,
        NEAR_ZERO_S_HALF_WIDTH,
        size=NEAR_TERMINAL_POINT_COUNT,
    ).astype(np.float64)

    tau = np.concatenate((test_tau, targeted_tau))
    s = np.concatenate((test_s, targeted_s))
    lower, upper = p4_derivative_bounds(
        tau,
        s,
        beta=equation.beta,
        lambda_f=equation.lambda_f,
    )
    sample_hash = hashlib.sha256(
        np.ascontiguousarray(np.column_stack((tau, s)), dtype=np.float64).tobytes()
    ).hexdigest()

    grid_records: list[dict[str, Any]] = []
    for n_space, n_steps in CONTAINMENT_LEVELS:
        tau_grid, r_grid, value = _monotone_half_line(
            beta=equation.beta,
            lambda_f=equation.lambda_f,
            horizon=equation.T,
            domain=12.0,
            n_space=n_space,
            n_steps=n_steps,
            n_output=513,
        )
        spacing = float(r_grid[1] - r_grid[0])
        radial_derivative = np.gradient(
            value, spacing, axis=1, edge_order=2
        )
        radial_derivative[:, 0] = 0.0
        radial_derivative[:, -1] = 1.0
        signed_derivative = np.sign(s) * _bilinear_uniform(
            radial_derivative,
            tau_grid,
            r_grid,
            tau,
            np.abs(s),
        )
        maximum, positive_count, q999 = _maximum_violation(
            signed_derivative, lower, upper
        )
        grid_records.append(
            {
                "n_space": n_space,
                "n_steps": n_steps,
                "ds": spacing,
                "maximum_violation": maximum,
                "positive_violation_count": positive_count,
                "violation_q999": q999,
            }
        )

    reference = load_norm_reference(str(reference_path.resolve()))
    reference_derivative = np.sign(s) * reference.derivative_spline.ev(
        tau, np.abs(s)
    )
    reference_maximum, reference_positive_count, reference_q999 = _maximum_violation(
        reference_derivative, lower, upper
    )
    cached_lower, cached_upper = p4_derivative_bounds_cached(tau, s)
    cached_reference_maximum, cached_positive_count, cached_q999 = _maximum_violation(
        reference_derivative, cached_lower, cached_upper
    )
    cached_endpoint_difference = float(
        max(
            np.max(np.abs(cached_lower - lower)),
            np.max(np.abs(cached_upper - upper)),
        )
    )
    violations = np.asarray(
        [record["maximum_violation"] for record in grid_records],
        dtype=np.float64,
    )
    strictly_decreasing = bool(np.all(np.diff(violations) < 0.0))
    finest_below_threshold = bool(violations[-1] < 1e-5)
    implemented_interval_below_threshold = bool(
        cached_reference_maximum < 1e-5
    )

    spot_checks = []
    for local_tau in (0.1, 0.25, 0.5):
        spot_lower, spot_upper = p4_derivative_bounds(local_tau, 0.0)
        spot_checks.append(
            {
                "tau": local_tau,
                "s": 0.0,
                "v_minus": float(spot_lower),
                "v_plus": float(spot_upper),
            }
        )

    return {
        "schema_version": 1,
        "base_seed": ROUND2_BASE_SEED,
        "dtype": "float64",
        "quadrature": {
            "rule": "Gauss-Hermite transformed for a standard normal",
            "nodes": GAUSS_HERMITE_ORDER,
        },
        "sample": {
            "test_distribution_points": CONTAINMENT_POINT_COUNT,
            "near_terminal_near_zero_points": NEAR_TERMINAL_POINT_COUNT,
            "near_terminal_tau_interval": [0.0, NEAR_TERMINAL_TAU_MAX],
            "near_zero_s_interval": [
                -NEAR_ZERO_S_HALF_WIDTH,
                NEAR_ZERO_S_HALF_WIDTH,
            ],
            "tau_s_sha256": sample_hash,
        },
        "successive_grid_results": grid_records,
        "strictly_decreasing_under_refinement": strictly_decreasing,
        "finest_grid_threshold": 1e-5,
        "finest_grid_below_threshold": finest_below_threshold,
        "accepted_reference_check": {
            "maximum_violation": reference_maximum,
            "positive_violation_count": reference_positive_count,
            "violation_q999": reference_q999,
            "reference_metadata": reference.metadata,
        },
        "implemented_cached_interval_check": {
            "maximum_endpoint_difference_from_direct_gh80": cached_endpoint_difference,
            "maximum_reference_violation": cached_reference_maximum,
            "positive_violation_count": cached_positive_count,
            "violation_q999": cached_q999,
            "reference_violation_below_1e-5": implemented_interval_below_threshold,
        },
        "spot_checks_at_s_zero": spot_checks,
        "proceed": (
            strictly_decreasing
            and finest_below_threshold
            and implemented_interval_below_threshold
        ),
    }
