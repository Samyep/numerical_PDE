"""Round-3 G1--G6 gates, including the production C3 reference gate."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from .equations import BASE_SEED, C3_REFERENCE, make_equation, make_points
from .reference_c3 import load_reference


def _skill(prediction: np.ndarray, truth: np.ndarray) -> float:
    return float(
        np.sqrt(np.mean((np.asarray(prediction) - np.asarray(truth)) ** 2))
        / np.std(truth, ddof=0)
    )


def time_only_source_gate(equation: Any, pde_id: str) -> dict[str, Any]:
    points = make_points(
        equation, pde_id=pde_id, n_points=2000, seed=BASE_SEED + 101
    )
    t, x = points["t"], points["x"]
    validation = points["is_validation"]
    exact = equation.exact_u(t, x)
    fzero = equation.fzero_reference(t, x)
    h = equation.T - t
    design = np.column_stack((h, h * h))
    coefficient, *_ = np.linalg.lstsq(
        design[validation], exact[validation] - fzero[validation], rcond=None
    )
    prediction = fzero + design @ coefficient
    test = ~validation
    return {
        "n_points": len(t),
        "fit_points": int(np.count_nonzero(validation)),
        "test_points": int(np.count_nonzero(test)),
        "coefficients_h_h2": coefficient.tolist(),
        "time_only_skill": _skill(prediction[test], exact[test]),
        "fzero_skill": _skill(fzero[test], exact[test]),
        "threshold": 0.15,
    }


def terminal_gradient_gate(equation: Any, pde_id: str) -> dict[str, Any]:
    points = make_points(
        equation, pde_id=pde_id, n_points=96, seed=BASE_SEED + 102
    )
    x = points["x"]
    terminal_error = float(
        np.max(np.abs(equation.exact_u(equation.T, x) - equation.terminal(x)))
    )
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, equation.d, 102])
    )
    indices = rng.choice(len(x), min(24, len(x)), replace=False)
    t = points["t"][indices]
    base = x[indices]
    exact = equation.exact_z(t, base) / equation.sigma
    coordinates = rng.choice(
        equation.d, min(8, equation.d), replace=False
    )
    step = 2.0e-4
    differences = []
    for coordinate in coordinates:
        plus = base.copy()
        minus = base.copy()
        plus[:, coordinate] += step
        minus[:, coordinate] -= step
        numerical = (
            equation.exact_u(t, plus) - equation.exact_u(t, minus)
        ) / (2.0 * step)
        differences.append(numerical - exact[:, coordinate])
    joined = np.concatenate(differences)
    return {
        "terminal_max_abs_error": terminal_error,
        "gradient_max_abs_error": float(np.max(np.abs(joined))),
        "gradient_rmse": float(np.sqrt(np.mean(joined**2))),
        "gradient_step": step,
        "gradient_coordinates_checked": coordinates.tolist(),
    }


def _ridge_fd_residual(equation: Any, pde_id: str) -> dict[str, Any]:
    points = make_points(
        equation, pde_id=pde_id, n_points=24, seed=BASE_SEED + 103
    )
    ht = min(5.0e-4, equation.T / 1000.0)
    hs = 2.0e-3
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, equation.d, 103]))
    t = rng.uniform(3.0 * ht, equation.T - 3.0 * ht, 24)
    x = points["x"][:24]
    w = equation.w
    u0 = equation.exact_u(t, x)
    u_t = (
        -equation.exact_u(t + 2.0 * ht, x)
        + 8.0 * equation.exact_u(t + ht, x)
        - 8.0 * equation.exact_u(t - ht, x)
        + equation.exact_u(t - 2.0 * ht, x)
    ) / (12.0 * ht)
    up1 = equation.exact_u(t, x + hs * w)
    um1 = equation.exact_u(t, x - hs * w)
    up2 = equation.exact_u(t, x + 2.0 * hs * w)
    um2 = equation.exact_u(t, x - 2.0 * hs * w)
    derivative = (-up2 + 8.0 * up1 - 8.0 * um1 + um2) / (12.0 * hs)
    second = (
        -up2 + 16.0 * up1 - 30.0 * u0 + 16.0 * um1 - um2
    ) / (12.0 * hs**2)
    z = equation.sigma * derivative[:, None] * w
    residual = u_t + second + equation.generator(u0, z)
    return {
        "kind": "independent fourth-order ridge finite-difference residual",
        "n_points": len(t),
        "time_step": ht,
        "spatial_step": hs,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
    }


def _coordinate_fd_residual(equation: Any, pde_id: str) -> dict[str, Any]:
    """Round-1-style fourth-order residual for the non-ridge LQG anchor."""

    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, equation.d, 107]))
    n_points = 12
    ht = min(2.0e-4, equation.T / 1000.0)
    hs = 2.0e-3
    t = rng.uniform(4.0 * ht, equation.T - 4.0 * ht, n_points)
    points = make_points(
        equation, pde_id=pde_id, n_points=n_points, seed=BASE_SEED + 107
    )
    x = points["x"]
    u0 = equation.exact_u(t, x)
    u_t = (
        -equation.exact_u(t + 2.0 * ht, x)
        + 8.0 * equation.exact_u(t + ht, x)
        - 8.0 * equation.exact_u(t - ht, x)
        + equation.exact_u(t - 2.0 * ht, x)
    ) / (12.0 * ht)
    gradient = np.empty((n_points, equation.d), dtype=np.float64)
    laplacian = np.zeros(n_points, dtype=np.float64)
    for coordinate in range(equation.d):
        plus1, minus1 = x.copy(), x.copy()
        plus2, minus2 = x.copy(), x.copy()
        plus1[:, coordinate] += hs
        minus1[:, coordinate] -= hs
        plus2[:, coordinate] += 2.0 * hs
        minus2[:, coordinate] -= 2.0 * hs
        up1 = equation.exact_u(t, plus1)
        um1 = equation.exact_u(t, minus1)
        up2 = equation.exact_u(t, plus2)
        um2 = equation.exact_u(t, minus2)
        gradient[:, coordinate] = (
            -up2 + 8.0 * up1 - 8.0 * um1 + um2
        ) / (12.0 * hs)
        laplacian += (
            -up2 + 16.0 * up1 - 30.0 * u0 + 16.0 * um1 - um2
        ) / (12.0 * hs**2)
    z = equation.sigma * gradient
    residual = u_t + laplacian + equation.generator(u0, z)
    return {
        "kind": "independent fourth-order coordinate finite-difference residual",
        "n_points": n_points,
        "time_step": ht,
        "spatial_step": hs,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
    }


def _analytic_residual(equation: Any, pde_id: str) -> dict[str, Any]:
    points = make_points(
        equation, pde_id=pde_id, n_points=512, seed=BASE_SEED + 104
    )
    if hasattr(equation, "analytic_residual"):
        residual = equation.analytic_residual(points["t"], points["x"])
    else:
        residual = equation.reference_residual(
            points["t"][:64], points["x"][:64]
        )
    return {
        "n_points": int(np.size(residual)),
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(np.asarray(residual) ** 2))),
    }


def _distances(
    equation: Any, state: np.ndarray, x: np.ndarray
) -> dict[str, np.ndarray]:
    z = np.asarray(state, dtype=np.float64)[..., 1:]
    if equation.family == "lqg":
        first, second = equation.dynamic_endpoints(x)
        low, high = np.minimum(first, second), np.maximum(first, second)
        squared_x = np.sum(x * x, axis=-1)
        coefficient = np.divide(
            np.sum(z * x, axis=-1), squared_x,
            out=np.zeros_like(squared_x), where=squared_x > 0.0,
        )
        parallel = coefficient[..., None] * x
        low_c = 2.0 * equation.sigma * equation.p_interval[0]
        high_c = 2.0 * equation.sigma * equation.p_interval[1]
        radius = high_c * np.sqrt(squared_x)
        return {
            "box": np.sqrt(np.sum(
                np.maximum(low - z, 0.0) ** 2
                + np.maximum(z - high, 0.0) ** 2,
                axis=-1,
            )),
            "ball": np.maximum(np.linalg.norm(z, axis=-1) - radius, 0.0),
            "span": np.linalg.norm(z - parallel, axis=-1),
            "segment": np.maximum.reduce((
                low_c - coefficient,
                coefficient - high_c,
                np.zeros_like(coefficient),
            )),
        }

    low = np.asarray(equation.box_low, dtype=np.float64)
    high = np.asarray(equation.box_high, dtype=np.float64)
    w = np.asarray(equation.w, dtype=np.float64)
    parallel = np.einsum("...d,d->...", z, w)[..., None] * w
    if equation.family == "ridge_lse":
        coefficient = -np.einsum("...d,d->...", z, w) / equation.sigma
    else:
        coefficient = np.einsum("...d,d->...", z, w) / equation.sigma
    lo, hi = equation.segment_coefficients
    return {
        "box": np.sqrt(np.sum(
            np.maximum(low - z, 0.0) ** 2
            + np.maximum(z - high, 0.0) ** 2,
            axis=-1,
        )),
        "ball": np.maximum(
            np.linalg.norm(z, axis=-1) - equation.z_ball_radius, 0.0
        ),
        "span": np.linalg.norm(z - parallel, axis=-1),
        "segment": np.maximum.reduce((
            lo - coefficient, coefficient - hi, np.zeros_like(coefficient)
        )),
    }


def containment_gate(
    equation: Any, pde_id: str, *, n_points: int = 100_000
) -> dict[str, Any]:
    points = make_points(
        equation, pde_id=pde_id, n_points=n_points, seed=BASE_SEED + 105
    )
    maxima: dict[str, float] = {}
    outside: dict[str, int] = {}
    for begin in range(0, n_points, 5000):
        end = min(begin + 5000, n_points)
        t, x = points["t"][begin:end], points["x"][begin:end]
        for name, distance in _distances(
            equation, equation.exact_state(t, x), x
        ).items():
            maxima[name] = max(maxima.get(name, 0.0), float(np.max(distance)))
            outside[name] = outside.get(name, 0) + int(
                np.count_nonzero(distance > 1.0e-10)
            )
    return {
        "n_points": n_points,
        "tolerance": 1.0e-10,
        "maximum_distance": maxima,
        "outside_counts": outside,
        "passed": all(count == 0 for count in outside.values()),
    }


def nonaffine_gate(equation: Any) -> dict[str, Any]:
    if equation.family == "lqg":
        direction = np.zeros(equation.d, dtype=np.float64)
        direction[0] = 1.0
        center = np.zeros(equation.d, dtype=np.float64)
        radius = 0.25
    else:
        direction = np.asarray(equation.w, dtype=np.float64)
        center = np.asarray(equation.z_center, dtype=np.float64)
        radius = min(float(equation.z_ball_radius) * 0.25, 0.25)
    values = np.asarray([
        equation.generator(np.array(0.0), (center - radius * direction)[None])[0],
        equation.generator(np.array(0.0), center[None])[0],
        equation.generator(np.array(0.0), (center + radius * direction)[None])[0],
    ])
    second = float(values[0] - 2.0 * values[1] + values[2])
    return {
        "probe_values": values.tolist(),
        "second_difference": second,
        "passed": abs(second) > 1.0e-10,
    }


def analytic_gates(pde_id: str, d: int) -> dict[str, Any]:
    equation = make_equation(pde_id, d)
    terminal_gradient = terminal_gradient_gate(equation, pde_id)
    time_only = time_only_source_gate(equation, pde_id)
    containment = containment_gate(equation, pde_id)
    nonaffine = nonaffine_gate(equation)
    analytic = _analytic_residual(equation, pde_id)
    finite_difference = None
    reference = None
    if pde_id == "C3":
        _, _, _, reference = load_reference(C3_REFERENCE)
        g1 = bool(
            reference["passed"]
            and terminal_gradient["terminal_max_abs_error"] <= 1.0e-10
            and terminal_gradient["gradient_max_abs_error"] <= 1.0e-6
        )
    elif equation.family == "lqg":
        finite_difference = _coordinate_fd_residual(equation, pde_id)
        g1 = bool(
            analytic["max_abs_residual"] <= 1.0e-6
            and finite_difference["max_abs_residual"] <= 1.0e-6
            and terminal_gradient["terminal_max_abs_error"] <= 1.0e-10
            and terminal_gradient["gradient_max_abs_error"] <= 1.0e-6
        )
    else:
        finite_difference = _ridge_fd_residual(equation, pde_id)
        g1 = bool(
            analytic["max_abs_residual"] <= 1.0e-6
            and finite_difference["max_abs_residual"] <= 1.0e-6
            and terminal_gradient["terminal_max_abs_error"] <= 1.0e-10
            and terminal_gradient["gradient_max_abs_error"] <= 1.0e-6
        )

    labels = equation.certificate_labels()
    label_ok = all(
        "PDE-derived" in label
        or "solution-informed" in label
        or "Riccati-derived" in label
        for label in labels.values()
    )
    return {
        "equation": equation.provenance(),
        "G1": {
            "passed": g1,
            "analytic": analytic,
            "finite_difference": finite_difference,
            "terminal_gradient": terminal_gradient,
            "production_reference": reference,
        },
        "G2": {"passed": time_only["time_only_skill"] >= 0.15, **time_only},
        "G5": {
            "passed": bool(label_ok and containment["passed"]),
            "labels": labels,
            "label_provenance_passed": label_ok,
            "containment": containment,
        },
        "G6": nonaffine,
        "picard_proxy": {"value": float(equation.picard_proxy), "guide": 1.5},
    }
