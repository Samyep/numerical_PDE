"""E0 analytic/reference gates for the mechanism benchmark suite."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .equations import BASE_SEED, make_points, stable_seed


def skill(prediction: np.ndarray, truth: np.ndarray) -> float:
    pred = np.asarray(prediction, dtype=np.float64).reshape(-1)
    ref = np.asarray(truth, dtype=np.float64).reshape(-1)
    return float(np.sqrt(np.mean((pred - ref) ** 2)) / np.std(ref, ddof=0))


def time_only_source_gate(equation: Any, *, n_points: int = 2000) -> dict[str, Any]:
    points = make_points(equation, n_points=n_points, seed=BASE_SEED + 101)
    t = points["t"]
    x = points["x"]
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
        "n_points": n_points,
        "fit_points": int(np.count_nonzero(validation)),
        "test_points": int(np.count_nonzero(test)),
        "coefficients_h_h2": coefficient.tolist(),
        "time_only_skill": skill(prediction[test], exact[test]),
        "fzero_skill": skill(fzero[test], exact[test]),
        "threshold": 0.15,
    }


def terminal_and_gradient_checks(equation: Any) -> dict[str, Any]:
    points = make_points(equation, n_points=96, seed=BASE_SEED + 102)
    x = points["x"]
    terminal_error = float(
        np.max(np.abs(equation.exact_u(equation.T, x) - equation.terminal(x)))
    )
    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, equation.d, stable_seed(equation.name), 102])
    )
    indices = rng.choice(len(x), size=min(24, len(x)), replace=False)
    t = points["t"][indices]
    base = x[indices]
    exact_gradient = equation.exact_z(t, base) / equation.sigma
    coordinate_count = min(8, equation.d)
    coordinates = rng.choice(equation.d, size=coordinate_count, replace=False)
    step = 2e-4
    maximum = 0.0
    rms_terms: list[np.ndarray] = []
    for coordinate in coordinates:
        plus = base.copy()
        minus = base.copy()
        plus[:, coordinate] += step
        minus[:, coordinate] -= step
        numerical = (
            equation.exact_u(t, plus) - equation.exact_u(t, minus)
        ) / (2.0 * step)
        difference = numerical - exact_gradient[:, coordinate]
        maximum = max(maximum, float(np.max(np.abs(difference))))
        rms_terms.append(difference)
    rms = float(np.sqrt(np.mean(np.concatenate(rms_terms) ** 2)))
    return {
        "terminal_max_abs_error": terminal_error,
        "gradient_max_abs_error": maximum,
        "gradient_rmse": rms,
        "gradient_coordinates_checked": coordinates.tolist(),
        "gradient_step": step,
    }


def finite_difference_residual(equation: Any, *, n_points: int = 12) -> dict[str, Any]:
    """Fourth-order finite-difference residual for analytic equations."""

    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, equation.d, stable_seed(equation.name), 103])
    )
    spatial_low = -0.5 if equation.family in {"vba", "burgers_fisher", "published_vb"} else -1.0
    spatial_high = -spatial_low
    x = rng.uniform(spatial_low, spatial_high, size=(n_points, equation.d))
    time_step = min(2e-4, equation.T / 1000.0)
    t = rng.uniform(4.0 * time_step, equation.T - 4.0 * time_step, size=n_points)
    spatial_step = 2e-3

    u_t = (
        -equation.exact_u(t + 2.0 * time_step, x)
        + 8.0 * equation.exact_u(t + time_step, x)
        - 8.0 * equation.exact_u(t - time_step, x)
        + equation.exact_u(t - 2.0 * time_step, x)
    ) / (12.0 * time_step)
    u0 = equation.exact_u(t, x)
    gradient = np.empty((n_points, equation.d), dtype=np.float64)
    laplacian = np.zeros(n_points, dtype=np.float64)
    for coordinate in range(equation.d):
        plus1 = x.copy()
        minus1 = x.copy()
        plus2 = x.copy()
        minus2 = x.copy()
        plus1[:, coordinate] += spatial_step
        minus1[:, coordinate] -= spatial_step
        plus2[:, coordinate] += 2.0 * spatial_step
        minus2[:, coordinate] -= 2.0 * spatial_step
        up1 = equation.exact_u(t, plus1)
        um1 = equation.exact_u(t, minus1)
        up2 = equation.exact_u(t, plus2)
        um2 = equation.exact_u(t, minus2)
        gradient[:, coordinate] = (-up2 + 8.0 * up1 - 8.0 * um1 + um2) / (
            12.0 * spatial_step
        )
        laplacian += (-up2 + 16.0 * up1 - 30.0 * u0 + 16.0 * um1 - um2) / (
            12.0 * spatial_step**2
        )
    drift = np.sum(np.asarray(equation.mu) * gradient, axis=-1)
    z = equation.sigma * gradient
    residual = u_t + drift + 0.5 * equation.sigma**2 * laplacian + equation.generator(u0, z)
    return {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
        "time_step": time_step,
        "spatial_step": spatial_step,
        "n_points": n_points,
    }


def numerical_reference_finite_difference_residual(
    equation: Any,
    *,
    n_points: int = 24,
    time_step: float = 5e-4,
    spatial_step: float = 2e-3,
) -> dict[str, Any]:
    """Independent fourth-order FD residual for the P4 ridge reference."""

    rng = np.random.default_rng(
        np.random.SeedSequence([BASE_SEED, equation.d, stable_seed(equation.name), 106])
    )
    t = rng.uniform(3.0 * time_step, equation.T - 3.0 * time_step, size=n_points)
    x = rng.uniform(-1.0, 1.0, size=(n_points, equation.d))
    u0 = equation.exact_u(t, x)
    u_t = (
        -equation.exact_u(t + 2.0 * time_step, x)
        + 8.0 * equation.exact_u(t + time_step, x)
        - 8.0 * equation.exact_u(t - time_step, x)
        + equation.exact_u(t - 2.0 * time_step, x)
    ) / (12.0 * time_step)
    plus1 = x + spatial_step * equation.w
    minus1 = x - spatial_step * equation.w
    plus2 = x + 2.0 * spatial_step * equation.w
    minus2 = x - 2.0 * spatial_step * equation.w
    up1 = equation.exact_u(t, plus1)
    um1 = equation.exact_u(t, minus1)
    up2 = equation.exact_u(t, plus2)
    um2 = equation.exact_u(t, minus2)
    directional_gradient = (-up2 + 8.0 * up1 - 8.0 * um1 + um2) / (
        12.0 * spatial_step
    )
    laplacian = (-up2 + 16.0 * up1 - 30.0 * u0 + 16.0 * um1 - um2) / (
        12.0 * spatial_step**2
    )
    z = equation.sigma * directional_gradient[:, None] * equation.w
    residual = u_t + laplacian + equation.generator(u0, z)
    return {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
        "time_step": time_step,
        "spatial_step": spatial_step,
        "n_points": n_points,
        "kind": "independent fourth-order finite-difference residual",
    }


def analytic_residual_check(equation: Any) -> dict[str, Any]:
    points = make_points(equation, n_points=512, seed=BASE_SEED + 104)
    if hasattr(equation, "analytic_residual"):
        residual = equation.analytic_residual(points["t"], points["x"])
        return {
            "max_abs_residual": float(np.max(np.abs(residual))),
            "rmse_residual": float(np.sqrt(np.mean(np.asarray(residual) ** 2))),
            "n_points": len(points["t"]),
        }
    residual = equation.reference_residual(points["t"][:64], points["x"][:64])
    return {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(np.asarray(residual) ** 2))),
        "n_points": 64,
        "kind": "spline-differential residual of numerical reference",
    }


def _certificate_distances(equation: Any, state: np.ndarray) -> dict[str, np.ndarray]:
    u = state[..., 0]
    z = state[..., 1:]
    low = np.asarray(equation.box_low, dtype=np.float64)
    high = np.asarray(equation.box_high, dtype=np.float64)
    distances: dict[str, np.ndarray] = {
        "box": np.sqrt(
            np.sum(np.maximum(low - z, 0.0) ** 2 + np.maximum(z - high, 0.0) ** 2, axis=-1)
        ),
        "ball": np.maximum(np.linalg.norm(z, axis=-1) - equation.z_ball_radius, 0.0),
    }
    if getattr(equation, "has_u_certificate", False):
        u_low, u_high = equation.u_interval
        distances["u_interval"] = np.maximum.reduce(
            (u_low - u, u - u_high, np.zeros_like(u))
        )
    if equation.family in {"vba", "burgers_fisher", "published_vb"}:
        distances["z_nonnegative"] = np.max(np.maximum(-z, 0.0), axis=-1)
    if equation.family in {"ridge_lse", "norm_hjb"}:
        w = equation.w
        parallel = np.einsum("...d,d->...", z, w)[..., None] * w
        distances["span"] = np.linalg.norm(z - parallel, axis=-1)
        if equation.family == "ridge_lse":
            coefficient = -np.einsum("...d,d->...", z, w) / equation.sigma
            lo, hi = equation.segment_coefficients
        else:
            coefficient = np.einsum("...d,d->...", z, w) / equation.sigma
            lo, hi = -1.0, 1.0
        distances["segment"] = np.maximum.reduce(
            (lo - coefficient, coefficient - hi, np.zeros_like(coefficient))
        )
    if equation.family == "multidirection_lse":
        # Span containment is checked through the row-space projector.
        vectors = equation.vectors
        q, _ = np.linalg.qr(vectors.T)
        projection = np.einsum("...d,dk,kj->...j", z, q, q.T)
        distances["span"] = np.linalg.norm(z - projection, axis=-1)
    return distances


def containment_gate(
    equation: Any,
    *,
    n_points: int = 100_000,
    chunk_size: int = 5_000,
) -> dict[str, Any]:
    # Use the actual fixed test distribution, including the 1/6 face-boundary
    # stratum of the VB family, rather than a convenient interior surrogate.
    points = make_points(equation, n_points=n_points, seed=BASE_SEED + 105)
    maxima: dict[str, float] = {}
    outside: dict[str, int] = {}
    for begin in range(0, n_points, chunk_size):
        end = min(begin + chunk_size, n_points)
        x = points["x"][begin:end]
        t = points["t"][begin:end]
        state = equation.exact_state(t, x)
        for name, distance in _certificate_distances(equation, state).items():
            maximum = float(np.max(distance))
            maxima[name] = max(maxima.get(name, 0.0), maximum)
            outside[name] = outside.get(name, 0) + int(np.count_nonzero(distance > 1e-10))
    return {
        "n_points": n_points,
        "boundary_points": int(np.count_nonzero(points["is_boundary"])),
        "tolerance": 1e-10,
        "maximum_distance": maxima,
        "outside_counts": outside,
        "passed": all(value == 0 for value in outside.values()),
    }


def nonaffine_gate(equation: Any) -> dict[str, Any]:
    if equation.family in {"ridge_lse", "norm_hjb", "multidirection_lse"}:
        z0 = np.asarray(equation.z_center, dtype=np.float64)
        direction = equation.w if hasattr(equation, "w") else equation.vectors[0] / np.linalg.norm(equation.vectors[0])
        radius = min(float(equation.z_ball_radius) * 0.25, 0.25)
        values = np.array(
            [
                equation.generator(np.array(0.0), (z0 - radius * direction)[None, :])[0],
                equation.generator(np.array(0.0), z0[None, :])[0],
                equation.generator(np.array(0.0), (z0 + radius * direction)[None, :])[0],
            ]
        )
    else:
        u = np.array([0.2, 0.5, 0.8], dtype=np.float64)
        # Follow a joint-state line through C.  Holding z fixed would make the
        # P2 product driver look affine in u even though it is bilinear on C.
        z = u[:, None] * np.full((1, equation.d), equation.z_upper, dtype=np.float64)
        values = equation.generator(u, z)
    second_difference = float(values[0] - 2.0 * values[1] + values[2])
    return {
        "probe_values": values.tolist(),
        "second_difference": second_difference,
        "passed": abs(second_difference) > 1e-10,
    }


def run_analytic_gates(equation: Any) -> dict[str, Any]:
    analytic = analytic_residual_check(equation)
    if equation.family == "norm_hjb":
        # Audit a stable band around the time/space grid scales rather than
        # relying on a single fortuitous finite-difference step.  Steps much
        # larger smear the even spline across s=0; steps below the stored grid
        # scale suffer subtraction cancellation in the second derivative.
        step_pairs = ((7.5e-4, 2.5e-3), (5e-4, 2e-3), (4e-4, 1.5e-3))
        step_audit = [
            numerical_reference_finite_difference_residual(
                equation, time_step=time_step, spatial_step=spatial_step
            )
            for time_step, spatial_step in step_pairs
        ]
        finite_difference = dict(step_audit[1])
        finite_difference["max_abs_residual"] = max(
            item["max_abs_residual"] for item in step_audit
        )
        finite_difference["rmse_residual"] = max(
            item["rmse_residual"] for item in step_audit
        )
        finite_difference["step_audit"] = step_audit
        finite_difference["kind"] = "independent fourth-order finite-difference residual; three-step stability audit"
    else:
        finite_difference = finite_difference_residual(equation)
    terminal_gradient = terminal_and_gradient_checks(equation)
    time_only = time_only_source_gate(equation)
    containment = containment_gate(equation)
    nonaffine = nonaffine_gate(equation)
    reference_refinement = None
    reference_domain = None
    reference_gradient = None
    if equation.family == "norm_hjb":
        reference_parameters = equation.certificate_parameters()
        reference_refinement = reference_parameters["reference_metadata"]
        reference_gradient = {
            "method": reference_parameters["gradient_evaluator"],
            "max_abs_difference_from_direct_spline_derivative": reference_parameters[
                "gradient_cache_max_difference"
            ],
            "threshold": reference_parameters["gradient_cache_threshold"],
            "passed": reference_parameters["gradient_cache_max_difference"] <= 1e-6,
        }
        domain_points = make_points(equation, n_points=1200, seed=BASE_SEED)
        maximum_abs_s = float(np.max(np.abs(domain_points["x"] @ equation.w)))
        half_domain = float(reference_refinement["domain"][1])
        reference_domain = {
            "half_domain": half_domain,
            "maximum_abs_s_on_fixed_test_points": maximum_abs_s,
            "required_buffer": 8.0,
            "actual_buffer": half_domain - maximum_abs_s,
            "passed": half_domain - maximum_abs_s >= 8.0,
        }
    g1 = (
        analytic["max_abs_residual"] <= 1e-6
        and finite_difference["max_abs_residual"] <= 1e-6
        and terminal_gradient["terminal_max_abs_error"] <= 1e-10
        and terminal_gradient["gradient_max_abs_error"] <= 1e-6
        and (
            reference_refinement is None
            or (
                reference_refinement["richardson_refinement_max_difference"] <= 1e-6
                and reference_refinement.get("storage_spline_passed", False)
                and reference_domain["passed"]
                and reference_parameters["gradient_cache_max_difference"] <= 1e-6
            )
        )
    )
    g2 = time_only["time_only_skill"] >= 0.15
    required_labels = {
        "ridge_lse": {"segment", "box", "ball", "span"},
        "norm_hjb": {"segment", "box", "ball", "span"},
        "vba": {"u_[0,1]", "z_nonnegative", "z_upper_sigma_over_4", "ball"},
        "burgers_fisher": {"u_[0,1]", "z_nonnegative", "z_upper_sigma_over_4", "ball"},
        "published_vb": {"u_[0,1]", "z_nonnegative", "z_upper_sigma_over_4", "ball"},
        "multidirection_lse": {"convex_hull", "box", "ball", "span"},
    }[equation.family]
    labels = equation.certificate_labels()
    missing_labels = sorted(required_labels - set(labels))
    g5 = not missing_labels and containment["passed"] and all(
        "PDE-derived" in label or "solution-informed" in label
        for label in labels.values()
    )
    g6 = bool(nonaffine["passed"])
    return {
        "equation": equation.provenance(),
        "G1": {"passed": g1, "threshold": "all residual/gradient errors <=1e-6", "analytic": analytic, "finite_difference": finite_difference, "terminal_gradient": terminal_gradient, "reference_refinement": reference_refinement, "reference_domain": reference_domain, "reference_gradient": reference_gradient},
        "G2": {"passed": g2, **time_only},
        "G5": {"passed": g5, "labels": labels, "required_labels": sorted(required_labels), "missing_labels": missing_labels, "containment": containment},
        "G6": nonaffine,
        "picard_proxy": {"value": float(equation.picard_proxy), "guide": 1.5},
    }
