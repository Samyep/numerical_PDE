"""P4 norm-driver ridge HJB and its independently refined 1-D reference."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import splu

from .equations import BASE_SEED, ExactEquation, _as_points


def _monotone_half_line(
    *,
    beta: float,
    lambda_f: float,
    horizon: float,
    domain: float,
    n_space: int,
    n_steps: int,
    n_output: int = 257,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Backward-Euler monotone upwind solve on r=|s| in [0,L].

    For the even convex terminal condition the derivative is nonnegative on
    r>=0, so ``|v_r|=v_r`` and the equation is the linear half-line problem
    ``v_tau=v_rr-lambda_f*v_r`` with Neumann slopes 0 and 1.  Backward Euler
    applied to the first-order upwind spatial operator is an M-matrix scheme.
    """

    if n_steps % (n_output - 1):
        raise ValueError("n_steps must be divisible by n_output-1")
    r = np.linspace(0.0, domain, n_space, dtype=np.float64)
    dx = float(r[1] - r[0])
    dt = horizon / n_steps

    lower = np.full(n_space - 1, 1.0 / dx**2 + lambda_f / dx, dtype=np.float64)
    diagonal = np.full(n_space, -2.0 / dx**2 - lambda_f / dx, dtype=np.float64)
    upper = np.full(n_space - 1, 1.0 / dx**2, dtype=np.float64)
    # Reflecting derivative at r=0.
    diagonal[0] = -2.0 / dx**2
    upper[0] = 2.0 / dx**2
    # Prescribed derivative v_r=1 at the remote right boundary.
    lower[-1] = 2.0 / dx**2
    diagonal[-1] = -2.0 / dx**2
    operator = diags((lower, diagonal, upper), offsets=(-1, 0, 1), format="csc")
    system = csc_matrix(diags(np.ones(n_space), 0, format="csc") - dt * operator)
    factor = splu(system)
    forcing = np.zeros(n_space, dtype=np.float64)
    forcing[-1] = 2.0 / dx - lambda_f

    value = np.logaddexp(beta * r, -beta * r) / beta - math.log(2.0) / beta
    output = np.empty((n_output, n_space), dtype=np.float64)
    output[0] = value
    stride = n_steps // (n_output - 1)
    output_index = 1
    for step in range(1, n_steps + 1):
        value = factor.solve(value + dt * forcing)
        if step % stride == 0:
            output[output_index] = value
            output_index += 1
    tau = np.linspace(0.0, horizon, n_output, dtype=np.float64)
    return tau, r, output


def _bilinear_uniform(
    values: np.ndarray,
    tau_grid: np.ndarray,
    r_grid: np.ndarray,
    tau: np.ndarray,
    r: np.ndarray,
) -> np.ndarray:
    tt = np.clip(np.asarray(tau, dtype=np.float64), tau_grid[0], tau_grid[-1])
    rr = np.clip(np.asarray(r, dtype=np.float64), r_grid[0], r_grid[-1])
    ti = (tt - tau_grid[0]) / (tau_grid[-1] - tau_grid[0]) * (len(tau_grid) - 1)
    ri = (rr - r_grid[0]) / (r_grid[-1] - r_grid[0]) * (len(r_grid) - 1)
    t0 = np.minimum(np.floor(ti).astype(int), len(tau_grid) - 2)
    r0 = np.minimum(np.floor(ri).astype(int), len(r_grid) - 2)
    at = ti - t0
    ar = ri - r0
    return (
        (1.0 - at) * (1.0 - ar) * values[t0, r0]
        + at * (1.0 - ar) * values[t0 + 1, r0]
        + (1.0 - at) * ar * values[t0, r0 + 1]
        + at * ar * values[t0 + 1, r0 + 1]
    )


def _build_norm_reference_interpolated_legacy(
    path: Path,
    *,
    beta: float = 2.0,
    lambda_f: float = 1.0,
    horizon: float = 0.25,
    domain: float = 12.0,
) -> dict[str, Any]:
    """Build five monotone refinements and save a Richardson reference."""

    path.parent.mkdir(parents=True, exist_ok=True)
    levels = (
        (2401, 512),
        (4801, 1024),
        (9601, 2048),
        (19201, 4096),
        (38401, 8192),
    )
    solutions = [
        _monotone_half_line(
            beta=beta,
            lambda_f=lambda_f,
            horizon=horizon,
            domain=domain,
            n_space=n_space,
            n_steps=n_steps,
        )
        for n_space, n_steps in levels
    ]
    tau = solutions[-1][0]
    coarse_r, coarse = solutions[0][1], solutions[0][2]
    fine_r, fine = solutions[1][1], solutions[1][2]
    ultra_r, ultra = solutions[2][1], solutions[2][2]
    super_r, super_fine = solutions[3][1], solutions[3][2]
    hyper_r, hyper_fine = solutions[4][1], solutions[4][2]

    coarse_on_fine = np.empty_like(fine)
    for index in range(len(tau)):
        coarse_on_fine[index] = np.interp(fine_r, coarse_r, coarse[index])
    extrapolated_cf = 2.0 * fine - coarse_on_fine

    fine_on_ultra = np.empty_like(ultra)
    for index in range(len(tau)):
        fine_on_ultra[index] = np.interp(ultra_r, fine_r, fine[index])
    extrapolated_fu = 2.0 * ultra - fine_on_ultra

    extrapolated_cf_on_ultra = np.empty_like(ultra)
    for index in range(len(tau)):
        extrapolated_cf_on_ultra[index] = np.interp(
            ultra_r, fine_r, extrapolated_cf[index]
        )
    second_order_cfu = (4.0 * extrapolated_fu - extrapolated_cf_on_ultra) / 3.0

    ultra_on_super = np.empty_like(super_fine)
    extrapolated_fu_on_super = np.empty_like(super_fine)
    for index in range(len(tau)):
        ultra_on_super[index] = np.interp(super_r, ultra_r, ultra[index])
        extrapolated_fu_on_super[index] = np.interp(
            super_r, ultra_r, extrapolated_fu[index]
        )
    extrapolated_us = 2.0 * super_fine - ultra_on_super
    second_order_fus = (4.0 * extrapolated_us - extrapolated_fu_on_super) / 3.0

    second_order_cfu_on_super = np.empty_like(super_fine)
    for index in range(len(tau)):
        second_order_cfu_on_super[index] = np.interp(
            super_r, ultra_r, second_order_cfu[index]
        )
    third_order_cfus = (8.0 * second_order_fus - second_order_cfu_on_super) / 7.0

    super_on_hyper = np.empty_like(hyper_fine)
    extrapolated_us_on_hyper = np.empty_like(hyper_fine)
    second_order_fus_on_hyper = np.empty_like(hyper_fine)
    for index in range(len(tau)):
        super_on_hyper[index] = np.interp(hyper_r, super_r, super_fine[index])
        extrapolated_us_on_hyper[index] = np.interp(
            hyper_r, super_r, extrapolated_us[index]
        )
        second_order_fus_on_hyper[index] = np.interp(
            hyper_r, super_r, second_order_fus[index]
        )
    extrapolated_sh = 2.0 * hyper_fine - super_on_hyper
    second_order_ush = (4.0 * extrapolated_sh - extrapolated_us_on_hyper) / 3.0
    third_order_fush = (8.0 * second_order_ush - second_order_fus_on_hyper) / 7.0

    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 4, 2]))
    test_tau = rng.uniform(0.0, horizon, size=10_000)
    test_r = rng.uniform(0.0, min(8.0, domain - 1.0), size=10_000)
    cfus = _bilinear_uniform(third_order_cfus, tau, super_r, test_tau, test_r)
    fush = _bilinear_uniform(third_order_fush, tau, hyper_r, test_tau, test_r)
    raw_super = _bilinear_uniform(super_fine, tau, super_r, test_tau, test_r)
    raw_hyper = _bilinear_uniform(hyper_fine, tau, hyper_r, test_tau, test_r)
    max_refinement = float(np.max(np.abs(fush - cfus)))
    raw_refinement = float(np.max(np.abs(raw_hyper - raw_super)))

    dr = float(hyper_r[1] - hyper_r[0])
    derivative = np.gradient(third_order_fush, dr, axis=1, edge_order=2)
    derivative[:, 0] = 0.0
    derivative[:, -1] = 1.0

    # The converged solve is performed on 38,401 points, but a cubic spline
    # does not need that density at evaluation time.  Store every fourth point
    # and audit the compression directly against the full converged spline.
    storage_stride = 4
    stored_r = hyper_r[::storage_stride]
    stored_value = third_order_fush[:, ::storage_stride]
    stored_derivative = derivative[:, ::storage_stride]
    full_spline = RectBivariateSpline(tau, hyper_r, third_order_fush, kx=3, ky=3, s=0.0)
    stored_spline = RectBivariateSpline(tau, stored_r, stored_value, kx=3, ky=3, s=0.0)
    compression_full = full_spline.ev(test_tau, test_r)
    compression_stored = stored_spline.ev(test_tau, test_r)
    compression_error = float(np.max(np.abs(compression_full - compression_stored)))
    # Keep the complete auditable reference-error budget below G1's 1e-6
    # tolerance: refinement disagreement plus storage interpolation error.
    combined_reference_error = max_refinement + compression_error
    metadata = {
        "schema_version": 1,
        "scheme": "monotone backward Euler; first-order upwind for |psi_s| on even half-line",
        "reference": "third Richardson extrapolation of the four finest monotone refinements",
        "beta": beta,
        "lambda_f": lambda_f,
        "T": horizon,
        "domain": [-domain, domain],
        "levels": [
            {"n_space": n_space, "n_steps": n_steps} for n_space, n_steps in levels
        ],
        "raw_superfine_hyperfine_max_difference": raw_refinement,
        "richardson_refinement_max_difference": max_refinement,
        "refinement_threshold": 1e-6,
        "refinement_passed": max_refinement <= 1e-6,
        "storage_stride": storage_stride,
        "stored_space_points": len(stored_r),
        "storage_spline_max_difference": compression_error,
        "combined_reference_error_bound": combined_reference_error,
        "storage_spline_threshold": max(1e-6 - max_refinement, 0.0),
        "storage_spline_passed": combined_reference_error <= 1e-6,
    }
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        tau=tau,
        r=stored_r,
        value=stored_value,
        derivative=stored_derivative,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    temporary.replace(path)
    return metadata


def build_norm_reference(
    path: Path,
    *,
    beta: float = 2.0,
    lambda_f: float = 1.0,
    horizon: float = 0.25,
    domain: float = 12.0,
) -> dict[str, Any]:
    """Build a smooth, independently refined monotone reference.

    Richardson combinations are formed only at nodes shared by every grid in
    a quartet.  Combining values after interpolating them to the finest grid
    creates an odd/even texture that is tiny in value but large in the second
    derivative; the common-node construction avoids that failure mode.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    levels = (
        (2401, 512),
        (4801, 1024),
        (9601, 2048),
        (19201, 4096),
        (38401, 8192),
        (76801, 16384),
    )
    solutions = [
        _monotone_half_line(
            beta=beta,
            lambda_f=lambda_f,
            horizon=horizon,
            domain=domain,
            n_space=n_space,
            n_steps=n_steps,
            n_output=513,
        )
        for n_space, n_steps in levels
    ]
    values = [solution[2] for solution in solutions]

    def third_order_common(start: int) -> np.ndarray:
        """Third Richardson estimate on ``values[start]``'s nodes."""

        u0, u1, u2, u3 = values[start : start + 4]
        r01 = 2.0 * u1[:, ::2] - u0
        r12 = 2.0 * u2[:, ::4] - u1[:, ::2]
        r23 = 2.0 * u3[:, ::8] - u2[:, ::4]
        s012 = (4.0 * r12 - r01) / 3.0
        s123 = (4.0 * r23 - r12) / 3.0
        return (8.0 * s123 - s012) / 7.0

    coarser_reference = third_order_common(1)
    reference = third_order_common(2)
    tau = solutions[2][0]
    stored_r = solutions[2][1]
    max_refinement = float(
        np.max(np.abs(reference[:, ::2] - coarser_reference))
    )
    raw_refinement = float(np.max(np.abs(values[5][:, ::2] - values[4])))
    metadata = {
        "schema_version": 2,
        "scheme": "monotone backward Euler; first-order upwind for |psi_s| on even half-line",
        "reference": "common-node third Richardson extrapolation of the four finest monotone refinements",
        "beta": beta,
        "lambda_f": lambda_f,
        "T": horizon,
        "domain": [-domain, domain],
        "levels": [
            {"n_space": n_space, "n_steps": n_steps} for n_space, n_steps in levels
        ],
        "raw_finest_pair_max_difference": raw_refinement,
        "richardson_refinement_max_difference": max_refinement,
        "refinement_threshold": 1e-6,
        "refinement_passed": max_refinement <= 1e-6,
        "stored_time_points": len(tau),
        "stored_space_points": len(stored_r),
        "storage_stride": 1,
        "storage_spline_max_difference": 0.0,
        "combined_reference_error_bound": max_refinement,
        "storage_spline_threshold": max(1e-6 - max_refinement, 0.0),
        "storage_spline_passed": max_refinement <= 1e-6,
    }
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        tau=tau,
        r=stored_r,
        value=reference,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    temporary.replace(path)
    return metadata


@dataclass
class _LoadedReference:
    tau: np.ndarray
    r: np.ndarray
    value: np.ndarray
    metadata: dict[str, Any]
    value_spline: RectBivariateSpline
    derivative_spline: RectBivariateSpline
    derivative_cache_max_difference: float


@lru_cache(maxsize=8)
def load_norm_reference(path_text: str) -> _LoadedReference:
    path = Path(path_text)
    with np.load(path, allow_pickle=False) as data:
        tau = data["tau"]
        r = data["r"]
        value = data["value"]
        metadata = json.loads(str(data["metadata_json"]))
    # The even extension encodes the Neumann condition at r=0 in the cubic
    # spline itself, avoiding one-sided endpoint artifacts in the gradient.
    signed_r = np.concatenate((-r[:0:-1], r))
    signed_value = np.concatenate((value[:, :0:-1], value), axis=1)
    value_spline = RectBivariateSpline(
        tau, signed_r, signed_value, kx=3, ky=3, s=0.0
    )
    # FITPACK's pointwise derivative evaluator has high fixed cost and is
    # called hundreds of times inside one MLP repetition.  Differentiate once
    # on the uniform reference grid, then interpolate that odd derivative as
    # an ordinary value spline.  This is the derivative of the same cubic
    # reference, not a second numerical PDE solve.
    radial_derivative = value_spline(tau, r, dx=0, dy=1, grid=True)
    radial_derivative[:, 0] = 0.0
    signed_derivative = np.concatenate(
        (-radial_derivative[:, :0:-1], radial_derivative), axis=1
    )
    derivative_spline = RectBivariateSpline(
        tau, signed_r, signed_derivative, kx=3, ky=3, s=0.0
    )
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 4, 3]))
    audit_tau = rng.uniform(tau[0], tau[-1], size=10_000)
    # Cover the full eight-unit reference buffer required by the protocol,
    # rather than only the narrow region occupied by root test points.
    audit_r = rng.uniform(0.0, min(8.0, r[-1]), size=10_000)
    direct_derivative = value_spline.ev(audit_tau, audit_r, dx=0, dy=1)
    cached_derivative = derivative_spline.ev(audit_tau, audit_r)
    derivative_cache_max_difference = float(
        np.max(np.abs(direct_derivative - cached_derivative))
    )
    return _LoadedReference(
        tau=tau,
        r=r,
        value=value,
        metadata=metadata,
        value_spline=value_spline,
        derivative_spline=derivative_spline,
        derivative_cache_max_difference=derivative_cache_max_difference,
    )


@dataclass(frozen=True)
class NormDriverHJB(ExactEquation):
    """P4 with a precomputed, independently refined one-dimensional solution."""

    d: int
    reference_path: str
    beta: float = 2.0
    lambda_f: float = 1.0
    direction_seed: int = BASE_SEED
    T: float = 0.25
    sigma: float = math.sqrt(2.0)
    name: str = "P4_norm_driver_hjb"
    family: str = "norm_hjb"

    def __post_init__(self) -> None:
        rng = np.random.default_rng(np.random.SeedSequence([self.direction_seed, self.d, 4]))
        w = rng.standard_normal(self.d, dtype=np.float64)
        w /= np.linalg.norm(w)
        object.__setattr__(self, "w", w)
        object.__setattr__(self, "mu", 0.0)
        metadata = load_norm_reference(self.reference_path).metadata
        for key, expected in (
            ("beta", self.beta),
            ("lambda_f", self.lambda_f),
            ("T", self.T),
        ):
            if not math.isclose(float(metadata[key]), expected, rel_tol=0.0, abs_tol=1e-14):
                raise ValueError(f"reference {key} does not match equation")

    def _tau_r_sign(self, t: np.ndarray | float, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        time, state = _as_points(t, x, self.d)
        s = np.einsum("...d,d->...", state, self.w)
        tau = np.broadcast_to(self.T - time, np.shape(s))
        return tau, np.abs(s), np.sign(s)

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        tau, r, _ = self._tau_r_sign(time, state)
        flat = load_norm_reference(self.reference_path).value_spline.ev(tau.reshape(-1), r.reshape(-1))
        result = flat.reshape(np.shape(tau))
        at_terminal = np.abs(tau) <= 8.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            terminal = np.logaddexp(self.beta * r, -self.beta * r) / self.beta - math.log(2.0) / self.beta
            result = np.where(at_terminal, terminal, result)
        return result

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        tau, r, sign = self._tau_r_sign(t, x)
        reference = load_norm_reference(self.reference_path)
        flat_tau = tau.reshape(-1)
        flat_r = r.reshape(-1)
        u = reference.value_spline.ev(flat_tau, flat_r).reshape(np.shape(tau))
        radial = reference.derivative_spline.ev(flat_tau, flat_r).reshape(np.shape(tau))
        at_terminal = np.abs(tau) <= 8.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            terminal = (
                np.logaddexp(self.beta * r, -self.beta * r) / self.beta
                - math.log(2.0) / self.beta
            )
            u = np.where(at_terminal, terminal, u)
            radial = np.where(at_terminal, np.tanh(self.beta * r), radial)
        z = self.sigma * (sign * radial)[..., None] * self.w
        return np.concatenate((u[..., None], z), axis=-1)

    def exact_grad(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        tau, r, sign = self._tau_r_sign(t, x)
        radial = load_norm_reference(self.reference_path).derivative_spline.ev(
            tau.reshape(-1), r.reshape(-1)
        ).reshape(np.shape(tau))
        at_terminal = np.abs(tau) <= 8.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            radial = np.where(at_terminal, np.tanh(self.beta * r), radial)
        return (sign * radial)[..., None] * self.w

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.exact_grad(t, x)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        s = np.einsum("...d,d->...", state, self.w)
        return np.logaddexp(self.beta * s, -self.beta * s) / self.beta - math.log(2.0) / self.beta

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        return -(self.lambda_f / self.sigma) * np.linalg.norm(
            np.asarray(z, dtype=np.float64), axis=-1
        )

    def fzero_reference(self, t: np.ndarray, x: np.ndarray, quadrature_order: int = 80) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        s = state @ self.w
        h = np.maximum(self.T - time, 0.0)
        shifted = s[:, None] + 2.0 * np.sqrt(h)[:, None] * nodes[None, :]
        terminal = (
            np.logaddexp(self.beta * shifted, -self.beta * shifted) / self.beta
            - math.log(2.0) / self.beta
        )
        return np.sum(terminal * weights[None, :], axis=1) / math.sqrt(math.pi)

    @property
    def segment_coefficients(self) -> tuple[float, float]:
        return (-1.0, 1.0)

    @property
    def segment_endpoints(self) -> tuple[np.ndarray, np.ndarray]:
        return -self.sigma * self.w, self.sigma * self.w

    @property
    def box_low(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.minimum(first, second)

    @property
    def box_high(self) -> np.ndarray:
        first, second = self.segment_endpoints
        return np.maximum(first, second)

    @property
    def z_ball_radius(self) -> float:
        return self.sigma

    @property
    def z_center(self) -> np.ndarray:
        return np.zeros(self.d, dtype=np.float64)

    @property
    def dose_scale(self) -> float:
        return self.sigma

    @property
    def picard_proxy(self) -> float:
        return self.T * self.lambda_f / self.sigma

    def certificate_labels(self) -> dict[str, str]:
        return {
            "span": "PDE-derived: translation invariance orthogonal to w",
            "segment": "PDE-derived: Lipschitz preservation |psi_s|<=sup|g'|=1",
            "box": "PDE-derived: coordinatewise hull of the certified segment",
            "ball": "PDE-derived: Euclidean hull of the certified segment",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        reference = load_norm_reference(self.reference_path)
        reference_path = Path(self.reference_path)
        try:
            anchor = reference_path.parts.index("invariant_region_mlp")
            portable_reference_path = Path(*reference_path.parts[anchor:]).as_posix()
        except ValueError:
            portable_reference_path = reference_path.as_posix()
        return {
            "w_seed": self.direction_seed,
            "w_sha256": hashlib.sha256(self.w.tobytes()).hexdigest(),
            "segment_coefficients": [-1.0, 1.0],
            "box_low": self.box_low.tolist(),
            "box_high": self.box_high.tolist(),
            "ball_radius": self.z_ball_radius,
            "reference_path": portable_reference_path,
            "reference_metadata": reference.metadata,
            "gradient_evaluator": "cubic interpolation of the value-spline derivative sampled on the reference grid",
            "gradient_cache_max_difference": reference.derivative_cache_max_difference,
            "gradient_cache_threshold": 1e-6,
        }

    def reference_residual(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        """PDE residual from independent derivatives of the reference spline."""

        tau, r, _ = self._tau_r_sign(t, x)
        spline = load_norm_reference(self.reference_path).value_spline
        flat_tau = tau.reshape(-1)
        flat_r = r.reshape(-1)
        psi_tau = spline.ev(flat_tau, flat_r, dx=1, dy=0)
        psi_s = spline.ev(flat_tau, flat_r, dx=0, dy=1)
        psi_ss = spline.ev(flat_tau, flat_r, dx=0, dy=2)
        residual = -psi_tau + psi_ss - self.lambda_f * np.abs(psi_s)
        return residual.reshape(np.shape(tau))

    def provenance(self) -> dict[str, Any]:
        return {**super().provenance(), "beta": self.beta, "lambda_f": self.lambda_f}
