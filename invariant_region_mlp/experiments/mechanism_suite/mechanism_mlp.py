"""Mechanism-suite adapters around the unchanged active-VB FullHistoryMLP."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any, Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
VB_DIR = HERE.parent / "active_vb_high_budget"
if str(VB_DIR) not in sys.path:
    sys.path.insert(0, str(VB_DIR))

# This class is imported and reused unchanged, as required by the protocol.
from vb_mlp_methods import FullHistoryMLP, MethodSpec, WorkDiagnostics  # noqa: E402


RAW = MethodSpec("raw", "raw")
BOX = MethodSpec("box", "box")
SEGMENT = MethodSpec("segment", "segment")
BALL = MethodSpec("ball", "ball")
SIGN_ONLY = MethodSpec("sign_only", "sign_only")
SPAN_ONLY = MethodSpec("span_only", "span_only")
BATCH_BOX = MethodSpec("batch_box", "batch_box")
ORACLE_STATE = MethodSpec("oracle_state", "oracle_state")
ORACLE_Z = MethodSpec("oracle_z", "oracle_z")
ORACLE_U = MethodSpec("oracle_u", "oracle_u")
Z_ZERO = MethodSpec("z_zero", "z_zero", 0.0)
F_ZERO = MethodSpec("f_zero", "f_zero", 0.0)
CENTRE = MethodSpec("centre", "centre", 0.0)
U_DOSE_SCALE = 0.5  # half-width of the certified P3 interval u in [0,1]


def dose_method(scale: float, *, projected: str | None = None) -> MethodSpec:
    if projected is None:
        return MethodSpec(f"dose_s{scale:g}", "dose", float(scale))
    return MethodSpec(f"dose_{projected}_s{scale:g}", f"dose_{projected}", float(scale))


def dose_u_method(scale: float, *, clipped: bool = False) -> MethodSpec:
    return MethodSpec(
        f"dose_u{'_box' if clipped else ''}_s{scale:g}",
        "dose_u_box" if clipped else "dose_u",
        float(scale),
    )


def shrink_method(factor: float) -> MethodSpec:
    return MethodSpec(f"shrink_c{factor:g}", "shrink", float(factor))


def illegal_method(factor: float) -> MethodSpec:
    return MethodSpec(f"illegal_a{factor:g}", "illegal", float(factor))


def looser_method(factor: float) -> MethodSpec:
    return MethodSpec(f"looser_a{factor:g}", "looser", float(factor))


def shrink_centre_method(factor: float) -> MethodSpec:
    return MethodSpec(f"shrink_centre_c{factor:g}", "shrink_centre", float(factor))


SHRINK_FACTORS = (0.1, 0.25, 0.5, 0.75)
ILLEGAL_FACTORS = (0.25, 0.5, 0.75, 0.9)
LOOSER_FACTORS = (1.25, 1.5, 2.0)
SHRINK_CENTRE_FACTORS = (0.25, 0.5, 0.75)


def tuning_methods() -> list[MethodSpec]:
    return (
        [shrink_method(value) for value in SHRINK_FACTORS]
        + [illegal_method(value) for value in ILLEGAL_FACTORS]
        + [looser_method(value) for value in LOOSER_FACTORS]
        + [shrink_centre_method(value) for value in SHRINK_CENTRE_FACTORS]
    )


def method_from_dict(payload: dict[str, Any]) -> MethodSpec:
    return MethodSpec(
        str(payload["name"]),
        str(payload["transform"]),
        float(payload.get("factor", 1.0)),
    )


@dataclass
class NoiseCorrelation:
    count: int = 0
    sum_u: float = 0.0
    sum_z: float = 0.0
    sum_u2: float = 0.0
    sum_z2: float = 0.0
    sum_uz: float = 0.0

    def add(self, u_error: np.ndarray, z_error: np.ndarray) -> None:
        u = np.asarray(u_error, dtype=np.float64).reshape(-1)
        # The scalar gradient noise channel is the mean coordinate error.  For
        # P2/P3 this is proportional to the error in sum(z), which is exactly
        # what enters the product driver.
        z = np.mean(np.asarray(z_error, dtype=np.float64), axis=-1).reshape(-1)
        finite = np.isfinite(u) & np.isfinite(z)
        if not np.any(finite):
            return
        u = u[finite]
        z = z[finite]
        self.count += len(u)
        self.sum_u += float(np.sum(u))
        self.sum_z += float(np.sum(z))
        self.sum_u2 += float(np.sum(u * u))
        self.sum_z2 += float(np.sum(z * z))
        self.sum_uz += float(np.sum(u * z))

    def merge(self, other: "NoiseCorrelation") -> None:
        for name in ("count", "sum_u", "sum_z", "sum_u2", "sum_z2", "sum_uz"):
            setattr(self, name, getattr(self, name) + getattr(other, name))

    def summary(self) -> dict[str, float | int | None]:
        if self.count < 2:
            return {"count": self.count, "correlation": None}
        count = float(self.count)
        covariance = self.sum_uz - self.sum_u * self.sum_z / count
        variance_u = self.sum_u2 - self.sum_u**2 / count
        variance_z = self.sum_z2 - self.sum_z**2 / count
        denominator = math.sqrt(max(variance_u, 0.0) * max(variance_z, 0.0))
        correlation = covariance / denominator if denominator > 0.0 else None
        return {
            "count": self.count,
            "mean_u_noise": self.sum_u / count,
            "mean_coordinate_z_noise": self.sum_z / count,
            "correlation": correlation,
        }


@dataclass
class ExtraDiagnostics:
    noise: NoiseCorrelation = field(default_factory=NoiseCorrelation)
    dose_normal_variates: int = 0
    s_bin_edges: np.ndarray = field(
        default_factory=lambda: np.array(
            [-np.inf, -2.0, -1.0, -0.25, 0.25, 1.0, 2.0, np.inf],
            dtype=np.float64,
        )
    )
    s_bin_count: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=np.int64))
    s_bin_error_sum: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=np.float64))
    s_bin_abs_error_sum: np.ndarray = field(default_factory=lambda: np.zeros(7, dtype=np.float64))

    def add_s_bins(self, s: np.ndarray, error: np.ndarray) -> None:
        scalar = np.asarray(s, dtype=np.float64).reshape(-1)
        difference = np.asarray(error, dtype=np.float64).reshape(-1)
        indices = np.digitize(scalar, self.s_bin_edges[1:-1], right=False)
        finite = np.isfinite(difference)
        for index in range(len(self.s_bin_count)):
            mask = (indices == index) & finite
            if np.any(mask):
                self.s_bin_count[index] += int(np.count_nonzero(mask))
                self.s_bin_error_sum[index] += float(np.sum(difference[mask]))
                self.s_bin_abs_error_sum[index] += float(np.sum(np.abs(difference[mask])))

    def merge(self, other: "ExtraDiagnostics") -> None:
        self.noise.merge(other.noise)
        self.dose_normal_variates += other.dose_normal_variates
        self.s_bin_count += other.s_bin_count
        self.s_bin_error_sum += other.s_bin_error_sum
        self.s_bin_abs_error_sum += other.s_bin_abs_error_sum

    def summary(self) -> dict[str, Any]:
        bins = []
        for index, count in enumerate(self.s_bin_count):
            bins.append(
                {
                    "low": float(self.s_bin_edges[index]),
                    "high": float(self.s_bin_edges[index + 1]),
                    "count": int(count),
                    "bias": float(self.s_bin_error_sum[index] / count) if count else None,
                    "mae": float(self.s_bin_abs_error_sum[index] / count) if count else None,
                }
            )
        return {
            "child_u_z_noise": self.noise.summary(),
            "dose_normal_variates": self.dose_normal_variates,
            "generator_bias_by_s_bin": bins,
        }


def _u_violation_energy(equation: Any, u: np.ndarray) -> np.ndarray:
    if not getattr(equation, "has_u_certificate", False):
        return np.zeros_like(u, dtype=np.float64)
    low, high = equation.u_interval
    return np.maximum(low - u, 0.0) ** 2 + np.maximum(u - high, 0.0) ** 2


def _box_energy(equation: Any, state: np.ndarray) -> np.ndarray:
    u = state[..., 0]
    z = state[..., 1:]
    low = np.asarray(equation.box_low, dtype=np.float64)
    high = np.asarray(equation.box_high, dtype=np.float64)
    return (
        _u_violation_energy(equation, u)
        + np.sum(np.maximum(low - z, 0.0) ** 2 + np.maximum(z - high, 0.0) ** 2, axis=-1)
    )


def _clip_u(equation: Any, u: np.ndarray) -> np.ndarray:
    if not getattr(equation, "has_u_certificate", False):
        return u
    low, high = equation.u_interval
    return np.clip(u, low, high)


def _project_box(equation: Any, state: np.ndarray, factor: float = 1.0) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    corrected[..., 0] = _clip_u(equation, corrected[..., 0])
    low = np.asarray(equation.box_low, dtype=np.float64)
    high = np.asarray(equation.box_high, dtype=np.float64)
    if factor != 1.0:
        if getattr(equation, "family", "") in {"vba", "burgers_fisher", "published_vb"}:
            high = factor * high
        else:
            centre = 0.5 * (low + high)
            low = centre + factor * (low - centre)
            high = centre + factor * (high - centre)
    corrected[..., 1:] = np.clip(corrected[..., 1:], low, high)
    return corrected


def _project_ball(equation: Any, state: np.ndarray) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    corrected[..., 0] = _clip_u(equation, corrected[..., 0])
    z = corrected[..., 1:]
    norm = np.linalg.norm(z, axis=-1)
    scale = np.minimum(1.0, equation.z_ball_radius / np.maximum(norm, 1e-300))
    corrected[..., 1:] = z * scale[..., None]
    return corrected


def _project_span(equation: Any, state: np.ndarray) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    w = np.asarray(equation.w, dtype=np.float64)
    z = corrected[..., 1:]
    corrected[..., 1:] = np.einsum("...d,d->...", z, w)[..., None] * w
    return corrected


def _project_segment(equation: Any, state: np.ndarray) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    corrected[..., 0] = _clip_u(equation, corrected[..., 0])
    w = np.asarray(equation.w, dtype=np.float64)
    z = corrected[..., 1:]
    if equation.family == "ridge_lse":
        coefficient = -np.einsum("...d,d->...", z, w) / equation.sigma
        low, high = equation.segment_coefficients
        coefficient = np.clip(coefficient, low, high)
        corrected[..., 1:] = -equation.sigma * coefficient[..., None] * w
    elif equation.family == "norm_hjb":
        coefficient = np.einsum("...d,d->...", z, w) / equation.sigma
        coefficient = np.clip(coefficient, -1.0, 1.0)
        corrected[..., 1:] = equation.sigma * coefficient[..., None] * w
    else:
        raise ValueError(f"segment projection is unavailable for {equation.family}")
    return corrected


class MechanismMLP(FullHistoryMLP):
    """PDE-generic correction layer over the unchanged stochastic recursion."""

    def __init__(
        self,
        equation: Any,
        M: int,
        method: MethodSpec,
        rng: np.random.Generator,
        *,
        dose_rng: np.random.Generator,
        time_beta_alpha: float = 0.5,
        trace_draws: bool = False,
    ) -> None:
        super().__init__(
            equation,
            M,
            method,
            rng,
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        self.dose_rng = dose_rng
        self.extra = ExtraDiagnostics()

    def _correct(self, state: np.ndarray) -> np.ndarray:  # pragma: no cover - guardrail
        raise RuntimeError("VB-specific _correct must not be used in mechanism_suite")

    def _dose_normal(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.dose_rng.standard_normal(shape, dtype=np.float64)
        self.extra.dose_normal_variates += int(np.prod(shape))
        return values

    def _candidate(self, state: np.ndarray, exact: np.ndarray) -> np.ndarray:
        transform = self.method.transform
        candidate = np.asarray(state, dtype=np.float64).copy()
        if transform == "oracle_state":
            candidate = exact.copy()
        elif transform == "oracle_z":
            candidate[..., 1:] = exact[..., 1:]
        elif transform == "oracle_u":
            candidate[..., 0] = exact[..., 0]
        elif transform in {"dose", "dose_box", "dose_segment"}:
            candidate = exact.copy()
            noise = self._dose_normal(candidate[..., 1:].shape)
            candidate[..., 1:] += self.method.factor * self.equation.dose_scale * noise
        elif transform in {"dose_u", "dose_u_box"}:
            candidate = exact.copy()
            noise = self._dose_normal(candidate[..., 0].shape)
            candidate[..., 0] += self.method.factor * U_DOSE_SCALE * noise
        return candidate

    def _transform(self, state: np.ndarray, exact: np.ndarray) -> np.ndarray:
        transform = self.method.transform
        candidate = self._candidate(state, exact)
        box_energy = _box_energy(self.equation, candidate)
        z = candidate[..., 1:]
        ball_energy = _u_violation_energy(self.equation, candidate[..., 0]) + np.maximum(
            np.linalg.norm(z, axis=-1) - self.equation.z_ball_radius, 0.0
        ) ** 2
        count = candidate[..., 0].size
        self.stats.correction_opportunities += count
        self.stats.pre_box_violations += int(np.count_nonzero(box_energy > 0.0))
        self.stats.pre_ball_violations += int(np.count_nonzero(ball_energy > 0.0))
        self.stats.box_overshoot_energy_sum += float(np.sum(box_energy))

        if transform in {"raw", "oracle_state", "oracle_z", "oracle_u", "dose", "dose_u"}:
            corrected = candidate
        elif transform in {"box", "dose_box"}:
            corrected = _project_box(self.equation, candidate)
        elif transform in {"segment", "dose_segment"}:
            corrected = _project_segment(self.equation, candidate)
        elif transform == "ball":
            corrected = _project_ball(self.equation, candidate)
        elif transform == "span_only":
            corrected = _project_span(self.equation, candidate)
        elif transform == "sign_only":
            corrected = candidate.copy()
            corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
            corrected[..., 1:] = np.maximum(corrected[..., 1:], 0.0)
        elif transform == "batch_box":
            corrected = candidate.copy()
            corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
            if self.equation.family in {"vba", "burgers_fisher", "published_vb"}:
                # Preserve the active-VB Batch-IR definition exactly: truncate
                # negative coordinates, then apply one contraction per parent
                # to every sibling so their largest coordinate is admissible.
                nonnegative = np.maximum(corrected[..., 1:], 0.0)
                maximum = np.max(nonnegative, axis=(1, 2))
                alpha = np.minimum(
                    1.0,
                    self.equation.z_upper / np.maximum(maximum, 1e-300),
                )
                corrected[..., 1:] = nonnegative * alpha[:, None, None]
            else:
                low = np.asarray(self.equation.box_low, dtype=np.float64)
                high = np.asarray(self.equation.box_high, dtype=np.float64)
                centre = 0.5 * (low + high)
                half = 0.5 * (high - low)
                deviation = corrected[..., 1:] - centre
                normalized = np.abs(deviation) / np.maximum(half, 1e-300)
                maximum = np.max(normalized, axis=(1, 2))
                alpha = np.minimum(1.0, 1.0 / np.maximum(maximum, 1e-300))
                corrected[..., 1:] = centre + alpha[:, None, None] * deviation
            self._record_alpha(alpha)
        elif transform == "z_zero":
            corrected = candidate.copy()
            corrected[..., 1:] = 0.0
        elif transform == "shrink":
            corrected = candidate.copy()
            corrected[..., 1:] = self.method.factor * corrected[..., 1:]
        elif transform == "illegal":
            corrected = _project_box(self.equation, candidate, self.method.factor)
        elif transform == "looser":
            corrected = _project_box(self.equation, candidate, self.method.factor)
        elif transform == "centre":
            corrected = candidate.copy()
            corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
            corrected[..., 1:] = self.equation.z_center
        elif transform == "shrink_centre":
            corrected = candidate.copy()
            corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
            centre = np.asarray(self.equation.z_center, dtype=np.float64)
            corrected[..., 1:] = centre + self.method.factor * (corrected[..., 1:] - centre)
        elif transform == "dose_u_box":
            corrected = candidate.copy()
            corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
        else:
            raise ValueError(f"unsupported mechanism transform {transform!r}")

        changed = np.any(corrected != candidate, axis=-1)
        method_energy = np.sum((corrected - candidate) ** 2, axis=-1)
        self.stats.method_geometry_violations += int(np.count_nonzero(method_energy > 0.0))
        self.stats.method_overshoot_energy_sum += float(np.sum(method_energy))
        self.stats.activated_states += int(np.count_nonzero(changed))
        return corrected

    def _generator(self, state: np.ndarray, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        b, k, width = state.shape
        if width != self.equation.d + 1:
            raise ValueError("invalid state width")
        exact = self.equation.exact_state(t, x)
        self.extra.noise.add(state[..., 0] - exact[..., 0], state[..., 1:] - exact[..., 1:])
        corrected = self._transform(state, exact)
        estimate = self.equation.generator(corrected[..., 0], corrected[..., 1:])
        # Reuse the exact state already needed by oracle transforms and noise
        # diagnostics.  Numerical ridge references make a second spline pass
        # here particularly expensive, and it is mathematically identical.
        truth = self.equation.generator(exact[..., 0], exact[..., 1:])
        error = estimate - truth
        count = b * k
        finite = np.isfinite(estimate) & np.isfinite(error)
        self.stats.f_evals += count
        self.stats.generator_count += count
        self.stats.nonfinite_generator_values += int(count - np.count_nonzero(finite))
        self.stats.generator_error_sum += float(np.sum(error))
        self.stats.generator_abs_error_sum += float(np.sum(np.abs(error)))
        self.stats.generator_squared_error_sum += float(np.sum(error**2))
        self.stats.generator_estimate_sum += float(np.sum(estimate))
        self.stats.generator_truth_sum += float(np.sum(truth))
        if hasattr(self.equation, "w"):
            s = np.einsum("...d,d->...", x, self.equation.w)
            self.extra.add_s_bins(s, error)
        return estimate


def error_metrics(
    prediction: np.ndarray,
    truth: np.ndarray,
    mask: np.ndarray,
) -> dict[str, float | int]:
    pred = np.asarray(prediction, dtype=np.float64)[mask]
    ref = np.asarray(truth, dtype=np.float64)[mask]
    value_error = pred[:, 0] - ref[:, 0]
    gradient_error = pred[:, 1:] - ref[:, 1:]
    finite_rows = np.all(np.isfinite(pred), axis=1)
    if not np.all(finite_rows):
        value_rmse = float("inf")
        value_skill = float("inf")
        value_relative = float("inf")
        gradient_relative = float("inf")
    else:
        value_rmse = float(np.sqrt(np.mean(value_error**2)))
        reference_std = float(np.std(ref[:, 0], ddof=0))
        value_skill = value_rmse / reference_std if reference_std > 0.0 else float("inf")
        value_relative = float(np.linalg.norm(value_error) / np.linalg.norm(ref[:, 0]))
        gradient_relative = float(
            np.linalg.norm(gradient_error) / np.linalg.norm(ref[:, 1:])
        )
    return {
        "n_points": len(pred),
        "skill": value_skill,
        "value_rmse": value_rmse,
        "value_relative_l2": value_relative,
        "value_bias": float(np.mean(value_error)),
        "value_mae": float(np.mean(np.abs(value_error))),
        "gradient_relative_l2": gradient_relative,
        "gradient_mae": float(np.mean(np.abs(gradient_error))),
        "nonfinite_state_count": int(np.size(pred) - np.count_nonzero(np.isfinite(pred))),
        "nonfinite_row_count": int(len(pred) - np.count_nonzero(finite_rows)),
    }


def run_single_repetition(
    *,
    pde_id: str,
    equation: Any,
    method: MethodSpec,
    n: int,
    M: int,
    repetition: int,
    t: np.ndarray,
    x: np.ndarray,
    is_validation: np.ndarray,
    base_seed: int,
    chunk_size: int,
    time_beta_alpha: float = 0.5,
    trace_draws: bool = False,
) -> dict[str, Any]:
    """Run one method-repetition with method-independent paired tree seeds."""

    t = np.asarray(t, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    is_validation = np.asarray(is_validation, dtype=bool)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    diagnostics = WorkDiagnostics()
    extras = ExtraDiagnostics()
    fingerprints: list[str] = []
    started = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(t), chunk_size)):
        end = min(begin + chunk_size, len(t))
        tree_seed = np.random.SeedSequence(
            [base_seed, equation.d, n, M, repetition, chunk_index]
        )
        dose_seed = np.random.SeedSequence(
            [base_seed, equation.d, n, M, repetition, chunk_index, 0xD05E]
        )
        solver = MechanismMLP(
            equation,
            M,
            method,
            np.random.default_rng(tree_seed),
            dose_rng=np.random.default_rng(dose_seed),
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        prediction[begin:end] = solver.solve(n, t[begin:end], x[begin:end])
        diagnostics.merge(solver.stats)
        extras.merge(solver.extra)
        if solver.draw_fingerprint is not None:
            fingerprints.append(solver.draw_fingerprint)
    elapsed = time.perf_counter() - started
    all_mask = np.ones(len(t), dtype=bool)
    test_mask = ~is_validation
    metadata = {
        "schema_version": 1,
        "implementation_revision": (
            "p4_cached_reference_derivative_v1"
            if equation.family == "norm_hjb"
            else "mechanism_suite_v1"
        ),
        "pde_id": pde_id,
        "equation_name": equation.name,
        "family": equation.family,
        "dimension": equation.d,
        "n": n,
        "M": M,
        "method": method.to_dict(),
        "repetition": repetition,
        "base_seed": base_seed,
        "chunk_size": chunk_size,
        "time_beta_alpha": time_beta_alpha,
        "dtype": "float64",
        "wall_clock_seconds": elapsed,
        "work": diagnostics.summary(),
        "extra_diagnostics": extras.summary(),
        "metrics": {
            "all": error_metrics(prediction, truth, all_mask),
            "validation": error_metrics(prediction, truth, is_validation),
            "test": error_metrics(prediction, truth, test_mask),
        },
        "draw_fingerprint": (
            hashlib.sha256("".join(fingerprints).encode("ascii")).hexdigest()
            if fingerprints
            else None
        ),
        "equation": equation.provenance(),
    }
    return {
        "prediction_u": prediction[:, 0],
        "truth_u": truth[:, 0],
        "is_validation": is_validation,
        "metadata": metadata,
    }


def save_repetition(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        prediction_u=np.asarray(result["prediction_u"], dtype=np.float64),
        truth_u=np.asarray(result["truth_u"], dtype=np.float64),
        is_validation=np.asarray(result["is_validation"], dtype=bool),
        metadata_json=np.asarray(json.dumps(result["metadata"], sort_keys=True)),
    )
    temporary.replace(path)


def load_repetition(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        return {
            "prediction_u": data["prediction_u"],
            "truth_u": data["truth_u"],
            "is_validation": data["is_validation"],
            "metadata": json.loads(str(data["metadata_json"])),
        }


def method_names(methods: Iterable[MethodSpec]) -> list[str]:
    return [method.name for method in methods]
