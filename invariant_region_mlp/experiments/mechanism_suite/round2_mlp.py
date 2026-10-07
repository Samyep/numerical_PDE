"""Round-2 certificate transforms over the unchanged FullHistoryMLP."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

from .mechanism_mlp import (
    ExtraDiagnostics,
    MechanismMLP,
    NoiseCorrelation,
    SEGMENT,
    WorkDiagnostics,
    _box_energy,
    _u_violation_energy,
    error_metrics,
)
from .round2_certificate import p4_derivative_bounds_cached

# MethodSpec is imported by mechanism_mlp from the unchanged implementation.
from .mechanism_mlp import MethodSpec  # noqa: E402


ROUND2_IMPLEMENTATION_REVISION = "mechanism_suite_r2_v2"
ROUND2_PREREGISTRATION_COMMIT = "f83ce81fbeaff6b7c12aeec90906fc6e6124a71e"

TIGHT_SEGMENT = MethodSpec("tight_segment", "tight_segment", 1.0)
TIGHT_CENTRE = MethodSpec("centre", "tight_centre", 0.0)


def tightness_method(theta: float) -> MethodSpec:
    theta = float(theta)
    if math.isclose(theta, 0.0, rel_tol=0.0, abs_tol=0.0):
        return SEGMENT
    if math.isclose(theta, 1.0, rel_tol=0.0, abs_tol=0.0):
        return TIGHT_SEGMENT
    return MethodSpec(
        f"tight_theta{theta:g}", "tight_segment", theta
    )


def tight_illegal_method(factor: float) -> MethodSpec:
    return MethodSpec(
        f"illegal_a{float(factor):g}", "tight_illegal", float(factor)
    )


@dataclass
class Round2ExtraDiagnostics(ExtraDiagnostics):
    """Add the z-noise denominator required by corrected localization."""

    s_bin_z_count: np.ndarray = field(
        default_factory=lambda: np.zeros(7, dtype=np.int64)
    )
    s_bin_z_squared_error_sum: np.ndarray = field(
        default_factory=lambda: np.zeros(7, dtype=np.float64)
    )

    def add_s_bins(
        self,
        s: np.ndarray,
        error: np.ndarray,
        *,
        z_error: np.ndarray | None = None,
    ) -> None:
        super().add_s_bins(s, error)
        if z_error is None:
            return
        scalar = np.asarray(s, dtype=np.float64).reshape(-1)
        squared_norm = np.sum(
            np.asarray(z_error, dtype=np.float64) ** 2, axis=-1
        ).reshape(-1)
        indices = np.digitize(scalar, self.s_bin_edges[1:-1], right=False)
        finite = np.isfinite(squared_norm)
        for index in range(len(self.s_bin_z_count)):
            mask = (indices == index) & finite
            if np.any(mask):
                self.s_bin_z_count[index] += int(np.count_nonzero(mask))
                self.s_bin_z_squared_error_sum[index] += float(
                    np.sum(squared_norm[mask])
                )

    def merge(self, other: "Round2ExtraDiagnostics") -> None:
        super().merge(other)
        self.s_bin_z_count += other.s_bin_z_count
        self.s_bin_z_squared_error_sum += other.s_bin_z_squared_error_sum

    def summary(self) -> dict[str, Any]:
        payload = super().summary()
        for index, item in enumerate(payload["generator_bias_by_s_bin"]):
            count = int(self.s_bin_z_count[index])
            item["z_error_count"] = count
            item["z_error_rms"] = (
                float(
                    math.sqrt(
                        self.s_bin_z_squared_error_sum[index] / count
                    )
                )
                if count
                else None
            )
            item["z_squared_error_sum"] = float(
                self.s_bin_z_squared_error_sum[index]
            )
        return payload


def _project_p4_tight(
    equation: Any,
    state: np.ndarray,
    t: np.ndarray,
    x: np.ndarray,
    *,
    theta: float,
    illegal_factor: float | None = None,
    centre_only: bool = False,
) -> np.ndarray:
    if equation.family != "norm_hjb":
        raise ValueError("tight P4 transforms require the norm_hjb equation")
    corrected = np.asarray(state, dtype=np.float64).copy()
    time = np.asarray(t, dtype=np.float64)
    location = np.einsum(
        "...d,d->...", np.asarray(x, dtype=np.float64), equation.w
    )
    tau = np.maximum(equation.T - time, 0.0)
    v_minus, v_plus = p4_derivative_bounds_cached(
        tau, location, lambda_f=equation.lambda_f
    )
    lower = (1.0 - theta) * -1.0 + theta * v_minus
    upper = (1.0 - theta) * 1.0 + theta * v_plus
    centre = 0.5 * (lower + upper)
    if illegal_factor is not None:
        half_width = 0.5 * (upper - lower) * illegal_factor
        lower = centre - half_width
        upper = centre + half_width
    coefficient = np.einsum(
        "...d,d->...", corrected[..., 1:], equation.w
    ) / equation.sigma
    if centre_only:
        coefficient = centre
    else:
        coefficient = np.clip(coefficient, lower, upper)
    corrected[..., 1:] = equation.sigma * coefficient[..., None] * equation.w
    return corrected


class Round2MechanismMLP(MechanismMLP):
    """MechanismMLP with state-dependent P4 certificates."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.extra = Round2ExtraDiagnostics()

    def _round2_transform(
        self,
        state: np.ndarray,
        exact: np.ndarray,
        t: np.ndarray,
        x: np.ndarray,
    ) -> np.ndarray:
        transform = self.method.transform
        if transform not in {"tight_segment", "tight_illegal", "tight_centre"}:
            return super()._transform(state, exact)

        candidate = self._candidate(state, exact)
        box_energy = _box_energy(self.equation, candidate)
        z = candidate[..., 1:]
        ball_energy = _u_violation_energy(
            self.equation, candidate[..., 0]
        ) + np.maximum(
            np.linalg.norm(z, axis=-1) - self.equation.z_ball_radius, 0.0
        ) ** 2
        count = candidate[..., 0].size
        self.stats.correction_opportunities += count
        self.stats.pre_box_violations += int(
            np.count_nonzero(box_energy > 0.0)
        )
        self.stats.pre_ball_violations += int(
            np.count_nonzero(ball_energy > 0.0)
        )
        self.stats.box_overshoot_energy_sum += float(np.sum(box_energy))

        if transform == "tight_segment":
            corrected = _project_p4_tight(
                self.equation,
                candidate,
                t,
                x,
                theta=float(self.method.factor),
            )
        elif transform == "tight_illegal":
            corrected = _project_p4_tight(
                self.equation,
                candidate,
                t,
                x,
                theta=1.0,
                illegal_factor=float(self.method.factor),
            )
        else:
            corrected = _project_p4_tight(
                self.equation,
                candidate,
                t,
                x,
                theta=1.0,
                centre_only=True,
            )

        changed = np.any(corrected != candidate, axis=-1)
        method_energy = np.sum((corrected - candidate) ** 2, axis=-1)
        self.stats.method_geometry_violations += int(
            np.count_nonzero(method_energy > 0.0)
        )
        self.stats.method_overshoot_energy_sum += float(np.sum(method_energy))
        self.stats.activated_states += int(np.count_nonzero(changed))
        return corrected

    def _generator(
        self, state: np.ndarray, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        b, k, width = state.shape
        if width != self.equation.d + 1:
            raise ValueError("invalid state width")
        exact = self.equation.exact_state(t, x)
        z_error = state[..., 1:] - exact[..., 1:]
        self.extra.noise.add(state[..., 0] - exact[..., 0], z_error)
        corrected = self._round2_transform(state, exact, t, x)
        estimate = self.equation.generator(
            corrected[..., 0], corrected[..., 1:]
        )
        truth = self.equation.generator(exact[..., 0], exact[..., 1:])
        error = estimate - truth
        count = b * k
        finite = np.isfinite(estimate) & np.isfinite(error)
        self.stats.f_evals += count
        self.stats.generator_count += count
        self.stats.nonfinite_generator_values += int(
            count - np.count_nonzero(finite)
        )
        self.stats.generator_error_sum += float(np.sum(error))
        self.stats.generator_abs_error_sum += float(np.sum(np.abs(error)))
        self.stats.generator_squared_error_sum += float(np.sum(error**2))
        self.stats.generator_estimate_sum += float(np.sum(estimate))
        self.stats.generator_truth_sum += float(np.sum(truth))
        if hasattr(self.equation, "w"):
            s = np.einsum("...d,d->...", x, self.equation.w)
            self.extra.add_s_bins(s, error, z_error=z_error)
        return estimate


def run_single_repetition_r2(
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
    study: str,
    time_beta_alpha: float = 0.5,
    trace_draws: bool = False,
) -> dict[str, Any]:
    """Run one fresh round-2 repetition with method-independent tree seeds."""

    t = np.asarray(t, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    is_validation = np.asarray(is_validation, dtype=bool)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    diagnostics = WorkDiagnostics()
    extras = Round2ExtraDiagnostics()
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
        solver = Round2MechanismMLP(
            equation,
            M,
            method,
            np.random.default_rng(tree_seed),
            dose_rng=np.random.default_rng(dose_seed),
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        prediction[begin:end] = solver.solve(
            n, t[begin:end], x[begin:end]
        )
        diagnostics.merge(solver.stats)
        extras.merge(solver.extra)
        if solver.draw_fingerprint is not None:
            fingerprints.append(solver.draw_fingerprint)
    elapsed = time.perf_counter() - started
    all_mask = np.ones(len(t), dtype=bool)
    test_mask = ~is_validation
    metadata = {
        "schema_version": 2,
        "round": 2,
        "study_origin": study,
        "implementation_revision": ROUND2_IMPLEMENTATION_REVISION,
        "preregistration_commit": ROUND2_PREREGISTRATION_COMMIT,
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
            "validation": error_metrics(
                prediction, truth, is_validation
            ),
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


def save_round2_repetition(path: Path, result: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        prediction_u=np.asarray(result["prediction_u"], dtype=np.float64),
        truth_u=np.asarray(result["truth_u"], dtype=np.float64),
        is_validation=np.asarray(result["is_validation"], dtype=bool),
        metadata_json=np.asarray(
            json.dumps(result["metadata"], sort_keys=True)
        ),
    )
    temporary.replace(path)
