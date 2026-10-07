"""Round-3 transforms over the unchanged ``FullHistoryMLP`` recursion."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import (
    ExtraDiagnostics,
    MethodSpec,
    WorkDiagnostics,
    _box_energy,
    _u_violation_energy,
    error_metrics,
)
from invariant_region_mlp.experiments.mechanism_suite.round2_mlp import (
    Round2ExtraDiagnostics,
    Round2MechanismMLP,
)


IMPLEMENTATION_REVISION = "mechanism_suite_r3_v1"
PREREGISTRATION_COMMIT = "5326a463a8336c2a56376541aad0ba01a1034dd2"


def _generic_segment(equation: Any, state: np.ndarray) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    w = np.asarray(equation.w, dtype=np.float64)
    low, high = equation.segment_coefficients
    coefficient = np.einsum("...d,d->...", corrected[..., 1:], w) / equation.sigma
    coefficient = np.clip(coefficient, low, high)
    corrected[..., 1:] = equation.sigma * coefficient[..., None] * w
    return corrected


def _lqg_geometry(
    equation: Any,
    state: np.ndarray,
    x: np.ndarray,
    transform: str,
) -> np.ndarray:
    corrected = np.asarray(state, dtype=np.float64).copy()
    location = np.asarray(x, dtype=np.float64)
    z = corrected[..., 1:]
    low_p, high_p = equation.p_interval
    low_coefficient = 2.0 * equation.sigma * low_p
    high_coefficient = 2.0 * equation.sigma * high_p
    squared_radius = np.sum(location * location, axis=-1)
    parallel_coefficient = np.divide(
        np.sum(z * location, axis=-1),
        squared_radius,
        out=np.zeros_like(squared_radius),
        where=squared_radius > 0.0,
    )
    if transform == "segment":
        coefficient = np.clip(
            parallel_coefficient, low_coefficient, high_coefficient
        )
        coefficient = np.where(squared_radius > 0.0, coefficient, 0.0)
        corrected[..., 1:] = coefficient[..., None] * location
    elif transform == "span_only":
        corrected[..., 1:] = parallel_coefficient[..., None] * location
    elif transform == "box":
        first = low_coefficient * location
        second = high_coefficient * location
        corrected[..., 1:] = np.clip(z, np.minimum(first, second), np.maximum(first, second))
    elif transform == "ball":
        radius = high_coefficient * np.sqrt(squared_radius)
        norm = np.linalg.norm(z, axis=-1)
        scale = np.minimum(1.0, radius / np.maximum(norm, 1e-300))
        corrected[..., 1:] = z * scale[..., None]
    elif transform == "centre":
        coefficient = 0.5 * (low_coefficient + high_coefficient)
        corrected[..., 1:] = coefficient * location
    else:  # pragma: no cover - caller guard
        raise ValueError(transform)
    return corrected


def _lqg_energies(
    equation: Any, state: np.ndarray, x: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    location = np.asarray(x, dtype=np.float64)
    z = np.asarray(state, dtype=np.float64)[..., 1:]
    first, second = equation.dynamic_endpoints(location)
    low = np.minimum(first, second)
    high = np.maximum(first, second)
    box = np.sum(
        np.maximum(low - z, 0.0) ** 2 + np.maximum(z - high, 0.0) ** 2,
        axis=-1,
    )
    radius = 2.0 * equation.sigma * equation.p_interval[1] * np.linalg.norm(
        location, axis=-1
    )
    ball = np.maximum(np.linalg.norm(z, axis=-1) - radius, 0.0) ** 2
    return box, ball


class Round3MechanismMLP(Round2MechanismMLP):
    """Add C3/LQG geometry and optional LQG terminal antithetics."""

    def _terminal_estimate(
        self, n: int, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        if self.method.name != "raw_antithetic":
            return super()._terminal_estimate(n, t, x)

        b = len(t)
        d = self.equation.d
        samples = self.M**n
        remaining = np.maximum(self.equation.T - t, 0.0)
        gx = self.equation.terminal(x)
        self.stats.terminal_g_evals += b
        output = np.zeros((b, d + 1), dtype=np.float64)
        at_terminal = remaining <= 4.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            output[at_terminal] = self.equation.exact_state(
                t[at_terminal], x[at_terminal]
            )
        active = ~at_terminal
        if not np.any(active):
            return output
        xa = x[active]
        ha = remaining[active]
        ga = gx[active]
        ba = len(ha)
        pair_count = samples // 2
        base = self._normal((ba, pair_count, d))
        if samples % 2:
            extra = self._normal((ba, 1, d))
            normal = np.concatenate((base, -base, extra), axis=1)
        else:
            normal = np.concatenate((base, -base), axis=1)
        xt = (
            xa[:, None, :]
            + self.equation.mu * ha[:, None, None]
            + self.equation.sigma * np.sqrt(ha)[:, None, None] * normal
        )
        gt = self.equation.terminal(xt)
        self.stats.terminal_g_evals += ba * samples
        self.stats.terminal_samples += ba * samples
        difference = gt - ga[:, None]
        output[active, 0] = ga + np.mean(difference, axis=1)
        weights = normal / np.sqrt(ha)[:, None, None]
        output[active, 1:] = np.mean(difference[:, :, None] * weights, axis=1)
        return output

    def _round2_transform(
        self,
        state: np.ndarray,
        exact: np.ndarray,
        t: np.ndarray,
        x: np.ndarray,
    ) -> np.ndarray:
        transform = self.method.transform
        is_cubic_segment = (
            self.equation.family == "cubic_hj" and transform == "segment"
        )
        is_lqg = self.equation.family == "lqg" and transform in {
            "segment", "box", "ball", "span_only", "centre"
        }
        if not is_cubic_segment and not is_lqg:
            return super()._round2_transform(state, exact, t, x)

        candidate = self._candidate(state, exact)
        if is_lqg:
            box_energy, ball_energy = _lqg_energies(
                self.equation, candidate, x
            )
            corrected = _lqg_geometry(
                self.equation, candidate, x, transform
            )
        else:
            box_energy = _box_energy(self.equation, candidate)
            z = candidate[..., 1:]
            ball_energy = _u_violation_energy(
                self.equation, candidate[..., 0]
            ) + np.maximum(
                np.linalg.norm(z, axis=-1) - self.equation.z_ball_radius,
                0.0,
            ) ** 2
            corrected = _generic_segment(self.equation, candidate)

        count = candidate[..., 0].size
        self.stats.correction_opportunities += count
        self.stats.pre_box_violations += int(np.count_nonzero(box_energy > 0.0))
        self.stats.pre_ball_violations += int(np.count_nonzero(ball_energy > 0.0))
        self.stats.box_overshoot_energy_sum += float(np.sum(box_energy))
        method_energy = np.sum((corrected - candidate) ** 2, axis=-1)
        changed = np.any(corrected != candidate, axis=-1)
        self.stats.method_geometry_violations += int(
            np.count_nonzero(method_energy > 0.0)
        )
        self.stats.method_overshoot_energy_sum += float(np.sum(method_energy))
        self.stats.activated_states += int(np.count_nonzero(changed))
        return corrected


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
    study: str,
    code_commit: str | None,
    trace_draws: bool = False,
) -> dict[str, Any]:
    """Run one paired float64 repetition with the registered seed tree."""

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
        solver = Round3MechanismMLP(
            equation,
            M,
            method,
            np.random.default_rng(tree_seed),
            dose_rng=np.random.default_rng(dose_seed),
            time_beta_alpha=0.5,
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
        "round": 3,
        "study_origin": study,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "code_commit_at_execution": code_commit,
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
        "time_beta_alpha": 0.5,
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
            if fingerprints else None
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

