"""Double-estimator recursion implemented as a new MechanismMLP subclass."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import subprocess
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import (
    ExtraDiagnostics,
    MechanismMLP,
    MethodSpec,
    WorkDiagnostics,
    _box_energy,
    _clip_u,
    _u_violation_energy,
)

from .equations import BASE_SEED, equation_kind, fixed_points, make_equation


IMPLEMENTATION_REVISION = "double_generator_v1"


@dataclass(frozen=True)
class StudyMethod:
    name: str
    transform: str = "raw"
    pathwise_terminal: bool = False
    double_generator: bool = False

    @property
    def method_spec(self) -> MethodSpec:
        return MethodSpec(self.name, self.transform, 1.0)


METHODS: dict[str, StudyMethod] = {
    "raw": StudyMethod("raw"),
    "path": StudyMethod("path", pathwise_terminal=True),
    "double": StudyMethod("double", double_generator=True),
    "double_path": StudyMethod(
        "double_path", pathwise_terminal=True, double_generator=True
    ),
    "oracle_state": StudyMethod("oracle_state", "oracle_state"),
    "centre": StudyMethod("centre", "centre"),
    "f_zero": StudyMethod("f_zero", "f_zero"),
    "box": StudyMethod("box", "box"),
    "sub_box": StudyMethod("sub_box", "sub_box"),
}


def methods_for_pde(pde_id: str) -> tuple[str, ...]:
    certificate = "sub_box" if pde_id == "MR" else "box"
    return (
        "raw",
        "path",
        "double",
        "double_path",
        "oracle_state",
        "centre",
        "f_zero",
        certificate,
    )


class DoubleEstimatorMLP(MechanismMLP):
    """Full-history MLP with optional pathwise terminal and double driver.

    The primary RNG is never consumed by an auxiliary recursive estimate.
    Consequently, its entire draw stream is identical to the corresponding
    raw/path run.  Auxiliary recursive trees use descendants of a separate
    SeedSequence and are independent conditional on their shared child point.
    """

    def __init__(
        self,
        equation: Any,
        M: int,
        study_method: StudyMethod,
        rng: np.random.Generator,
        *,
        aux_seed_sequence: np.random.SeedSequence,
        dose_rng: np.random.Generator,
        time_beta_alpha: float = 0.5,
        trace_draws: bool = False,
    ) -> None:
        super().__init__(
            equation,
            M,
            study_method.method_spec,
            rng,
            dose_rng=dose_rng,
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        self.study_method = study_method
        self._aux_seed_sequence = aux_seed_sequence
        self._trace_aux_draws = trace_draws
        self._aux_draw_hash = hashlib.sha256() if trace_draws else None
        self.recursive_invocations = 0
        self.auxiliary_recursive_invocations = 0

    @property
    def auxiliary_draw_fingerprint(self) -> str | None:
        return None if self._aux_draw_hash is None else self._aux_draw_hash.hexdigest()

    @property
    def full_draw_fingerprint(self) -> str | None:
        if self.draw_fingerprint is None:
            return None
        payload = self.draw_fingerprint + (self.auxiliary_draw_fingerprint or "")
        return hashlib.sha256(payload.encode("ascii")).hexdigest()

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        if not self.study_method.pathwise_terminal:
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
        normal = self._normal((len(xa), samples, d))
        terminal_points = (
            xa[:, None, :]
            + self.equation.mu * ha[:, None, None]
            + self.equation.sigma * np.sqrt(ha)[:, None, None] * normal
        )
        terminal_values = self.equation.terminal(terminal_points)
        self.stats.terminal_g_evals += len(xa) * samples
        self.stats.terminal_samples += len(xa) * samples
        output[active, 0] = ga + np.mean(terminal_values - ga[:, None], axis=1)
        terminal_times = np.full(terminal_points.shape[:-1], self.equation.T)
        terminal_z = self.equation.exact_z(terminal_times, terminal_points)
        output[active, 1:] = np.mean(terminal_z, axis=1)
        return output

    def _transform(self, state: np.ndarray, exact: np.ndarray) -> np.ndarray:
        if self.method.transform != "sub_box":
            return super()._transform(state, exact)
        if self.equation.family != "multiridge":
            raise ValueError("sub_box is only defined for the multi-ridge equation")

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
        self.stats.pre_box_violations += int(np.count_nonzero(box_energy > 0.0))
        self.stats.pre_ball_violations += int(np.count_nonzero(ball_energy > 0.0))
        self.stats.box_overshoot_energy_sum += float(np.sum(box_energy))

        corrected = candidate.copy()
        corrected[..., 0] = _clip_u(self.equation, corrected[..., 0])
        coordinates = corrected[..., 1:] @ self.equation.Qh
        coordinates = np.clip(
            coordinates, self.equation.sub_lo, self.equation.sub_hi
        )
        corrected[..., 1:] = coordinates @ self.equation.Qh.T
        method_energy = np.sum((corrected - candidate) ** 2, axis=-1)
        changed = np.any(corrected != candidate, axis=-1)
        self.stats.method_geometry_violations += int(
            np.count_nonzero(method_energy > 0.0)
        )
        self.stats.method_overshoot_energy_sum += float(np.sum(method_energy))
        self.stats.activated_states += int(np.count_nonzero(changed))
        return corrected

    def _double_value(self, first: np.ndarray, second: np.ndarray) -> np.ndarray:
        z1 = np.asarray(first[..., 1:], dtype=np.float64)
        z2 = np.asarray(second[..., 1:], dtype=np.float64)
        kind = equation_kind(self.equation)
        if kind == "quadratic":
            return -0.5 * np.sum(z1 * z2, axis=-1)
        if kind == "game":
            split = int(self.equation.dA)
            return (
                -0.5
                * float(self.equation.a)
                * np.sum(z1[..., :split] * z2[..., :split], axis=-1)
                + 0.5
                * float(self.equation.b)
                * np.sum(z1[..., split:] * z2[..., split:], axis=-1)
            )
        if kind == "norm":
            norm = np.linalg.norm(z1, axis=-1)
            direction = np.zeros_like(z1)
            nonzero = norm > 0.0
            direction[nonzero] = z1[nonzero] / norm[nonzero, None]
            # In this repository z=sigma*grad(u), hence the equation's
            # generator coefficient is lambda_f/sigma.
            coefficient = float(self.equation.lambda_f / self.equation.sigma)
            return -coefficient * np.sum(direction * z2, axis=-1)
        raise AssertionError(kind)

    def _record_generator(
        self, estimate: np.ndarray, truth: np.ndarray
    ) -> None:
        error = estimate - truth
        count = int(estimate.size)
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

    def _generator_double(
        self,
        first: np.ndarray,
        second: np.ndarray,
        t: np.ndarray,
        x: np.ndarray,
    ) -> np.ndarray:
        if first.shape != second.shape or first.shape[-1] != self.equation.d + 1:
            raise ValueError("double states have incompatible shapes")
        exact = self.equation.exact_state(t, x)
        average = 0.5 * (first + second)
        # The average is the state used wherever the ordinary recursion would
        # inspect U; only f itself receives the two separate estimates.
        self.extra.noise.add(
            average[..., 0] - exact[..., 0],
            average[..., 1:] - exact[..., 1:],
        )
        self._transform(average, exact)  # identity for registered double methods; diagnostics only
        estimate = self._double_value(first, second)
        truth = self.equation.generator(exact[..., 0], exact[..., 1:])
        self._record_generator(estimate, truth)
        return estimate

    def _independent_solve(
        self, n: int, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        bundle = self._aux_seed_sequence.spawn(1)[0]
        primary_seed, auxiliary_seed, dose_seed = bundle.spawn(3)
        child = DoubleEstimatorMLP(
            self.equation,
            self.M,
            self.study_method,
            np.random.default_rng(primary_seed),
            aux_seed_sequence=auxiliary_seed,
            dose_rng=np.random.default_rng(dose_seed),
            time_beta_alpha=self.time_beta_alpha,
            trace_draws=self._trace_aux_draws,
        )
        result = child.solve(n, t, x)
        self.stats.merge(child.stats)
        self.extra.merge(child.extra)
        self.recursive_invocations += child.recursive_invocations
        self.auxiliary_recursive_invocations += child.recursive_invocations
        if self._aux_draw_hash is not None and child.full_draw_fingerprint is not None:
            self._aux_draw_hash.update(child.full_draw_fingerprint.encode("ascii"))
        return result

    def solve(
        self,
        n: int,
        t: np.ndarray,
        x: np.ndarray,
        *,
        collect_root: bool = False,
    ) -> np.ndarray:
        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.equation.d)
        b = len(t)
        self.recursive_invocations += 1
        self.stats.recursively_evaluated_states += b
        if n <= 0:
            return np.zeros((b, self.equation.d + 1), dtype=np.float64)

        output = self._terminal_estimate(n, t, x)
        if collect_root:
            self.root_terminal = output.copy()
            self.root_level_u_corrections = np.zeros(
                (max(n - 1, 0), b), dtype=np.float64
            )
        if self.method.transform == "f_zero":
            return output

        remaining = self.equation.T - t
        alpha = self.time_beta_alpha
        for level in range(1, n):
            siblings = self.M ** (n - level)
            r = self._power_time((b, siblings))
            dt = remaining[:, None] * r
            normal = self._normal((b, siblings, self.equation.d))
            self.stats.transition_samples += b * siblings
            xr = (
                x[:, None, :]
                + self.equation.mu * dt[:, :, None]
                + self.equation.sigma * np.sqrt(dt)[:, :, None] * normal
            )
            tr = t[:, None] + dt
            flat_t = tr.reshape(-1)
            flat_x = xr.reshape(-1, self.equation.d)

            high_first = self.solve(level, flat_t, flat_x).reshape(
                b, siblings, self.equation.d + 1
            )
            if self.study_method.double_generator:
                high_second = self._independent_solve(
                    level, flat_t, flat_x
                ).reshape(b, siblings, self.equation.d + 1)
                difference = self._generator_double(
                    high_first, high_second, tr, xr
                )
            else:
                difference = self._generator(high_first, tr, xr)

            if level > 1:
                low_first = self.solve(level - 1, flat_t, flat_x).reshape(
                    b, siblings, self.equation.d + 1
                )
                if self.study_method.double_generator:
                    low_second = self._independent_solve(
                        level - 1, flat_t, flat_x
                    ).reshape(b, siblings, self.equation.d + 1)
                    difference = difference - self._generator_double(
                        low_first, low_second, tr, xr
                    )
                else:
                    difference = difference - self._generator(low_first, tr, xr)

            importance = remaining[:, None] * r ** (1.0 - alpha) / alpha
            value_increment = np.mean(importance * difference, axis=1)
            ebl = normal / np.sqrt(
                np.maximum(dt, np.finfo(np.float64).tiny)
            )[:, :, None]
            gradient_increment = np.mean(
                importance[:, :, None]
                * difference[:, :, None]
                * ebl,
                axis=1,
            )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
            if collect_root and self.root_level_u_corrections is not None:
                self.root_level_u_corrections[level - 1] = value_increment
        return output


def _metrics(
    prediction: np.ndarray, truth: np.ndarray, mask: np.ndarray
) -> dict[str, float | int]:
    pred = np.asarray(prediction, dtype=np.float64)[mask]
    ref = np.asarray(truth, dtype=np.float64)[mask]
    finite_rows = np.all(np.isfinite(pred), axis=1)
    nonfinite_count = int(pred.size - np.count_nonzero(np.isfinite(pred)))
    if len(pred) == 0 or not np.all(finite_rows):
        return {
            "n_points": int(len(pred)),
            "skill": float("inf"),
            "value_rmse": float("inf"),
            "value_relative_l2": float("inf"),
            "value_bias": float("nan"),
            "gradient_relative_l2": float("inf"),
            "nonfinite_state_count": nonfinite_count,
        }
    value_error = pred[:, 0] - ref[:, 0]
    reference_std = float(np.std(ref[:, 0], ddof=0))
    value_rmse = float(np.sqrt(np.mean(value_error**2)))
    return {
        "n_points": int(len(pred)),
        "skill": value_rmse / reference_std if reference_std > 0.0 else float("inf"),
        "value_rmse": value_rmse,
        "value_relative_l2": float(
            np.linalg.norm(value_error) / np.linalg.norm(ref[:, 0])
        ),
        "value_bias": float(np.mean(value_error)),
        "gradient_relative_l2": float(
            np.linalg.norm(pred[:, 1:] - ref[:, 1:])
            / np.linalg.norm(ref[:, 1:])
        ),
        "nonfinite_state_count": nonfinite_count,
    }


@lru_cache(maxsize=1)
def code_commit() -> str:
    root = Path(__file__).resolve().parents[3]
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def run_task(task: dict[str, Any], *, trace_draws: bool = False) -> dict[str, Any]:
    """Run one complete 1,200-point method/repetition task."""

    pde_id = str(task["pde_id"])
    d = int(task["d"])
    n = int(task["n"])
    M = int(task["M"])
    repetition = int(task["repetition"])
    chunk_size = int(task["chunk_size"])
    method_name = str(task["method"])
    study_method = METHODS[method_name]
    equation = make_equation(pde_id, d)
    t, x, is_validation = fixed_points(pde_id, d)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    diagnostics = WorkDiagnostics()
    extras = ExtraDiagnostics()
    recursive_invocations = 0
    auxiliary_invocations = 0
    primary_fingerprints: list[str] = []
    auxiliary_fingerprints: list[str] = []

    started = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(t), chunk_size)):
        end = min(begin + chunk_size, len(t))
        seed_parts = [BASE_SEED, d, n, M, repetition, chunk_index]
        primary_seed = np.random.SeedSequence(seed_parts)
        auxiliary_seed = np.random.SeedSequence(seed_parts + [0xD0B1E])
        dose_seed = np.random.SeedSequence(seed_parts + [0xD05E])
        solver = DoubleEstimatorMLP(
            equation,
            M,
            study_method,
            np.random.default_rng(primary_seed),
            aux_seed_sequence=auxiliary_seed,
            dose_rng=np.random.default_rng(dose_seed),
            trace_draws=trace_draws,
        )
        prediction[begin:end] = solver.solve(n, t[begin:end], x[begin:end])
        diagnostics.merge(solver.stats)
        extras.merge(solver.extra)
        recursive_invocations += solver.recursive_invocations
        auxiliary_invocations += solver.auxiliary_recursive_invocations
        if solver.draw_fingerprint is not None:
            primary_fingerprints.append(solver.draw_fingerprint)
        if solver.auxiliary_draw_fingerprint is not None:
            auxiliary_fingerprints.append(solver.auxiliary_draw_fingerprint)
    elapsed = time.perf_counter() - started

    work = diagnostics.summary()
    generator = work.get("generator") or {}
    test_metrics = _metrics(prediction, truth, ~is_validation)
    validation_metrics = _metrics(prediction, truth, is_validation)
    all_metrics = _metrics(prediction, truth, np.ones(len(t), dtype=bool))
    point_hash = hashlib.sha256()
    point_hash.update(np.ascontiguousarray(t).tobytes())
    point_hash.update(np.ascontiguousarray(x).tobytes())
    point_hash.update(np.ascontiguousarray(is_validation).tobytes())
    return {
        "schema_version": 1,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "code_commit": code_commit(),
        "pde_id": pde_id,
        "family": equation.family,
        "dimension": d,
        "n": n,
        "M": M,
        "method": method_name,
        "repetition": repetition,
        "base_seed": BASE_SEED,
        "seed_scheme": "SeedSequence([20261201,d,n,M,rep,chunk_index])",
        "aux_seed_scheme": "primary seed parts + [0xD0B1E], spawned recursively",
        "chunk_size": chunk_size,
        "point_count": len(t),
        "point_fingerprint": point_hash.hexdigest(),
        "dtype": "float64",
        "skill": test_metrics["skill"],
        "value_rmse": test_metrics["value_rmse"],
        "value_relative_l2": test_metrics["value_relative_l2"],
        "value_bias": test_metrics["value_bias"],
        "gradient_relative_l2": test_metrics["gradient_relative_l2"],
        "validation_skill": validation_metrics["skill"],
        "all_skill": all_metrics["skill"],
        "mean_generator_bias": generator.get("bias"),
        "generator_rmse": (
            math.sqrt(float(generator["mse"])) if generator else None
        ),
        "generator_calls": int(work["f_evals"]),
        "recursive_calls": int(recursive_invocations),
        "auxiliary_recursive_calls": int(auxiliary_invocations),
        "recursively_evaluated_states": int(work["recursively_evaluated_states"]),
        "terminal_samples": int(work["terminal_samples"]),
        "transition_samples": int(work["transition_samples"]),
        "wall_time_seconds": elapsed,
        "nonfinite_state_count": int(test_metrics["nonfinite_state_count"]),
        "nonfinite_generator_count": int(work["nonfinite_generator_values"]),
        "primary_draw_fingerprint": (
            hashlib.sha256("".join(primary_fingerprints).encode("ascii")).hexdigest()
            if primary_fingerprints
            else None
        ),
        "auxiliary_draw_fingerprint": (
            hashlib.sha256("".join(auxiliary_fingerprints).encode("ascii")).hexdigest()
            if auxiliary_fingerprints
            else None
        ),
        "equation": equation.provenance(),
        "extra_diagnostics": extras.summary(),
    }

