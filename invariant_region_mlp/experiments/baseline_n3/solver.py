"""Correction-level baseline as a strict subclass of the frozen solver.

The only estimator change is in the level-1 gradient correction.  The value
correction and all deeper gradient corrections are inherited mathematically
unchanged.  The implementation is kept here so the frozen double-estimator
module remains byte-for-byte untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import subprocess
import time
from pathlib import Path
from typing import Any, Literal

import numpy as np

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    fixed_points,
    make_equation,
)
from invariant_region_mlp.experiments.double_estimator.solver import (
    DoubleEstimatorMLP,
    _metrics,
)
from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import (
    ExtraDiagnostics,
    MethodSpec,
    WorkDiagnostics,
)


IMPLEMENTATION_REVISION = "correction_baseline_v1"
BaselineKind = Literal["none", "terminal", "oracle"]


@dataclass(frozen=True)
class BaselineMethod:
    name: str
    pathwise_terminal: bool = False
    double_generator: bool = False
    baseline_kind: BaselineKind = "none"
    transform: str = "raw"

    @property
    def method_spec(self) -> MethodSpec:
        return MethodSpec(self.name, self.transform, 1.0)


METHODS: dict[str, BaselineMethod] = {
    "raw": BaselineMethod("raw"),
    "path": BaselineMethod("path", pathwise_terminal=True),
    "double": BaselineMethod("double", double_generator=True),
    "double_path": BaselineMethod(
        "double_path", pathwise_terminal=True, double_generator=True
    ),
    "path_base": BaselineMethod(
        "path_base", pathwise_terminal=True, baseline_kind="terminal"
    ),
    "double_path_base": BaselineMethod(
        "double_path_base",
        pathwise_terminal=True,
        double_generator=True,
        baseline_kind="terminal",
    ),
    "path_base_oracle": BaselineMethod(
        "path_base_oracle", pathwise_terminal=True, baseline_kind="oracle"
    ),
}


class DrawTape:
    """Ordered NumPy draw tape shared by primary and auxiliary recursion.

    In record mode the NumPy solver generates its normal/uniform arrays and
    appends copies.  In replay mode the same arrays are returned without
    consuming an RNG.  The GPU solver consumes the public ``take`` interface.
    """

    def __init__(
        self,
        mode: Literal["record", "replay"] = "record",
        entries: list[tuple[str, np.ndarray]] | None = None,
    ) -> None:
        self.mode = mode
        self.entries: list[tuple[str, np.ndarray]] = [] if entries is None else entries
        self.cursor = 0

    def put(self, kind: str, values: np.ndarray) -> None:
        if self.mode != "record":
            raise RuntimeError("cannot append to a replay tape")
        self.entries.append((kind, np.asarray(values, dtype=np.float64).copy()))

    def take(self, kind: str, shape: tuple[int, ...]) -> np.ndarray:
        if self.mode != "replay":
            raise RuntimeError("cannot consume a record tape")
        if self.cursor >= len(self.entries):
            raise IndexError(f"draw tape exhausted while requesting {kind} {shape}")
        observed_kind, values = self.entries[self.cursor]
        self.cursor += 1
        if observed_kind != kind or values.shape != shape:
            raise ValueError(
                f"draw tape mismatch at {self.cursor - 1}: "
                f"expected {(kind, shape)}, observed {(observed_kind, values.shape)}"
            )
        return values.copy()

    def replay(self) -> "DrawTape":
        return DrawTape("replay", self.entries)

    @property
    def exhausted(self) -> bool:
        return self.cursor == len(self.entries)


class CorrectionBaselineMLP(DoubleEstimatorMLP):
    """Frozen double-estimator recursion plus a level-1 control variate."""

    def __init__(
        self,
        equation: Any,
        M: int,
        study_method: BaselineMethod,
        rng: np.random.Generator,
        *,
        aux_seed_sequence: np.random.SeedSequence,
        dose_rng: np.random.Generator,
        time_beta_alpha: float = 0.5,
        trace_draws: bool = False,
        draw_tape: DrawTape | None = None,
    ) -> None:
        super().__init__(
            equation,
            M,
            study_method,
            rng,
            aux_seed_sequence=aux_seed_sequence,
            dose_rng=dose_rng,
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        self.baseline_method = study_method
        self.draw_tape = draw_tape
        self.root_level1_integrand_variance: np.ndarray | None = None
        self.root_level1_gradient_increment: np.ndarray | None = None
        self.max_abs_level1_gradient_integrand = 0.0

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        if not self.baseline_method.pathwise_terminal:
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
        mean_z = np.mean(terminal_z, axis=1)
        # A mean of bitwise-identical FP64 values can depend on the sample
        # count through summation roundoff.  Preserve the exact constant in
        # this special case; this makes the preregistered linear-terminal
        # identity exact without affecting nonconstant terminal gradients.
        constant = np.all(terminal_z == terminal_z[:, :1, :], axis=(1, 2))
        if np.any(constant):
            mean_z[constant] = terminal_z[constant, 0, :]
        output[active, 1:] = mean_z
        return output

    def _normal(self, shape: tuple[int, ...]) -> np.ndarray:
        if self.draw_tape is None:
            return super()._normal(shape)
        if self.draw_tape.mode == "record":
            values = super()._normal(shape)
            self.draw_tape.put("normal", values)
            return values
        values = self.draw_tape.take("normal", shape)
        self.stats.standard_normal_variates += int(np.prod(shape))
        self._record_draw(values)
        return values

    def _power_time(self, shape: tuple[int, ...]) -> np.ndarray:
        if self.draw_tape is None:
            return super()._power_time(shape)
        if self.draw_tape.mode == "record":
            uniform = self.rng.random(shape, dtype=np.float64)
            self.draw_tape.put("uniform", uniform)
        else:
            uniform = self.draw_tape.take("uniform", shape)
        self.stats.time_uniform_variates += int(np.prod(shape))
        self._record_draw(uniform)
        return np.maximum(uniform, np.finfo(np.float64).tiny) ** (
            1.0 / self.time_beta_alpha
        )

    def _independent_solve(
        self, n: int, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        """Frozen auxiliary seeding, with the new subclass propagated."""

        bundle = self._aux_seed_sequence.spawn(1)[0]
        primary_seed, auxiliary_seed, dose_seed = bundle.spawn(3)
        child = type(self)(
            self.equation,
            self.M,
            self.baseline_method,
            np.random.default_rng(primary_seed),
            aux_seed_sequence=auxiliary_seed,
            dose_rng=np.random.default_rng(dose_seed),
            time_beta_alpha=self.time_beta_alpha,
            trace_draws=self._trace_aux_draws,
            draw_tape=self.draw_tape,
        )
        result = child.solve(n, t, x)
        self.stats.merge(child.stats)
        self.extra.merge(child.extra)
        self.recursive_invocations += child.recursive_invocations
        self.auxiliary_recursive_invocations += child.recursive_invocations
        self.max_abs_level1_gradient_integrand = max(
            self.max_abs_level1_gradient_integrand,
            child.max_abs_level1_gradient_integrand,
        )
        if self._aux_draw_hash is not None and child.full_draw_fingerprint is not None:
            self._aux_draw_hash.update(child.full_draw_fingerprint.encode("ascii"))
        return result

    def _baseline_value(
        self, terminal: np.ndarray, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray | None:
        kind = self.baseline_method.baseline_kind
        if kind == "none":
            return None
        state = terminal if kind == "terminal" else self.equation.exact_state(t, x)
        # _generator supplies the required plain-f accounting and diagnostics.
        return self._generator(state[:, None, :], t[:, None], x[:, None, :])[:, 0]

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
            self.root_level1_integrand_variance = np.full(b, np.nan, dtype=np.float64)
            self.root_level1_gradient_increment = np.full(
                (b, self.equation.d), np.nan, dtype=np.float64
            )
        if self.method.transform == "f_zero":
            return output

        baseline = self._baseline_value(output, t, x) if n >= 2 else None
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
            if self.baseline_method.double_generator:
                high_second = self._independent_solve(
                    level, flat_t, flat_x
                ).reshape(b, siblings, self.equation.d + 1)
                difference = self._generator_double(high_first, high_second, tr, xr)
            else:
                difference = self._generator(high_first, tr, xr)

            if level > 1:
                low_first = self.solve(level - 1, flat_t, flat_x).reshape(
                    b, siblings, self.equation.d + 1
                )
                if self.baseline_method.double_generator:
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
            gradient_difference = difference
            if level == 1 and baseline is not None:
                gradient_difference = difference - baseline[:, None]
            ebl = normal / np.sqrt(
                np.maximum(dt, np.finfo(np.float64).tiny)
            )[:, :, None]
            gradient_integrand = (
                importance[:, :, None]
                * gradient_difference[:, :, None]
                * ebl
            )
            gradient_increment = np.mean(gradient_integrand, axis=1)
            if level == 1 and gradient_integrand.size:
                self.max_abs_level1_gradient_integrand = max(
                    self.max_abs_level1_gradient_integrand,
                    float(np.max(np.abs(gradient_integrand))),
                )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
            if collect_root and self.root_level_u_corrections is not None:
                self.root_level_u_corrections[level - 1] = value_increment
            if collect_root and level == 1:
                ddof = 1 if siblings > 1 else 0
                assert self.root_level1_integrand_variance is not None
                assert self.root_level1_gradient_increment is not None
                self.root_level1_integrand_variance[:] = np.mean(
                    np.var(gradient_integrand, axis=1, ddof=ddof), axis=1
                )
                self.root_level1_gradient_increment[:] = gradient_increment
        return output


def _git_commit() -> str:
    root = Path(__file__).resolve().parents[3]
    return subprocess.run(
        [
            "git",
            "log",
            "-1",
            "--format=%H",
            "--",
            "invariant_region_mlp/experiments/baseline_n3/solver.py",
        ],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def run_task(task: dict[str, Any], *, trace_draws: bool = False) -> dict[str, Any]:
    """Run one complete 1,200-point Study-F row."""

    pde_id = str(task["pde_id"])
    d = int(task["d"])
    n = int(task["n"])
    M = int(task["M"])
    repetition = int(task["repetition"])
    chunk_size = int(task["chunk_size"])
    method_name = str(task["method"])
    method = METHODS[method_name]
    equation = make_equation(pde_id, d)
    t, x, is_validation = fixed_points(pde_id, d)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    root_variance = np.full(len(t), np.nan, dtype=np.float64)
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
        solver = CorrectionBaselineMLP(
            equation,
            M,
            method,
            np.random.default_rng(np.random.SeedSequence(seed_parts)),
            aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
            dose_rng=np.random.default_rng(
                np.random.SeedSequence(seed_parts + [0xD05E])
            ),
            trace_draws=trace_draws,
        )
        prediction[begin:end] = solver.solve(
            n, t[begin:end], x[begin:end], collect_root=True
        )
        if solver.root_level1_integrand_variance is not None:
            root_variance[begin:end] = solver.root_level1_integrand_variance
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
    finite_root_variance = root_variance[(~is_validation) & np.isfinite(root_variance)]
    point_hash = hashlib.sha256()
    point_hash.update(np.ascontiguousarray(t).tobytes())
    point_hash.update(np.ascontiguousarray(x).tobytes())
    point_hash.update(np.ascontiguousarray(is_validation).tobytes())
    gaussian_vectors = int(work["terminal_samples"] + work["transition_samples"])
    return {
        "schema_version": 1,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "code_commit": _git_commit(),
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
        "root_gradient_skill": test_metrics["gradient_relative_l2"],
        "validation_skill": validation_metrics["skill"],
        "all_skill": all_metrics["skill"],
        "mean_generator_bias": generator.get("bias"),
        "generator_rmse": math.sqrt(float(generator["mse"])) if generator else None,
        "generator_calls": int(work["f_evals"]),
        "recursive_calls": int(recursive_invocations),
        "auxiliary_recursive_calls": int(auxiliary_invocations),
        "recursively_evaluated_states": int(work["recursively_evaluated_states"]),
        "terminal_samples": int(work["terminal_samples"]),
        "transition_samples": int(work["transition_samples"]),
        "gaussian_vectors_drawn": gaussian_vectors,
        "gaussian_coordinates_drawn": gaussian_vectors * d,
        "root_level1_gradient_integrand_variance": (
            float(np.mean(finite_root_variance)) if len(finite_root_variance) else None
        ),
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


__all__ = [
    "BaselineMethod",
    "CorrectionBaselineMLP",
    "DrawTape",
    "IMPLEMENTATION_REVISION",
    "METHODS",
    "run_task",
]
