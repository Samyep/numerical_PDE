"""Corrected full-history MLP variants for the active VB benchmark.

All methods traverse the same stochastic tree for a fixed seed.  Constraint
corrections are applied only to completed child states, immediately before
the child is reused by the nonlinear generator.  The returned root state is
never clipped.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
import platform
import time
from typing import Any, Iterable

import numpy as np

from vb_equation import ViscousBurgersEquation


@dataclass(frozen=True)
class MethodSpec:
    name: str
    transform: str
    factor: float = 1.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


RAW = MethodSpec("raw", "raw")
SAMPLE_BOX = MethodSpec("sample_box", "sample_box")
Z_ONLY_BOX = MethodSpec("z_only_box", "z_only_box")
SAMPLE_BALL = MethodSpec("sample_ball", "sample_ball")
BATCH_BOX = MethodSpec("batch_box", "batch_box")
Z_ZERO = MethodSpec("z_zero", "z_zero", 0.0)
F_ZERO = MethodSpec("f_zero", "f_zero", 0.0)


def shrink_method(c: float) -> MethodSpec:
    return MethodSpec(f"shrink_c{c:g}", "shrink", float(c))


def box_factor_method(factor: float) -> MethodSpec:
    return MethodSpec(f"box_factor_a{factor:g}", "box_factor", float(factor))


def core_methods() -> list[MethodSpec]:
    return [RAW, SAMPLE_BOX, Z_ONLY_BOX, SAMPLE_BALL, BATCH_BOX, Z_ZERO, F_ZERO]


def tuning_methods() -> list[MethodSpec]:
    # c=0 is exactly z_zero and c=1 is exactly raw.  a=0 is also z_zero for
    # this generator, while a=1 is sample_box.  Avoid recomputing aliases.
    return [
        shrink_method(0.1),
        shrink_method(0.25),
        shrink_method(0.5),
        shrink_method(0.75),
        box_factor_method(0.25),
        box_factor_method(0.5),
        box_factor_method(0.75),
        box_factor_method(1.25),
        box_factor_method(1.5),
        box_factor_method(2.0),
    ]


def pilot_methods() -> list[MethodSpec]:
    return core_methods() + tuning_methods()


def parse_method(name: str) -> MethodSpec:
    fixed = {method.name: method for method in core_methods()}
    if name in fixed:
        return fixed[name]
    if name.startswith("shrink_c"):
        return shrink_method(float(name.removeprefix("shrink_c")))
    if name.startswith("box_factor_a"):
        return box_factor_method(float(name.removeprefix("box_factor_a")))
    raise ValueError(f"unknown method {name!r}")


def terminal_ebl_weight(standard_normal: np.ndarray, remaining_time: np.ndarray) -> np.ndarray:
    """Correct EBL weight for z=sigma*grad(u): xi/sqrt(T-t)."""

    remaining_time = np.asarray(remaining_time, dtype=np.float64)
    return np.asarray(standard_normal, dtype=np.float64) / np.sqrt(remaining_time)[..., None]


@dataclass
class WorkDiagnostics:
    terminal_g_evals: int = 0
    terminal_samples: int = 0
    f_evals: int = 0
    recursively_evaluated_states: int = 0
    transition_samples: int = 0
    standard_normal_variates: int = 0
    time_uniform_variates: int = 0
    correction_opportunities: int = 0
    pre_box_violations: int = 0
    pre_ball_violations: int = 0
    method_geometry_violations: int = 0
    activated_states: int = 0
    box_overshoot_energy_sum: float = 0.0
    method_overshoot_energy_sum: float = 0.0
    generator_count: int = 0
    generator_error_sum: float = 0.0
    generator_abs_error_sum: float = 0.0
    generator_squared_error_sum: float = 0.0
    generator_estimate_sum: float = 0.0
    generator_truth_sum: float = 0.0
    nonfinite_generator_values: int = 0
    batch_alpha_count: int = 0
    batch_alpha_sum: float = 0.0
    batch_alpha_active: int = 0
    batch_alpha_hist: np.ndarray = field(
        default_factory=lambda: np.zeros(1000, dtype=np.int64)
    )

    def merge(self, other: "WorkDiagnostics") -> None:
        for name in (
            "terminal_g_evals",
            "terminal_samples",
            "f_evals",
            "recursively_evaluated_states",
            "transition_samples",
            "standard_normal_variates",
            "time_uniform_variates",
            "correction_opportunities",
            "pre_box_violations",
            "pre_ball_violations",
            "method_geometry_violations",
            "activated_states",
            "generator_count",
            "nonfinite_generator_values",
            "batch_alpha_count",
            "batch_alpha_active",
        ):
            setattr(self, name, getattr(self, name) + getattr(other, name))
        for name in (
            "box_overshoot_energy_sum",
            "method_overshoot_energy_sum",
            "generator_error_sum",
            "generator_abs_error_sum",
            "generator_squared_error_sum",
            "generator_estimate_sum",
            "generator_truth_sum",
            "batch_alpha_sum",
        ):
            setattr(self, name, float(getattr(self, name) + getattr(other, name)))
        self.batch_alpha_hist += other.batch_alpha_hist

    def _alpha_quantile(self, probability: float) -> float | None:
        if self.batch_alpha_count == 0:
            return None
        target = probability * max(self.batch_alpha_count - 1, 0)
        index = int(np.searchsorted(np.cumsum(self.batch_alpha_hist), target, side="right"))
        index = min(index, len(self.batch_alpha_hist) - 1)
        return float((index + 0.5) / len(self.batch_alpha_hist))

    def summary(self) -> dict[str, Any]:
        opportunities = max(self.correction_opportunities, 1)
        gc = max(self.generator_count, 1)
        payload: dict[str, Any] = {
            "terminal_g_evals": self.terminal_g_evals,
            "terminal_samples": self.terminal_samples,
            "f_evals": self.f_evals,
            "recursively_evaluated_states": self.recursively_evaluated_states,
            "transition_samples": self.transition_samples,
            "total_stochastic_samples": self.terminal_samples + self.transition_samples,
            "standard_normal_variates": self.standard_normal_variates,
            "time_uniform_variates": self.time_uniform_variates,
            "correction_opportunities": self.correction_opportunities,
            "pre_box_violation_rate": self.pre_box_violations / opportunities,
            "pre_ball_violation_rate": self.pre_ball_violations / opportunities,
            "method_geometry_violation_rate": self.method_geometry_violations / opportunities,
            "projection_activation_rate": self.activated_states / opportunities,
            "mean_box_overshoot_energy": self.box_overshoot_energy_sum / opportunities,
            "mean_method_overshoot_energy": self.method_overshoot_energy_sum / opportunities,
            "nonfinite_generator_values": self.nonfinite_generator_values,
        }
        if self.generator_count:
            bias = self.generator_error_sum / gc
            payload["generator"] = {
                "count": self.generator_count,
                "mse": self.generator_squared_error_sum / gc,
                "bias": bias,
                "absolute_bias": abs(bias),
                "mae": self.generator_abs_error_sum / gc,
                "mean_estimate": self.generator_estimate_sum / gc,
                "mean_truth": self.generator_truth_sum / gc,
            }
        else:
            payload["generator"] = None
        if self.batch_alpha_count:
            payload["batch_alpha"] = {
                "count": self.batch_alpha_count,
                "mean": self.batch_alpha_sum / self.batch_alpha_count,
                "median": self._alpha_quantile(0.5),
                "p10": self._alpha_quantile(0.1),
                "p90": self._alpha_quantile(0.9),
                "activation_rate": self.batch_alpha_active / self.batch_alpha_count,
                "histogram_bins": len(self.batch_alpha_hist),
            }
        else:
            payload["batch_alpha"] = None
        return payload


class FullHistoryMLP:
    """Vectorized, float64, full-history MLP with corrected EBL scaling."""

    def __init__(
        self,
        equation: ViscousBurgersEquation,
        M: int,
        method: MethodSpec,
        rng: np.random.Generator,
        *,
        time_beta_alpha: float = 0.5,
        trace_draws: bool = False,
    ) -> None:
        if M < 1:
            raise ValueError("M must be at least one")
        if not 0.0 < time_beta_alpha <= 1.0:
            raise ValueError("time_beta_alpha must lie in (0,1]")
        self.equation = equation
        self.M = int(M)
        self.method = method
        self.rng = rng
        self.time_beta_alpha = float(time_beta_alpha)
        self.stats = WorkDiagnostics()
        self._draw_hash = hashlib.sha256() if trace_draws else None
        self.root_terminal: np.ndarray | None = None
        self.root_level_u_corrections: np.ndarray | None = None

    @property
    def draw_fingerprint(self) -> str | None:
        return None if self._draw_hash is None else self._draw_hash.hexdigest()

    def _record_draw(self, values: np.ndarray) -> None:
        if self._draw_hash is not None:
            array = np.ascontiguousarray(values, dtype=np.float64)
            self._draw_hash.update(np.asarray(array.shape, dtype=np.int64).tobytes())
            self._draw_hash.update(array.tobytes())

    def _normal(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.rng.standard_normal(shape, dtype=np.float64)
        self.stats.standard_normal_variates += int(np.prod(shape))
        self._record_draw(values)
        return values

    def _power_time(self, shape: tuple[int, ...]) -> np.ndarray:
        uniform = self.rng.random(shape, dtype=np.float64)
        self.stats.time_uniform_variates += int(np.prod(shape))
        self._record_draw(uniform)
        # Inverse CDF of Beta(alpha,1), with a positive floor for EBL weights.
        return np.maximum(uniform, np.finfo(np.float64).tiny) ** (1.0 / self.time_beta_alpha)

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        b = len(t)
        d = self.equation.d
        samples = self.M**n
        remaining = np.maximum(self.equation.T - t, 0.0)
        gx = self.equation.terminal(x)
        self.stats.terminal_g_evals += b
        output = np.zeros((b, d + 1), dtype=np.float64)
        at_terminal = remaining <= 4.0 * np.finfo(np.float64).eps
        if np.any(at_terminal):
            output[at_terminal] = self.equation.exact_state(t[at_terminal], x[at_terminal])
        active = ~at_terminal
        if not np.any(active):
            return output
        ta = t[active]
        xa = x[active]
        ha = remaining[active]
        ga = gx[active]
        ba = len(ta)
        normal = self._normal((ba, samples, d))
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
        weights = terminal_ebl_weight(normal, ha[:, None])
        output[active, 1:] = np.mean(difference[:, :, None] * weights, axis=1)
        return output

    def _record_alpha(self, alpha: np.ndarray) -> None:
        alpha = np.asarray(alpha, dtype=np.float64).reshape(-1)
        self.stats.batch_alpha_count += len(alpha)
        self.stats.batch_alpha_sum += float(np.sum(alpha))
        self.stats.batch_alpha_active += int(np.count_nonzero(alpha < 1.0 - 1e-14))
        hist, _ = np.histogram(alpha, bins=1000, range=(0.0, 1.0))
        self.stats.batch_alpha_hist += hist

    def _correct(self, state: np.ndarray) -> np.ndarray:
        """Correct a B x K sibling array and record pre-correction diagnostics."""

        state = np.asarray(state, dtype=np.float64)
        u = state[..., 0]
        z = state[..., 1:]
        nstates = u.size
        eq = self.equation
        u_lo = np.maximum(-u, 0.0)
        u_hi = np.maximum(u - 1.0, 0.0)
        z_lo = np.maximum(-z, 0.0)
        z_hi = np.maximum(z - eq.z_upper, 0.0)
        box_energy_by_state = u_lo**2 + u_hi**2 + np.sum(z_lo**2 + z_hi**2, axis=-1)
        z_norm = np.linalg.norm(z, axis=-1)
        ball_energy_by_state = u_lo**2 + u_hi**2 + np.maximum(z_norm - eq.z_ball_radius, 0.0) ** 2
        self.stats.correction_opportunities += nstates
        self.stats.pre_box_violations += int(np.count_nonzero(box_energy_by_state > 0.0))
        self.stats.pre_ball_violations += int(np.count_nonzero(ball_energy_by_state > 0.0))
        self.stats.box_overshoot_energy_sum += float(np.sum(box_energy_by_state))

        transform = self.method.transform
        corrected = state.copy()
        if transform == "raw":
            method_energy = box_energy_by_state
        elif transform == "sample_box":
            method_energy = box_energy_by_state
            corrected[..., 0] = np.clip(u, 0.0, 1.0)
            corrected[..., 1:] = np.clip(z, 0.0, eq.z_upper)
        elif transform == "z_only_box":
            method_energy = np.sum(z_lo**2 + z_hi**2, axis=-1)
            corrected[..., 1:] = np.clip(z, 0.0, eq.z_upper)
        elif transform == "box_factor":
            upper = self.method.factor * eq.z_upper
            factor_hi = np.maximum(z - upper, 0.0)
            method_energy = u_lo**2 + u_hi**2 + np.sum(z_lo**2 + factor_hi**2, axis=-1)
            corrected[..., 0] = np.clip(u, 0.0, 1.0)
            corrected[..., 1:] = np.clip(z, 0.0, upper)
        elif transform == "sample_ball":
            method_energy = ball_energy_by_state
            corrected[..., 0] = np.clip(u, 0.0, 1.0)
            scale = np.minimum(1.0, eq.z_ball_radius / np.maximum(z_norm, 1e-300))
            corrected[..., 1:] = z * scale[..., None]
        elif transform == "batch_box":
            # Sign feasibility is enforced coordinatewise.  A single common
            # alpha per parent then contracts every nonnegative sibling so the
            # maximum coordinate is at most sigma/4.  This is identity on an
            # already feasible sibling batch and is not misrepresented as the
            # coordinatewise Euclidean box projection.
            method_energy = box_energy_by_state
            corrected[..., 0] = np.clip(u, 0.0, 1.0)
            nonnegative = np.maximum(z, 0.0)
            max_coordinate = np.max(nonnegative, axis=(1, 2))
            alpha = np.minimum(1.0, eq.z_upper / np.maximum(max_coordinate, 1e-300))
            corrected[..., 1:] = nonnegative * alpha[:, None, None]
            self._record_alpha(alpha)
        elif transform == "z_zero":
            method_energy = np.sum(z**2, axis=-1)
            corrected[..., 1:] = 0.0
        elif transform == "shrink":
            method_energy = np.sum(((1.0 - self.method.factor) * z) ** 2, axis=-1)
            corrected[..., 1:] = self.method.factor * z
        else:
            raise ValueError(f"correction {transform!r} cannot be used before f")

        self.stats.method_geometry_violations += int(np.count_nonzero(method_energy > 0.0))
        self.stats.method_overshoot_energy_sum += float(np.sum(method_energy))
        changed = np.any(corrected != state, axis=-1)
        self.stats.activated_states += int(np.count_nonzero(changed))
        return corrected

    def _generator(self, state: np.ndarray, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        b, k, width = state.shape
        if width != self.equation.d + 1:
            raise ValueError("invalid state width")
        corrected = self._correct(state)
        estimate = self.equation.generator(corrected[..., 0], corrected[..., 1:])
        truth = self.equation.exact_generator(t, x)
        error = estimate - truth
        count = b * k
        finite = np.isfinite(estimate) & np.isfinite(error)
        self.stats.f_evals += count
        self.stats.generator_count += count
        self.stats.nonfinite_generator_values += int(count - np.count_nonzero(finite))
        # Preserve infinities in the aggregate diagnostics; they are scientific
        # failures, not values to hide with nan_to_num or clipping.
        self.stats.generator_error_sum += float(np.sum(error))
        self.stats.generator_abs_error_sum += float(np.sum(np.abs(error)))
        self.stats.generator_squared_error_sum += float(np.sum(error**2))
        self.stats.generator_estimate_sum += float(np.sum(estimate))
        self.stats.generator_truth_sum += float(np.sum(truth))
        return estimate

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
        self.stats.recursively_evaluated_states += b
        if n <= 0:
            return np.zeros((b, self.equation.d + 1), dtype=np.float64)

        output = self._terminal_estimate(n, t, x)
        if collect_root:
            self.root_terminal = output.copy()
            self.root_level_u_corrections = np.zeros((max(n - 1, 0), b), dtype=np.float64)
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
            brownian = np.sqrt(dt)[:, :, None] * normal
            xr = (
                x[:, None, :]
                + self.equation.mu * dt[:, :, None]
                + self.equation.sigma * brownian
            )
            tr = t[:, None] + dt
            high = self.solve(level, tr.reshape(-1), xr.reshape(-1, self.equation.d))
            high = high.reshape(b, siblings, self.equation.d + 1)
            difference = self._generator(high, tr, xr)
            if level > 1:
                low = self.solve(level - 1, tr.reshape(-1), xr.reshape(-1, self.equation.d))
                low = low.reshape(b, siblings, self.equation.d + 1)
                difference = difference - self._generator(low, tr, xr)

            importance = remaining[:, None] * r ** (1.0 - alpha) / alpha
            value_increment = np.mean(importance * difference, axis=1)
            ebl = normal / np.sqrt(np.maximum(dt, np.finfo(np.float64).tiny))[:, :, None]
            gradient_increment = np.mean(
                importance[:, :, None] * difference[:, :, None] * ebl,
                axis=1,
            )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
            if collect_root and self.root_level_u_corrections is not None:
                self.root_level_u_corrections[level - 1] = value_increment
        return output


def _error_metrics(
    prediction: np.ndarray,
    truth: np.ndarray,
    mask: np.ndarray,
) -> dict[str, float]:
    pred = prediction[mask]
    ref = truth[mask]
    if len(pred) == 0:
        return {
            "n_points": 0,
            "value_relative_l2": float("nan"),
            "value_mae": float("nan"),
            "value_bias": float("nan"),
            "gradient_relative_l2": float("nan"),
            "gradient_mae": float("nan"),
            "full_state_relative_l2": float("nan"),
            "nonfinite_state_fraction": float("nan"),
        }
    value_error = pred[:, 0] - ref[:, 0]
    gradient_error = pred[:, 1:] - ref[:, 1:]
    full_error = pred - ref
    return {
        "n_points": int(np.count_nonzero(mask)),
        "value_relative_l2": float(np.linalg.norm(value_error) / np.linalg.norm(ref[:, 0])),
        "value_mae": float(np.mean(np.abs(value_error))),
        "value_bias": float(np.mean(value_error)),
        "gradient_relative_l2": float(np.linalg.norm(gradient_error) / np.linalg.norm(ref[:, 1:])),
        "gradient_mae": float(np.mean(np.abs(gradient_error))),
        "full_state_relative_l2": float(np.linalg.norm(full_error) / np.linalg.norm(ref)),
        "nonfinite_state_fraction": float(1.0 - np.mean(np.isfinite(pred))),
    }


def run_single_repetition(
    *,
    equation: ViscousBurgersEquation,
    method: MethodSpec,
    n: int,
    M: int,
    repetition: int,
    t: np.ndarray,
    x: np.ndarray,
    is_validation: np.ndarray,
    base_seed: int = 20261006,
    chunk_size: int = 16,
    time_beta_alpha: float = 0.5,
    trace_draws: bool = False,
) -> dict[str, Any]:
    """Run one paired repetition and return arrays plus JSON-safe metadata."""

    t = np.asarray(t, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    is_validation = np.asarray(is_validation, dtype=bool)
    truth = equation.exact_state(t, x)
    prediction = np.empty_like(truth)
    root_terminal = np.empty_like(truth)
    nonlinear_u = np.empty(len(t), dtype=np.float64)
    diagnostics = WorkDiagnostics()
    chunk_fingerprints: list[str] = []
    start = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(t), chunk_size)):
        end = min(begin + chunk_size, len(t))
        seed_sequence = np.random.SeedSequence(
            [base_seed, equation.d, n, M, repetition, chunk_index]
        )
        solver = FullHistoryMLP(
            equation,
            M,
            method,
            np.random.default_rng(seed_sequence),
            time_beta_alpha=time_beta_alpha,
            trace_draws=trace_draws,
        )
        value = solver.solve(n, t[begin:end], x[begin:end], collect_root=True)
        prediction[begin:end] = value
        if solver.root_terminal is None:
            raise RuntimeError("root terminal state was not recorded")
        root_terminal[begin:end] = solver.root_terminal
        nonlinear_u[begin:end] = value[:, 0] - solver.root_terminal[:, 0]
        diagnostics.merge(solver.stats)
        if solver.draw_fingerprint is not None:
            chunk_fingerprints.append(solver.draw_fingerprint)
    elapsed = time.perf_counter() - start
    all_mask = np.ones(len(t), dtype=bool)
    test_mask = ~is_validation
    metadata = {
        "schema_version": 1,
        "method": method.to_dict(),
        "dimension": equation.d,
        "n": n,
        "M": M,
        "repetition": repetition,
        "paired_seed": base_seed,
        "chunk_size": chunk_size,
        "time_beta_alpha": time_beta_alpha,
        "dtype": "float64",
        "corrected_terminal_ebl": "standard_normal / sqrt(T-t)",
        "zero_level_generator_summand_elided": True,
        "wall_clock_seconds": elapsed,
        "work": diagnostics.summary(),
        "metrics": {
            "all": _error_metrics(prediction, truth, all_mask),
            "validation": _error_metrics(prediction, truth, is_validation),
            "test": _error_metrics(prediction, truth, test_mask),
        },
        "draw_fingerprint": (
            hashlib.sha256("".join(chunk_fingerprints).encode("ascii")).hexdigest()
            if chunk_fingerprints
            else None
        ),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "equation": equation.provenance(),
    }
    return {
        "prediction": prediction,
        "truth": truth,
        "root_terminal": root_terminal,
        "nonlinear_u_correction": nonlinear_u,
        "metadata": metadata,
    }


def save_repetition(path: str, result: dict[str, Any]) -> None:
    """Write an immutable compressed repetition artifact."""

    np.savez_compressed(
        path,
        prediction=result["prediction"],
        truth=result["truth"],
        root_terminal=result["root_terminal"],
        nonlinear_u_correction=result["nonlinear_u_correction"],
        metadata_json=np.asarray(json.dumps(result["metadata"], sort_keys=True)),
    )


def load_repetition(path: str) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        return {
            "prediction": data["prediction"],
            "truth": data["truth"],
            "root_terminal": data["root_terminal"],
            "nonlinear_u_correction": data["nonlinear_u_correction"],
            "metadata": json.loads(str(data["metadata_json"])),
        }


def method_names(methods: Iterable[MethodSpec]) -> list[str]:
    return [method.name for method in methods]
