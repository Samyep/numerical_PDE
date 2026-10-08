"""Study-A adapters around the unchanged full-history MLP recursion."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np

from lqg import BASE_SEED, PROBLEM_LQG, STUDY_A, LQGEquation
from network import FrozenMLP, NetworkCounters


HERE = Path(__file__).resolve().parent
VB_DIR = HERE.parent / "active_vb_high_budget"
if str(VB_DIR) not in sys.path:
    sys.path.insert(0, str(VB_DIR))

# Imported and reused without modification, as required by the frozen protocol.
from vb_mlp_methods import FullHistoryMLP, MethodSpec, WorkDiagnostics  # noqa: E402


METHODS: dict[str, MethodSpec] = {
    "mlp": MethodSpec("mlp", "plain"),
    "mlp_clip": MethodSpec("mlp_clip", "plain_clip", 10.0),
    "scasml": MethodSpec("scasml", "defect_bismut_clip", 0.1),
    "scasml_noclip": MethodSpec("scasml_noclip", "defect_bismut"),
    "path": MethodSpec("path", "defect_path"),
    "path_clip": MethodSpec("path_clip", "defect_path_clip", 0.1),
    "oracle_state": MethodSpec("oracle_state", "defect_oracle"),
}


def _counter_difference(after: dict[str, int], before: dict[str, int]) -> dict[str, int]:
    return {key: int(after[key] - before.get(key, 0)) for key in after}


class StudyAMLP(FullHistoryMLP):
    """Plain or surrogate-defect MLP using the frozen base recursion."""

    def __init__(
        self,
        equation: LQGEquation,
        surrogate: FrozenMLP | None,
        M: int,
        method: MethodSpec,
        rng: np.random.Generator,
        *,
        trace_draws: bool = False,
    ) -> None:
        super().__init__(
            equation,
            M,
            method,
            rng,
            time_beta_alpha=0.5,
            trace_draws=trace_draws,
        )
        self.base_equation = equation
        self.surrogate = surrogate
        self.is_defect = method.transform.startswith("defect")
        if self.is_defect and surrogate is None:
            raise ValueError("defect methods require a surrogate")

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        eq = self.base_equation
        b = len(t)
        samples = self.M**n
        remaining = np.maximum(eq.T - t, 0.0)
        output = np.zeros((b, eq.d + 1), dtype=np.float64)
        at_terminal = remaining <= 4.0 * np.finfo(np.float64).eps
        active = ~at_terminal

        def terminal_fields(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            flat = np.asarray(points, dtype=np.float64).reshape(-1, eq.d)
            if not self.is_defect:
                value = eq.terminal(flat)
                z = eq.sigma_grad_terminal(flat)
            else:
                assert self.surrogate is not None
                value, z = self.surrogate.terminal_defect(eq, flat)
            return value.reshape(points.shape[:-1]), z.reshape(*points.shape[:-1], eq.d)

        if np.any(at_terminal):
            value, z = terminal_fields(x[at_terminal])
            output[at_terminal, 0] = value
            output[at_terminal, 1:] = z
        if not np.any(active):
            return output
        xa = x[active]
        ha = remaining[active]
        ba = len(xa)
        normal = self._normal((ba, samples, eq.d))
        xt = xa[:, None, :] + eq.sigma * np.sqrt(ha)[:, None, None] * normal
        gt, zt = terminal_fields(xt)
        gx, _ = terminal_fields(xa)
        difference = gt - gx[:, None]
        output[active, 0] = gx + np.mean(difference, axis=1)
        if self.method.transform in {"defect_path", "defect_path_clip"}:
            output[active, 1:] = np.mean(zt, axis=1)
        else:
            weights = normal / np.sqrt(np.maximum(ha, np.finfo(np.float64).tiny))[:, None, None]
            output[active, 1:] = np.mean(difference[:, :, None] * weights, axis=1)
        self.stats.terminal_g_evals += ba * samples + ba
        self.stats.terminal_samples += ba * samples
        return output

    def _clip_state(self, state: np.ndarray) -> np.ndarray:
        transform = self.method.transform
        if transform == "plain_clip":
            return np.clip(state, -10.0, 10.0)
        if transform in {"defect_bismut_clip", "defect_path_clip"}:
            return np.clip(state, -0.1, 0.1)
        return state

    def _generator(self, state: np.ndarray, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        eq = self.base_equation
        b, k, width = state.shape
        if width != eq.d + 1:
            raise ValueError("invalid state width")
        count = b * k
        candidate = np.asarray(state, dtype=np.float64)
        corrected = self._clip_state(candidate)
        if corrected is not candidate:
            changed = np.any(corrected != candidate, axis=-1)
            self.stats.correction_opportunities += count
            self.stats.method_geometry_violations += int(np.count_nonzero(changed))
            self.stats.activated_states += int(np.count_nonzero(changed))
            self.stats.method_overshoot_energy_sum += float(np.sum((corrected - candidate) ** 2))

        flat_t = np.asarray(t, dtype=np.float64).reshape(-1)
        flat_x = np.asarray(x, dtype=np.float64).reshape(-1, eq.d)
        if not self.is_defect:
            estimate = eq.generator(corrected[..., 0], corrected[..., 1:])
            self.stats.f_evals += count
            self.stats.nonfinite_generator_values += int(count - np.count_nonzero(np.isfinite(estimate)))
            return estimate

        assert self.surrogate is not None
        residual, surrogate_value, z_surrogate = self.surrogate.residual(flat_t, flat_x)
        residual = residual.reshape(b, k)
        z_surrogate = z_surrogate.reshape(b, k, eq.d)
        if self.method.transform == "defect_oracle":
            exact_u, exact_z = eq.hopf_cole(flat_t, flat_x)
            corrected = np.concatenate(
                [
                    (exact_u - surrogate_value)[:, None],
                    exact_z - z_surrogate.reshape(-1, eq.d),
                ],
                axis=1,
            ).reshape(b, k, eq.d + 1)
        estimate = (
            eq.generator(corrected[..., 0], z_surrogate + corrected[..., 1:])
            - eq.generator(np.zeros((b, k)), z_surrogate)
            + residual
        )
        # The residual and f(z_surrogate) cancel in the mechanism comparison.
        # Use the deterministic Hopf--Cole evaluator only for the true z state.
        _, exact_z = eq.hopf_cole(flat_t, flat_x)
        exact_z = exact_z.reshape(b, k, eq.d)
        truth = (
            eq.generator(np.zeros((b, k)), exact_z)
            - eq.generator(np.zeros((b, k)), z_surrogate)
            + residual
        )
        error = estimate - truth
        finite = np.isfinite(estimate) & np.isfinite(error)
        self.stats.f_evals += count
        self.stats.generator_count += count
        self.stats.nonfinite_generator_values += int(count - np.count_nonzero(finite))
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
        output = super().solve(n, t, x, collect_root=collect_root)
        # The public SCaSML implementation clips the completed state returned
        # by each recursive level, including the root.  These two registered
        # controls reproduce that behavior; no-clipping methods remain raw.
        return self._clip_state(output)


def error_metrics(prediction: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    pred = np.asarray(prediction, dtype=np.float64)[mask]
    ref = np.asarray(truth, dtype=np.float64)[mask]
    finite_rows = np.all(np.isfinite(pred), axis=1)
    if not np.all(finite_rows):
        return {
            "n_points": len(pred),
            "value_relative_l2": float("inf"),
            "skill": float("inf"),
            "value_rmse": float("inf"),
            "value_mae": float("inf"),
            "value_bias": float("nan"),
            "gradient_relative_l2": float("inf"),
            "gradient_mae": float("inf"),
            "nonfinite_state_count": int(np.size(pred) - np.count_nonzero(np.isfinite(pred))),
            "nonfinite_row_count": int(len(pred) - np.count_nonzero(finite_rows)),
        }
    value_error = pred[:, 0] - ref[:, 0]
    gradient_error = pred[:, 1:] - ref[:, 1:]
    denominator = np.linalg.norm(ref[:, 0])
    gradient_denominator = np.linalg.norm(ref[:, 1:])
    rmse = float(np.sqrt(np.mean(value_error**2)))
    reference_std = float(np.std(ref[:, 0], ddof=0))
    return {
        "n_points": len(pred),
        "value_relative_l2": (
            float(np.linalg.norm(value_error) / denominator) if denominator > 0.0 else float("inf")
        ),
        "skill": rmse / reference_std if reference_std > 0.0 else float("inf"),
        "value_rmse": rmse,
        "value_mae": float(np.mean(np.abs(value_error))),
        "value_bias": float(np.mean(value_error)),
        "gradient_relative_l2": (
            float(np.linalg.norm(gradient_error) / gradient_denominator)
            if gradient_denominator > 0.0
            else float("inf")
        ),
        "gradient_mae": float(np.mean(np.abs(gradient_error))),
        "nonfinite_state_count": 0,
        "nonfinite_row_count": 0,
    }


def run_repetition(
    *,
    equation: LQGEquation,
    surrogate: FrozenMLP | None,
    method_name: str,
    checkpoint: int,
    network_seed: int,
    n: int,
    M: int,
    repetition: int,
    t: np.ndarray,
    x: np.ndarray,
    is_validation: np.ndarray,
    reference_state: np.ndarray,
    chunk_size: int = 8,
    trace_draws: bool = False,
) -> tuple[dict[str, Any], np.ndarray]:
    """Run one paired MC repetition and return a CSV-ready row plus prediction."""

    if method_name not in METHODS:
        raise ValueError(f"unknown method {method_name!r}")
    method = METHODS[method_name]
    prediction = np.empty_like(reference_state)
    diagnostics = WorkDiagnostics()
    fingerprints: list[str] = []
    before = NetworkCounters().as_dict() if surrogate is None else surrogate.counters.as_dict()
    started = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(t), chunk_size)):
        end = min(begin + chunk_size, len(t))
        tree_seed = np.random.SeedSequence(
            [
                BASE_SEED,
                STUDY_A,
                PROBLEM_LQG,
                equation.d,
                network_seed,
                50,
                checkpoint,
                n,
                M,
                repetition,
                chunk_index,
            ]
        )
        solver = StudyAMLP(
            equation,
            surrogate,
            M,
            method,
            np.random.default_rng(tree_seed),
            trace_draws=trace_draws,
        )
        estimate = solver.solve(n, t[begin:end], x[begin:end])
        if method.transform.startswith("defect"):
            assert surrogate is not None
            surrogate_u, surrogate_z = surrogate.value_z(t[begin:end], x[begin:end])
            estimate[:, 0] += surrogate_u
            estimate[:, 1:] += surrogate_z
        prediction[begin:end] = estimate
        diagnostics.merge(solver.stats)
        if solver.draw_fingerprint is not None:
            fingerprints.append(solver.draw_fingerprint)
    wall = time.perf_counter() - started
    after = NetworkCounters().as_dict() if surrogate is None else surrogate.counters.as_dict()
    network_work = _counter_difference(after, before)
    test_mask = ~is_validation
    metrics_test = error_metrics(prediction, reference_state, test_mask)
    metrics_validation = error_metrics(prediction, reference_state, is_validation)
    work = diagnostics.summary()
    generator = work.get("generator")
    row: dict[str, Any] = {
        "study": "A",
        "dimension": equation.d,
        "network_seed": network_seed,
        "checkpoint": checkpoint,
        "n": n,
        "M": M,
        "repetition": repetition,
        "method": method_name,
        "laplacian": (
            "none"
            if surrogate is None
            else ("exact" if surrogate.exact_laplacian else f"hutchinson_{len(surrogate.probes)}")
        ),
        "value_relative_l2": metrics_test["value_relative_l2"],
        "skill": metrics_test["skill"],
        "value_rmse": metrics_test["value_rmse"],
        "value_mae": metrics_test["value_mae"],
        "value_bias": metrics_test["value_bias"],
        "gradient_relative_l2": metrics_test["gradient_relative_l2"],
        "gradient_mae": metrics_test["gradient_mae"],
        "validation_value_relative_l2": metrics_validation["value_relative_l2"],
        "wall_clock_seconds": wall,
        "device": "numpy+" + ("none" if surrogate is None else str(surrogate.device)),
        "network_forward_calls": network_work["forward_calls"],
        "network_forward_points": network_work["forward_points"],
        "network_backward_calls": network_work["backward_calls"],
        "network_backward_points": network_work["backward_points"],
        "laplacian_probe_points": network_work["laplacian_probe_points"],
        "generator_calls": work["f_evals"],
        "terminal_samples": work["terminal_samples"],
        "transition_samples": work["transition_samples"],
        "nonfinite_generator_values": work["nonfinite_generator_values"],
        "nonfinite_state_count": metrics_test["nonfinite_state_count"],
        "generator_bias": None if generator is None else generator["bias"],
        "generator_absolute_bias": None if generator is None else generator["absolute_bias"],
        "generator_mae": None if generator is None else generator["mae"],
        "draw_fingerprint": (
            hashlib.sha256("".join(fingerprints).encode("ascii")).hexdigest()
            if fingerprints
            else None
        ),
    }
    return row, prediction


def baseline_row(
    *,
    equation: LQGEquation,
    surrogate: FrozenMLP | None,
    method_name: str,
    checkpoint: int,
    network_seed: int,
    n: int,
    M: int,
    repetition: int,
    t: np.ndarray,
    x: np.ndarray,
    is_validation: np.ndarray,
    reference_state: np.ndarray,
    f_zero_u: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray]:
    started = time.perf_counter()
    before = NetworkCounters().as_dict() if surrogate is None else surrogate.counters.as_dict()
    if method_name == "surrogate":
        if surrogate is None:
            raise ValueError("surrogate baseline requires network")
        u, z = surrogate.value_z(t, x)
        prediction = np.concatenate([u[:, None], z], axis=1)
    elif method_name == "f_zero":
        # The preregistered MC reference supplies the value.  Its z channel is
        # not needed for any criterion, so retain NaNs and report z as N/A.
        prediction = np.concatenate(
            [f_zero_u[:, None], np.full((len(t), equation.d), np.nan)], axis=1
        )
    else:
        raise ValueError(method_name)
    elapsed = time.perf_counter() - started
    after = NetworkCounters().as_dict() if surrogate is None else surrogate.counters.as_dict()
    work = _counter_difference(after, before)
    test = ~is_validation
    value_error = prediction[test, 0] - reference_state[test, 0]
    rel = float(np.linalg.norm(value_error) / np.linalg.norm(reference_state[test, 0]))
    rmse = float(np.sqrt(np.mean(value_error**2)))
    row: dict[str, Any] = {
        "study": "A",
        "dimension": equation.d,
        "network_seed": network_seed,
        "checkpoint": checkpoint,
        "n": n,
        "M": M,
        "repetition": repetition,
        "method": method_name,
        "laplacian": "none",
        "value_relative_l2": rel,
        "skill": rmse / float(np.std(reference_state[test, 0], ddof=0)),
        "value_rmse": rmse,
        "value_mae": float(np.mean(np.abs(value_error))),
        "value_bias": float(np.mean(value_error)),
        "gradient_relative_l2": (
            None
            if method_name == "f_zero"
            else float(
                np.linalg.norm(prediction[test, 1:] - reference_state[test, 1:])
                / np.linalg.norm(reference_state[test, 1:])
            )
        ),
        "gradient_mae": (
            None
            if method_name == "f_zero"
            else float(np.mean(np.abs(prediction[test, 1:] - reference_state[test, 1:])))
        ),
        "validation_value_relative_l2": float(
            np.linalg.norm(prediction[is_validation, 0] - reference_state[is_validation, 0])
            / np.linalg.norm(reference_state[is_validation, 0])
        ),
        "wall_clock_seconds": elapsed,
        "device": "reference" if method_name == "f_zero" else "numpy+" + str(surrogate.device),
        "network_forward_calls": work["forward_calls"],
        "network_forward_points": work["forward_points"],
        "network_backward_calls": work["backward_calls"],
        "network_backward_points": work["backward_points"],
        "laplacian_probe_points": work["laplacian_probe_points"],
        "generator_calls": 0,
        "terminal_samples": 0,
        "transition_samples": 0,
        "nonfinite_generator_values": 0,
        "nonfinite_state_count": 0,
        "generator_bias": None,
        "generator_absolute_bias": None,
        "generator_mae": None,
        "draw_fingerprint": None,
    }
    return row, prediction
