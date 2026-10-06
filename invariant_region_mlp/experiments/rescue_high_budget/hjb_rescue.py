"""Adaptive high-budget rescue runner for the Rosenbrock HJB benchmark."""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import jax
import jax.random as jr
import numpy as np
from scipy.special import logsumexp, roots_laguerre

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "rescue_high_budget"
RAW_ROOT = RESULT_ROOT / "raw" / "hjb"
POINT_ROOT = RESULT_ROOT / "test_points"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common_work_accounting import (  # noqa: E402
    DrawRecorder,
    WorkCounters,
    atomic_json,
    combine_fingerprints,
    environment_summary,
    git_output,
    metadata_array,
    save_npz_atomic,
    utc_now,
)


CERTIFIED_RADIUS = math.sqrt(15.0)
DEFAULT_BASE_SEED = 20261007
OLD_RAW_REL_L2_D100 = 2.2799487622119137
OLD_SAMPLEWISE_REL_L2_D100 = 0.6207713181675171


@dataclass(frozen=True)
class HJBMethod:
    name: str
    transform: str
    factor: float = 1.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


RAW = HJBMethod("raw", "raw")
SAMPLEWISE = HJBMethod("samplewise", "samplewise")
F_ZERO = HJBMethod("f_zero", "f_zero", 0.0)


def parse_method(name: str) -> HJBMethod:
    fixed = {item.name: item for item in (RAW, SAMPLEWISE, F_ZERO)}
    if name in fixed:
        return fixed[name]
    if name.startswith("tight_a"):
        return HJBMethod(name, "radius_factor", float(name.removeprefix("tight_a")))
    raise ValueError(f"unknown HJB method {name!r}")


class RosenbrockHJB:
    T = 1.0
    sigma = math.sqrt(2.0)

    def __init__(self, d: int) -> None:
        self.d = int(d)
        # This exactly preserves the existing validated headline construction:
        # JAX default float32 uniforms are cast to float64 afterwards.
        self.c1 = np.asarray(
            jr.uniform(jr.PRNGKey(0), (self.d - 1,), minval=0.5, maxval=1.5),
            dtype=np.float64,
        )
        self.c2 = np.asarray(
            jr.uniform(jr.PRNGKey(1), (self.d - 1,), minval=0.5, maxval=1.5),
            dtype=np.float64,
        )
        matrix = np.zeros((self.d, self.d), dtype=np.float64)
        index = np.arange(self.d - 1)
        matrix[index, index] += self.c1
        matrix[index + 1, index + 1] += self.c1 + self.c2
        matrix[index, index + 1] -= self.c1
        matrix[index + 1, index] -= self.c1
        self.matrix = matrix
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(matrix)
        self.trace = float(np.trace(matrix))
        self.lambda_max = float(self.eigenvalues[-1])

    def quadratic(self, x: np.ndarray) -> np.ndarray:
        difference = x[..., :-1] - x[..., 1:]
        return np.sum(
            self.c1 * difference**2 + self.c2 * x[..., 1:] ** 2,
            axis=-1,
        )

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return np.log((1.0 + self.quadratic(x)) / 2.0)

    @staticmethod
    def generator(z: np.ndarray) -> np.ndarray:
        return -0.5 * np.sum(z**2, axis=-1)

    def exact_state(
        self,
        t: np.ndarray,
        x: np.ndarray,
        *,
        n_quad: int = 64,
        chunk_size: int = 128,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Stable Hopf--Cole/Gauss--Laguerre value and z=sqrt(2) grad(u)."""

        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.d)
        nodes, weights = roots_laguerre(n_quad)
        log_weights = np.log(weights)
        values = np.empty(len(t), dtype=np.float64)
        gradients = np.empty((len(t), self.d), dtype=np.float64)
        for start in range(0, len(t), chunk_size):
            end = min(len(t), start + chunk_size)
            tc = t[start:end]
            xc = x[start:end]
            y = xc @ self.eigenvectors
            tau = self.T - tc
            q = self.quadratic(xc)
            c = 1.0 + q + 2.0 * tau * self.trace
            scaled_nodes = nodes[None, :, None] / c[:, None, None]
            denominator = (
                1.0
                + 4.0
                * tau[:, None, None]
                * scaled_nodes
                * self.eigenvalues[None, None, :]
            )
            log_j = -0.5 * np.sum(np.log(denominator), axis=2) - scaled_nodes[
                :, :, 0
            ] * np.sum(
                self.eigenvalues[None, None, :]
                * y[:, None, :] ** 2
                / denominator,
                axis=2,
            )
            log_terms = (
                log_weights[None, :]
                - np.log(c)[:, None]
                + nodes[None, :] * (1.0 - 1.0 / c)[:, None]
                + log_j
            )
            log_integral = logsumexp(log_terms, axis=1)
            values[start:end] = -math.log(2.0) - log_integral
            alpha = np.exp(log_terms - log_integral[:, None])
            numerator = np.sum(
                alpha[:, :, None] * (scaled_nodes / denominator), axis=1
            )
            gradient_y = 2.0 * (self.eigenvalues[None, :] * y) * numerator
            gradients[start:end] = math.sqrt(2.0) * (
                gradient_y @ self.eigenvectors.T
            )
        return values, gradients

    def provenance(self) -> dict[str, Any]:
        coefficient_hash = hashlib.sha256(
            np.ascontiguousarray(np.concatenate([self.c1, self.c2])).tobytes()
        ).hexdigest()
        return {
            "problem": "Rosenbrock HJB",
            "d": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "driver": "-0.5*||z||^2",
            "certified_radius": CERTIFIED_RADIUS,
            "matrix_construction": "jax.random.PRNGKey(0/1), uniform[0.5,1.5), float32 then cast float64",
            "jax_version": jax.__version__,
            "coefficient_sha256": coefficient_hash,
            "trace_A": self.trace,
            "lambda_max": self.lambda_max,
            "reference": "stable scaled Hopf--Cole / Gauss--Laguerre quadrature",
            "time_distribution": "Uniform(0,1), matching validated headline",
        }


def make_test_points(
    d: int, n_interior: int, n_boundary: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed + d)
    interior = rng.normal(size=(n_interior, d))
    interior /= np.linalg.norm(interior, axis=1, keepdims=True)
    interior *= rng.random(n_interior)[:, None] ** (1.0 / d)
    boundary = rng.normal(size=(n_boundary, d))
    boundary /= np.linalg.norm(boundary, axis=1, keepdims=True)
    x = np.vstack([interior, boundary]).astype(np.float64)
    t = rng.random(n_interior + n_boundary).astype(np.float64)
    return t, x


class HJBFullHistoryMLP:
    """Float64 full-history MLP; projection acts only before nonlinear reuse."""

    def __init__(
        self,
        equation: RosenbrockHJB,
        *,
        M: int,
        method: HJBMethod,
        seed: int | np.random.SeedSequence,
        diagnostic_limit: int = 0,
        trace_draws: bool = True,
    ) -> None:
        self.equation = equation
        self.M = int(M)
        self.method = method
        self.draws = DrawRecorder(np.random.default_rng(seed), enabled=trace_draws)
        self.work = WorkCounters()
        self.diagnostic_limit = int(diagnostic_limit)
        self.diag_t: list[np.ndarray] = []
        self.diag_x: list[np.ndarray] = []
        self.diag_f: list[np.ndarray] = []
        self.diag_count = 0
        self.root_terminal: np.ndarray | None = None

    def _normal(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.draws.normal(shape)
        self.work.standard_normal_variates += int(np.prod(shape))
        return values

    def _uniform(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.draws.uniform(shape)
        self.work.time_uniform_variates += int(np.prod(shape))
        return values

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        b = len(t)
        samples = self.M**n
        remaining = np.maximum(
            self.equation.T - t, np.finfo(np.float64).tiny
        )
        gx = self.equation.terminal(x)
        normal = self._normal((b, samples, self.equation.d))
        xt = (
            x[:, None, :]
            + self.equation.sigma
            * np.sqrt(remaining)[:, None, None]
            * normal
        )
        gt = self.equation.terminal(xt)
        output = np.empty((b, self.equation.d + 1), dtype=np.float64)
        output[:, 0] = np.mean(gt, axis=1)
        output[:, 1:] = np.mean(
            (gt - gx[:, None])[:, :, None]
            * normal
            / np.sqrt(remaining)[:, None, None],
            axis=1,
        )
        self.work.terminal_g_evals += b * (samples + 1)
        self.work.terminal_samples += b * samples
        return output

    def _correct(self, state: np.ndarray) -> np.ndarray:
        corrected = np.asarray(state, dtype=np.float64).copy()
        z = corrected[..., 1:]
        norm = np.linalg.norm(z, axis=-1)
        certified_overshoot = np.maximum(norm - CERTIFIED_RADIUS, 0.0)
        self.work.correction_opportunities += norm.size
        self.work.constraint_violations += int(
            np.count_nonzero(certified_overshoot > 0.0)
        )
        self.work.overshoot_energy_sum += float(np.sum(certified_overshoot**2))
        if self.method.transform == "raw":
            return corrected
        if self.method.transform == "samplewise":
            radius = CERTIFIED_RADIUS
        elif self.method.transform == "radius_factor":
            radius = self.method.factor * CERTIFIED_RADIUS
        else:
            raise ValueError(f"unsupported HJB transform {self.method.transform!r}")
        scale = np.minimum(1.0, radius / np.maximum(norm, 1e-300))
        corrected[..., 1:] = z * scale[..., None]
        self.work.activated_states += int(np.count_nonzero(scale < 1.0))
        return corrected

    def _capture_diagnostic(
        self,
        t: np.ndarray,
        x: np.ndarray,
        estimate: np.ndarray,
    ) -> None:
        remaining = self.diagnostic_limit - self.diag_count
        if remaining <= 0:
            return
        tf = np.asarray(t).reshape(-1)
        xf = np.asarray(x).reshape(-1, self.equation.d)
        ff = np.asarray(estimate).reshape(-1)
        take = min(remaining, len(tf))
        self.diag_t.append(tf[:take].copy())
        self.diag_x.append(xf[:take].copy())
        self.diag_f.append(ff[:take].copy())
        self.diag_count += take

    def _generator(self, state: np.ndarray, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        corrected = self._correct(state)
        estimate = self.equation.generator(corrected[..., 1:])
        self.work.f_evals += estimate.size
        self.work.nonfinite_generators += int(
            estimate.size - np.count_nonzero(np.isfinite(estimate))
        )
        self._capture_diagnostic(t, x, estimate)
        return estimate

    def diagnostic_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if not self.diag_t:
            return (
                np.empty(0),
                np.empty((0, self.equation.d)),
                np.empty(0),
            )
        return (
            np.concatenate(self.diag_t),
            np.concatenate(self.diag_x),
            np.concatenate(self.diag_f),
        )

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
        self.work.recursively_evaluated_states += b
        if n <= 0:
            return np.zeros((b, self.equation.d + 1), dtype=np.float64)
        output = self._terminal_estimate(n, t, x)
        if collect_root:
            self.root_terminal = output.copy()
        if self.method.transform == "f_zero":
            return output

        remaining = self.equation.T - t
        for level in range(1, n):
            siblings = self.M ** (n - level)
            r = self._uniform((b, siblings))
            dt = remaining[:, None] * r
            normal = self._normal((b, siblings, self.equation.d))
            xr = (
                x[:, None, :]
                + self.equation.sigma * np.sqrt(dt)[:, :, None] * normal
            )
            tr = t[:, None] + dt
            self.work.transition_samples += b * siblings
            high = self.solve(level, tr.reshape(-1), xr.reshape(-1, self.equation.d))
            high = high.reshape(b, siblings, self.equation.d + 1)
            difference = self._generator(high, tr, xr)
            if level > 1:
                low = self.solve(
                    level - 1, tr.reshape(-1), xr.reshape(-1, self.equation.d)
                )
                low = low.reshape(b, siblings, self.equation.d + 1)
                difference -= self._generator(low, tr, xr)
            value_increment = remaining * np.mean(difference, axis=1)
            ebl = normal / np.sqrt(
                np.maximum(dt, np.finfo(np.float64).tiny)
            )[:, :, None]
            gradient_increment = remaining[:, None] * np.mean(
                difference[:, :, None] * ebl, axis=1
            )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment

        self.work.nonfinite_states += int(
            output.size - np.count_nonzero(np.isfinite(output))
        )
        return output


def point_path(stage: str, d: int, n_interior: int, n_boundary: int, seed: int) -> Path:
    return POINT_ROOT / (
        f"hjb_{stage}_d{d:03d}_i{n_interior}_b{n_boundary}_seed{seed}.npz"
    )


def prepare_points(
    stage: str,
    d: int,
    n_interior: int,
    n_boundary: int,
    seed: int,
    *,
    n_quad: int = 64,
) -> Path:
    path = point_path(stage, d, n_interior, n_boundary, seed)
    if path.exists():
        return path
    equation = RosenbrockHJB(d)
    t, x = make_test_points(d, n_interior, n_boundary, seed)
    truth_u, truth_z = equation.exact_state(t, x, n_quad=n_quad)
    save_npz_atomic(
        path,
        t=t,
        x=x,
        truth_u=truth_u,
        truth_z=truth_z,
        is_boundary=np.concatenate(
            [np.zeros(n_interior, dtype=bool), np.ones(n_boundary, dtype=bool)]
        ),
        provenance_json=np.asarray(
            json.dumps(
                {
                    "stage": stage,
                    "dimension": d,
                    "n_interior": n_interior,
                    "n_boundary": n_boundary,
                    "seed": seed,
                    "n_quad": n_quad,
                    "equation": equation.provenance(),
                },
                sort_keys=True,
            )
        ),
    )
    return path


def dynamic_chunk_size(d: int, n: int, M: int, requested: int) -> int:
    root_samples = M**n
    memory_limited = max(1, int(120_000_000 / max(16 * d * root_samples, 1)))
    return max(1, min(requested, memory_limited))


def task_path(stage: str, d: int, n: int, M: int, method: HJBMethod, rep: int) -> Path:
    return (
        RAW_ROOT
        / stage
        / f"d{d:03d}_n{n:02d}_M{M:03d}_{method.name.replace('.', 'p')}_rep{rep:02d}.npz"
    )


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    output_path = Path(task["output"])
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    d, n, M, rep = (int(task[key]) for key in ("d", "n", "M", "rep"))
    method = HJBMethod(**task["method"])
    equation = RosenbrockHJB(d)
    with np.load(task["point_path"], allow_pickle=False) as points:
        t = points["t"]
        x = points["x"]
        truth_u = points["truth_u"]
        truth_z = points["truth_z"]

    chunk_size = dynamic_chunk_size(d, n, M, int(task["chunk_size"]))
    prediction = np.empty((len(t), d + 1), dtype=np.float64)
    root_terminal_u = np.empty(len(t), dtype=np.float64)
    combined_work = WorkCounters()
    fingerprints: list[str] = []
    diag_t: list[np.ndarray] = []
    diag_x: list[np.ndarray] = []
    diag_f: list[np.ndarray] = []
    diag_remaining = int(task["diagnostic_limit"])
    started = time.perf_counter()
    for start in range(0, len(t), chunk_size):
        end = min(len(t), start + chunk_size)
        seed = np.random.SeedSequence(
            [int(task["base_seed"]), d, n, M, rep, start]
        )
        solver = HJBFullHistoryMLP(
            equation,
            M=M,
            method=method,
            seed=seed,
            diagnostic_limit=diag_remaining,
        )
        state = solver.solve(n, t[start:end], x[start:end], collect_root=True)
        prediction[start:end] = state
        assert solver.root_terminal is not None
        root_terminal_u[start:end] = solver.root_terminal[:, 0]
        combined_work.merge(solver.work)
        if solver.draws.fingerprint is not None:
            fingerprints.append(solver.draws.fingerprint)
        td, xd, fd = solver.diagnostic_arrays()
        if len(td):
            diag_t.append(td)
            diag_x.append(xd)
            diag_f.append(fd)
            diag_remaining -= len(td)
    elapsed = time.perf_counter() - started

    value_error = prediction[:, 0] - truth_u
    gradient_error = prediction[:, 1:] - truth_z
    gradient_error_sq = np.sum(gradient_error**2, axis=1)
    truth_gradient_sq = np.sum(truth_z**2, axis=1)
    generator_diagnostic: dict[str, Any] | None = None
    if diag_t:
        td = np.concatenate(diag_t)
        xd = np.concatenate(diag_x)
        estimate = np.concatenate(diag_f)
        _, z_true_diag = equation.exact_state(td, xd, n_quad=32)
        truth_f = equation.generator(z_true_diag)
        error = estimate - truth_f
        generator_diagnostic = {
            "count": len(estimate),
            "mse": float(np.mean(error**2)),
            "bias": float(np.mean(error)),
            "absolute_bias": float(abs(np.mean(error))),
            "mae": float(np.mean(np.abs(error))),
            "mean_estimate": float(np.mean(estimate)),
            "mean_truth": float(np.mean(truth_f)),
            "mean_truth_squared": float(np.mean(truth_f**2)),
        }

    metadata = {
        "schema_version": 1,
        "problem": "hjb",
        "stage": task["stage"],
        "dimension": d,
        "n": n,
        "M": M,
        "repetition": rep,
        "method": method.to_dict(),
        "base_seed": task["base_seed"],
        "point_path": str(Path(task["point_path"]).relative_to(RESULT_ROOT)).replace(
            "\\", "/"
        ),
        "n_points": len(t),
        "chunk_size": chunk_size,
        "dtype": "float64",
        "corrected_terminal_ebl": "standard_normal / sqrt(T-t)",
        "zero_level_generator_summand_elided": True,
        "heuristic_clipping": False,
        "wall_clock_seconds": elapsed,
        "draw_fingerprint": combine_fingerprints(fingerprints),
        "work": combined_work.summary(),
        "metrics": {
            "value_relative_l2": float(
                np.linalg.norm(value_error) / np.linalg.norm(truth_u)
            ),
            "value_mae": float(np.mean(np.abs(value_error))),
            "value_bias": float(np.mean(value_error)),
            "gradient_relative_l2": float(
                np.sqrt(np.sum(gradient_error_sq) / np.sum(truth_gradient_sq))
            ),
            "gradient_mae": float(np.mean(np.sqrt(gradient_error_sq))),
        },
        "generator_diagnostic": generator_diagnostic,
        "equation": equation.provenance(),
        "environment": environment_summary(),
    }
    save_npz_atomic(
        output_path,
        prediction_u=prediction[:, 0],
        truth_u=truth_u,
        gradient_error_squared=gradient_error_sq,
        truth_gradient_squared=truth_gradient_sq,
        root_terminal_u=root_terminal_u,
        nonlinear_u_correction=prediction[:, 0] - root_terminal_u,
        metadata_json=metadata_array(metadata),
    )
    return {
        "path": str(output_path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "d": d,
        "n": n,
        "M": M,
        "method": method.name,
        "rep": rep,
        "value_relative_l2": metadata["metrics"]["value_relative_l2"],
        "wall_clock_seconds": elapsed,
    }


def read_task_summary(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
    return {
        "path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "d": metadata["dimension"],
        "n": metadata["n"],
        "M": metadata["M"],
        "method": metadata["method"]["name"],
        "rep": metadata["repetition"],
        "value_relative_l2": metadata["metrics"]["value_relative_l2"],
        "wall_clock_seconds": metadata["wall_clock_seconds"],
    }


def parse_configs(value: str) -> list[tuple[int, int]]:
    result = []
    for item in value.split(","):
        n_text, m_text = item.split(":", maxsplit=1)
        result.append((int(n_text), int(m_text)))
    return result


def run_stage(args: argparse.Namespace) -> Path:
    dimensions = [int(item) for item in args.dimensions.split(",") if item]
    configs = parse_configs(args.configs)
    methods = [parse_method(item) for item in args.methods.split(",") if item]
    point_paths = {
        d: prepare_points(
            args.stage,
            d,
            args.n_interior,
            args.n_boundary,
            args.point_seed,
        )
        for d in dimensions
    }
    tasks = []
    for d in dimensions:
        for n, M in configs:
            for rep in range(args.repetitions):
                for method in methods:
                    output = task_path(args.stage, d, n, M, method, rep)
                    tasks.append(
                        {
                            "stage": args.stage,
                            "d": d,
                            "n": n,
                            "M": M,
                            "rep": rep,
                            "method": method.to_dict(),
                            "point_path": str(point_paths[d]),
                            "base_seed": args.base_seed,
                            "chunk_size": args.chunk_size,
                            "diagnostic_limit": args.diagnostic_limit,
                            "output": str(output),
                        }
                    )

    completed = []
    pending = []
    for task in tasks:
        path = Path(task["output"])
        if path.exists() and args.resume:
            completed.append(read_task_summary(path))
        elif path.exists():
            raise FileExistsError(f"use --resume to reuse {path}")
        else:
            pending.append(task)

    manifest_path = RESULT_ROOT / f"hjb_{args.stage}_manifest.json"
    manifest = {
        "schema_version": 1,
        "problem": "hjb",
        "stage": args.stage,
        "status": "running",
        "created_or_resumed_utc": utc_now(),
        "branch": git_output(REPOSITORY_ROOT, "branch", "--show-current"),
        "experiment_commit_at_start": git_output(REPOSITORY_ROOT, "rev-parse", "HEAD"),
        "dimensions": dimensions,
        "configs": configs,
        "methods": [method.to_dict() for method in methods],
        "repetitions": args.repetitions,
        "point_geometry": {
            "n_interior": args.n_interior,
            "n_boundary": args.n_boundary,
            "point_seed": args.point_seed,
        },
        "base_seed": args.base_seed,
        "tasks_total": len(tasks),
        "tasks_completed": len(completed),
        "results": sorted(completed, key=lambda item: item["path"]),
    }
    atomic_json(manifest_path, manifest)
    print(
        f"HJB {args.stage}: total={len(tasks)} existing={len(completed)} "
        f"pending={len(pending)} workers={args.workers}",
        flush=True,
    )
    started = time.perf_counter()
    if pending:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_task, task): task for task in pending}
            for index, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                completed.append(result)
                print(
                    f"HJB progress {index}/{len(pending)} d={result['d']} n={result['n']} "
                    f"M={result['M']} {result['method']} rep={result['rep']} "
                    f"relL2={result['value_relative_l2']:.6g} "
                    f"sec={result['wall_clock_seconds']:.3f}",
                    flush=True,
                )
                manifest["tasks_completed"] = len(completed)
                manifest["results"] = sorted(completed, key=lambda item: item["path"])
                manifest["elapsed_this_invocation_seconds"] = time.perf_counter() - started
                atomic_json(manifest_path, manifest)
    manifest["status"] = "complete"
    manifest["completed_utc"] = utc_now()
    manifest["tasks_completed"] = len(completed)
    manifest["results"] = sorted(completed, key=lambda item: item["path"])
    manifest["elapsed_this_invocation_seconds"] = time.perf_counter() - started
    atomic_json(manifest_path, manifest)
    return manifest_path


def run_legacy_reproduction(args: argparse.Namespace) -> Path:
    d, n, M = 100, 2, 10
    equation = RosenbrockHJB(d)
    t, x = make_test_points(d, 1000, 200, 20261005)
    truth_u, _ = equation.exact_state(t, x, n_quad=32)
    rows = []
    predictions: dict[str, list[np.ndarray]] = {"raw": [], "samplewise": []}
    started = time.perf_counter()
    for rep in range(10):
        for method in (RAW, SAMPLEWISE):
            solver = HJBFullHistoryMLP(
                equation,
                M=M,
                method=method,
                seed=910000 + d * 100 + rep,
                diagnostic_limit=0,
            )
            values = np.empty(len(t), dtype=np.float64)
            for start in range(0, len(t), 96):
                end = min(len(t), start + 96)
                values[start:end] = solver.solve(n, t[start:end], x[start:end])[:, 0]
            predictions[method.name].append(values)
    for method in (RAW, SAMPLEWISE):
        pred = np.stack(predictions[method.name])
        relative = np.linalg.norm(pred - truth_u[None, :], axis=1) / np.linalg.norm(
            truth_u
        )
        expected = (
            OLD_RAW_REL_L2_D100 if method is RAW else OLD_SAMPLEWISE_REL_L2_D100
        )
        rows.append(
            {
                "method": method.name,
                "relative_l2_mean": float(np.mean(relative)),
                "relative_l2_values": relative.tolist(),
                "archived_mean": expected,
                "absolute_difference": float(abs(np.mean(relative) - expected)),
            }
        )
    payload = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "protocol": {
            "d": d,
            "n": n,
            "M": M,
            "repetitions": 10,
            "points": 1200,
            "chunk_size": 96,
            "seed_schedule": "910000 + d*100 + rep",
            "quadrature_nodes": 32,
        },
        "rows": rows,
        "max_absolute_difference": max(row["absolute_difference"] for row in rows),
        "elapsed_seconds": time.perf_counter() - started,
        "equation": equation.provenance(),
    }
    output = RESULT_ROOT / "hjb_legacy_reproduction.json"
    atomic_json(output, payload)
    print(json.dumps(payload, indent=2), flush=True)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    stage = subparsers.add_parser("stage")
    stage.add_argument("--stage", required=True)
    stage.add_argument("--dimensions", default="100")
    stage.add_argument("--configs", required=True)
    stage.add_argument("--methods", default="raw,samplewise,f_zero")
    stage.add_argument("--repetitions", type=int, default=3)
    stage.add_argument("--n-interior", type=int, default=300)
    stage.add_argument("--n-boundary", type=int, default=60)
    stage.add_argument("--point-seed", type=int, default=20261005)
    stage.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    stage.add_argument("--chunk-size", type=int, default=32)
    stage.add_argument("--diagnostic-limit", type=int, default=512)
    stage.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    stage.add_argument("--resume", action="store_true")
    subparsers.add_parser("legacy")
    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    if arguments.command == "legacy":
        run_legacy_reproduction(arguments)
    else:
        run_stage(arguments)
