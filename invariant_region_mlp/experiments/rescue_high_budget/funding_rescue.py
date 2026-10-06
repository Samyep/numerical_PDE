"""Resumable high-budget rescue runner for the 100D nonlinear Funding PDE."""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
RESULT_ROOT = PROJECT_ROOT / "results" / "rescue_high_budget"
RAW_ROOT = RESULT_ROOT / "raw" / "funding"
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common_work_accounting import (  # noqa: E402
    DrawRecorder,
    WorkCounters,
    atomic_json,
    environment_summary,
    git_output,
    metadata_array,
    save_npz_atomic,
    utc_now,
)


REFERENCE_VALUE = 21.299
DEFAULT_BASE_SEED = 20261005


@dataclass(frozen=True)
class FundingMethod:
    name: str
    transform: str
    factor: float = 1.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


RAW = FundingMethod("raw", "raw")
SAMPLEWISE = FundingMethod("samplewise", "samplewise")
Z_ZERO = FundingMethod("z_zero", "z_zero", 0.0)
F_ZERO = FundingMethod("f_zero", "f_zero", 0.0)


def parse_method(name: str) -> FundingMethod:
    fixed = {item.name: item for item in (RAW, SAMPLEWISE, Z_ZERO, F_ZERO)}
    if name in fixed:
        return fixed[name]
    if name.startswith("shrink_c"):
        factor = float(name.removeprefix("shrink_c"))
        return FundingMethod(name, "shrink", factor)
    if name.startswith("tight_a"):
        factor = float(name.removeprefix("tight_a"))
        return FundingMethod(name, "radius_factor", factor)
    raise ValueError(f"unknown Funding method {name!r}")


class FundingEquation:
    d = 100
    T = 0.5
    sigma = 0.2
    mu = 0.06
    R_l = 0.04
    R_b = 0.06

    def terminal(self, x: np.ndarray) -> np.ndarray:
        maximum = np.max(x, axis=-1)
        return np.maximum(maximum - 120.0, 0.0) - 2.0 * np.maximum(
            maximum - 150.0, 0.0
        )

    def generator(self, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        z_sum = np.sum(z, axis=-1)
        return (
            -self.R_l * y
            - ((self.mu - self.R_l) / self.sigma) * z_sum
            + (self.R_b - self.R_l)
            * np.maximum(z_sum / self.sigma - y, 0.0)
        )

    def delta_radius(self, t: np.ndarray) -> np.ndarray:
        return np.exp(0.5 * self.sigma**2 * (self.T - t))

    def provenance(self) -> dict[str, Any]:
        return {
            "problem": "100D nonlinear Funding",
            "d": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "mu": self.mu,
            "R_l": self.R_l,
            "R_b": self.R_b,
            "reference_value": REFERENCE_VALUE,
            "time_distribution": "Beta(1/2,1)",
            "certificate": "||z/(sigma*x)||_2 <= exp(sigma^2*(T-t)/2)",
        }


class FundingFullHistoryMLP:
    """Float64 full-history MLP with corrected EBL and pre-generator projection."""

    def __init__(
        self,
        *,
        M: int,
        method: FundingMethod,
        seed: int,
        time_beta_alpha: float = 0.5,
        trace_draws: bool = True,
    ) -> None:
        self.M = int(M)
        self.method = method
        self.equation = FundingEquation()
        self.alpha = float(time_beta_alpha)
        if self.M < 1 or not 0.0 < self.alpha <= 1.0:
            raise ValueError("invalid M or time beta alpha")
        self.draws = DrawRecorder(np.random.default_rng(seed), enabled=trace_draws)
        self.work = WorkCounters()
        self.root_terminal: np.ndarray | None = None
        self.root_level_u_corrections: np.ndarray | None = None

    def _normal(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.draws.normal(shape)
        self.work.standard_normal_variates += int(np.prod(shape))
        return values

    def _power_time(self, shape: tuple[int, ...]) -> np.ndarray:
        values = self.draws.power(self.alpha, shape)
        self.work.time_uniform_variates += int(np.prod(shape))
        return values

    def _terminal_estimate(self, n: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        b = len(t)
        samples = self.M**n
        remaining = np.maximum(self.equation.T - t, np.finfo(np.float64).tiny)
        gx = self.equation.terminal(x)
        normal = self._normal((b, samples, self.equation.d))
        brownian = np.sqrt(remaining)[:, None, None] * normal
        xt = x[:, None, :] * np.exp(
            (self.equation.mu - 0.5 * self.equation.sigma**2)
            * remaining[:, None, None]
            + self.equation.sigma * brownian
        )
        gt = self.equation.terminal(xt)
        difference = gt - gx[:, None]
        output = np.empty((b, self.equation.d + 1), dtype=np.float64)
        output[:, 0] = gx + np.mean(difference, axis=1)
        output[:, 1:] = np.mean(
            difference[:, :, None]
            * normal
            / np.sqrt(remaining)[:, None, None],
            axis=1,
        )
        self.work.terminal_g_evals += b * (samples + 1)
        self.work.terminal_samples += b * samples
        return output

    def _correct(
        self, state: np.ndarray, t: np.ndarray, x: np.ndarray
    ) -> np.ndarray:
        corrected = np.asarray(state, dtype=np.float64).copy()
        z = corrected[..., 1:]
        denominator = self.equation.sigma * np.maximum(x, 1e-300)
        delta = z / denominator
        norm = np.linalg.norm(delta, axis=-1)
        certified_radius = self.equation.delta_radius(t)
        overshoot = np.maximum(norm - certified_radius, 0.0)
        self.work.correction_opportunities += norm.size
        self.work.constraint_violations += int(np.count_nonzero(overshoot > 0.0))
        self.work.overshoot_energy_sum += float(np.sum(overshoot**2))

        transform = self.method.transform
        if transform == "raw":
            return corrected
        if transform == "samplewise":
            target_radius = certified_radius
        elif transform == "radius_factor":
            target_radius = self.method.factor * certified_radius
        elif transform == "z_zero":
            target_radius = np.zeros_like(certified_radius)
        elif transform == "shrink":
            corrected[..., 1:] = self.method.factor * z
            self.work.activated_states += int(
                np.count_nonzero(np.any(corrected[..., 1:] != z, axis=-1))
            )
            return corrected
        else:
            raise ValueError(f"unsupported transform {transform!r}")

        scale = np.minimum(1.0, target_radius / np.maximum(norm, 1e-300))
        corrected[..., 1:] = denominator * (delta * scale[..., None])
        self.work.activated_states += int(np.count_nonzero(scale < 1.0))
        return corrected

    def _generator(self, state: np.ndarray, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        corrected = self._correct(state, t, x)
        estimate = self.equation.generator(corrected[..., 0], corrected[..., 1:])
        self.work.f_evals += estimate.size
        self.work.nonfinite_generators += int(
            estimate.size - np.count_nonzero(np.isfinite(estimate))
        )
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
        self.work.recursively_evaluated_states += b
        if n <= 0:
            return np.zeros((b, self.equation.d + 1), dtype=np.float64)

        output = self._terminal_estimate(n, t, x)
        if collect_root:
            self.root_terminal = output.copy()
            self.root_level_u_corrections = np.zeros((max(n - 1, 0), b))
        if self.method.transform == "f_zero":
            return output

        remaining = self.equation.T - t
        for level in range(1, n):
            siblings = self.M ** (n - level)
            r = self._power_time((b, siblings))
            dt = remaining[:, None] * r
            normal = self._normal((b, siblings, self.equation.d))
            brownian = np.sqrt(dt)[:, :, None] * normal
            xr = x[:, None, :] * np.exp(
                (self.equation.mu - 0.5 * self.equation.sigma**2)
                * dt[:, :, None]
                + self.equation.sigma * brownian
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

            importance = remaining[:, None] * r ** (1.0 - self.alpha) / self.alpha
            value_increment = np.mean(importance * difference, axis=1)
            ebl = normal / np.sqrt(
                np.maximum(dt, np.finfo(np.float64).tiny)
            )[:, :, None]
            gradient_increment = np.mean(
                importance[:, :, None] * difference[:, :, None] * ebl,
                axis=1,
            )
            output[:, 0] += value_increment
            output[:, 1:] += gradient_increment
            if collect_root and self.root_level_u_corrections is not None:
                self.root_level_u_corrections[level - 1] = value_increment

        self.work.nonfinite_states += int(
            output.size - np.count_nonzero(np.isfinite(output))
        )
        return output


def block_path(
    *, n: int, M: int, method: FundingMethod, start: int, base_seed: int
) -> Path:
    return (
        RAW_ROOT
        / f"seed{base_seed}"
        / f"n{n:02d}_M{M:03d}"
        / method.name.replace(".", "p")
        / f"block_{start:03d}.npz"
    )


def run_block(task: dict[str, Any]) -> dict[str, Any]:
    n = int(task["n"])
    M = int(task["M"])
    start = int(task["start"])
    count = int(task["count"])
    base_seed = int(task["base_seed"])
    method = FundingMethod(**task["method"])
    output_path = Path(task["output"])
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    seed = base_seed + 100000 * n + 1000 * M + start
    solver = FundingFullHistoryMLP(M=M, method=method, seed=seed)
    t = np.zeros(count, dtype=np.float64)
    x = np.full((count, FundingEquation.d), 100.0, dtype=np.float64)
    started = time.perf_counter()
    prediction = solver.solve(n, t, x, collect_root=True)
    elapsed = time.perf_counter() - started
    assert solver.root_terminal is not None
    metadata = {
        "schema_version": 1,
        "problem": "funding",
        "n": n,
        "M": M,
        "method": method.to_dict(),
        "root_start": start,
        "root_count": count,
        "base_seed": base_seed,
        "solver_seed": seed,
        "dtype": "float64",
        "corrected_terminal_ebl": "standard_normal / sqrt(T-t)",
        "zero_level_generator_summand_elided": True,
        "wall_clock_seconds": elapsed,
        "draw_fingerprint": solver.draws.fingerprint,
        "work": solver.work.summary(),
        "equation": FundingEquation().provenance(),
        "environment": environment_summary(),
    }
    save_npz_atomic(
        output_path,
        prediction_state=prediction,
        root_terminal_state=solver.root_terminal,
        nonlinear_u_correction=prediction[:, 0] - solver.root_terminal[:, 0],
        metadata_json=metadata_array(metadata),
    )
    errors = np.abs(prediction[:, 0] - REFERENCE_VALUE)
    return {
        "path": str(output_path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "n": n,
        "M": M,
        "method": method.name,
        "root_start": start,
        "root_count": count,
        "mae": float(np.mean(errors)),
        "wall_clock_seconds": elapsed,
    }


def parse_configs(value: str) -> list[tuple[int, int]]:
    result = []
    for item in value.split(","):
        n_text, m_text = item.split(":", maxsplit=1)
        result.append((int(n_text), int(m_text)))
    return result


def read_block_summary(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        prediction = data["prediction_state"][:, 0]
    return {
        "path": str(path.relative_to(RESULT_ROOT)).replace("\\", "/"),
        "n": metadata["n"],
        "M": metadata["M"],
        "method": metadata["method"]["name"],
        "root_start": metadata["root_start"],
        "root_count": metadata["root_count"],
        "mae": float(np.mean(np.abs(prediction - REFERENCE_VALUE))),
        "wall_clock_seconds": metadata["wall_clock_seconds"],
    }


def run_stage(args: argparse.Namespace) -> Path:
    configs = parse_configs(args.configs)
    methods = [parse_method(item) for item in args.methods.split(",") if item]
    tasks = []
    for n, M in configs:
        for method in methods:
            for start in range(0, args.roots, args.block_size):
                count = min(args.block_size, args.roots - start)
                output = block_path(
                    n=n,
                    M=M,
                    method=method,
                    start=start,
                    base_seed=args.base_seed,
                )
                tasks.append(
                    {
                        "n": n,
                        "M": M,
                        "method": method.to_dict(),
                        "start": start,
                        "count": count,
                        "base_seed": args.base_seed,
                        "output": str(output),
                    }
                )

    completed = []
    pending = []
    for task in tasks:
        path = Path(task["output"])
        if path.exists() and args.resume:
            completed.append(read_block_summary(path))
        elif path.exists():
            raise FileExistsError(f"use --resume to reuse {path}")
        else:
            pending.append(task)

    manifest_path = RESULT_ROOT / f"funding_{args.stage}_manifest.json"
    manifest = {
        "schema_version": 1,
        "problem": "funding",
        "stage": args.stage,
        "status": "running",
        "created_or_resumed_utc": utc_now(),
        "branch": git_output(REPOSITORY_ROOT, "branch", "--show-current"),
        "experiment_commit_at_start": git_output(REPOSITORY_ROOT, "rev-parse", "HEAD"),
        "configs": configs,
        "methods": [method.to_dict() for method in methods],
        "roots": args.roots,
        "block_size": args.block_size,
        "base_seed": args.base_seed,
        "tasks_total": len(tasks),
        "tasks_completed": len(completed),
        "results": sorted(completed, key=lambda item: item["path"]),
    }
    atomic_json(manifest_path, manifest)
    print(
        f"Funding {args.stage}: total={len(tasks)} existing={len(completed)} "
        f"pending={len(pending)} workers={args.workers}",
        flush=True,
    )
    started = time.perf_counter()
    if pending:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_block, task): task for task in pending}
            for index, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                completed.append(result)
                print(
                    f"Funding progress {index}/{len(pending)} n={result['n']} "
                    f"M={result['M']} {result['method']} block={result['root_start']} "
                    f"MAE={result['mae']:.6g} sec={result['wall_clock_seconds']:.3f}",
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--configs", required=True, help="comma-separated n:M pairs")
    parser.add_argument("--methods", default="raw,samplewise")
    parser.add_argument("--roots", type=int, default=30)
    parser.add_argument("--block-size", type=int, default=10)
    parser.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--resume", action="store_true")
    return parser


if __name__ == "__main__":
    run_stage(build_parser().parse_args())
