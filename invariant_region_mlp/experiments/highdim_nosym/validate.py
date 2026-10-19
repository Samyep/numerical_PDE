"""Pre-production validation for the frozen B6 benchmark and adapters."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import DrawTape, METHODS
from invariant_region_mlp.experiments.double_estimator.protocol import chunk_size
from invariant_region_mlp.experiments.gpu_mlp.solver import TapeTorchDraws
from invariant_region_mlp.experiments.nonridge.solver import (
    NonRidgeCorrectionBaselineMLP,
)

from .equations import array_sha256, fixed_points
from .protocol import BASE_SEED, CALIBRATION_FREEZE_PATH, RESULTS_ROOT, ensure_directories
from .solver import B6TorchEquation, load_b6_selection, make_equation, make_torch_solver


TERMINAL_TOLERANCE = 1e-12
AUTODIFF_RESIDUAL_TOLERANCE = 1e-10
FINITE_DIFFERENCE_TOLERANCE = 1e-6
MONTE_CARLO_STANDARD_ERRORS = 4.0
GPU_ADAPTER_TOLERANCE = 1e-12
DRAW_REPLAY_TOLERANCE = 1e-10
NONRIDGE_SOURCE_COMMIT = "d384f4e8"


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def _relative_linf(observed: np.ndarray, expected: np.ndarray) -> float:
    denominator = max(float(np.max(np.abs(expected))), np.finfo(float).tiny)
    return float(np.max(np.abs(observed - expected))) / denominator


def validate_calibration_freeze() -> list[dict[str, Any]]:
    payload = json.loads(CALIBRATION_FREEZE_PATH.read_text(encoding="utf-8"))
    rows = []
    for d in (100, 400):
        selected = payload["selections"][str(d)]
        row = selected["selected_row"]
        beta, gamma, _ = load_b6_selection(d)
        checks = {
            "qualified_or_limited": selected["status"] in {
                "qualified", "calibration-limited"
            },
            "selection_matches_row": beta == row["beta"] and gamma == row["gamma"],
            "unequal_norms": row["unequal_norm_range"] > 0.0,
            "solver_not_touched": payload["production_solver_touched_B6"] is False,
        }
        rows.append(
            {
                "stage": "calibration-freeze",
                "pde_id": "B6",
                "dimension": d,
                "beta": beta,
                "gamma": gamma,
                "checks": json.dumps(checks, sort_keys=True),
                "passed": all(checks.values()),
            }
        )
    return rows


def validate_terminal_identity() -> list[dict[str, Any]]:
    records = []
    for d in (100, 400):
        equation = make_equation("B6", d)
        _t, x, _mask = fixed_points("B6", d)
        x = x[:200]
        state = equation.exact_state(np.full(len(x), equation.T), x)
        expected_u, expected_z = equation.terminal_and_z(x)
        u_error = _relative_linf(state[:, 0], expected_u)
        z_error = _relative_linf(state[:, 1:], expected_z)
        maximum = max(u_error, z_error)
        records.append(
            {
                "stage": "terminal-identity",
                "pde_id": "B6",
                "dimension": d,
                "points": len(x),
                "u_relative_linf_error": u_error,
                "z_relative_linf_error": z_error,
                "max_relative_error": maximum,
                "threshold": TERMINAL_TOLERANCE,
                "passed": maximum <= TERMINAL_TOLERANCE,
            }
        )
    return records


def _autodiff_residual_one(
    adapter: B6TorchEquation, t: float, x: np.ndarray
) -> tuple[float, float]:
    point = torch.cat(
        (
            torch.tensor([t], dtype=torch.float64, device=adapter.device),
            torch.tensor(x, dtype=torch.float64, device=adapter.device),
        )
    )

    def scalar(combined: torch.Tensor) -> torch.Tensor:
        time = combined[0]
        state = combined[1:]
        logits = torch.matmul(state, adapter.b.T) + (adapter.T - time) * adapter.n2
        return -(torch.logsumexp(logits, dim=0) - adapter.log_count)

    gradient = torch.func.jacrev(scalar)(point)
    hessian = torch.func.hessian(scalar)(point)
    time_derivative = gradient[0]
    laplacian = torch.diagonal(hessian[1:, 1:]).sum()
    gradient_square = gradient[1:].square().sum()
    residual = time_derivative + laplacian - gradient_square
    scale = time_derivative.abs() + laplacian.abs() + gradient_square.abs()
    return float(residual.detach().cpu()), float(scale.detach().cpu())


def validate_autodiff_residual(device: torch.device) -> list[dict[str, Any]]:
    records = []
    for d in (100, 400):
        equation = make_equation("B6", d)
        adapter = B6TorchEquation(equation, device)
        t, x, is_validation = fixed_points("B6", d)
        selected = np.flatnonzero(~is_validation)[:200]
        residuals = np.empty(len(selected), dtype=np.float64)
        scales = np.empty(len(selected), dtype=np.float64)
        started = time.perf_counter()
        for offset, index in enumerate(selected):
            residuals[offset], scales[offset] = _autodiff_residual_one(
                adapter, float(t[index]), x[index]
            )
        relative = float(
            np.max(np.abs(residuals))
            / max(float(np.max(scales)), np.finfo(float).tiny)
        )
        records.append(
            {
                "stage": "autodiff-pde-residual",
                "pde_id": "B6",
                "dimension": d,
                "points": len(selected),
                "max_absolute_residual": float(np.max(np.abs(residuals))),
                "max_scale": float(np.max(scales)),
                "relative_residual": relative,
                "threshold": AUTODIFF_RESIDUAL_TOLERANCE,
                "elapsed_seconds": time.perf_counter() - started,
                "passed": relative <= AUTODIFF_RESIDUAL_TOLERANCE,
            }
        )
    return records


def validate_finite_difference() -> list[dict[str, Any]]:
    records = []
    step = 1e-4
    for d in (100, 400):
        equation = make_equation("B6", d)
        t, x, is_validation = fixed_points("B6", d)
        selected = np.flatnonzero(~is_validation)[:100]
        tt, xx = t[selected], x[selected]
        rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, d, 66, 3]))
        directions = rng.standard_normal((len(xx), 5, d))
        directions /= np.linalg.norm(directions, axis=-1, keepdims=True)
        expanded_t = np.broadcast_to(tt[:, None], (len(xx), 5))
        finite_difference = (
            equation.exact_u(expanded_t, xx[:, None, :] + step * directions)
            - equation.exact_u(expanded_t, xx[:, None, :] - step * directions)
        ) / (2.0 * step)
        z = equation.exact_z(tt, xx)
        expected = np.einsum("bd,bkd->bk", z / equation.sigma, directions)
        error = _relative_linf(finite_difference, expected)
        records.append(
            {
                "stage": "z-finite-difference",
                "pde_id": "B6",
                "dimension": d,
                "points": len(xx),
                "directions_per_point": 5,
                "finite_difference_step": step,
                "max_relative_error": error,
                "threshold": FINITE_DIFFERENCE_TOLERANCE,
                "passed": error <= FINITE_DIFFERENCE_TOLERANCE,
            }
        )
    return records


def _hopf_cole_component_mc(
    equation: Any,
    t: np.ndarray,
    x: np.ndarray,
    *,
    samples: int = 1_000_000,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Unbiased Hopf--Cole MC over Brownian noise and mixture component.

    Sampling the component index outside the expectation is exactly equivalent
    to evaluating the finite mixture in exp(-g).  Antithetic Gaussian pairs
    are the independent units used for the reported standard errors.
    """

    if samples % 2:
        raise ValueError("antithetic sample count must be even")
    pairs = samples // 2
    estimates = np.empty(len(t), dtype=np.float64)
    errors = np.empty(len(t), dtype=np.float64)
    started = time.perf_counter()
    for point in range(len(t)):
        rng = np.random.default_rng(
            np.random.SeedSequence([BASE_SEED, equation.d, 66, 4, point])
        )
        total = 0.0
        total_square = 0.0
        count = 0
        tau = max(equation.T - float(t[point]), 0.0)
        for begin in range(0, pairs, 100_000):
            size = min(100_000, pairs - begin)
            component = rng.integers(0, equation.K, size=size)
            normal = rng.standard_normal(size)
            vectors = equation.b[component]
            center = vectors @ x[point]
            spread = math.sqrt(2.0 * tau) * np.sqrt(equation.n2[component]) * normal
            pair = 0.5 * (np.exp(center + spread) + np.exp(center - spread))
            total += float(np.sum(pair))
            total_square += float(np.sum(pair * pair))
            count += size
        mean = total / count
        variance = max((total_square - count * mean * mean) / (count - 1), 0.0)
        estimates[point] = -math.log(mean)
        errors[point] = math.sqrt(variance / count) / mean
    return estimates, errors, time.perf_counter() - started


def validate_hopf_cole() -> list[dict[str, Any]]:
    records = []
    for d in (100, 400):
        equation = make_equation("B6", d)
        t, x, is_validation = fixed_points("B6", d)
        selected = np.flatnonzero(~is_validation)[:20]
        tt, xx = t[selected], x[selected]
        expected = equation.exact_u(tt, xx)
        observed, standard_error, elapsed = _hopf_cole_component_mc(
            equation, tt, xx
        )
        standardized = np.abs(observed - expected) / np.maximum(
            standard_error, np.finfo(float).tiny
        )
        maximum = float(np.max(standardized))
        records.append(
            {
                "stage": "hopf-cole-monte-carlo",
                "pde_id": "B6",
                "dimension": d,
                "points": len(xx),
                "samples_per_point": 1_000_000,
                "independent_antithetic_pairs": 500_000,
                "max_standardized_error": maximum,
                "threshold": MONTE_CARLO_STANDARD_ERRORS,
                "elapsed_seconds": elapsed,
                "passed": maximum <= MONTE_CARLO_STANDARD_ERRORS,
            }
        )
    return records


def validate_gpu_adapter(device: torch.device) -> list[dict[str, Any]]:
    records = []
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 66, 5]))
    for d in (100, 400):
        equation = make_equation("B6", d)
        adapter = B6TorchEquation(equation, device)
        t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), 32)
        x = rng.uniform(-1.0, 1.0, (32, d))
        first = rng.standard_normal((32, d + 1))
        second = rng.standard_normal((32, d + 1))
        tt = torch.tensor(t, dtype=torch.float64, device=device)
        xx = torch.tensor(x, dtype=torch.float64, device=device)
        one = torch.tensor(first, dtype=torch.float64, device=device)
        two = torch.tensor(second, dtype=torch.float64, device=device)
        with torch.no_grad():
            gpu = {
                "terminal": adapter.terminal(xx).cpu().numpy(),
                "terminal_gradient": adapter.terminal_gradient_z(xx).cpu().numpy(),
                "exact_state": adapter.exact_state(tt, xx).cpu().numpy(),
                "generator": adapter.generator(one).cpu().numpy(),
                "double_generator": adapter.double_generator(one, two).cpu().numpy(),
            }
        expected = {
            "terminal": equation.terminal(x),
            "terminal_gradient": equation.terminal_gradient_z(x),
            "exact_state": equation.exact_state(t, x),
            "generator": equation.generator(first[:, 0], first[:, 1:]),
            "double_generator": equation.double_generator(
                first[:, 1:], second[:, 1:]
            ),
        }
        errors = {
            key: _relative_linf(gpu[key], expected[key]) for key in gpu
        }
        maximum = max(errors.values())
        records.append(
            {
                "stage": "gpu-adapter",
                "pde_id": "B6",
                "dimension": d,
                **{f"{key}_relative_error": value for key, value in errors.items()},
                "max_relative_error": maximum,
                "threshold": GPU_ADAPTER_TOLERANCE,
                "passed": maximum <= GPU_ADAPTER_TOLERANCE,
            }
        )
    return records


def _numpy_tape_run(
    equation: Any, n: int, M: int, method: str
) -> tuple[np.ndarray, DrawTape, np.ndarray, np.ndarray]:
    d, repetition, chunk_index = equation.d, 0, 0
    t, x, _mask = fixed_points("B6", d)
    size = chunk_size(d)
    t, x = t[:size], x[:size]
    seed_parts = [BASE_SEED, d, n, M, repetition, chunk_index]
    tape = DrawTape("record")
    solver = NonRidgeCorrectionBaselineMLP(
        equation,
        M,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        draw_tape=tape,
    )
    return solver.solve(n, t, x), tape, t, x


def validate_draw_replay(device: torch.device) -> list[dict[str, Any]]:
    records = []
    equation = make_equation("B6", 100)
    for n, M in ((2, 8), (3, 4)):
        for method in ("path", "double_path_base"):
            expected, tape, t, x = _numpy_tape_run(equation, n, M, method)
            replay = tape.replay()
            seed_parts = [BASE_SEED, 100, n, M, 0, 0]
            solver = make_torch_solver(
                equation,
                M,
                method,
                device,
                seed_sequence=np.random.SeedSequence(seed_parts),
                aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
                draw_source=TapeTorchDraws(device, replay),
            )
            with torch.no_grad():
                observed = solver.solve(
                    n,
                    torch.tensor(t.copy(), dtype=torch.float64, device=device),
                    torch.tensor(x.copy(), dtype=torch.float64, device=device),
                ).cpu().numpy()
            u_error = _relative_linf(observed[:, 0], expected[:, 0])
            z_error = _relative_linf(observed[:, 1:], expected[:, 1:])
            maximum = max(u_error, z_error)
            records.append(
                {
                    "stage": "draw-tape-replay",
                    "pde_id": "B6",
                    "dimension": 100,
                    "n": n,
                    "M": M,
                    "method": method,
                    "u_relative_linf_error": u_error,
                    "z_relative_linf_error": z_error,
                    "max_relative_error": maximum,
                    "threshold": DRAW_REPLAY_TOLERANCE,
                    "draw_entries": len(tape.entries),
                    "tape_exhausted": replay.exhausted,
                    "passed": maximum <= DRAW_REPLAY_TOLERANCE and replay.exhausted,
                }
            )
    return records


def validate_s05_bank_hash() -> list[dict[str, Any]]:
    records = []
    for d in (100, 400):
        t, x, mask = fixed_points("S05", d)
        digest = hashlib.sha256()
        digest.update(np.ascontiguousarray(t).tobytes())
        digest.update(np.ascontiguousarray(x).tobytes())
        digest.update(np.ascontiguousarray(mask).tobytes())
        observed = digest.hexdigest()
        source_path = (
            "invariant_region_mlp/results/nonridge/raw/stage_a/"
            f"S05/d{d}/n2_M64/path/rep00.json"
        )
        source = json.loads(_git("show", f"{NONRIDGE_SOURCE_COMMIT}:{source_path}"))
        expected = source["point_fingerprint"]
        records.append(
            {
                "stage": "s05-bank-reuse",
                "pde_id": "S05",
                "dimension": d,
                "point_count": len(t),
                "test_count": int(np.count_nonzero(~mask)),
                "array_sha256": array_sha256(t, x, mask),
                "point_fingerprint": observed,
                "source_point_fingerprint": expected,
                "source_commit": NONRIDGE_SOURCE_COMMIT,
                "passed": observed == expected,
            }
        )
    return records


def validate(device_name: str = "cuda") -> dict[str, Any]:
    ensure_directories()
    device = torch.device(device_name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    stages = (
        ("calibration-freeze", validate_calibration_freeze),
        ("terminal-identity", validate_terminal_identity),
        ("autodiff-pde-residual", lambda: validate_autodiff_residual(device)),
        ("z-finite-difference", validate_finite_difference),
        ("hopf-cole-monte-carlo", validate_hopf_cole),
        ("gpu-adapter", lambda: validate_gpu_adapter(device)),
        ("draw-tape-replay", lambda: validate_draw_replay(device)),
        ("s05-bank-reuse", validate_s05_bank_hash),
    )
    records: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for name, function in stages:
        started = time.perf_counter()
        current = function()
        records.extend(current)
        passed = bool(current) and all(bool(row["passed"]) for row in current)
        summary[name] = {
            "checks": len(current),
            "failures": sum(not bool(row["passed"]) for row in current),
            "verdict": "PASS" if passed else "FAIL",
            "elapsed_seconds": time.perf_counter() - started,
        }
        print(json.dumps({"stage": name, **summary[name]}), flush=True)
        if not passed:
            break
    fields = sorted({key for row in records for key in row})
    with (RESULTS_ROOT / "validation.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    payload = {
        "all_passed": len(summary) == len(stages)
        and all(item["verdict"] == "PASS" for item in summary.values()),
        "stages": summary,
        "device": str(device),
        "calibration_freeze_commit": _git(
            "log", "--diff-filter=A", "-1", "--format=%H", "--",
            str(CALIBRATION_FREEZE_PATH.relative_to(_repository_root())),
        ).strip(),
        "code_commit": _git("rev-parse", "HEAD").strip(),
    }
    from .protocol import write_json_atomic

    write_json_atomic(RESULTS_ROOT / "validation_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = validate(args.device)
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
