"""Pre-MLP equation validation V-1 through V-6."""

from __future__ import annotations

import csv
import json
import math
import platform
import subprocess
import sys
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import DrawTape, METHODS
from invariant_region_mlp.experiments.double_estimator.equations import BASE_SEED
from invariant_region_mlp.experiments.double_estimator.protocol import chunk_size
from invariant_region_mlp.experiments.gpu_mlp.solver import TapeTorchDraws, torch_seed
from invariant_region_mlp.experiments.mechanism_suite.equations import stable_seed

from .equations import (
    FAMILIES,
    PROFILE_SLOPES,
    SpectralSeparableHJB,
    fixed_points,
    make_equation,
)
from .protocol import RESULTS_ROOT, REVISION, write_json_atomic
from .solver import (
    NonRidgeCorrectionBaselineMLP,
    NonRidgeTorchMLP,
    SpectralTorchEquation,
    current_commit,
)


V1_TOLERANCE = 1e-11
V2_TOLERANCE = 1e-11
V3_TOLERANCE = 1e-6
V4_STANDARD_ERRORS = 4.0
V5_TOLERANCE = 1e-12
V6_TOLERANCE = 1e-10


def _device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _relative_linf(observed: np.ndarray, expected: np.ndarray) -> float:
    numerator = float(np.max(np.abs(observed - expected)))
    denominator = max(float(np.max(np.abs(expected))), np.finfo(np.float64).tiny)
    return numerator / denominator


def record_environment() -> dict[str, Any]:
    device = _device()
    payload: dict[str, Any] = {
        "nonridge_revision": REVISION,
        "code_commit": current_commit(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "device": str(device),
        "dtype": "float64",
        "dct_implementation": "Makhoul FFT type-II DCT/IDCT in PyTorch",
    }
    if torch.cuda.is_available():
        properties = torch.cuda.get_device_properties(torch.cuda.current_device())
        payload.update(
            {
                "gpu_name": properties.name,
                "gpu_total_memory_bytes": properties.total_memory,
                "gpu_compute_capability": f"{properties.major}.{properties.minor}",
                "fp64_throughput_class": "consumer GPU; reduced-rate FP64",
            }
        )
        try:
            completed = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=driver_version",
                    "--format=csv,noheader",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            payload["nvidia_driver"] = completed.stdout.strip().splitlines()[0]
        except (OSError, subprocess.CalledProcessError, IndexError):
            payload["nvidia_driver"] = None
    write_json_atomic(RESULTS_ROOT / "environment.json", payload)
    return payload


def validate_v1() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pde_id in FAMILIES:
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            t, x, _mask = fixed_points(pde_id, d)
            reference = equation.exact_state_order(t[:200], x[:200], order=192)
            observed = equation.exact_state_order(t[:200], x[:200], order=96)
            u_error = _relative_linf(observed[:, 0], reference[:, 0])
            z_error = _relative_linf(observed[:, 1:], reference[:, 1:])
            maximum = max(u_error, z_error)
            records.append(
                {
                    "stage": "V-1",
                    "pde_id": pde_id,
                    "dimension": d,
                    "points": 200,
                    "u_relative_linf_error": u_error,
                    "z_relative_linf_error": z_error,
                    "max_relative_error": maximum,
                    "threshold": V1_TOLERANCE,
                    "passed": maximum <= V1_TOLERANCE,
                }
            )
    return records


def validate_v2() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    slopes = PROFILE_SLOPES
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 2]))
    for d in (100, 400):
        equation = make_equation("S0", d)
        t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), 200)
        x = rng.uniform(-1.0, 1.0, (200, d))
        y = equation.rotate(x)
        coefficients = np.zeros(d, dtype=np.float64)
        coefficients[0] = 1.0
        v, derivative = equation.one_dimensional_state(
            coefficients[None, :],
            (equation.T - t)[:, None],
            y,
            order=96,
        )
        observed_u = np.sum(v, axis=1)
        observed_z = equation.sigma * equation.inverse_rotate(derivative)
        logits = (
            y[:, :1] * slopes
            + (equation.T - t)[:, None] * slopes**2
        )
        normalization = np.logaddexp.reduce(logits, axis=1)
        expected_u = -(normalization - math.log(3.0))
        weights = np.exp(logits - normalization[:, None])
        coefficient = -np.sum(weights * slopes, axis=1)
        gradient_y = np.zeros_like(y)
        gradient_y[:, 0] = coefficient
        expected_z = equation.sigma * equation.inverse_rotate(gradient_y)
        u_error = _relative_linf(observed_u, expected_u)
        z_error = _relative_linf(observed_z, expected_z)
        maximum = max(u_error, z_error)
        records.append(
            {
                "stage": "V-2",
                "pde_id": "P1-special",
                "dimension": d,
                "points": 200,
                "u_relative_linf_error": u_error,
                "z_relative_linf_error": z_error,
                "max_relative_error": maximum,
                "threshold": V2_TOLERANCE,
                "passed": maximum <= V2_TOLERANCE,
            }
        )
    return records


def validate_v3() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 3]))
    step = 1e-4
    available = np.concatenate(
        [make_equation(pde_id, d).c for pde_id in FAMILIES for d in (100, 400)]
    )
    c = rng.choice(available, size=100, replace=True)
    tau = rng.uniform(0.01, 0.24, 100)
    y = rng.uniform(-2.0, 2.0, 100)
    value, derivative = SpectralSeparableHJB.one_dimensional_state(c, tau, y)
    tau_plus, _ = SpectralSeparableHJB.one_dimensional_state(c, tau + step, y)
    tau_minus, _ = SpectralSeparableHJB.one_dimensional_state(c, tau - step, y)
    y_plus, _ = SpectralSeparableHJB.one_dimensional_state(c, tau, y + step)
    y_minus, _ = SpectralSeparableHJB.one_dimensional_state(c, tau, y - step)
    d_tau = (tau_plus - tau_minus) / (2.0 * step)
    d_yy = (y_plus - 2.0 * value + y_minus) / (step * step)
    residual = d_tau - d_yy + derivative**2
    residual_max = float(np.max(np.abs(residual)))
    records.append(
        {
            "stage": "V-3",
            "pde_id": "one-dimensional-v",
            "dimension": 1,
            "points": 100,
            "finite_difference_step": step,
            "max_absolute_error": residual_max,
            "threshold": V3_TOLERANCE,
            "passed": residual_max <= V3_TOLERANCE,
        }
    )

    for pde_id in FAMILIES:
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            t, x, is_validation = fixed_points(pde_id, d)
            selected = np.flatnonzero(~is_validation)[:100]
            tt, xx = t[selected], x[selected]
            z = equation.exact_z(tt, xx)
            directions = rng.normal(size=(len(xx), 5, d))
            directions /= np.linalg.norm(directions, axis=-1, keepdims=True)
            plus = xx[:, None, :] + step * directions
            minus = xx[:, None, :] - step * directions
            expanded_t = np.broadcast_to(tt[:, None], (len(xx), 5))
            finite_difference = (
                equation.exact_u(expanded_t, plus)
                - equation.exact_u(expanded_t, minus)
            ) / (2.0 * step)
            expected = np.einsum("bd,bkd->bk", z / equation.sigma, directions)
            error = _relative_linf(finite_difference, expected)
            records.append(
                {
                    "stage": "V-3",
                    "pde_id": pde_id,
                    "dimension": d,
                    "points": len(xx),
                    "directions_per_point": 5,
                    "finite_difference_step": step,
                    "max_relative_error": error,
                    "threshold": V3_TOLERANCE,
                    "passed": error <= V3_TOLERANCE,
                }
            )
    return records


def _monte_carlo_hopf_cole(
    equation: SpectralSeparableHJB,
    t: np.ndarray,
    x: np.ndarray,
    *,
    samples: int = 1_000_000,
    batch_size: int = 4096,
) -> tuple[np.ndarray, np.ndarray, float]:
    device = _device()
    adapter = SpectralTorchEquation(equation, device)
    y = adapter.rotate(torch.as_tensor(x, dtype=torch.float64, device=device))
    scale = torch.sqrt(
        2.0
        * torch.clamp(
            equation.T - torch.as_tensor(t, dtype=torch.float64, device=device),
            min=0.0,
        )
    )
    generator = torch.Generator(device=device)
    seed = np.random.SeedSequence(
        [BASE_SEED, equation.d, stable_seed(equation.pde_id), 4]
    )
    generator.manual_seed(torch_seed(seed))
    total = torch.zeros(len(x), dtype=torch.float64, device=device)
    total_square = torch.zeros_like(total)
    started = time.perf_counter()
    with torch.no_grad():
        for begin in range(0, samples, batch_size):
            count = min(batch_size, samples - begin)
            normal = torch.randn(
                (count, equation.d),
                dtype=torch.float64,
                device=device,
                generator=generator,
            )
            shifted = y[:, None, :] + scale[:, None, None] * normal[None, :, :]
            terminal = torch.sum(adapter.c * adapter.profile(shifted), dim=-1)
            weights = torch.exp(-terminal)
            total += torch.sum(weights, dim=1)
            total_square += torch.sum(weights.square(), dim=1)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - started
    mean = total / samples
    variance = torch.clamp(
        (total_square - samples * mean.square()) / (samples - 1), min=0.0
    )
    estimate = -torch.log(mean)
    standard_error = torch.sqrt(variance / samples) / mean
    return estimate.cpu().numpy(), standard_error.cpu().numpy(), elapsed


def validate_v4() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pde_id in FAMILIES:
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            t, x, is_validation = fixed_points(pde_id, d)
            selected = np.flatnonzero(~is_validation)[:20]
            tt, xx = t[selected], x[selected]
            expected = equation.exact_u(tt, xx)
            observed, standard_error, elapsed = _monte_carlo_hopf_cole(
                equation, tt, xx
            )
            standardized = np.abs(observed - expected) / np.maximum(
                standard_error, np.finfo(np.float64).tiny
            )
            maximum = float(np.max(standardized))
            records.append(
                {
                    "stage": "V-4",
                    "pde_id": pde_id,
                    "dimension": d,
                    "points": len(xx),
                    "monte_carlo_samples": 1_000_000,
                    "max_standardized_error": maximum,
                    "threshold": V4_STANDARD_ERRORS,
                    "elapsed_seconds": elapsed,
                    "passed": maximum <= V4_STANDARD_ERRORS,
                }
            )
    return records


def validate_v5() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    device = _device()
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, 5]))
    for pde_id in FAMILIES:
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            adapter = SpectralTorchEquation(equation, device)
            x = rng.uniform(-1.0, 1.0, (32, d))
            first = rng.normal(size=(32, d + 1))
            second = rng.normal(size=(32, d + 1))
            xx = torch.as_tensor(x, dtype=torch.float64, device=device)
            one = torch.as_tensor(first, dtype=torch.float64, device=device)
            two = torch.as_tensor(second, dtype=torch.float64, device=device)
            with torch.no_grad():
                gpu_terminal = adapter.terminal(xx).cpu().numpy()
                gpu_gradient = adapter.terminal_gradient_z(xx).cpu().numpy()
                gpu_generator = adapter.generator(one).cpu().numpy()
                gpu_double = adapter.double_generator(one, two).cpu().numpy()
            errors = {
                "terminal": _relative_linf(gpu_terminal, equation.terminal(x)),
                "terminal_gradient": _relative_linf(
                    gpu_gradient, equation.terminal_gradient_z(x)
                ),
                "generator": _relative_linf(
                    gpu_generator,
                    equation.generator(first[:, 0], first[:, 1:]),
                ),
                "double_generator": _relative_linf(
                    gpu_double,
                    equation.double_generator(first[:, 1:], second[:, 1:]),
                ),
            }
            benchmark = torch.as_tensor(
                rng.normal(size=(1024, d)), dtype=torch.float64, device=device
            )
            with torch.no_grad():
                adapter.rotate(benchmark)
                if device.type == "cuda":
                    torch.cuda.synchronize(device)
                started = time.perf_counter()
                for _ in range(20):
                    adapter.rotate(benchmark)
                if device.type == "cuda":
                    torch.cuda.synchronize(device)
            milliseconds = 1000.0 * (time.perf_counter() - started) / 20.0
            maximum = max(errors.values())
            records.append(
                {
                    "stage": "V-5",
                    "pde_id": pde_id,
                    "dimension": d,
                    **{f"{key}_relative_error": value for key, value in errors.items()},
                    "max_relative_error": maximum,
                    "threshold": V5_TOLERANCE,
                    "dct_batch_points": 1024,
                    "dct_milliseconds_per_call": milliseconds,
                    "dct_implementation": "Makhoul FFT",
                    "device": str(device),
                    "passed": maximum <= V5_TOLERANCE,
                }
            )
    return records


def _numpy_tape_run(
    pde_id: str,
    n: int,
    M: int,
    method: str,
) -> tuple[np.ndarray, DrawTape, np.ndarray, np.ndarray]:
    d, repetition, chunk_index = 100, 0, 0
    equation = make_equation(pde_id, d)
    t, x, _ = fixed_points(pde_id, d)
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


def validate_v6() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    device = _device()
    for pde_id in ("S0", "S05"):
        equation = make_equation(pde_id, 100)
        for n, M in ((2, 8), (3, 4)):
            for method in ("path", "double_path_base"):
                expected, tape, t, x = _numpy_tape_run(pde_id, n, M, method)
                replay = tape.replay()
                seed_parts = [BASE_SEED, 100, n, M, 0, 0]
                solver = NonRidgeTorchMLP(
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
                        "stage": "V-6",
                        "pde_id": pde_id,
                        "dimension": 100,
                        "n": n,
                        "M": M,
                        "method": method,
                        "u_relative_linf_error": u_error,
                        "z_relative_linf_error": z_error,
                        "max_relative_error": maximum,
                        "threshold": V6_TOLERANCE,
                        "draw_entries": len(tape.entries),
                        "tape_exhausted": replay.exhausted,
                        "passed": maximum <= V6_TOLERANCE and replay.exhausted,
                    }
                )
    return records


def validate() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    environment = record_environment()
    records: list[dict[str, Any]] = []
    stages = (
        ("V-1", validate_v1),
        ("V-2", validate_v2),
        ("V-3", validate_v3),
        ("V-4", validate_v4),
        ("V-5", validate_v5),
        ("V-6", validate_v6),
    )
    for label, function in stages:
        started = time.perf_counter()
        stage_records = function()
        records.extend(stage_records)
        print(
            json.dumps(
                {
                    "stage": label,
                    "checks": len(stage_records),
                    "passed": sum(bool(row["passed"]) for row in stage_records),
                    "elapsed_seconds": time.perf_counter() - started,
                }
            ),
            flush=True,
        )
        if not all(bool(row["passed"]) for row in stage_records):
            break

    fields = sorted({key for row in records for key in row})
    with (RESULTS_ROOT / "validation.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    summary: dict[str, Any] = {
        "nonridge_revision": REVISION,
        "code_commit": current_commit(),
        "device": environment["device"],
        "all_passed": True,
    }
    for label, _function in stages:
        rows = [row for row in records if row["stage"] == label]
        verdict = "PASS" if rows and all(bool(row["passed"]) for row in rows) else (
            "FAIL" if rows else "NOT EVALUATED (EARLIER VALIDATION FAILURE)"
        )
        summary[label] = {
            "verdict": verdict,
            "checks": len(rows),
            "failures": sum(not bool(row["passed"]) for row in rows),
        }
        if verdict != "PASS":
            summary["all_passed"] = False
    write_json_atomic(RESULTS_ROOT / "validation_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    result = validate()
    raise SystemExit(0 if result["all_passed"] else 1)
