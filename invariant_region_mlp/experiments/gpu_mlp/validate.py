"""G-1 equation checks and G-2 draw-tape exactness validation."""

from __future__ import annotations

import csv
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import (
    CorrectionBaselineMLP,
    DrawTape,
    METHODS,
)
from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    fixed_points,
    make_equation,
)
from invariant_region_mlp.experiments.double_estimator.protocol import chunk_size

from .protocol import RESULTS_ROOT, REVISION, write_json_atomic
from .solver import GPU_METHODS, TapeTorchDraws, TorchEquation, TorchMLP


G1_TOLERANCE = 1e-12
G2_TOLERANCE = 1e-10
G2_CELLS = ((2, 8), (3, 4), (3, 6))


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _current_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _relative_linf(observed: np.ndarray, expected: np.ndarray) -> float:
    numerator = float(np.max(np.abs(observed - expected)))
    denominator = max(float(np.max(np.abs(expected))), np.finfo(np.float64).tiny)
    return numerator / denominator


def _device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def record_environment() -> dict[str, Any]:
    device = _device()
    payload: dict[str, Any] = {
        "gpu_mlp_revision": REVISION,
        "code_commit": _current_commit(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "device": str(device),
        "dtype": "float64",
    }
    if torch.cuda.is_available():
        index = torch.cuda.current_device()
        properties = torch.cuda.get_device_properties(index)
        payload.update(
            {
                "gpu_name": properties.name,
                "gpu_total_memory_bytes": properties.total_memory,
                "gpu_compute_capability": f"{properties.major}.{properties.minor}",
                "fp64_throughput_class": "consumer GPU; reduced-rate FP64",
            }
        )
        try:
            smi = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=driver_version",
                    "--format=csv,noheader",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            payload["nvidia_driver"] = smi.stdout.strip().splitlines()[0]
        except (OSError, subprocess.CalledProcessError, IndexError):
            payload["nvidia_driver"] = None
    write_json_atomic(RESULTS_ROOT / "environment.json", payload)
    return payload


def run_g1() -> list[dict[str, Any]]:
    device = _device()
    rng = np.random.default_rng(20261201)
    records: list[dict[str, Any]] = []
    for pde_id in ("P1", "MR"):
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            adapter = TorchEquation(equation, device)
            t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), 32).astype(np.float64)
            x = rng.uniform(-1.0, 1.0, (32, d)).astype(np.float64)
            state = rng.normal(0.0, 1.0, (32, d + 1)).astype(np.float64)
            tt = torch.as_tensor(t, dtype=torch.float64, device=device)
            xx = torch.as_tensor(x, dtype=torch.float64, device=device)
            ss = torch.as_tensor(state, dtype=torch.float64, device=device)
            with torch.no_grad():
                torch_terminal = adapter.terminal(xx).cpu().numpy()
                torch_exact = adapter.exact_state(tt, xx).cpu().numpy()
                torch_generator = adapter.generator(ss).cpu().numpy()
            numpy_terminal = equation.terminal(x)
            numpy_exact = equation.exact_state(t, x)
            numpy_generator = equation.generator(state[:, 0], state[:, 1:])
            for quantity, observed, expected in (
                ("terminal", torch_terminal, numpy_terminal),
                ("exact_state", torch_exact, numpy_exact),
                ("generator", torch_generator, numpy_generator),
            ):
                error = _relative_linf(observed, expected)
                records.append(
                    {
                        "stage": "G-1",
                        "pde_id": pde_id,
                        "dimension": d,
                        "quantity": quantity,
                        "relative_linf_error": error,
                        "threshold": G1_TOLERANCE,
                        "passed": error <= G1_TOLERANCE,
                        "device": str(device),
                    }
                )
    return records


def _numpy_tape_run(
    pde_id: str, n: int, M: int, method: str
) -> tuple[np.ndarray, DrawTape, np.ndarray, np.ndarray]:
    d, repetition, chunk_index = 100, 0, 0
    equation = make_equation(pde_id, d)
    t, x, _ = fixed_points(pde_id, d)
    size = chunk_size(d)
    t, x = t[:size], x[:size]
    seed_parts = [BASE_SEED, d, n, M, repetition, chunk_index]
    tape = DrawTape("record")
    solver = CorrectionBaselineMLP(
        equation,
        M,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        draw_tape=tape,
    )
    return solver.solve(n, t, x), tape, t, x


def run_g2() -> list[dict[str, Any]]:
    device = _device()
    records: list[dict[str, Any]] = []
    for pde_id in ("P1", "MR"):
        equation = make_equation(pde_id, 100)
        for n, M in G2_CELLS:
            for method in GPU_METHODS:
                numpy_output, record_tape, t, x = _numpy_tape_run(
                    pde_id, n, M, method
                )
                replay = record_tape.replay()
                seed_parts = [BASE_SEED, 100, n, M, 0, 0]
                source = TapeTorchDraws(device, replay)
                solver = TorchMLP(
                    equation,
                    M,
                    method,
                    device,
                    seed_sequence=np.random.SeedSequence(seed_parts),
                    aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
                    draw_source=source,
                )
                with torch.no_grad():
                    torch_output = solver.solve(
                        n,
                        torch.tensor(t, dtype=torch.float64, device=device),
                        torch.tensor(x, dtype=torch.float64, device=device),
                    ).cpu().numpy()
                u_error = _relative_linf(torch_output[:, 0], numpy_output[:, 0])
                z_error = _relative_linf(torch_output[:, 1:], numpy_output[:, 1:])
                maximum = max(u_error, z_error)
                records.append(
                    {
                        "stage": "G-2",
                        "pde_id": pde_id,
                        "dimension": 100,
                        "n": n,
                        "M": M,
                        "method": method,
                        "u_relative_linf_error": u_error,
                        "z_relative_linf_error": z_error,
                        "max_relative_error": maximum,
                        "threshold": G2_TOLERANCE,
                        "draw_entries": len(record_tape.entries),
                        "tape_exhausted": replay.exhausted,
                        "passed": maximum <= G2_TOLERANCE and replay.exhausted,
                        "device": str(device),
                    }
                )
    return records


def validate() -> dict[str, Any]:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    environment = record_environment()
    records = run_g1() + run_g2()
    fields = sorted({key for row in records for key in row})
    with (RESULTS_ROOT / "validation.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    g1 = [row for row in records if row["stage"] == "G-1"]
    g2 = [row for row in records if row["stage"] == "G-2"]
    payload = {
        "gpu_mlp_revision": REVISION,
        "code_commit": _current_commit(),
        "device": environment["device"],
        "G-1": {
            "verdict": "PASS" if g1 and all(row["passed"] for row in g1) else "FAIL",
            "checks": len(g1),
            "failures": sum(not row["passed"] for row in g1),
            "max_relative_error": max(row["relative_linf_error"] for row in g1),
            "threshold": G1_TOLERANCE,
        },
        "G-2": {
            "verdict": "PASS" if len(g2) == 36 and all(row["passed"] for row in g2) else "FAIL",
            "checks": len(g2),
            "expected_checks": 36,
            "failures": sum(not row["passed"] for row in g2),
            "max_relative_error": max(row["max_relative_error"] for row in g2),
            "threshold": G2_TOLERANCE,
        },
    }
    write_json_atomic(RESULTS_ROOT / "validation_summary.json", payload)
    print(json.dumps(payload, indent=2), flush=True)
    return payload


if __name__ == "__main__":
    validate()
