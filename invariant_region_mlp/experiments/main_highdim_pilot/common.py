"""Shared loading, hashing, and device helpers for the pilot."""

from __future__ import annotations

import platform
import subprocess
from typing import Any

import numpy as np
import torch

from .equations import (
    CoupledLQGEquation,
    DenseFourierDriver,
    LQGActorBenchmark,
)
from .protocol import (
    DENSE_PARAMETER_SEED,
    INPUT_ROOT,
    LQG_COEFFICIENT_SEED,
    TERMINAL_TIME,
)


def load_dense() -> DenseFourierDriver:
    with np.load(INPUT_ROOT / "dense_fourier_parameters.npz") as data:
        return DenseFourierDriver(
            data["omega"],
            data["signs"],
            data["phases"],
            T=TERMINAL_TIME,
            parameter_seed=DENSE_PARAMETER_SEED,
        )


def load_lqg() -> CoupledLQGEquation:
    with np.load(INPUT_ROOT / "coupled_lqg_parameters.npz") as data:
        return CoupledLQGEquation(
            data["c1"],
            data["c2"],
            T=TERMINAL_TIME,
            coefficient_seed=LQG_COEFFICIENT_SEED,
        )


def load_lqg_benchmark() -> LQGActorBenchmark:
    return LQGActorBenchmark(load_lqg())


def load_bank(pde: str) -> dict[str, np.ndarray]:
    stem = {
        "dense_fourier": "dense_fourier_banks.npz",
        "coupled_lqg": "coupled_lqg_banks.npz",
    }[pde]
    with np.load(INPUT_ROOT / stem) as data:
        return {key: data[key] for key in data.files}


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def environment_record() -> dict[str, Any]:
    gpu: dict[str, Any] | None = None
    if torch.cuda.is_available():
        properties = torch.cuda.get_device_properties(0)
        gpu = {
            "name": properties.name,
            "total_memory_bytes": int(properties.total_memory),
            "capability": list(torch.cuda.get_device_capability(0)),
            "cuda_runtime": torch.version.cuda,
        }
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "git_commit_at_recording": commit,
        "gpu": gpu,
    }


__all__ = [
    "environment_record",
    "load_bank",
    "load_dense",
    "load_lqg",
    "load_lqg_benchmark",
    "synchronize",
]

