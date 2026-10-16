"""Frozen protocol, paths, and small utilities for MCV training."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.common import (
    atomic_json,
    atomic_torch_save,
    driver_for,
    metrics,
    stable_seed,
    state_fingerprint,
    synchronize,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_mcv_training"
RAW_DIR = RESULTS_DIR / "raw"
LABEL_DIR = RESULTS_DIR / "labels"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
STRONG_RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_strong_baselines"
HIGH_BUDGET_RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_high_budget"

BENCHMARKS = ("nonseparable_d100_k2", "joint_concave_d400_k2")
SEEDS = (0, 1)


def ensure_dirs() -> None:
    for directory in (
        RESULTS_DIR, RAW_DIR, LABEL_DIR, CHECKPOINT_DIR, PREDICTION_DIR,
        FIGURE_DIR, SLURM_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def a3_checkpoint(benchmark: str, seed: int) -> Path:
    path = (
        STRONG_RESULTS_DIR / "checkpoints"
        / f"tangent_actor_{benchmark}_s{seed}_a3.pt"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def late_checkpoint(benchmark: str, seed: int) -> tuple[Path, int]:
    iteration = 6 if benchmark == "nonseparable_d100_k2" else (7 if seed == 0 else 8)
    path = (
        HIGH_BUDGET_RESULTS_DIR / "checkpoints"
        / f"tangent_actor_{benchmark}_s{seed}_a{iteration}_points8192.pt"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    return path, iteration


def validation_bank(
    benchmark: str, seed: int, *, count: int = 128
) -> tuple[np.ndarray, np.ndarray, str]:
    """Reuse the fixed held-out high-budget bank without opening Ground Truth."""
    path = (
        HIGH_BUDGET_RESULTS_DIR / "cache"
        / f"nested_points_{benchmark}_s{seed}_e1.npz"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path) as archive:
        t = archive["valid_t"][:count].copy()
        x = archive["valid_x"][:count].copy()
    return t, x, state_fingerprint(t, x)


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode())
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows for {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


__all__ = [
    "BENCHMARKS", "CHECKPOINT_DIR", "FIGURE_DIR", "HIGH_BUDGET_RESULTS_DIR",
    "LABEL_DIR", "PACKAGE_DIR", "PREDICTION_DIR", "PROJECT_DIR", "RAW_DIR",
    "RESULTS_DIR", "SEEDS", "SLURM_DIR", "STRONG_RESULTS_DIR", "a3_checkpoint",
    "array_sha256", "atomic_json", "atomic_torch_save", "driver_for", "ensure_dirs",
    "late_checkpoint", "metrics", "read_json", "stable_seed", "state_fingerprint",
    "synchronize", "validation_bank", "write_csv",
]
