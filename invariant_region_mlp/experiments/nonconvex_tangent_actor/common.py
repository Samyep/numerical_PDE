"""Paths, deterministic banks, and atomic output helpers."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import time
from typing import Any, Iterable

import numpy as np
import torch


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "nonconvex_tangent_actor"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
REFERENCE_DIR = RESULTS_DIR / "references"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_LOG_DIR = RESULTS_DIR / "slurm"


def ensure_dirs() -> None:
    for path in (
        RESULTS_DIR,
        RAW_DIR,
        CACHE_DIR,
        CHECKPOINT_DIR,
        PREDICTION_DIR,
        REFERENCE_DIR,
        FIGURE_DIR,
        SLURM_LOG_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "little") % (2**32)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def state_bank(
    benchmark: str,
    d: int,
    split: str,
    count: int,
    *,
    T: float,
    seed: int,
    include_terminal: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a full-ambient bank; active coordinates are never constructed."""
    rng = np.random.default_rng(
        stable_seed("nonconvex_tangent_actor_v1", benchmark, split, count, d, seed)
    )
    high = T if include_terminal else np.nextafter(T, 0.0)
    t = rng.uniform(0.0, high, count).astype(np.float64)
    x = rng.standard_normal((count, d)).astype(np.float64)
    return t, x


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def atomic_torch_save(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    torch.save(payload, temporary)
    temporary.replace(path)


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    materialized = list(rows)
    if not materialized:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in materialized:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(materialized)
    temporary.replace(path)


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def environment_payload(stage: str) -> dict[str, Any]:
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        commit = "unknown"
    return {
        "stage": stage,
        "host": platform.node(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "git_commit_at_start": commit,
        "started_unix": time.time(),
    }


def metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    predicted = np.asarray(prediction, dtype=np.float64)
    exact = np.asarray(truth, dtype=np.float64)
    error = predicted - exact
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(
            np.linalg.norm(error) / max(np.linalg.norm(exact), np.finfo(float).tiny)
        ),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "max_abs": float(np.max(np.abs(error))),
    }


__all__ = [
    "CACHE_DIR",
    "CHECKPOINT_DIR",
    "FIGURE_DIR",
    "PACKAGE_DIR",
    "PREDICTION_DIR",
    "PROJECT_DIR",
    "RAW_DIR",
    "REFERENCE_DIR",
    "RESULTS_DIR",
    "SLURM_LOG_DIR",
    "array_sha256",
    "atomic_json",
    "atomic_npz",
    "atomic_torch_save",
    "ensure_dirs",
    "environment_payload",
    "metrics",
    "seed_everything",
    "stable_seed",
    "state_bank",
    "synchronize",
    "write_csv",
]
