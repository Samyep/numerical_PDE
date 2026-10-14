"""Paths, deterministic banks, atomic output, and scalar metrics."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import random
from typing import Any, Sequence

import numpy as np
import torch


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "convex_concave_semilinear_pilot"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"


def ensure_dirs() -> None:
    for path in (RESULTS_DIR, CHECKPOINT_DIR, CACHE_DIR, PREDICTION_DIR, FIGURE_DIR):
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


def fixed_points(
    name: str, d: int, split: str, count: int, *, T: float, include_time: bool
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(stable_seed("tangent_actor_v1", name, split, count, d))
    t = (
        rng.uniform(0.0, np.nextafter(T, 0.0), count)
        if include_time else np.zeros(count)
    )
    return t.astype(np.float64), rng.standard_normal((count, d)).astype(np.float64)


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.shape).encode())
        digest.update(str(value.dtype).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as handle:
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


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("empty CSV")
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def model_metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    truth = np.asarray(truth, dtype=np.float64).reshape(-1)
    error = prediction - truth
    rmse = float(np.sqrt(np.mean(error**2)))
    std = float(np.std(truth))
    return {
        "rmse": rmse,
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "bias": float(np.mean(error)),
        "truth_mean": float(np.mean(truth)),
        "truth_std": std,
        "rmse_over_truth_std": rmse / std,
    }
