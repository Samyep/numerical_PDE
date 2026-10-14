"""Deterministic data banks, atomic I/O, timing, and metric helpers."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import time
from typing import Any, Iterator, Sequence

import numpy as np
import torch


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "hjb_compute_frontier"
DOCS_DIR = PROJECT_DIR / "docs" / "hjb_compute_frontier"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
CACHE_DIR = RESULTS_DIR / "cache"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
PROTOCOL_PATH = PACKAGE_DIR / "protocol_frozen.json"


def ensure_dirs() -> None:
    for path in (
        RESULTS_DIR,
        DOCS_DIR,
        CHECKPOINT_DIR,
        CACHE_DIR,
        PREDICTION_DIR,
        FIGURE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**32)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.shape).encode())
        digest.update(str(value.dtype).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


def atomic_torch_save(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    torch.save(payload, temporary)
    temporary.replace(path)


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def load_protocol() -> dict[str, Any]:
    with PROTOCOL_PATH.open() as handle:
        return json.load(handle)


def fixed_points(
    equation_name: str,
    d: int,
    split: str,
    count: int,
    *,
    include_time: bool,
    T: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a sealed deterministic bank, disjoint by split namespace."""

    if split not in {"train", "label_valid", "policy_valid", "test", "audit"}:
        raise ValueError(f"unknown split {split}")
    rng = np.random.default_rng(
        stable_seed("hjb_compute_frontier_v1", equation_name, d, split, count)
    )
    if include_time:
        t = rng.uniform(0.0, np.nextafter(T, 0.0), count).astype(np.float64)
    else:
        t = np.zeros(count, dtype=np.float64)
    x = rng.standard_normal((count, d), dtype=np.float64)
    t.setflags(write=False)
    x.setflags(write=False)
    return t, x


def named_points(
    equation_name: str,
    d: int,
    namespace: str,
    count: int,
    *,
    T: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic training/validation points for a named algorithm phase."""

    rng = np.random.default_rng(
        stable_seed("hjb_compute_frontier_v1", equation_name, namespace, d, count)
    )
    t = rng.uniform(0.0, np.nextafter(T, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, d), dtype=np.float64)
    return t, x


def model_metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    truth = np.asarray(truth, dtype=np.float64).reshape(-1)
    error = prediction - truth
    rmse = float(np.sqrt(np.mean(error**2)))
    relative_l2 = float(
        np.linalg.norm(error) / max(np.linalg.norm(truth), np.finfo(np.float64).tiny)
    )
    truth_std = float(np.std(truth, ddof=0))
    return {
        "relative_l2": relative_l2,
        "rmse": rmse,
        "bias": float(np.mean(error)),
        "rmse_over_truth_std": rmse / max(truth_std, np.finfo(np.float64).tiny),
        "truth_mean": float(np.mean(truth)),
        "truth_std": truth_std,
    }


def vector_relative_error(prediction: np.ndarray, truth: np.ndarray) -> float:
    prediction = np.asarray(prediction, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    return float(
        np.linalg.norm(prediction - truth)
        / max(np.linalg.norm(truth), np.finfo(np.float64).tiny)
    )


@contextlib.contextmanager
def measured_cuda(device: torch.device) -> Iterator[dict[str, Any]]:
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    record: dict[str, Any] = {}
    try:
        yield record
    finally:
        if device.type == "cuda":
            torch.cuda.synchronize(device)
            record["peak_gpu_bytes"] = int(torch.cuda.max_memory_allocated(device))
        else:
            record["peak_gpu_bytes"] = None
        record["wall_seconds"] = time.perf_counter() - started


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty CSV")
    columns = list(fields or rows[0].keys())
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def require_cuda() -> torch.device:
    if not torch.cuda.is_available():
        raise RuntimeError("this experiment requires an allocated CUDA device")
    return torch.device("cuda")
