"""Shared deterministic data, serialization, and runtime helpers."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
from pathlib import Path
import random
import time
from typing import Any, Iterator

import numpy as np
import torch


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "main_b1_sixmethods"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
FIGURES_DIR = RESULTS_DIR / "figures"
CONFIG_PATH = PACKAGE_DIR / "config.json"


def load_config() -> dict[str, Any]:
    with CONFIG_PATH.open() as handle:
        config = json.load(handle)
    dpi_profile_path = os.environ.get("B1_DPI_PROFILE_PATH")
    if dpi_profile_path:
        with Path(dpi_profile_path).open() as handle:
            profile = json.load(handle)
        config["dpi"].update(profile["resource_adapted_final_configuration"])
    return config


def ensure_result_dirs() -> None:
    for directory in (RESULTS_DIR, CHECKPOINT_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**32)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def fixed_test_points(d: int, k: int, count: int = 960) -> np.ndarray:
    rng = np.random.default_rng(stable_seed("B1", "heldout", d, k, count))
    points = rng.standard_normal((count, d), dtype=np.float64)
    points.setflags(write=False)
    return points


def fixed_reference_points(d: int, k: int, count: int = 1024) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(stable_seed("B1", "reference_calibration", d, k, count))
    t = rng.uniform(0.0, np.nextafter(0.25, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, d), dtype=np.float64)
    t.setflags(write=False)
    x.setflags(write=False)
    return t, x


def training_points(
    d: int,
    k: int,
    method: str,
    seed: int,
    count: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(stable_seed("B1", "train", method, d, k, seed, count))
    t = rng.uniform(0.0, np.nextafter(0.25, 0.0), count).astype(np.float64)
    x = rng.standard_normal((count, d), dtype=np.float64)
    return t, x


def array_sha256(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in arrays:
        data = np.ascontiguousarray(value)
        digest.update(str(data.shape).encode())
        digest.update(str(data.dtype).encode())
        digest.update(data.tobytes())
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


@contextlib.contextmanager
def measured_cuda(device: torch.device) -> Iterator[dict[str, Any]]:
    if device.type == "cuda":
        torch.cuda.synchronize(device)
        torch.cuda.reset_peak_memory_stats(device)
    start = time.perf_counter()
    record: dict[str, Any] = {}
    try:
        yield record
    finally:
        if device.type == "cuda":
            torch.cuda.synchronize(device)
            record["peak_gpu_bytes"] = int(torch.cuda.max_memory_allocated(device))
        else:
            record["peak_gpu_bytes"] = None
        record["wall_seconds"] = time.perf_counter() - start
