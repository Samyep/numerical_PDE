"""Shared paths, deterministic seeds, and atomic artifact helpers."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import random
import tempfile
from typing import Any, Iterable

import numpy as np
import torch


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "original_pde_benchmark_comparison"
RAW_DIR = RESULTS_DIR / "raw"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
LABEL_DIR = RESULTS_DIR / "labels"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
OFFICIAL_WORK_DIR = RESULTS_DIR / "official_work"


def ensure_dirs() -> None:
    for path in (
        RESULTS_DIR,
        RAW_DIR,
        CHECKPOINT_DIR,
        LABEL_DIR,
        PREDICTION_DIR,
        FIGURE_DIR,
        SLURM_DIR,
        OFFICIAL_WORK_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**31 - 1)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, prefix=path.name + ".", suffix=".tmp", delete=False
    ) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_torch_save(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    os.replace(temporary, path)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = [
    "CHECKPOINT_DIR",
    "FIGURE_DIR",
    "LABEL_DIR",
    "OFFICIAL_WORK_DIR",
    "PACKAGE_DIR",
    "PREDICTION_DIR",
    "PROJECT_DIR",
    "RAW_DIR",
    "RESULTS_DIR",
    "SLURM_DIR",
    "atomic_json",
    "atomic_torch_save",
    "ensure_dirs",
    "seed_everything",
    "sha256_file",
    "stable_seed",
    "synchronize",
    "write_csv",
]
