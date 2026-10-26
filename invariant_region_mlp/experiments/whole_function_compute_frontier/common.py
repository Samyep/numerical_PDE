"""Shared deterministic I/O, metrics, and stopping helpers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import random
from typing import Any

import numpy as np
import torch

from .protocol import (
    PLATEAU_MIN_SECONDS,
    PLATEAU_RELATIVE_IMPROVEMENT,
    PLATEAU_WINDOW,
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def save_npz_atomic(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def value_metrics(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    truth = np.asarray(truth, dtype=np.float64).reshape(-1)
    if prediction.shape != truth.shape or prediction.size == 0:
        raise ValueError((prediction.shape, truth.shape))
    if not np.isfinite(prediction).all() or not np.isfinite(truth).all():
        raise FloatingPointError("nonfinite prediction or truth")
    error = prediction - truth
    truth_l1 = float(np.abs(truth).sum())
    truth_l2 = float(np.linalg.norm(truth))
    if truth_l1 <= 0.0 or truth_l2 <= 0.0:
        raise ValueError("zero truth norm")
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relmae": float(np.abs(error).sum() / truth_l1),
        "relative_l2": float(np.linalg.norm(error) / truth_l2),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
    }


def plateau_reached(training_seconds: float, metrics: list[float]) -> bool:
    if training_seconds < PLATEAU_MIN_SECONDS or len(metrics) < PLATEAU_WINDOW:
        return False
    recent = np.asarray(metrics[-PLATEAU_WINDOW:], dtype=np.float64)
    if not np.isfinite(recent).all():
        return False
    baseline = max(abs(float(recent[0])), 1e-12)
    improvement = (float(recent[0]) - float(np.min(recent[1:]))) / baseline
    return improvement < PLATEAU_RELATIVE_IMPROVEMENT


def rng_state() -> dict[str, Any]:
    result: dict[str, Any] = {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch": torch.random.get_rng_state(),
    }
    if torch.cuda.is_available():
        result["cuda"] = torch.cuda.get_rng_state_all()
    return result


def restore_rng(state: dict[str, Any]) -> None:
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.random.set_rng_state(state["torch"])
    if torch.cuda.is_available() and "cuda" in state:
        torch.cuda.set_rng_state_all(state["cuda"])


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)

