"""Common protocol, fixed data, metrics, and artifact helpers.

The final 960-point test split is never passed to a training or checkpoint
selection function.  A separate 240-point validation split is used for all
method-specific choices.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import tempfile
from typing import Any, Iterable

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
    make_driver,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    PeriodicReference2D,
    solve_periodic_reference_2d,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_strong_baselines"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
CALIBRATION_DIR = RESULTS_DIR / "calibration"

BUDGETS = (60.0, 120.0, 240.0, 480.0, 960.0)
TRAINING_SEEDS = (0, 1)

BENCHMARKS: dict[str, dict[str, Any]] = {
    "nonseparable_d100_k2": {"family": "nonseparable", "d": 100, "k": 2},
    "joint_concave_d400_k2": {"family": "joint_concave", "d": 400, "k": 2},
}

OFFICIAL_SOURCES = {
    "deep_bsde": {
        "url": "https://github.com/frankhan91/DeepBSDE",
        "commit": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
    },
    "sdgd_pinn": {
        "url": "https://github.com/zheyuanhu01/SDGD_PINN",
        "commit": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
    },
    "deep_picard": {
        "url": "https://github.com/DeepOptimalControl/DeepPicardIteration",
        "commit": "e8616c2f531c9f95299b246a13c1180b8f02822a",
    },
    "scasml": {
        "url": "https://github.com/Francis-Fan-create/SCaSML",
        "commit": "c2412b942b1720084accf9ac9651c6085e073919",
    },
    "full_history_mlp": {
        "url": "https://arxiv.org/abs/2310.12545",
        "implementation_url": "https://github.com/Francis-Fan-create/SCaSML/blob/main/solvers/MLP_full_history.py",
        "implementation_commit": "c2412b942b1720084accf9ac9651c6085e073919",
    },
}


def ensure_dirs() -> None:
    for directory in (
        RESULTS_DIR,
        CHECKPOINT_DIR,
        RAW_DIR,
        PREDICTION_DIR,
        FIGURE_DIR,
        SLURM_DIR,
        CALIBRATION_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("::".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "little") % (2**31 - 1)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def driver_for(benchmark: str) -> ConcaveSemilinearDriver:
    spec = BENCHMARKS[benchmark]
    return make_driver(spec["family"], d=spec["d"], k=spec["k"], T=0.5)


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_torch_save(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    torch.save(payload, temporary)
    os.replace(temporary, path)


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialized:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in materialized:
        for key in row:
            if key not in fields:
                fields.append(key)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(materialized)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def synchronize(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def state_fingerprint(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.dtype.str.encode())
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


def _reference_cache_path(benchmark: str) -> Path:
    return RESULTS_DIR / f"reference_{benchmark}.npz"


def load_or_build_reference(benchmark: str) -> PeriodicReference2D:
    """Load a compact independent 2D spectral reference or build it once."""
    ensure_dirs()
    path = _reference_cache_path(benchmark)
    if path.exists():
        archive = np.load(path)
        return PeriodicReference2D(
            name=str(archive["name"].item()),
            times=archive["times"],
            space=archive["space"],
            values=archive["values"],
            derivative_1=archive["derivative_1"],
            derivative_2=archive["derivative_2"],
            dt=float(archive["dt"].item()),
            store_dt=float(archive["store_dt"].item()),
            n=int(archive["n"].item()),
        )
    driver = driver_for(benchmark)
    reference = solve_periodic_reference_2d(driver, n=128, dt=0.00025, store_dt=0.002)
    np.savez_compressed(
        path,
        name=np.asarray(reference.name),
        times=reference.times,
        space=reference.space,
        values=reference.values,
        derivative_1=reference.derivative_1,
        derivative_2=reference.derivative_2,
        dt=np.asarray(reference.dt),
        store_dt=np.asarray(reference.store_dt),
        n=np.asarray(reference.n),
    )
    return reference


def evaluation_data(benchmark: str) -> dict[str, np.ndarray]:
    """Fixed 240 validation and 960 final-test initial states at t=0."""
    ensure_dirs()
    path = RESULTS_DIR / f"evaluation_points_{benchmark}.npz"
    if path.exists():
        with np.load(path) as archive:
            return {name: archive[name] for name in archive.files}
    driver = driver_for(benchmark)
    reference = load_or_build_reference(benchmark)
    rng = np.random.default_rng(stable_seed("strong_baseline_points", benchmark))
    x = rng.standard_normal((1200, driver.d)).astype(np.float64)
    t = np.zeros(1200, dtype=np.float64)
    working = reference.value(t, x @ driver.B)
    original = driver.working_to_original(working)
    validation = np.arange(240, dtype=np.int64)
    test = np.arange(240, 1200, dtype=np.int64)
    payload = {
        "t": t,
        "x": x,
        "reference_working": working,
        "reference_original": original,
        "validation_indices": validation,
        "test_indices": test,
    }
    np.savez_compressed(path, **payload)
    atomic_json(
        RESULTS_DIR / f"evaluation_points_{benchmark}.json",
        {
            "benchmark": benchmark,
            "point_count": 1200,
            "validation_count": 240,
            "test_count": 960,
            "distribution": "x~N(0,I_d), t=0",
            "fingerprint": state_fingerprint(t, x, validation, test),
            "reference": "independent 2D periodic Fourier/Strang solver, n=128, dt=0.00025",
            "test_used_for_training_or_selection": False,
        },
    )
    return payload


def metrics(prediction: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=np.float64).reshape(-1)
    reference = np.asarray(reference, dtype=np.float64).reshape(-1)
    error = prediction - reference
    rmse = float(np.sqrt(np.mean(error**2)))
    denominator = float(np.sqrt(np.mean(reference**2)))
    reference_std = float(np.std(reference, ddof=1))
    return {
        "rmse": rmse,
        "relative_l2": rmse / denominator,
        "rmse_over_reference_std": rmse / reference_std,
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
        "max_abs": float(np.max(np.abs(error))),
        "reference_rms": denominator,
        "reference_std": reference_std,
        "nonfinite_count": int(np.size(prediction) - np.isfinite(prediction).sum()),
    }


def closest_budget(actual_seconds: float) -> tuple[float, float]:
    budget = min(BUDGETS, key=lambda target: abs(target - actual_seconds))
    return budget, 100.0 * (actual_seconds - budget) / budget


def parameter_count(module: torch.nn.Module) -> int:
    return sum(parameter.numel() for parameter in module.parameters())


__all__ = [
    "BENCHMARKS", "BUDGETS", "CALIBRATION_DIR", "CHECKPOINT_DIR", "FIGURE_DIR",
    "OFFICIAL_SOURCES", "PACKAGE_DIR", "PREDICTION_DIR", "RAW_DIR", "RESULTS_DIR",
    "SLURM_DIR", "TRAINING_SEEDS", "atomic_json", "atomic_torch_save",
    "closest_budget", "driver_for", "ensure_dirs", "evaluation_data",
    "load_or_build_reference", "metrics", "parameter_count", "seed_everything",
    "stable_seed", "state_fingerprint", "synchronize", "write_csv",
]
