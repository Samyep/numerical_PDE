"""Paths, deterministic nested sample banks, and atomic artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.common import (
    array_sha256,
    atomic_json,
    atomic_npz,
    atomic_torch_save,
    seed_everything,
    stable_seed,
    write_csv,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_scaling_diagnosis"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"


def ensure_dirs() -> None:
    for path in (
        RESULTS_DIR,
        CACHE_DIR,
        CHECKPOINT_DIR,
        PREDICTION_DIR,
        FIGURE_DIR,
        SLURM_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def nested_points(
    name: str,
    d: int,
    split: str,
    count: int,
    *,
    T: float,
    include_time: bool,
    training_seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a count-independent stream, so smaller banks are exact prefixes."""
    base_seed = stable_seed("tangent_actor_scaling_v1", name, d, split, training_seed)
    time_rng = np.random.default_rng(stable_seed(base_seed, "time"))
    state_rng = np.random.default_rng(stable_seed(base_seed, "state"))
    t = (
        time_rng.uniform(0.0, np.nextafter(T, 0.0), count)
        if include_time
        else np.zeros(count)
    )
    x = state_rng.standard_normal((count, d))
    return t.astype(np.float64), x.astype(np.float64)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


__all__: Sequence[str] = (
    "RESULTS_DIR",
    "CACHE_DIR",
    "CHECKPOINT_DIR",
    "PREDICTION_DIR",
    "FIGURE_DIR",
    "SLURM_DIR",
    "ensure_dirs",
    "nested_points",
    "array_sha256",
    "atomic_json",
    "atomic_npz",
    "atomic_torch_save",
    "seed_everything",
    "stable_seed",
    "write_csv",
    "json_safe",
)

