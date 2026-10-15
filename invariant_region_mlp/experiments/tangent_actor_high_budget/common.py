"""Paths and frozen protocol shared by the high-budget study."""

from __future__ import annotations

from pathlib import Path

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.common import (
    BENCHMARKS,
    BUDGETS,
    atomic_json,
    atomic_torch_save,
    driver_for,
    ensure_dirs as ensure_strong_dirs,
    evaluation_data,
    metrics,
    parameter_count,
    seed_everything,
    stable_seed,
    state_fingerprint,
    synchronize,
    write_csv,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_high_budget"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
CACHE_DIR = RESULTS_DIR / "cache"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
STRONG_RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_strong_baselines"


def ensure_dirs() -> None:
    ensure_strong_dirs()
    for directory in (
        RESULTS_DIR,
        CHECKPOINT_DIR,
        RAW_DIR,
        PREDICTION_DIR,
        CACHE_DIR,
        FIGURE_DIR,
        SLURM_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


__all__ = [
    "BENCHMARKS", "BUDGETS", "CACHE_DIR", "CHECKPOINT_DIR", "FIGURE_DIR",
    "PACKAGE_DIR", "PREDICTION_DIR", "PROJECT_DIR", "RAW_DIR", "RESULTS_DIR",
    "SLURM_DIR", "STRONG_RESULTS_DIR", "atomic_json", "atomic_torch_save",
    "driver_for", "ensure_dirs", "evaluation_data", "metrics", "parameter_count",
    "seed_everything", "stable_seed", "state_fingerprint", "synchronize", "write_csv",
]

