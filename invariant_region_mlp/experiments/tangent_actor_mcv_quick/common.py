"""Frozen paths, checkpoints, and protocol for the MCV feasibility study."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.common import (
    atomic_json,
    driver_for,
    metrics,
    stable_seed,
    state_fingerprint,
    synchronize,
    write_csv,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_mcv_quick"
RAW_DIR = RESULTS_DIR / "raw"
PREDICTION_DIR = RESULTS_DIR / "predictions"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_DIR = RESULTS_DIR / "slurm"
STRONG_RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_strong_baselines"
HIGH_BUDGET_RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_high_budget"

BENCHMARKS = ("nonseparable_d100_k2", "joint_concave_d400_k2")
INFERENCE_PATHS = (64, 128, 512)
TIME_STEPS = 32
MC_REPETITIONS = 8
STATE_COUNT = 128
CALIBRATION_PATHS = 512


def ensure_dirs() -> None:
    for directory in (RESULTS_DIR, RAW_DIR, PREDICTION_DIR, FIGURE_DIR, SLURM_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def checkpoint_specs(benchmark: str, seed: int) -> list[dict[str, Any]]:
    if benchmark not in BENCHMARKS or seed not in (0, 1):
        raise KeyError((benchmark, seed))
    strong = STRONG_RESULTS_DIR / "checkpoints"
    high = HIGH_BUDGET_RESULTS_DIR / "checkpoints"
    if benchmark == "nonseparable_d100_k2":
        late_iteration = 6
    else:
        late_iteration = 7 if seed == 0 else 8
    specs = [
        {
            "checkpoint_key": "a1",
            "iteration": 1,
            "path": strong / f"tangent_actor_{benchmark}_s{seed}_a1.pt",
            "source_study": "tangent_actor_strong_baselines",
        },
        {
            "checkpoint_key": "a3",
            "iteration": 3,
            "path": strong / f"tangent_actor_{benchmark}_s{seed}_a3.pt",
            "source_study": "tangent_actor_strong_baselines",
        },
        {
            "checkpoint_key": "late",
            "iteration": late_iteration,
            "path": high / (
                f"tangent_actor_{benchmark}_s{seed}_a{late_iteration}_points8192.pt"
            ),
            "source_study": "tangent_actor_high_budget",
        },
    ]
    missing = [str(item["path"]) for item in specs if not item["path"].exists()]
    if missing:
        raise FileNotFoundError("missing Actor checkpoints: " + ", ".join(missing))
    return specs


def state_bank(benchmark: str) -> tuple[np.ndarray, np.ndarray, str]:
    """Load states and indices only; numerical reference stays unopened."""
    path = STRONG_RESULTS_DIR / f"evaluation_points_{benchmark}.npz"
    with np.load(path) as archive:
        indices = archive["validation_indices"][:STATE_COUNT].copy()
        x = archive["x"][indices].copy()
    return x, indices, state_fingerprint(x, indices)


def posthoc_reference(benchmark: str, indices: np.ndarray) -> np.ndarray:
    path = STRONG_RESULTS_DIR / f"evaluation_points_{benchmark}.npz"
    with np.load(path) as archive:
        return archive["reference_original"][indices].copy()


__all__ = [
    "BENCHMARKS", "CALIBRATION_PATHS", "FIGURE_DIR", "HIGH_BUDGET_RESULTS_DIR",
    "INFERENCE_PATHS", "MC_REPETITIONS", "PACKAGE_DIR", "PREDICTION_DIR",
    "PROJECT_DIR", "RAW_DIR", "RESULTS_DIR", "SLURM_DIR", "STATE_COUNT",
    "STRONG_RESULTS_DIR", "TIME_STEPS", "atomic_json", "checkpoint_specs",
    "driver_for", "ensure_dirs", "metrics", "posthoc_reference", "stable_seed",
    "state_bank", "synchronize", "write_csv",
]
