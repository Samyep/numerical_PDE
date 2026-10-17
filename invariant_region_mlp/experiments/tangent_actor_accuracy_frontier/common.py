"""Paths and provenance helpers for the accuracy-frontier follow-up."""

from __future__ import annotations

from pathlib import Path

from invariant_region_mlp.experiments.nonconvex_tangent_actor.common import (
    array_sha256,
    atomic_json,
    atomic_npz,
    atomic_torch_save,
    environment_payload,
    metrics,
    seed_everything,
    stable_seed,
    synchronize,
    write_csv,
)


PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parents[1]
RESULTS_DIR = PROJECT_DIR / "results" / "tangent_actor_accuracy_frontier"
RAW_DIR = RESULTS_DIR / "raw"
CACHE_DIR = RESULTS_DIR / "cache"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
PREDICTION_DIR = RESULTS_DIR / "predictions"
REFERENCE_DIR = RESULTS_DIR / "references"
FIGURE_DIR = RESULTS_DIR / "figures"
SLURM_LOG_DIR = RESULTS_DIR / "slurm"

OLD_RESULTS_DIR = PROJECT_DIR / "results" / "nonconvex_tangent_actor"
OLD_CHECKPOINT_DIR = OLD_RESULTS_DIR / "checkpoints"
OLD_PREDICTION_DIR = OLD_RESULTS_DIR / "predictions"
OLD_REFERENCE_DIR = OLD_RESULTS_DIR / "references"


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


def state_bank(
    benchmark: str,
    d: int,
    split: str,
    count: int,
    *,
    T: float,
    seed: int,
    include_terminal: bool = False,
):
    """Use the frozen v1 state-bank rule so comparisons share distributions."""
    from invariant_region_mlp.experiments.nonconvex_tangent_actor.common import (
        state_bank as frozen_state_bank,
    )

    return frozen_state_bank(
        benchmark,
        d,
        split,
        count,
        T=T,
        seed=seed,
        include_terminal=include_terminal,
    )


__all__ = [
    "CACHE_DIR",
    "CHECKPOINT_DIR",
    "FIGURE_DIR",
    "OLD_CHECKPOINT_DIR",
    "OLD_PREDICTION_DIR",
    "OLD_REFERENCE_DIR",
    "OLD_RESULTS_DIR",
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

