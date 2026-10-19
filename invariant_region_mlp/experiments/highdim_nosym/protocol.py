"""Frozen constants, task grids, and artifact routing for the study."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Any, Final

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    PACKAGE_ROOT,
    POINT_COUNT,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor import (
    MATURE,
)


REVISION: Final = "highdim_nosym_v1"
DIMENSIONS: Final = (100, 400)
PDES: Final = ("S05", "B6")
TRAINING_SEEDS: Final = (0,)
OPTIONAL_TRAINING_SEEDS: Final = (1, 2)
TA_UPDATES: Final = 4
TA_PRIMARY_STEPS: Final = 64
TA_AUDIT_STEPS: Final = (128, 256)
TA_PATHS: Final = 256
TA_REPLICATES: Final = 4
FH_REPETITIONS: Final = 3
FH_RAW_PATH_GRID: Final = (
    (2, 16), (2, 64), (2, 256),
    (3, 8), (3, 16), (3, 24),
    (4, 6), (4, 10),
)
FH_REPAIRED_GRID: Final = (
    (3, 16), (3, 24), (3, 32),
    (4, 6), (4, 10), (4, 16),
)
B6_BETAS: Final = (1.0, 1.5, 2.0, 3.0, 4.0)
B6_GAMMAS: Final = (0.0, 1.0)
B6_HEAT_SAMPLES: Final = 4096

RESULTS_ROOT = PACKAGE_ROOT / "results" / "highdim_nosym"
RAW_ROOT = RESULTS_ROOT / "raw"
CHECKPOINT_ROOT = RESULTS_ROOT / "checkpoints"
CACHE_ROOT = RESULTS_ROOT / "cache"
PREDICTION_ROOT = RESULTS_ROOT / "predictions"
FIGURES_ROOT = RESULTS_ROOT / "figures"
REPORT_PATH = PACKAGE_ROOT / "docs" / "HIGHDIM_NOSYM_REPORT.md"
CALIBRATION_FREEZE_PATH = RESULTS_ROOT / "B6_CALIBRATION_FREEZE.json"


def ensure_directories() -> None:
    for path in (
        RESULTS_ROOT,
        RAW_ROOT,
        CHECKPOINT_ROOT,
        CACHE_ROOT,
        PREDICTION_ROOT,
        FIGURES_ROOT,
    ):
        path.mkdir(parents=True, exist_ok=True)


def actor_configuration() -> dict[str, Any]:
    """Return the imported final-suite configuration without mutating it."""

    return asdict(MATURE)


def fh_tasks() -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for pde_id in PDES:
        for dimension in DIMENSIONS:
            for method, grid in (
                ("raw", FH_RAW_PATH_GRID),
                ("path", FH_RAW_PATH_GRID),
                ("double_path_base", FH_REPAIRED_GRID),
            ):
                for n, M in grid:
                    for repetition in range(FH_REPETITIONS):
                        tasks.append(
                            {
                                "pde_id": pde_id,
                                "d": dimension,
                                "n": n,
                                "M": M,
                                "method": method,
                                "repetition": repetition,
                            }
                        )
    return tasks


def fh_task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / "fh"
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / str(task["method"])
        / f"rep{int(task['repetition']):02d}.json"
    )


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


__all__ = [
    "BASE_SEED",
    "B6_BETAS",
    "B6_GAMMAS",
    "B6_HEAT_SAMPLES",
    "CACHE_ROOT",
    "CALIBRATION_FREEZE_PATH",
    "CHECKPOINT_ROOT",
    "DIMENSIONS",
    "FH_RAW_PATH_GRID",
    "FH_REPAIRED_GRID",
    "FH_REPETITIONS",
    "FIGURES_ROOT",
    "OPTIONAL_TRAINING_SEEDS",
    "PDES",
    "POINT_COUNT",
    "PREDICTION_ROOT",
    "RAW_ROOT",
    "REPORT_PATH",
    "RESULTS_ROOT",
    "REVISION",
    "TA_AUDIT_STEPS",
    "TA_PATHS",
    "TA_PRIMARY_STEPS",
    "TA_REPLICATES",
    "TA_UPDATES",
    "TRAINING_SEEDS",
    "actor_configuration",
    "ensure_directories",
    "fh_task_path",
    "fh_tasks",
    "write_json_atomic",
]

