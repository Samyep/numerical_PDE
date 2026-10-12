"""Frozen task matrix and resumable artifact routing for the non-ridge study."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

from invariant_region_mlp.experiments.double_estimator.equations import PACKAGE_ROOT

from .solver import IMPLEMENTATION_REVISION, gpu_chunk_size


RESULTS_ROOT = PACKAGE_ROOT / "results" / "nonridge"
RAW_ROOT = RESULTS_ROOT / "raw"
REVISION = "nonridge_v1"
REPETITIONS = 10
DIMENSIONS = (100, 400)
STAGE_A_FAMILIES = ("S0", "S05")
STAGE_B_FAMILIES = ("S1",)

FULL_GRID: tuple[tuple[int, int, tuple[str, ...]], ...] = (
    (2, 16, ("path",)),
    (2, 64, ("path", "raw")),
    (2, 256, ("path",)),
    (3, 16, ("raw", "path", "path_base", "double_path", "double_path_base")),
    (3, 24, ("double_path", "double_path_base")),
)


def make_task(
    stage: str,
    pde_id: str,
    d: int,
    n: int,
    M: int,
    method: str,
    repetition: int,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "pde_id": pde_id,
        "d": int(d),
        "n": int(n),
        "M": int(M),
        "method": method,
        "repetition": int(repetition),
        "chunk_size": gpu_chunk_size(d, n, M),
    }


def _full_grid_tasks(stage: str, families: Iterable[str]) -> list[dict[str, Any]]:
    return [
        make_task(stage, family, d, n, M, method, repetition)
        for family in families
        for d in DIMENSIONS
        for n, M, methods in FULL_GRID
        for method in methods
        for repetition in range(REPETITIONS)
    ]


def stage_a_tasks() -> list[dict[str, Any]]:
    return sorted(_full_grid_tasks("stage_a", STAGE_A_FAMILIES), key=task_sort_key)


def stage_b_s1_tasks() -> list[dict[str, Any]]:
    return sorted(_full_grid_tasks("stage_b_s1", STAGE_B_FAMILIES), key=task_sort_key)


def stage_b_extension_tasks() -> list[dict[str, Any]]:
    tasks = [
        make_task("stage_b_extension", family, d, 3, 32, method, repetition)
        for family in ("S0", "S05", "S1")
        for d in DIMENSIONS
        for method in ("double_path", "double_path_base")
        for repetition in range(REPETITIONS)
    ]
    tasks.extend(
        make_task("stage_b_extension", family, d, 2, 512, "path", repetition)
        for family in ("S0", "S05", "S1")
        for d in DIMENSIONS
        for repetition in range(REPETITIONS)
    )
    return sorted(tasks, key=task_sort_key)


def stage_b_n4_tasks() -> list[dict[str, Any]]:
    return sorted(
        [
            make_task(
                "stage_b_n4",
                family,
                400,
                4,
                M,
                "double_path_base",
                repetition,
            )
            for family in ("S0", "S05", "S1")
            for M in (6, 10)
            for repetition in range(REPETITIONS)
        ],
        key=task_sort_key,
    )


def stage_b_tasks() -> list[dict[str, Any]]:
    return stage_b_s1_tasks() + stage_b_extension_tasks() + stage_b_n4_tasks()


def all_tasks() -> list[dict[str, Any]]:
    return stage_a_tasks() + stage_b_tasks()


def identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(task["pde_id"]),
        int(task["d"]),
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["repetition"]),
    )


def task_sort_key(task: dict[str, Any]) -> tuple[Any, ...]:
    method_order = {
        "raw": 0,
        "path": 1,
        "path_base": 2,
        "double_path": 3,
        "double_path_base": 4,
    }
    estimated_work = int(task["M"]) ** int(task["n"])
    if str(task["method"]).startswith("double"):
        estimated_work *= 2
    return (
        -int(task["d"]),
        estimated_work,
        str(task["pde_id"]),
        int(task["n"]),
        int(task["M"]),
        method_order[str(task["method"])],
        int(task["repetition"]),
    )


def task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / str(task["stage"])
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / str(task["method"])
        / f"rep{int(task['repetition']):02d}.json"
    )


def task_is_complete(task: dict[str, Any]) -> bool:
    path = task_path(task)
    if not path.is_file():
        return False
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (
        row.get("nonridge_revision") == REVISION
        and row.get("implementation_revision") == IMPLEMENTATION_REVISION
        and row.get("pde_id") == task["pde_id"]
        and int(row.get("dimension", -1)) == int(task["d"])
        and int(row.get("n", -1)) == int(task["n"])
        and int(row.get("M", -1)) == int(task["M"])
        and row.get("method") == task["method"]
        and int(row.get("repetition", -1)) == int(task["repetition"])
        and int(row.get("chunk_size", -1)) == int(task["chunk_size"])
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
    "DIMENSIONS",
    "FULL_GRID",
    "RAW_ROOT",
    "REPETITIONS",
    "RESULTS_ROOT",
    "REVISION",
    "all_tasks",
    "identity",
    "stage_a_tasks",
    "stage_b_extension_tasks",
    "stage_b_n4_tasks",
    "stage_b_s1_tasks",
    "stage_b_tasks",
    "task_is_complete",
    "task_path",
    "write_json_atomic",
]
