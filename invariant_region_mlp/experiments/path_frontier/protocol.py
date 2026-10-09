"""Frozen task matrix and artifact routing for the path-frontier study."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.double_estimator.equations import PACKAGE_ROOT
from invariant_region_mlp.experiments.double_estimator.protocol import (
    EC_CELLS,
    REPETITIONS,
    chunk_size,
    make_task,
    task_is_complete as double_task_is_complete,
    task_path as double_task_path,
    write_json_atomic,
)
from invariant_region_mlp.experiments.double_estimator.solver import (
    IMPLEMENTATION_REVISION as DOUBLE_IMPLEMENTATION_REVISION,
)


RESULTS_ROOT = PACKAGE_ROOT / "results" / "path_frontier"
RAW_ROOT = RESULTS_ROOT / "raw"
DOUBLE_RESULTS_ROOT = PACKAGE_ROOT / "results" / "double_estimator"
PATH_FRONTIER_REVISION = "path_frontier_v1"
PROBLEMS = (("P1", 100), ("P1", 400), ("MR", 100), ("MR", 400))
METHODS = ("raw", "path", "double", "double_path")
EXPECTED_PATH_ROWS = len(PROBLEMS) * len(EC_CELLS) * REPETITIONS
EXPECTED_COMBINED_ROWS = EXPECTED_PATH_ROWS * len(METHODS)


def path_tasks() -> list[dict[str, Any]]:
    """Return the 480 registered path tasks in a deterministic heavy-first order."""

    tasks = [
        make_task("path_frontier", pde_id, d, n, M, "path", repetition)
        for pde_id, d in PROBLEMS
        for n, M in EC_CELLS
        for repetition in range(REPETITIONS)
    ]
    return sorted(
        tasks,
        key=lambda task: (
            -(int(task["d"]) * int(task["M"]) ** int(task["n"])),
            str(task["pde_id"]),
            -int(task["n"]),
            -int(task["M"]),
            int(task["repetition"]),
        ),
    )


def identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(task["pde_id"]),
        int(task["d"]),
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["repetition"]),
    )


def new_task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / f"rep{int(task['repetition']):02d}.json"
    )


def existing_double_path_task(task: dict[str, Any]) -> dict[str, Any]:
    candidate = task.copy()
    candidate["stage"] = "core"
    candidate["method"] = "path"
    return candidate


def has_reusable_path(task: dict[str, Any]) -> bool:
    return double_task_is_complete(existing_double_path_task(task))


def new_task_is_complete(task: dict[str, Any]) -> bool:
    path = new_task_path(task)
    if not path.is_file():
        return False
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (
        row.get("path_frontier_revision") == PATH_FRONTIER_REVISION
        and row.get("implementation_revision") == DOUBLE_IMPLEMENTATION_REVISION
        and row.get("pde_id") == task["pde_id"]
        and int(row.get("dimension", -1)) == int(task["d"])
        and int(row.get("n", -1)) == int(task["n"])
        and int(row.get("M", -1)) == int(task["M"])
        and row.get("method") == "path"
        and int(row.get("repetition", -1)) == int(task["repetition"])
        and int(row.get("chunk_size", -1)) == chunk_size(int(task["d"]))
    )


def source_path(task: dict[str, Any], method: str) -> tuple[Path, str]:
    candidate = task.copy()
    candidate["method"] = method
    if method != "path" or has_reusable_path(candidate):
        return double_task_path(candidate), "double_estimator"
    return new_task_path(candidate), "path_frontier"


__all__ = [
    "DOUBLE_RESULTS_ROOT",
    "EC_CELLS",
    "EXPECTED_COMBINED_ROWS",
    "EXPECTED_PATH_ROWS",
    "METHODS",
    "PATH_FRONTIER_REVISION",
    "PROBLEMS",
    "RAW_ROOT",
    "REPETITIONS",
    "RESULTS_ROOT",
    "has_reusable_path",
    "identity",
    "new_task_is_complete",
    "new_task_path",
    "path_tasks",
    "source_path",
    "write_json_atomic",
]
