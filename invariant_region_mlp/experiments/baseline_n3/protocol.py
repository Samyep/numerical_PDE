"""Frozen Study-F task matrix, reuse rules, and artifact routing."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.double_estimator.equations import PACKAGE_ROOT
from invariant_region_mlp.experiments.double_estimator.protocol import (
    chunk_size,
    make_task,
    task_is_complete as double_task_is_complete,
    task_path as double_task_path,
    write_json_atomic,
)
from invariant_region_mlp.experiments.path_frontier.protocol import (
    source_path as frontier_source_path,
)

from .solver import IMPLEMENTATION_REVISION


RESULTS_ROOT = PACKAGE_ROOT / "results" / "baseline_n3"
RAW_ROOT = RESULTS_ROOT / "raw"
REVISION = "baseline_n3_v1"
PROBLEMS = (("P1", 100), ("P1", 400), ("MR", 100), ("MR", 400))
REPETITIONS = 10
N3_CELLS = ((3, 16), (3, 24), (3, 32))
N3_METHODS = (
    "path",
    "double_path",
    "path_base",
    "double_path_base",
    "path_base_oracle",
)
N2_CELLS = tuple((2, value) for value in (8, 16, 32, 64, 128, 256, 512))
EXPECTED_ROWS = len(PROBLEMS) * REPETITIONS * (
    len(N3_CELLS) * len(N3_METHODS) + len(N2_CELLS)
)
EXPECTED_NEW_ROWS = 600


def identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(task["pde_id"]),
        int(task["d"]),
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["repetition"]),
    )


def all_tasks() -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for pde_id, d in PROBLEMS:
        for repetition in range(REPETITIONS):
            for n, M in N2_CELLS:
                tasks.append(make_task("floor", pde_id, d, n, M, "path", repetition))
            for n, M in N3_CELLS:
                for method in N3_METHODS:
                    tasks.append(make_task("n3", pde_id, d, n, M, method, repetition))
    return sorted(tasks, key=task_sort_key)


def task_sort_key(task: dict[str, Any]) -> tuple[Any, ...]:
    n = int(task["n"])
    M = int(task["M"])
    method = str(task["method"])
    if n == 3 and M in (16, 24):
        priority = 0
    elif n == 3 and M == 32:
        priority = 1
    elif n == 2 and M == 256:
        priority = 2
    elif n == 2 and M == 512:
        priority = 3
    else:
        priority = 4
    method_order = {
        "path": 0,
        "path_base": 1,
        "path_base_oracle": 2,
        "double_path": 3,
        "double_path_base": 4,
    }
    return (
        priority,
        -(int(task["d"]) * M**n),
        str(task["pde_id"]),
        method_order[method],
        int(task["repetition"]),
    )


def tasks_for_stage(stage: str) -> list[dict[str, Any]]:
    tasks = all_tasks()
    if stage == "all":
        return tasks
    if stage == "priority1":
        return [task for task in tasks if int(task["n"]) == 3 and int(task["M"]) in (16, 24)]
    if stage == "priority2":
        return [task for task in tasks if int(task["n"]) == 3 and int(task["M"]) == 32]
    if stage == "floor256":
        return [task for task in tasks if int(task["n"]) == 2 and int(task["M"]) == 256]
    if stage == "floor512":
        return [task for task in tasks if int(task["n"]) == 2 and int(task["M"]) == 512]
    raise ValueError(stage)


def reusable_source(task: dict[str, Any]) -> tuple[Path, str] | None:
    n, M, method = int(task["n"]), int(task["M"]), str(task["method"])
    if n == 2 and method == "path" and M <= 128:
        path, source = frontier_source_path(task, "path")
        return path, source
    if n == 3 and M == 16 and method == "path":
        path, source = frontier_source_path(task, "path")
        return path, source
    if n == 3 and M == 16 and method == "double_path":
        candidate = task.copy()
        candidate["stage"] = "core"
        return double_task_path(candidate), "double_estimator"
    return None


def reusable_is_complete(task: dict[str, Any]) -> bool:
    source = reusable_source(task)
    if source is None:
        return False
    path, label = source
    candidate = task.copy()
    if label == "double_estimator":
        return double_task_is_complete(candidate)
    if not path.is_file():
        return False
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (
        row.get("pde_id") == task["pde_id"]
        and int(row.get("dimension", -1)) == int(task["d"])
        and int(row.get("n", -1)) == int(task["n"])
        and int(row.get("M", -1)) == int(task["M"])
        and row.get("method") == task["method"]
        and int(row.get("repetition", -1)) == int(task["repetition"])
    )


def task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / str(task["method"])
        / f"rep{int(task['repetition']):02d}.json"
    )


def task_is_complete(task: dict[str, Any]) -> bool:
    if reusable_source(task) is not None:
        return reusable_is_complete(task)
    path = task_path(task)
    if not path.is_file():
        return False
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (
        row.get("baseline_n3_revision") == REVISION
        and row.get("implementation_revision") == IMPLEMENTATION_REVISION
        and row.get("pde_id") == task["pde_id"]
        and int(row.get("dimension", -1)) == int(task["d"])
        and int(row.get("n", -1)) == int(task["n"])
        and int(row.get("M", -1)) == int(task["M"])
        and row.get("method") == task["method"]
        and int(row.get("repetition", -1)) == int(task["repetition"])
        and int(row.get("chunk_size", -1)) == chunk_size(int(task["d"]))
    )


def source_path(task: dict[str, Any]) -> tuple[Path, str]:
    reusable = reusable_source(task)
    return reusable if reusable is not None else (task_path(task), "baseline_n3")


__all__ = [
    "EXPECTED_NEW_ROWS",
    "EXPECTED_ROWS",
    "N2_CELLS",
    "N3_CELLS",
    "N3_METHODS",
    "PROBLEMS",
    "RAW_ROOT",
    "REPETITIONS",
    "RESULTS_ROOT",
    "REVISION",
    "all_tasks",
    "identity",
    "reusable_is_complete",
    "reusable_source",
    "source_path",
    "task_is_complete",
    "task_path",
    "tasks_for_stage",
    "write_json_atomic",
]
