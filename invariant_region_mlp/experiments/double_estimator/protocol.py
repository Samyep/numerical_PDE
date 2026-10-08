"""Frozen task matrix and resumable artifact layout."""

from __future__ import annotations

from collections.abc import Iterable
import json
from pathlib import Path
from typing import Any

from .equations import PACKAGE_ROOT
from .solver import IMPLEMENTATION_REVISION, methods_for_pde


RESULTS_ROOT = PACKAGE_ROOT / "results" / "double_estimator"
RAW_ROOT = RESULTS_ROOT / "raw"
CORE_CELLS = ((2, 32), (3, 6), (3, 10), (4, 6))
EC_CELLS = (
    (2, 8),
    (2, 16),
    (2, 32),
    (2, 64),
    (2, 128),
    (3, 4),
    (3, 6),
    (3, 10),
    (3, 16),
    (4, 3),
    (4, 4),
    (4, 6),
)
REPETITIONS = 10


def chunk_size(d: int) -> int:
    if d <= 20:
        return 8
    if d <= 100:
        return 4
    return 1


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
        "chunk_size": chunk_size(d),
    }


def cross_tasks(
    stage: str,
    problems: Iterable[tuple[str, int]],
    cells: Iterable[tuple[int, int]],
    methods: Any,
) -> list[dict[str, Any]]:
    return [
        make_task(stage, pde_id, d, n, M, method, repetition)
        for pde_id, d in problems
        for n, M in cells
        for method in (methods(pde_id) if callable(methods) else methods)
        for repetition in range(REPETITIONS)
    ]


def core_tasks() -> list[dict[str, Any]]:
    problems = [
        ("P1", 20),
        ("P1", 100),
        ("P1", 400),
        ("MR", 100),
        ("MR", 400),
    ]
    return cross_tasks("core", problems, CORE_CELLS, methods_for_pde)


def c2_tasks() -> list[dict[str, Any]]:
    problems = [
        (pde_id, d)
        for pde_id in ("C2-convex", "C2-cancel", "C2-flip")
        for d in (20, 100, 400)
    ]
    return cross_tasks("c2", problems, CORE_CELLS, methods_for_pde)


def ec_tasks() -> list[dict[str, Any]]:
    problems = [
        ("P1", 100),
        ("P1", 400),
        ("MR", 100),
        ("MR", 400),
    ]
    return cross_tasks(
        "ec", problems, EC_CELLS, ("raw", "double", "double_path")
    )


def p4_tasks() -> list[dict[str, Any]]:
    return cross_tasks(
        "p4",
        (("P4", 20), ("P4", 100)),
        CORE_CELLS,
        methods_for_pde,
    )


def task_identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(task["pde_id"]),
        int(task["d"]),
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["repetition"]),
    )


def all_tasks() -> list[dict[str, Any]]:
    by_identity: dict[tuple[Any, ...], dict[str, Any]] = {}
    stages: dict[tuple[Any, ...], set[str]] = {}
    for task in core_tasks() + c2_tasks() + ec_tasks() + p4_tasks():
        identity = task_identity(task)
        stages.setdefault(identity, set()).add(str(task["stage"]))
        by_identity.setdefault(identity, task.copy())
    for identity, task in by_identity.items():
        task["stages"] = sorted(stages[identity])
    return list(by_identity.values())


def tasks_for_stage(stage: str) -> list[dict[str, Any]]:
    if stage == "core":
        tasks = core_tasks()
    elif stage == "c2":
        tasks = c2_tasks()
    elif stage == "ec":
        tasks = ec_tasks()
    elif stage == "p4":
        tasks = p4_tasks()
    elif stage == "all":
        return sorted(all_tasks(), key=task_sort_key)
    else:
        raise ValueError(stage)
    identities = {task_identity(task) for task in tasks}
    canonical = {
        task_identity(task): task for task in all_tasks() if task_identity(task) in identities
    }
    return sorted(canonical.values(), key=task_sort_key)


def task_sort_key(task: dict[str, Any]) -> tuple[Any, ...]:
    method_order = {
        "raw": 0,
        "double": 1,
        "path": 2,
        "double_path": 3,
        "box": 4,
        "sub_box": 4,
        "oracle_state": 5,
        "centre": 6,
        "f_zero": 7,
    }
    return (
        int(task["n"]),
        int(task["M"]),
        int(task["d"]),
        str(task["pde_id"]),
        method_order[str(task["method"])],
        int(task["repetition"]),
    )


def safe_name(value: str) -> str:
    return value.replace("-", "m").replace(".", "p").replace("+", "")


def task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / safe_name(str(task["pde_id"]))
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / safe_name(str(task["method"]))
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
        row.get("implementation_revision") == IMPLEMENTATION_REVISION
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

