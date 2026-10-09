"""Registered Study-G tasks and resumable artifact paths."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from invariant_region_mlp.experiments.double_estimator.equations import PACKAGE_ROOT
from invariant_region_mlp.experiments.double_estimator.protocol import write_json_atomic

from .solver import IMPLEMENTATION_REVISION


RESULTS_ROOT = PACKAGE_ROOT / "results" / "gpu_mlp"
RAW_ROOT = RESULTS_ROOT / "raw"
REVISION = "gpu_mlp_v1"
PROBLEMS = (("P1", 100), ("P1", 400), ("MR", 100), ("MR", 400))
G3_CELLS = ((2, 32), (3, 10))
G3_METHODS = ("path", "double_path", "path_base", "double_path_base")
G3_REPETITIONS = tuple(range(100, 130))
G4_CELLS = ((3, 16), (3, 32), (4, 6), (4, 10), (4, 16))
G5_CELLS = ((4, 6), (4, 10), (4, 16))
G5_METHODS = ("path", "double_path", "double_path_base")
G5_REPETITIONS = tuple(range(200, 210))


def identity(task: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(task["backend"]),
        str(task["pde_id"]),
        int(task["d"]),
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["repetition"]),
    )


def g3_tasks(backend: str | None = None) -> list[dict[str, Any]]:
    backends = (backend,) if backend else ("numpy_cpu", "torch_gpu")
    tasks = [
        {
            "stage": "g3",
            "backend": item_backend,
            "pde_id": pde,
            "d": d,
            "n": n,
            "M": M,
            "method": method,
            "repetition": rep,
        }
        for item_backend in backends
        for pde, d in PROBLEMS
        for n, M in G3_CELLS
        for method in G3_METHODS
        for rep in G3_REPETITIONS
    ]
    return sorted(
        tasks,
        key=lambda task: (
            -(int(task["d"]) * int(task["M"]) ** int(task["n"])),
            str(task["backend"]),
            str(task["pde_id"]),
            str(task["method"]),
            int(task["repetition"]),
        ),
    )


def g5_tasks() -> list[dict[str, Any]]:
    return sorted(
        [
            {
                "stage": "g5",
                "backend": "torch_gpu",
                "pde_id": pde,
                "d": d,
                "n": n,
                "M": M,
                "method": method,
                "repetition": rep,
            }
            for pde, d in PROBLEMS
            for n, M in G5_CELLS
            for method in G5_METHODS
            for rep in G5_REPETITIONS
        ],
        key=lambda task: (
            -(int(task["d"]) * int(task["M"]) ** int(task["n"])),
            str(task["pde_id"]),
            str(task["method"]),
            int(task["repetition"]),
        ),
    )


def task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / str(task["stage"])
        / str(task["backend"])
        / str(task["pde_id"])
        / f"d{int(task['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / str(task["method"])
        / f"rep{int(task['repetition']):03d}.json"
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
        row.get("gpu_mlp_revision") == REVISION
        and row.get("pde_id") == task["pde_id"]
        and int(row.get("dimension", -1)) == int(task["d"])
        and int(row.get("n", -1)) == int(task["n"])
        and int(row.get("M", -1)) == int(task["M"])
        and row.get("method") == task["method"]
        and row.get("backend") == task["backend"]
        and int(row.get("repetition", -1)) == int(task["repetition"])
        and (
            task["backend"] == "numpy_cpu"
            or row.get("implementation_revision") == IMPLEMENTATION_REVISION
        )
    )


__all__ = [
    "G3_CELLS",
    "G3_METHODS",
    "G3_REPETITIONS",
    "G4_CELLS",
    "G5_CELLS",
    "G5_METHODS",
    "G5_REPETITIONS",
    "PROBLEMS",
    "RAW_ROOT",
    "RESULTS_ROOT",
    "REVISION",
    "g3_tasks",
    "g5_tasks",
    "identity",
    "task_is_complete",
    "task_path",
    "write_json_atomic",
]
