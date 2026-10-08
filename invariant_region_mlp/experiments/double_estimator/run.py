"""Resumable local runner for the pre-registered study."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any

import numpy as np

from .protocol import (
    RESULTS_ROOT,
    all_tasks,
    task_identity,
    task_is_complete,
    task_path,
    tasks_for_stage,
    write_json_atomic,
)
from .solver import run_task


CSV_FIELDS = (
    "pde_id",
    "family",
    "dimension",
    "n",
    "M",
    "method",
    "repetition",
    "skill",
    "value_rmse",
    "value_relative_l2",
    "value_bias",
    "gradient_relative_l2",
    "validation_skill",
    "all_skill",
    "mean_generator_bias",
    "generator_rmse",
    "generator_calls",
    "recursive_calls",
    "auxiliary_recursive_calls",
    "recursively_evaluated_states",
    "terminal_samples",
    "transition_samples",
    "wall_time_seconds",
    "nonfinite_state_count",
    "nonfinite_generator_count",
    "base_seed",
    "chunk_size",
    "point_count",
    "point_fingerprint",
    "seed_scheme",
    "aux_seed_scheme",
    "implementation_revision",
    "code_commit",
)


def _worker(task: dict[str, Any]) -> tuple[dict[str, Any], str]:
    row = run_task(task)
    path = task_path(task)
    write_json_atomic(path, row)
    return row, str(path)


def collect() -> int:
    rows: list[dict[str, Any]] = []
    missing: list[tuple[Any, ...]] = []
    for task in all_tasks():
        if not task_is_complete(task):
            missing.append(task_identity(task))
            continue
        rows.append(json.loads(task_path(task).read_text(encoding="utf-8")))
    rows.sort(
        key=lambda row: (
            row["pde_id"],
            row["dimension"],
            row["n"],
            row["M"],
            row["method"],
            row["repetition"],
        )
    )
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    destination = RESULTS_ROOT / "rows.csv"
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in CSV_FIELDS})
    progress = {
        "complete": len(rows),
        "total_unique": len(all_tasks()),
        "missing": len(missing),
        "missing_preview": missing[:50],
    }
    write_json_atomic(RESULTS_ROOT / "progress.json", progress)
    print(json.dumps(progress, indent=2), flush=True)
    return 0


def record_environment() -> None:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parents[3],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    payload = {
        "code_commit": commit,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "logical_cpus": os.cpu_count(),
        "dtype": "float64",
    }
    write_json_atomic(RESULTS_ROOT / "environment.json", payload)


def run_stage(stage: str, workers: int, limit: int | None) -> int:
    tasks = [task for task in tasks_for_stage(stage) if not task_is_complete(task)]
    if limit is not None:
        tasks = tasks[:limit]
    record_environment()
    print(
        json.dumps(
            {
                "stage": stage,
                "pending": len(tasks),
                "workers": workers,
                "limit": limit,
            }
        ),
        flush=True,
    )
    if not tasks:
        return collect()
    started = time.perf_counter()
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_worker, task): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            row, path = future.result()
            completed += 1
            print(
                json.dumps(
                    {
                        "done": completed,
                        "total": len(tasks),
                        "pde": row["pde_id"],
                        "d": row["dimension"],
                        "cell": [row["n"], row["M"]],
                        "method": row["method"],
                        "rep": row["repetition"],
                        "skill": row["skill"],
                        "bias": row["mean_generator_bias"],
                        "seconds": row["wall_time_seconds"],
                        "artifact": path,
                    },
                    allow_nan=True,
                ),
                flush=True,
            )
    print(
        json.dumps(
            {
                "stage_complete": stage,
                "tasks": completed,
                "elapsed_seconds": time.perf_counter() - started,
            }
        ),
        flush=True,
    )
    return collect()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "stage", choices=("core", "c2", "ec", "p4", "all", "collect")
    )
    parser.add_argument("--workers", type=int, default=max(1, min(6, os.cpu_count() or 1)))
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.stage == "collect":
        return collect()
    return run_stage(args.stage, max(1, args.workers), args.limit)


if __name__ == "__main__":
    raise SystemExit(main())

