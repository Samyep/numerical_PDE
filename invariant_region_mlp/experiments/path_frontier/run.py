"""Resumable runner, collector, and paired-RNG audit."""

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

from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    fixed_points,
    make_equation,
)
from invariant_region_mlp.experiments.double_estimator.protocol import (
    chunk_size,
    task_is_complete as double_task_is_complete,
)
from invariant_region_mlp.experiments.double_estimator.run import (
    CSV_FIELDS as DOUBLE_CSV_FIELDS,
)
from invariant_region_mlp.experiments.double_estimator.solver import (
    DoubleEstimatorMLP,
    METHODS as DOUBLE_METHODS,
    run_task,
)

from .protocol import (
    EXPECTED_COMBINED_ROWS,
    EXPECTED_PATH_ROWS,
    METHODS,
    PATH_FRONTIER_REVISION,
    RESULTS_ROOT,
    REPETITIONS,
    has_reusable_path,
    identity,
    new_task_is_complete,
    new_task_path,
    path_tasks,
    source_path,
    write_json_atomic,
)


CSV_FIELDS = tuple(
    dict.fromkeys(
        DOUBLE_CSV_FIELDS
        + (
            "dtype",
            "primary_draw_fingerprint",
            "auxiliary_draw_fingerprint",
            "source_result",
            "source_artifact",
            "path_frontier_revision",
        )
    )
)


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _current_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _worker(task: dict[str, Any]) -> tuple[dict[str, Any], str]:
    # The separate audit traces one complete registered chunk for every task
    # identity. Avoid hashing tens of gigabytes in each 1,200-point estimator.
    row = run_task(task, trace_draws=False)
    row["path_frontier_revision"] = PATH_FRONTIER_REVISION
    row["source_result"] = "path_frontier"
    path = new_task_path(task)
    write_json_atomic(path, row)
    return row, str(path)


def _load_flat_row(path: Path, source: str) -> dict[str, Any]:
    row = json.loads(path.read_text(encoding="utf-8"))
    flat = {key: value for key, value in row.items() if not isinstance(value, dict)}
    flat["source_result"] = source
    try:
        flat["source_artifact"] = path.relative_to(_repository_root()).as_posix()
    except ValueError:
        flat["source_artifact"] = str(path)
    flat.setdefault("path_frontier_revision", None)
    return flat


def collect() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    missing: list[tuple[Any, ...]] = []
    reused_path = 0
    new_path = 0
    for task in path_tasks():
        for method in METHODS:
            candidate = task.copy()
            candidate["method"] = method
            path, source = source_path(candidate, method)
            complete = (
                double_task_is_complete(candidate)
                if source == "double_estimator"
                else new_task_is_complete(candidate)
            )
            if not complete:
                missing.append(identity(candidate))
                continue
            row = _load_flat_row(path, source)
            expected = identity(candidate)
            observed = (
                row.get("pde_id"),
                int(row.get("dimension", -1)),
                int(row.get("n", -1)),
                int(row.get("M", -1)),
                row.get("method"),
                int(row.get("repetition", -1)),
            )
            if observed != expected:
                raise ValueError(f"artifact identity mismatch: {path}: {observed} != {expected}")
            rows.append(row)
            if method == "path":
                if source == "double_estimator":
                    reused_path += 1
                else:
                    new_path += 1

    rows.sort(
        key=lambda row: (
            str(row["pde_id"]),
            int(row["dimension"]),
            int(row["n"]),
            int(row["M"]),
            str(row["method"]),
            int(row["repetition"]),
        )
    )
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    with (RESULTS_ROOT / "rows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in CSV_FIELDS})

    progress = {
        "combined_complete": len(rows),
        "combined_expected": EXPECTED_COMBINED_ROWS,
        "missing": len(missing),
        "missing_preview": [list(item) for item in missing[:50]],
        "path_rows_expected": EXPECTED_PATH_ROWS,
        "path_rows_reused": reused_path,
        "path_rows_new": new_path,
    }
    write_json_atomic(RESULTS_ROOT / "progress.json", progress)
    print(json.dumps(progress, indent=2), flush=True)
    return progress


def record_environment(workers: int) -> None:
    payload = {
        "code_commit": _current_commit(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "logical_cpus": os.cpu_count(),
        "workers": workers,
        "dtype": "float64",
    }
    write_json_atomic(RESULTS_ROOT / "environment.json", payload)


def run_missing(workers: int, limit: int | None) -> int:
    registered = path_tasks()
    reused = [task for task in registered if has_reusable_path(task)]
    pending = [
        task
        for task in registered
        if not has_reusable_path(task) and not new_task_is_complete(task)
    ]
    if limit is not None:
        pending = pending[:limit]
    record_environment(workers)
    print(
        json.dumps(
            {
                "registered_path_rows": len(registered),
                "reused_path_rows": len(reused),
                "pending_new_path_rows": len(pending),
                "workers": workers,
                "limit": limit,
            }
        ),
        flush=True,
    )
    if not pending:
        collect()
        return 0

    started = time.perf_counter()
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_worker, task): task for task in pending}
        for future in as_completed(futures):
            task = futures[future]
            row, artifact = future.result()
            completed += 1
            print(
                json.dumps(
                    {
                        "done": completed,
                        "total": len(pending),
                        "pde": row["pde_id"],
                        "d": row["dimension"],
                        "cell": [row["n"], row["M"]],
                        "rep": row["repetition"],
                        "skill": row["skill"],
                        "seconds": row["wall_time_seconds"],
                        "primary_draw_fingerprint": row["primary_draw_fingerprint"],
                        "artifact": artifact,
                        "task": list(identity(task)),
                    },
                    allow_nan=True,
                ),
                flush=True,
            )
    print(
        json.dumps(
            {
                "run_complete": completed,
                "elapsed_seconds": time.perf_counter() - started,
            }
        ),
        flush=True,
    )
    collect()
    return 0


def _traced_solver(
    task: dict[str, Any], method: str
) -> tuple[DoubleEstimatorMLP, np.ndarray, np.ndarray]:
    d = int(task["d"])
    n = int(task["n"])
    M = int(task["M"])
    repetition = int(task["repetition"])
    pde_id = str(task["pde_id"])
    equation = make_equation(pde_id, d)
    t, x, _ = fixed_points(pde_id, d)
    size = chunk_size(d)
    seed_parts = [BASE_SEED, d, n, M, repetition, 0]
    solver = DoubleEstimatorMLP(
        equation,
        M,
        DOUBLE_METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        trace_draws=True,
    )
    return solver, t[:size], x[:size]


def _fingerprint_worker(task: dict[str, Any]) -> dict[str, Any]:
    raw, t, x = _traced_solver(task, "raw")
    path, _, _ = _traced_solver(task, "path")
    raw.solve(int(task["n"]), t, x)
    path.solve(int(task["n"]), t, x)
    return {
        "identity": list(identity(task)),
        "chunk_index": 0,
        "chunk_size": len(t),
        "raw_fingerprint": raw.draw_fingerprint,
        "path_fingerprint": path.draw_fingerprint,
        "match": raw.draw_fingerprint == path.draw_fingerprint,
        "raw_terminal_samples": raw.stats.terminal_samples,
        "path_terminal_samples": path.stats.terminal_samples,
        "raw_transition_samples": raw.stats.transition_samples,
        "path_transition_samples": path.stats.transition_samples,
    }


def audit_fingerprints(workers: int) -> int:
    """Trace one complete registered chunk for every paired task identity.

    This is deliberately not a rerun of any 1,200-point result row. All chunks
    for a task have the same shape; the chunk index only changes its seed.
    """

    tasks = path_tasks()
    started = time.perf_counter()
    records: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_fingerprint_worker, task) for task in tasks]
        for index, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            records.append(record)
            if index % 25 == 0 or index == len(tasks):
                print(
                    json.dumps(
                        {
                            "fingerprint_pairs_done": index,
                            "total": len(tasks),
                            "mismatches": sum(not item["match"] for item in records),
                        }
                    ),
                    flush=True,
                )
    records.sort(key=lambda item: tuple(item["identity"]))
    mismatches = [item for item in records if not item["match"]]
    count_mismatches = [
        item
        for item in records
        if item["raw_terminal_samples"] != item["path_terminal_samples"]
        or item["raw_transition_samples"] != item["path_transition_samples"]
    ]
    payload = {
        "audit_revision": PATH_FRONTIER_REVISION,
        "code_commit": _current_commit(),
        "scope": (
            "one complete registered chunk for every (PDE,d,n,M,repetition); "
            "no 1200-point existing result row was rerun"
        ),
        "pairs": len(records),
        "fingerprint_mismatches": len(mismatches),
        "draw_count_mismatches": len(count_mismatches),
        "passed": not mismatches and not count_mismatches,
        "elapsed_seconds": time.perf_counter() - started,
        "mismatch_details": mismatches[:20],
    }
    write_json_atomic(RESULTS_ROOT / "fingerprint_audit.json", payload)
    print(json.dumps(payload, indent=2), flush=True)
    return 0 if payload["passed"] else 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("run", "collect", "audit"))
    parser.add_argument(
        "--workers", type=int, default=max(1, min(8, os.cpu_count() or 1))
    )
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.action == "collect":
        collect()
        return 0
    if args.action == "audit":
        return audit_fingerprints(max(1, args.workers))
    return run_missing(max(1, args.workers), args.limit)


if __name__ == "__main__":
    raise SystemExit(main())
