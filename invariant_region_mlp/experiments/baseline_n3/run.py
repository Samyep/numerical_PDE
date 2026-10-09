"""Resumable Study-F runner, collector, and paired-draw audit."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
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
from invariant_region_mlp.experiments.double_estimator.protocol import chunk_size
from invariant_region_mlp.experiments.double_estimator.run import (
    CSV_FIELDS as DOUBLE_CSV_FIELDS,
)

from .protocol import (
    EXPECTED_NEW_ROWS,
    EXPECTED_ROWS,
    RESULTS_ROOT,
    REVISION,
    all_tasks,
    identity,
    reusable_source,
    source_path,
    task_is_complete,
    task_path,
    tasks_for_stage,
    write_json_atomic,
)
from .solver import CorrectionBaselineMLP, METHODS, run_task


CSV_FIELDS = tuple(
    dict.fromkeys(
        DOUBLE_CSV_FIELDS
        + (
            "dtype",
            "root_gradient_skill",
            "gaussian_vectors_drawn",
            "gaussian_coordinates_drawn",
            "root_level1_gradient_integrand_variance",
            "primary_draw_fingerprint",
            "auxiliary_draw_fingerprint",
            "source_result",
            "source_artifact",
            "baseline_n3_revision",
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
    row = run_task(task, trace_draws=False)
    row["baseline_n3_revision"] = REVISION
    row["source_result"] = "baseline_n3"
    path = task_path(task)
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
    d = int(flat["dimension"])
    terminal = int(flat.get("terminal_samples") or 0)
    transition = int(flat.get("transition_samples") or 0)
    vectors = terminal + transition
    flat.setdefault("root_gradient_skill", flat.get("gradient_relative_l2"))
    flat.setdefault("gaussian_vectors_drawn", vectors)
    flat.setdefault("gaussian_coordinates_drawn", vectors * d)
    flat.setdefault("root_level1_gradient_integrand_variance", None)
    flat.setdefault("baseline_n3_revision", None)
    return flat


def collect() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    missing: list[tuple[Any, ...]] = []
    reused = 0
    newly_computed = 0
    for task in all_tasks():
        if not task_is_complete(task):
            missing.append(identity(task))
            continue
        path, source = source_path(task)
        row = _load_flat_row(path, source)
        expected = identity(task)
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
        if source == "baseline_n3":
            newly_computed += 1
        else:
            reused += 1

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
    payload = {
        "complete": len(rows),
        "expected": EXPECTED_ROWS,
        "missing": len(missing),
        "missing_preview": [list(item) for item in missing[:50]],
        "reused": reused,
        "new": newly_computed,
        "expected_new": EXPECTED_NEW_ROWS,
    }
    write_json_atomic(RESULTS_ROOT / "progress.json", payload)
    print(json.dumps(payload, indent=2), flush=True)
    return payload


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


def run_missing(stage: str, workers: int, limit: int | None) -> int:
    registered = tasks_for_stage(stage)
    pending = [
        task
        for task in registered
        if reusable_source(task) is None and not task_is_complete(task)
    ]
    if limit is not None:
        pending = pending[:limit]
    record_environment(workers)
    print(
        json.dumps(
            {
                "stage": stage,
                "registered": len(registered),
                "pending_new": len(pending),
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
                        "method": row["method"],
                        "rep": row["repetition"],
                        "skill": row["skill"],
                        "gradient_skill": row["root_gradient_skill"],
                        "seconds": row["wall_time_seconds"],
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
                "stage_complete": stage,
                "tasks": completed,
                "elapsed_seconds": time.perf_counter() - started,
            }
        ),
        flush=True,
    )
    collect()
    return 0


def _traced_solver(task: dict[str, Any], method_name: str) -> tuple[CorrectionBaselineMLP, np.ndarray, np.ndarray]:
    d, n, M = int(task["d"]), int(task["n"]), int(task["M"])
    repetition, pde_id = int(task["repetition"]), str(task["pde_id"])
    equation = make_equation(pde_id, d)
    t, x, _ = fixed_points(pde_id, d)
    size = chunk_size(d)
    seed_parts = [BASE_SEED, d, n, M, repetition, 0]
    solver = CorrectionBaselineMLP(
        equation,
        M,
        METHODS[method_name],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
        trace_draws=True,
    )
    return solver, t[:size], x[:size]


def _audit_pair(task: dict[str, Any], baseline_method: str, counterpart: str) -> dict[str, Any]:
    first, t, x = _traced_solver(task, baseline_method)
    second, _, _ = _traced_solver(task, counterpart)
    first.solve(int(task["n"]), t, x)
    second.solve(int(task["n"]), t, x)
    primary_match = first.draw_fingerprint == second.draw_fingerprint
    auxiliary_match = first.auxiliary_draw_fingerprint == second.auxiliary_draw_fingerprint
    return {
        "identity": list(identity({**task, "method": baseline_method})),
        "counterpart": counterpart,
        "chunk_index": 0,
        "chunk_size": len(t),
        "baseline_primary": first.draw_fingerprint,
        "counterpart_primary": second.draw_fingerprint,
        "baseline_auxiliary": first.auxiliary_draw_fingerprint,
        "counterpart_auxiliary": second.auxiliary_draw_fingerprint,
        "primary_match": primary_match,
        "auxiliary_match": auxiliary_match,
        "draw_counts_match": (
            first.stats.standard_normal_variates == second.stats.standard_normal_variates
            and first.stats.time_uniform_variates == second.stats.time_uniform_variates
        ),
    }


def audit_fingerprints(workers: int) -> int:
    templates = [
        task
        for task in all_tasks()
        if int(task["n"]) == 3 and str(task["method"]) == "path"
    ]
    specs = [
        (task, baseline, counterpart)
        for task in templates
        for baseline, counterpart in (
            ("path_base", "path"),
            ("path_base_oracle", "path"),
            ("double_path_base", "double_path"),
        )
    ]
    records: list[dict[str, Any]] = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_audit_pair, *spec) for spec in specs]
        for index, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if index % 20 == 0 or index == len(futures):
                print(
                    json.dumps(
                        {
                            "pairs_done": index,
                            "total": len(futures),
                            "mismatches": sum(
                                not (
                                    row["primary_match"]
                                    and row["auxiliary_match"]
                                    and row["draw_counts_match"]
                                )
                                for row in records
                            ),
                        }
                    ),
                    flush=True,
                )
    records.sort(key=lambda row: (tuple(row["identity"]), row["counterpart"]))
    failures = [
        row
        for row in records
        if not (row["primary_match"] and row["auxiliary_match"] and row["draw_counts_match"])
    ]
    payload = {
        "revision": REVISION,
        "code_commit": _current_commit(),
        "scope": "one complete registered chunk for every n=3 cell/repetition and each baseline pairing",
        "pairs": len(records),
        "mismatches": len(failures),
        "passed": not failures,
        "elapsed_seconds": time.perf_counter() - started,
        "records_sha256": hashlib.sha256(
            json.dumps(records, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "failure_details": failures[:20],
    }
    write_json_atomic(RESULTS_ROOT / "fingerprint_audit.json", payload)
    print(json.dumps(payload, indent=2), flush=True)
    return 0 if payload["passed"] else 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=("priority1", "priority2", "floor256", "floor512", "all", "collect", "audit"),
    )
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.action == "collect":
        collect()
        return 0
    if args.action == "audit":
        return audit_fingerprints(max(1, args.workers))
    return run_missing(args.action, max(1, args.workers), args.limit)


if __name__ == "__main__":
    raise SystemExit(main())
