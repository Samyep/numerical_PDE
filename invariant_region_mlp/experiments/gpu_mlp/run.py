"""Resumable G-3 agreement, G-4 speed, and gated G-5 runners."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.baseline_n3.solver import (
    CorrectionBaselineMLP,
    METHODS,
    run_task as run_cpu_task,
)
from invariant_region_mlp.experiments.double_estimator.equations import (
    BASE_SEED,
    fixed_points,
    make_equation,
)
from invariant_region_mlp.experiments.double_estimator.protocol import chunk_size

from .protocol import (
    G4_CELLS,
    RESULTS_ROOT,
    REVISION,
    g3_tasks,
    g5_tasks,
    identity,
    task_is_complete,
    task_path,
    write_json_atomic,
)
from .solver import run_gpu_task
from .validate import validate


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


def _cpu_worker(task: dict[str, Any]) -> tuple[dict[str, Any], str]:
    cpu_task = {
        **task,
        "chunk_size": chunk_size(int(task["d"])),
    }
    row = run_cpu_task(cpu_task, trace_draws=False)
    row.update(
        {
            "backend": "numpy_cpu",
            "paired_with_gpu": False,
            "gpu_mlp_revision": REVISION,
            "gpu_study_code_commit": _current_commit(),
        }
    )
    vectors = int(row["terminal_samples"] + row["transition_samples"])
    row.setdefault("gaussian_vectors_drawn", vectors)
    row.setdefault("gaussian_coordinates_drawn", vectors * int(task["d"]))
    path = task_path(task)
    write_json_atomic(path, row)
    return row, str(path)


def _gpu_worker(task: dict[str, Any]) -> tuple[dict[str, Any], str]:
    row = run_gpu_task(task, device="cuda" if torch.cuda.is_available() else "cpu")
    row.update(
        {
            "backend": task["backend"],
            "gpu_mlp_revision": REVISION,
            "code_commit": _current_commit(),
        }
    )
    path = task_path(task)
    write_json_atomic(path, row)
    return row, str(path)


def _progress(stage: str, tasks: list[dict[str, Any]]) -> dict[str, Any]:
    complete = sum(task_is_complete(task) for task in tasks)
    payload = {
        "stage": stage,
        "complete": complete,
        "expected": len(tasks),
        "missing": len(tasks) - complete,
    }
    write_json_atomic(RESULTS_ROOT / f"{stage}_progress.json", payload)
    return payload


def run_g3_cpu(workers: int, limit: int | None = None) -> int:
    tasks = [task for task in g3_tasks("numpy_cpu") if not task_is_complete(task)]
    if limit is not None:
        tasks = tasks[:limit]
    print(json.dumps({"stage": "g3_cpu", "pending": len(tasks), "workers": workers}), flush=True)
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_cpu_worker, task): task for task in tasks}
        for index, future in enumerate(as_completed(futures), start=1):
            row, artifact = future.result()
            print(
                json.dumps(
                    {
                        "done": index,
                        "total": len(tasks),
                        "identity": list(identity(futures[future])),
                        "skill": row["skill"],
                        "seconds": row["wall_time_seconds"],
                        "artifact": artifact,
                    },
                    allow_nan=True,
                ),
                flush=True,
            )
    print(json.dumps({"elapsed_seconds": time.perf_counter() - started}), flush=True)
    _progress("g3_cpu", g3_tasks("numpy_cpu"))
    return 0


def run_gpu_rows(stage: str, tasks: list[dict[str, Any]], limit: int | None = None) -> int:
    pending = [task for task in tasks if not task_is_complete(task)]
    if limit is not None:
        pending = pending[:limit]
    print(
        json.dumps(
            {
                "stage": stage,
                "pending": len(pending),
                "device": "cuda" if torch.cuda.is_available() else "cpu",
            }
        ),
        flush=True,
    )
    started = time.perf_counter()
    for index, task in enumerate(pending, start=1):
        row, artifact = _gpu_worker(task)
        print(
            json.dumps(
                {
                    "done": index,
                    "total": len(pending),
                    "identity": list(identity(task)),
                    "skill": row["skill"],
                    "seconds": row["wall_time_seconds"],
                    "peak_gpu_memory_bytes": row["peak_gpu_memory_bytes"],
                    "artifact": artifact,
                },
                allow_nan=True,
            ),
            flush=True,
        )
    print(json.dumps({"elapsed_seconds": time.perf_counter() - started}), flush=True)
    _progress(stage, tasks)
    return 0


def _one_cpu_chunk_time(task: dict[str, Any]) -> float:
    d, n, M = int(task["d"]), int(task["n"]), int(task["M"])
    pde_id, method = str(task["pde_id"]), str(task["method"])
    repetition = int(task["repetition"])
    equation = make_equation(pde_id, d)
    t, x, _ = fixed_points(pde_id, d)
    size = chunk_size(d)
    seed_parts = [BASE_SEED, d, n, M, repetition, 0]
    solver = CorrectionBaselineMLP(
        equation,
        M,
        METHODS[method],
        np.random.default_rng(np.random.SeedSequence(seed_parts)),
        aux_seed_sequence=np.random.SeedSequence(seed_parts + [0xD0B1E]),
        dose_rng=np.random.default_rng(np.random.SeedSequence(seed_parts + [0xD05E])),
    )
    started = time.perf_counter()
    solver.solve(n, t[:size], x[:size])
    return time.perf_counter() - started


def run_g4() -> int:
    if not torch.cuda.is_available():
        write_json_atomic(
            RESULTS_ROOT / "g4_status.json",
            {"verdict": "NOT EVALUATED (NO GPU)", "rows": 0},
        )
        return 0
    rows: list[dict[str, Any]] = []
    for index, (n, M) in enumerate(G4_CELLS):
        task = {
            "stage": "g4",
            "backend": "torch_gpu",
            "pde_id": "P1",
            "d": 400,
            "n": n,
            "M": M,
            "method": "double_path_base",
            "repetition": 300 + index,
        }
        gpu_row = run_gpu_task(task, device="cuda")
        one_chunk = _one_cpu_chunk_time(task)
        chunks = math.ceil(1200 / chunk_size(400))
        estimate = one_chunk * chunks
        if estimate <= 7200.0:
            full = run_cpu_task({**task, "chunk_size": chunk_size(400)})
            cpu_seconds = float(full["wall_time_seconds"])
            extrapolated = False
        else:
            cpu_seconds = estimate
            extrapolated = True
        rows.append(
            {
                "pde_id": "P1",
                "dimension": 400,
                "n": n,
                "M": M,
                "method": "double_path_base",
                "gpu_wall_time_seconds": gpu_row["wall_time_seconds"],
                "cpu_wall_time_seconds": cpu_seconds,
                "cpu_extrapolated": extrapolated,
                "cpu_one_chunk_seconds": one_chunk,
                "cpu_chunk_size": chunk_size(400),
                "cpu_logical_cores": os.cpu_count(),
                "speedup_cpu_over_gpu": cpu_seconds / gpu_row["wall_time_seconds"],
                "peak_gpu_memory_bytes": gpu_row["peak_gpu_memory_bytes"],
                "gpu_chunk_size": gpu_row["chunk_size"],
                "gpu_skill": gpu_row["skill"],
            }
        )
        print(json.dumps(rows[-1], allow_nan=True), flush=True)
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    with (RESULTS_ROOT / "speed.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json_atomic(
        RESULTS_ROOT / "g4_status.json",
        {"verdict": "REPORTED", "rows": len(rows), "code_commit": _current_commit()},
    )
    return 0


def _gates_pass() -> tuple[bool, dict[str, Any]]:
    validation_path = RESULTS_ROOT / "validation_summary.json"
    agreement_path = RESULTS_ROOT / "agreement_summary.json"
    validation = (
        json.loads(validation_path.read_text(encoding="utf-8"))
        if validation_path.is_file()
        else {}
    )
    agreement = (
        json.loads(agreement_path.read_text(encoding="utf-8"))
        if agreement_path.is_file()
        else {}
    )
    passed = (
        validation.get("G-2", {}).get("verdict") == "PASS"
        and agreement.get("verdict") == "PASS"
    )
    return passed, {"validation": validation.get("G-2"), "agreement": agreement}


def run_g5(limit: int | None = None) -> int:
    passed, evidence = _gates_pass()
    if not passed:
        write_json_atomic(
            RESULTS_ROOT / "g5_status.json",
            {"verdict": "NOT RUN (G-2/G-3 GATE)", "evidence": evidence},
        )
        print(json.dumps({"g5": "blocked by gate", "evidence": evidence}, indent=2), flush=True)
        return 2
    result = run_gpu_rows("g5", g5_tasks(), limit)
    progress = _progress("g5", g5_tasks())
    write_json_atomic(
        RESULTS_ROOT / "g5_status.json",
        {
            "verdict": "REPORTED" if progress["missing"] == 0 else "NOT EVALUATED (COMPUTE)",
            **progress,
        },
    )
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=("validate", "g3_cpu", "g3_gpu", "g3", "g4", "g5", "progress"),
    )
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.action == "validate":
        validate()
        return 0
    if args.action == "g3_cpu":
        return run_g3_cpu(max(1, args.workers), args.limit)
    if args.action == "g3_gpu":
        return run_gpu_rows("g3_gpu", g3_tasks("torch_gpu"), args.limit)
    if args.action == "g3":
        result = run_g3_cpu(max(1, args.workers), args.limit)
        if result:
            return result
        return run_gpu_rows("g3_gpu", g3_tasks("torch_gpu"), args.limit)
    if args.action == "g4":
        return run_g4()
    if args.action == "g5":
        return run_g5(args.limit)
    print(json.dumps({
        "g3_cpu": _progress("g3_cpu", g3_tasks("numpy_cpu")),
        "g3_gpu": _progress("g3_gpu", g3_tasks("torch_gpu")),
        "g5": _progress("g5", g5_tasks()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
