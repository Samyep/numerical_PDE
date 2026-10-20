"""Resumable Full-History grid runner with audited S05 reuse."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.nonridge.solver import (
    IMPLEMENTATION_REVISION as NONRIDGE_IMPLEMENTATION_REVISION,
    gpu_chunk_size,
)

from .equations import fixed_points
from .protocol import (
    REVISION,
    ensure_directories,
    fh_task_path,
    fh_tasks,
    write_json_atomic,
)
from .solver import run_fh_task


NONRIDGE_SOURCE_COMMIT = "d384f4e8"


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _git_bytes(*args: str) -> bytes:
    return subprocess.run(
        ["git", *args],
        cwd=_repository_root(),
        check=True,
        capture_output=True,
    ).stdout


def _point_fingerprint(
    t: np.ndarray, x: np.ndarray, is_validation: np.ndarray
) -> str:
    digest = hashlib.sha256()
    for value in (t, x, is_validation):
        digest.update(np.ascontiguousarray(value).tobytes())
    return digest.hexdigest()


def _source_stage(task: dict[str, Any]) -> str | None:
    n, M, method, d = (
        int(task["n"]),
        int(task["M"]),
        str(task["method"]),
        int(task["d"]),
    )
    stage_a = {
        (2, 16, "path"),
        (2, 64, "raw"),
        (2, 64, "path"),
        (2, 256, "path"),
        (3, 16, "raw"),
        (3, 16, "path"),
        (3, 16, "double_path_base"),
        (3, 24, "double_path_base"),
    }
    if (n, M, method) in stage_a:
        return "stage_a"
    if (n, M, method) == (3, 32, "double_path_base"):
        return "stage_b_extension"
    if d == 400 and n == 4 and M in (6, 10) and method == "double_path_base":
        return "stage_b_n4"
    return None


def _source_path(task: dict[str, Any], stage: str) -> str:
    return (
        "invariant_region_mlp/results/nonridge/raw/"
        f"{stage}/S05/d{int(task['d'])}/"
        f"n{int(task['n'])}_M{int(task['M'])}/"
        f"{task['method']}/rep{int(task['repetition']):02d}.json"
    )


def reuse_s05_task(task: dict[str, Any]) -> dict[str, Any] | None:
    if task["pde_id"] != "S05" or int(task["repetition"]) >= 3:
        return None
    stage = _source_stage(task)
    if stage is None:
        return None
    source_path = _source_path(task, stage)
    try:
        raw = _git_bytes("show", f"{NONRIDGE_SOURCE_COMMIT}:{source_path}")
    except subprocess.CalledProcessError:
        return None
    source = json.loads(raw)
    t, x, mask = fixed_points("S05", int(task["d"]))
    expected = {
        "pde_id": "S05",
        "dimension": int(task["d"]),
        "n": int(task["n"]),
        "M": int(task["M"]),
        "method": str(task["method"]),
        "repetition": int(task["repetition"]),
        "implementation_revision": NONRIDGE_IMPLEMENTATION_REVISION,
        "chunk_size": gpu_chunk_size(
            int(task["d"]), int(task["n"]), int(task["M"])
        ),
        "point_fingerprint": _point_fingerprint(t, x, mask),
    }
    mismatches = {
        key: {"expected": value, "observed": source.get(key)}
        for key, value in expected.items()
        if source.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"S05 reuse audit failed for {source_path}: {mismatches}")
    derived = dict(source)
    derived.update(
        {
            "study_revision": REVISION,
            "reused": True,
            "reuse_source_commit": NONRIDGE_SOURCE_COMMIT,
            "reuse_source_path": source_path,
            "reuse_source_sha256": hashlib.sha256(raw).hexdigest(),
            "reuse_audit": expected,
            "original_code_commit": source.get("code_commit"),
            "prediction_path": None,
            "mc_rmse_contribution": None,
            "non_mc_rmse_proxy": None,
        }
    )
    write_json_atomic(fh_task_path(task), derived)
    return derived


def task_complete(task: dict[str, Any]) -> bool:
    path = fh_task_path(task)
    if not path.is_file():
        return False
    try:
        row = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return all(
        (
            row.get("pde_id") == task["pde_id"],
            int(row.get("dimension", -1)) == int(task["d"]),
            int(row.get("n", -1)) == int(task["n"]),
            int(row.get("M", -1)) == int(task["M"]),
            row.get("method") == task["method"],
            int(row.get("repetition", -1)) == int(task["repetition"]),
        )
    )


def filtered_tasks(
    *,
    pde_id: str | None = None,
    dimension: int | None = None,
    method: str | None = None,
    n: int | None = None,
    M: int | None = None,
    repetition: int | None = None,
) -> list[dict[str, Any]]:
    result = []
    for task in fh_tasks():
        if pde_id is not None and task["pde_id"] != pde_id:
            continue
        if dimension is not None and task["d"] != dimension:
            continue
        if method is not None and task["method"] != method:
            continue
        if n is not None and task["n"] != n:
            continue
        if M is not None and task["M"] != M:
            continue
        if repetition is not None and task["repetition"] != repetition:
            continue
        result.append(task)
    return result


def run_tasks(
    tasks: list[dict[str, Any]],
    *,
    device: str = "cuda",
    reuse: bool = True,
    reuse_only: bool = False,
    limit: int | None = None,
) -> dict[str, Any]:
    ensure_directories()
    completed = reused = executed = failed = 0
    selected = tasks if limit is None else tasks[:limit]
    for index, task in enumerate(selected, start=1):
        path = fh_task_path(task)
        if task_complete(task):
            completed += 1
            continue
        if int(task["repetition"]) > 0:
            first = dict(task)
            first["repetition"] = 0
            first_path = fh_task_path(first)
            if first_path.is_file():
                first_row = json.loads(first_path.read_text(encoding="utf-8"))
                if float(first_row.get("task_wall_seconds", 0.0)) > 1800.0:
                    write_json_atomic(
                        path.with_suffix(".curtailed.json"),
                        {
                            "study_revision": REVISION,
                            "task": task,
                            "reason": "rep0 exceeded the pre-registered 30-minute per-seed threshold",
                            "rep0_seconds": first_row["task_wall_seconds"],
                        },
                    )
                    print(
                        json.dumps(
                            {
                                "progress": f"{index}/{len(selected)}",
                                "status": "curtailed_after_rep0",
                                **task,
                            }
                        ),
                        flush=True,
                    )
                    continue
        if reuse:
            row = reuse_s05_task(task)
            if row is not None:
                reused += 1
                print(
                    json.dumps(
                        {"progress": f"{index}/{len(selected)}", "status": "reused", **task}
                    ),
                    flush=True,
                )
                continue
        if reuse_only:
            continue
        started = time.perf_counter()
        try:
            row = run_fh_task(task, device=device)
            row["task_wall_seconds"] = time.perf_counter() - started
            row["over_30_minute_threshold"] = row["task_wall_seconds"] > 1800.0
            write_json_atomic(path, row)
            executed += 1
            print(
                json.dumps(
                    {
                        "progress": f"{index}/{len(selected)}",
                        "status": "executed",
                        "rmse": row["value_rmse"],
                        "seconds": row["task_wall_seconds"],
                        **task,
                    }
                ),
                flush=True,
            )
        except Exception as error:
            failed += 1
            write_json_atomic(
                path.with_suffix(".failure.json"),
                {
                    "study_revision": REVISION,
                    "task": task,
                    "error": repr(error),
                    "elapsed_seconds": time.perf_counter() - started,
                },
            )
            print(
                json.dumps(
                    {
                        "progress": f"{index}/{len(selected)}",
                        "status": "failed",
                        "error": repr(error),
                        **task,
                    }
                ),
                flush=True,
            )
    return {
        "requested": len(selected),
        "already_complete": completed,
        "reused": reused,
        "executed": executed,
        "failed": failed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pde", choices=("S05", "B6"))
    parser.add_argument("--dimension", type=int, choices=(100, 400))
    parser.add_argument("--method", choices=("raw", "path", "double_path_base"))
    parser.add_argument("--n", type=int)
    parser.add_argument("--M", type=int)
    parser.add_argument("--repetition", type=int)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--no-reuse", action="store_true")
    parser.add_argument("--reuse-only", action="store_true")
    args = parser.parse_args()
    tasks = filtered_tasks(
        pde_id=args.pde,
        dimension=args.dimension,
        method=args.method,
        n=args.n,
        M=args.M,
        repetition=args.repetition,
    )
    result = run_tasks(
        tasks,
        device=args.device,
        reuse=not args.no_reuse,
        reuse_only=args.reuse_only,
        limit=args.limit,
    )
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
