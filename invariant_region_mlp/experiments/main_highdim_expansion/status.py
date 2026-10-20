"""Attempt ledger and cumulative six-GH200-hour cell budget."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
from typing import Any

from .protocol import CELL_GPU_SECONDS_CAP, DIMENSIONS, METHODS, PDES, SOURCE_ROOT, STATUS_DIR


TERMINAL_STATUSES = ("COMPLETE", "OUT_OF_TIME", "OUT_OF_MEMORY", "FAILED")
INTERMEDIATE_STATUSES = ("RUNNING", "CHECKPOINTED")


def key(method: str, pde: str, dimension: int) -> str:
    return f"{method}__{pde}__d{dimension}"


def validate(method: str, pde: str, dimension: int) -> None:
    if method not in METHODS:
        raise ValueError(method)
    if pde not in PDES:
        raise ValueError(pde)
    if int(dimension) not in DIMENSIONS:
        raise ValueError(dimension)


def record(args: argparse.Namespace) -> Path:
    validate(args.method, args.pde, args.dimension)
    if args.status not in TERMINAL_STATUSES + INTERMEDIATE_STATUSES:
        raise ValueError(args.status)
    STATUS_DIR.mkdir(parents=True, exist_ok=True)
    job_id = str(
        args.job_id
        or os.environ.get("SLURM_ARRAY_JOB_ID")
        or os.environ.get("SLURM_JOB_ID")
        or "local"
    )
    task_id = str(args.task_id or os.environ.get("SLURM_ARRAY_TASK_ID") or "0")
    path = STATUS_DIR / f"{key(args.method, args.pde, args.dimension)}__job{job_id}_{task_id}.json"
    previous: dict[str, Any] = {}
    if path.is_file():
        previous = json.loads(path.read_text(encoding="utf-8"))
    git_head = previous.get("git_head")
    if not git_head:
        git_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
        ).strip()
    payload = {
        **previous,
        "method": args.method,
        "pde": args.pde,
        "dimension": int(args.dimension),
        "seed": 0,
        "status": args.status,
        "job_id": job_id,
        "task_id": task_id,
        "started_utc": args.started_utc or previous.get("started_utc"),
        "finished_utc": args.finished_utc,
        "actual_gpu_seconds": float(args.actual_gpu_seconds or 0.0),
        "exit_code": args.exit_code,
        "stdout": args.stdout,
        "stderr": args.stderr,
        "git_head": git_head,
        "cell_gpu_hour_cap": CELL_GPU_SECONDS_CAP / 3600.0,
        "notes": args.notes,
    }
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def attempts() -> list[dict[str, Any]]:
    if not STATUS_DIR.exists():
        return []
    result = []
    for path in sorted(STATUS_DIR.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["status_file"] = str(path)
            result.append(payload)
        except (OSError, json.JSONDecodeError):
            continue
    return result


def cell_budget(method: str, pde: str, dimension: int) -> dict[str, Any]:
    validate(method, pde, dimension)
    rows = [
        row
        for row in attempts()
        if row.get("method") == method
        and row.get("pde") == pde
        and int(row.get("dimension", -1)) == int(dimension)
    ]
    finished = [row for row in rows if row.get("status") != "RUNNING"]
    used = sum(float(row.get("actual_gpu_seconds") or 0.0) for row in finished)
    ordered = sorted(
        finished,
        key=lambda row: (str(row.get("finished_utc") or ""), str(row.get("status_file") or "")),
    )
    if ordered and ordered[-1].get("status") in TERMINAL_STATUSES:
        status = str(ordered[-1]["status"])
    elif any(row.get("status") == "CHECKPOINTED" for row in rows):
        status = "CHECKPOINTED"
    elif any(row.get("status") == "RUNNING" for row in rows):
        status = "RUNNING"
    else:
        status = "PENDING"
    if used >= CELL_GPU_SECONDS_CAP and status not in {"COMPLETE", "OUT_OF_MEMORY", "FAILED"}:
        status = "OUT_OF_TIME"
    return {
        "status": status,
        "used_gpu_seconds": used,
        "remaining_gpu_seconds": max(0.0, CELL_GPU_SECONDS_CAP - used),
        "attempt_count": len(finished),
    }


def grid() -> list[dict[str, Any]]:
    return [
        {"method": method, "pde": pde, "dimension": d, **cell_budget(method, pde, d)}
        for pde in PDES
        for d in DIMENSIONS
        for method in METHODS
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    write = sub.add_parser("record")
    for flag in ("method", "pde", "status"):
        write.add_argument(f"--{flag}", required=True)
    write.add_argument("--dimension", required=True, type=int)
    write.add_argument("--job-id")
    write.add_argument("--task-id")
    write.add_argument("--started-utc")
    write.add_argument("--finished-utc")
    write.add_argument("--actual-gpu-seconds", type=float)
    write.add_argument("--exit-code", type=int)
    write.add_argument("--stdout")
    write.add_argument("--stderr")
    write.add_argument("--notes")
    budget = sub.add_parser("budget")
    budget.add_argument("--method", required=True)
    budget.add_argument("--pde", required=True)
    budget.add_argument("--dimension", type=int, required=True)
    sub.add_parser("grid")
    args = parser.parse_args()
    if args.command == "record":
        print(record(args))
    elif args.command == "grid":
        print(json.dumps(grid(), indent=2, sort_keys=True))
    else:
        result = cell_budget(args.method, args.pde, args.dimension)
        print(
            result["status"],
            int(result["remaining_gpu_seconds"]),
            int(result["used_gpu_seconds"]),
            result["attempt_count"],
        )


if __name__ == "__main__":
    main()
