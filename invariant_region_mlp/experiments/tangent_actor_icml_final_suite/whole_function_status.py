"""Status ledger for the 140-cell whole-function dimension benchmark."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
from typing import Any

from .protocol import PDES, RESULTS_DIR


DIMENSIONS = (20, 50, 100, 200, 400)
METHODS = ("tangent_actor", "dpi", "pinn_hte", "sdgd_pinn")
FINAL_STATUSES = ("COMPLETE", "OUT_OF_TIME", "OUT_OF_MEMORY", "FAILED", "PENDING")
INTERMEDIATE_STATUSES = ("RUNNING", "CHECKPOINTED")
STATUS_DIR = RESULTS_DIR / "whole_function_cell_status"


def key(method: str, pde: str, dimension: int) -> str:
    return f"{method}__{pde}__d{dimension}"


def validate(method: str, pde: str, dimension: int) -> None:
    if method not in METHODS:
        raise ValueError(f"unknown whole-function method {method!r}")
    if pde not in {item.key for item in PDES}:
        raise ValueError(f"unknown PDE {pde!r}")
    if dimension not in DIMENSIONS:
        raise ValueError(f"dimension {dimension} is outside the frozen whole-function grid")


def record(args: argparse.Namespace) -> Path:
    validate(args.method, args.pde, args.dimension)
    if args.status not in FINAL_STATUSES and args.status not in INTERMEDIATE_STATUSES:
        raise ValueError(args.status)
    STATUS_DIR.mkdir(parents=True, exist_ok=True)
    job_id = str(args.job_id or os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID") or "local")
    task_id = str(args.task_id or os.environ.get("SLURM_ARRAY_TASK_ID") or "0")
    path = STATUS_DIR / f"{key(args.method, args.pde, args.dimension)}__job{job_id}_{task_id}.json"
    previous: dict[str, Any] = {}
    if path.exists():
        previous = json.loads(path.read_text(encoding="utf-8"))
    try:
        git_head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[3], text=True,
        ).strip()
    except Exception:
        git_head = None
    payload = {
        **previous,
        "method": args.method,
        "pde": args.pde,
        "dimension": args.dimension,
        "seed": 0,
        "status": args.status,
        "job_id": job_id,
        "array_task_id": task_id,
        "started_utc": args.started_utc or previous.get("started_utc"),
        "finished_utc": args.finished_utc,
        "actual_gpu_seconds": args.actual_gpu_seconds,
        "exit_code": args.exit_code,
        "stdout": args.stdout,
        "stderr": args.stderr,
        "git_head": git_head,
        "cell_gpu_hour_cap": 6.0,
        "notes": args.notes,
    }
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def attempts() -> list[dict[str, Any]]:
    if not STATUS_DIR.exists():
        return []
    records = []
    for path in sorted(STATUS_DIR.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["status_file"] = str(path)
            records.append(payload)
        except (json.JSONDecodeError, OSError):
            continue
    return records


def cell_budget(method: str, pde: str, dimension: int) -> dict[str, Any]:
    """Summarize cumulative attempts without treating an allocation boundary as OOT."""
    validate(method, pde, dimension)
    rows = [
        item for item in attempts()
        if item.get("method") == method
        and item.get("pde") == pde
        and int(item.get("dimension", -1)) == dimension
    ]
    used = sum(
        float(item.get("actual_gpu_seconds") or 0.0)
        for item in rows if item.get("status") != "RUNNING"
    )
    finals = [item for item in rows if item.get("status") in FINAL_STATUSES]
    finals.sort(key=lambda item: (str(item.get("finished_utc") or ""), str(item.get("status_file") or "")))
    if finals:
        status = str(finals[-1]["status"])
    elif any(item.get("status") == "CHECKPOINTED" for item in rows):
        status = "CHECKPOINTED"
    elif any(item.get("status") == "RUNNING" for item in rows):
        status = "RUNNING"
    else:
        status = "PENDING"
    return {
        "status": status,
        "used_gpu_seconds": used,
        "remaining_gpu_seconds": max(0.0, 21600.0 - used),
        "attempt_count": len([item for item in rows if item.get("status") != "RUNNING"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    write = subparsers.add_parser("record")
    write.add_argument("--method", required=True)
    write.add_argument("--pde", required=True)
    write.add_argument("--dimension", required=True, type=int)
    write.add_argument("--status", required=True)
    write.add_argument("--job-id")
    write.add_argument("--task-id")
    write.add_argument("--started-utc")
    write.add_argument("--finished-utc")
    write.add_argument("--actual-gpu-seconds", type=float)
    write.add_argument("--exit-code", type=int)
    write.add_argument("--stdout")
    write.add_argument("--stderr")
    write.add_argument("--notes")
    show = subparsers.add_parser("list")
    budget = subparsers.add_parser("budget")
    budget.add_argument("--method", required=True)
    budget.add_argument("--pde", required=True)
    budget.add_argument("--dimension", required=True, type=int)
    args = parser.parse_args()
    if args.command == "record":
        print(record(args))
    elif args.command == "list":
        print(json.dumps(attempts(), indent=2, sort_keys=True))
    else:
        payload = cell_budget(args.method, args.pde, args.dimension)
        print(
            payload["status"],
            int(payload["remaining_gpu_seconds"]),
            int(payload["used_gpu_seconds"]),
            payload["attempt_count"],
        )


if __name__ == "__main__":
    main()
