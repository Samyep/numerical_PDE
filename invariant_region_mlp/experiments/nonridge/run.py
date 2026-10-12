"""Resumable GPU execution for the pre-registered non-ridge task matrix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import time
from typing import Any

import torch

from .calibrate import calibrate
from .protocol import (
    RESULTS_ROOT,
    REVISION,
    all_tasks,
    identity,
    stage_a_tasks,
    stage_b_extension_tasks,
    stage_b_n4_tasks,
    stage_b_s1_tasks,
    stage_b_tasks,
    task_is_complete,
    task_path,
    write_json_atomic,
)
from .solver import run_gpu_task
from .validate import validate


CURTAIL_PATH = RESULTS_ROOT / "curtailed_cells.json"
INTERIM_MARKER = RESULTS_ROOT / "interim_push_verified.json"


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _validation_passed() -> bool:
    path = RESULTS_ROOT / "validation_summary.json"
    if not path.is_file():
        return False
    return bool(json.loads(path.read_text(encoding="utf-8")).get("all_passed"))


def _cell_key(task: dict[str, Any]) -> str:
    return "|".join(
        str(task[key]) for key in ("stage", "pde_id", "d", "n", "M", "method")
    )


def _load_curtailments() -> dict[str, Any]:
    if not CURTAIL_PATH.is_file():
        return {}
    return json.loads(CURTAIL_PATH.read_text(encoding="utf-8"))


def _is_required(task: dict[str, Any], curtailments: dict[str, Any]) -> bool:
    maximum = int(curtailments.get(_cell_key(task), {}).get("max_repetitions", 10))
    return int(task["repetition"]) < maximum


def effective_tasks(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    curtailments = _load_curtailments()
    return [task for task in tasks if _is_required(task, curtailments)]


def progress(label: str, tasks: list[dict[str, Any]]) -> dict[str, Any]:
    registered = len(tasks)
    required = effective_tasks(tasks)
    complete = sum(task_is_complete(task) for task in required)
    payload = {
        "stage": label,
        "complete": complete,
        "expected": len(required),
        "registered": registered,
        "curtailed_rows": registered - len(required),
        "missing": len(required) - complete,
    }
    write_json_atomic(RESULTS_ROOT / f"{label}_progress.json", payload)
    return payload


def run_stage(
    label: str,
    tasks: list[dict[str, Any]],
    *,
    limit: int | None = None,
) -> int:
    if not _validation_passed():
        raise RuntimeError("V-1 through V-6 must all pass before any MLP run")
    if not torch.cuda.is_available():
        raise RuntimeError("the preregistered verdict rows require a CUDA GPU")
    curtailments = _load_curtailments()
    pending = [
        task
        for task in tasks
        if _is_required(task, curtailments) and not task_is_complete(task)
    ]
    if limit is not None:
        pending = pending[:limit]
    started = time.perf_counter()
    for index, task in enumerate(pending, start=1):
        if not _is_required(task, curtailments):
            continue
        row = run_gpu_task(task, device="cuda", chunk_override=int(task["chunk_size"]))
        row["nonridge_revision"] = REVISION
        row["stage"] = task["stage"]
        path = task_path(task)
        write_json_atomic(path, row)
        if (
            task["stage"] == "stage_a"
            and float(row["wall_time_seconds"]) > 1800.0
        ):
            key = _cell_key(task)
            curtailments[key] = {
                "max_repetitions": 5,
                "trigger_repetition": int(task["repetition"]),
                "trigger_wall_time_seconds": float(row["wall_time_seconds"]),
                "rule": "one Stage-A repetition exceeded 30 minutes",
            }
            write_json_atomic(CURTAIL_PATH, curtailments)
        print(
            json.dumps(
                {
                    "done": index,
                    "pending_at_start": len(pending),
                    "identity": list(identity(task)),
                    "skill": row["skill"],
                    "seconds": row["wall_time_seconds"],
                    "chunk_size": row["chunk_size"],
                    "peak_gpu_memory_bytes": row["peak_gpu_memory_bytes"],
                    "artifact": str(path),
                },
                allow_nan=True,
            ),
            flush=True,
        )
        torch.cuda.empty_cache()
    result = progress(label, tasks)
    result["elapsed_seconds_this_invocation"] = time.perf_counter() - started
    print(json.dumps(result, indent=2), flush=True)
    return 0


def collect_rows() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for task in effective_tasks(all_tasks()):
        if not task_is_complete(task):
            continue
        records.append(json.loads(task_path(task).read_text(encoding="utf-8")))
    if records:
        fields = sorted(
            {
                key
                for row in records
                for key, value in row.items()
                if not isinstance(value, (dict, list))
            }
        )
        with (RESULTS_ROOT / "rows.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(
                {key: value for key, value in row.items() if key in fields}
                for row in records
            )
    return records


def verify_interim_push() -> dict[str, Any]:
    if INTERIM_MARKER.is_file():
        return json.loads(INTERIM_MARKER.read_text(encoding="utf-8"))
    stage_status = progress("stage_a", stage_a_tasks())
    report = _repository_root() / "invariant_region_mlp" / "docs" / "NONRIDGE_REPORT.md"
    if stage_status["missing"] or not report.is_file():
        raise RuntimeError("Stage A and its interim report must be complete before Stage B")
    root = _repository_root()
    porcelain = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if porcelain:
        raise RuntimeError("commit Stage A results and interim report before Stage B")
    local = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()
    try:
        upstream = subprocess.run(
            ["git", "rev-parse", "@{u}"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except subprocess.CalledProcessError as error:
        raise RuntimeError("push the Stage A commit and set its upstream before Stage B") from error
    if local != upstream:
        raise RuntimeError("local Stage A commit is not the pushed upstream commit")
    payload = {
        "verified": True,
        "stage_a_commit": local,
        "upstream_commit": upstream,
        "verified_unix_time": time.time(),
    }
    write_json_atomic(INTERIM_MARKER, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=(
            "validate",
            "calibrate",
            "stage_a",
            "stage_b_s1",
            "stage_b_extension",
            "stage_b_n4",
            "stage_b",
            "all",
            "collect",
            "progress",
        ),
    )
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.action == "validate":
        result = validate()
        return 0 if result["all_passed"] else 1
    if args.action == "calibrate":
        if not _validation_passed():
            raise RuntimeError("validation must pass before calibration")
        calibrate()
        return 0
    if args.action == "stage_a":
        return run_stage("stage_a", stage_a_tasks(), limit=args.limit)
    if args.action in {"stage_b_s1", "stage_b_extension", "stage_b_n4", "stage_b"}:
        verify_interim_push()
        if args.action == "stage_b_s1":
            return run_stage("stage_b_s1", stage_b_s1_tasks(), limit=args.limit)
        if args.action == "stage_b_extension":
            return run_stage(
                "stage_b_extension", stage_b_extension_tasks(), limit=args.limit
            )
        if args.action == "stage_b_n4":
            return run_stage("stage_b_n4", stage_b_n4_tasks(), limit=args.limit)
        for label, tasks in (
            ("stage_b_s1", stage_b_s1_tasks()),
            ("stage_b_extension", stage_b_extension_tasks()),
            ("stage_b_n4", stage_b_n4_tasks()),
        ):
            result = run_stage(label, tasks, limit=args.limit)
            if result:
                return result
        return 0
    if args.action == "all":
        result = run_stage("stage_a", stage_a_tasks(), limit=args.limit)
        if result:
            return result
        raise RuntimeError("Stage A must be analyzed, committed, and pushed before Stage B")
    if args.action == "collect":
        rows = collect_rows()
        print(json.dumps({"rows": len(rows)}, indent=2))
        return 0
    print(
        json.dumps(
            {
                "stage_a": progress("stage_a", stage_a_tasks()),
                "stage_b_s1": progress("stage_b_s1", stage_b_s1_tasks()),
                "stage_b_extension": progress(
                    "stage_b_extension", stage_b_extension_tasks()
                ),
                "stage_b_n4": progress("stage_b_n4", stage_b_n4_tasks()),
                "stage_b_total": progress("stage_b_total", stage_b_tasks()),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
