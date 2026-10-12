"""Create transparent GH200 accounting from completed top-level Slurm jobs."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import re
import subprocess
from typing import Any

from .common import RESULTS_DIR, atomic_json


TERMINAL_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "TIMEOUT",
}


def _gpu_count(alloc_tres: str) -> int:
    match = re.search(r"(?:^|,)gres/gpu=(\d+)(?:,|$)", alloc_tres)
    return int(match.group(1)) if match else 0


def collect(job_specs: list[str]) -> dict[str, Any]:
    roles: dict[str, str] = {}
    for item in job_specs:
        role, separator, job_id = item.partition("=")
        if not separator or not role or not job_id.isdigit():
            raise ValueError(f"expected ROLE=JOB_ID, received {item!r}")
        roles[job_id] = role

    command = [
        "sacct",
        "-j",
        ",".join(roles),
        "--starttime",
        "2026-10-12",
        "--noheader",
        "--parsable2",
        "--format=JobIDRaw,JobName,State,ElapsedRaw,ExitCode,AllocTRES,NodeList",
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    reader = csv.reader(io.StringIO(completed.stdout), delimiter="|")
    records: dict[str, dict[str, Any]] = {}
    for row in reader:
        if len(row) < 7 or row[0] not in roles:
            continue
        job_id, name, state, elapsed, exit_code, alloc_tres, node = row[:7]
        state = state.split()[0].rstrip("+")
        gpu_count = _gpu_count(alloc_tres)
        elapsed_seconds = int(elapsed)
        records[job_id] = {
            "job_id": int(job_id),
            "role": roles[job_id],
            "job_name": name,
            "state": state,
            "elapsed_seconds": elapsed_seconds,
            "exit_code": exit_code,
            "gpu_count": gpu_count,
            "gpu_seconds": elapsed_seconds * gpu_count,
            "alloc_tres": alloc_tres,
            "node": node,
        }

    missing = sorted(set(roles) - set(records))
    if missing:
        raise RuntimeError(f"sacct returned no top-level records for jobs {missing}")
    active = [
        record for record in records.values() if record["state"] not in TERMINAL_STATES
    ]
    if active:
        raise RuntimeError(
            "accounting requires terminal jobs: "
            + ", ".join(f"{row['job_id']}={row['state']}" for row in active)
        )

    jobs = [records[job_id] for job_id in roles]
    total_gpu_seconds = sum(row["gpu_seconds"] for row in jobs)
    return {
        "schema_version": 1,
        "source": "Slurm sacct top-level allocation records",
        "includes_failed_and_superseded_allocations": True,
        "jobs": jobs,
        "total_actual_gpu_seconds": total_gpu_seconds,
        "total_actual_gpu_hours": total_gpu_seconds / 3600.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--job",
        action="append",
        required=True,
        help="role and numeric job ID formatted as ROLE=JOB_ID",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS_DIR / "execution_accounting.json",
    )
    args = parser.parse_args()
    payload = collect(args.job)
    atomic_json(args.output, payload)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
