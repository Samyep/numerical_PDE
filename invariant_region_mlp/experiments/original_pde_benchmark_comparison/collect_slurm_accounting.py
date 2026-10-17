"""Capture auditable Slurm allocation records for this experiment."""

from __future__ import annotations

import argparse
import getpass
import json
import re
import subprocess
from typing import Any

from .common import RESULTS_DIR, atomic_json


def _sacct(*, start: str, allocations_only: bool) -> list[dict[str, str]]:
    fields = [
        "JobIDRaw",
        "JobName%60",
        "Partition",
        "State",
        "ElapsedRaw",
        "ExitCode",
        "AllocTRES%100",
        "Submit",
        "Start",
        "End",
        "MaxRSS",
    ]
    command = [
        "sacct",
        "-S",
        start,
        "-u",
        getpass.getuser(),
        "-n",
        "-P",
        f"--format={','.join(fields)}",
    ]
    if allocations_only:
        command.insert(1, "-X")
    output = subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout
    names = [field.split("%", 1)[0] for field in fields]
    rows: list[dict[str, str]] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        values = line.split("|")
        if len(values) < len(names):
            values.extend([""] * (len(names) - len(values)))
        rows.append(dict(zip(names, values, strict=False)))
    return rows


def collect(start: str) -> dict[str, Any]:
    allocations = [
        row for row in _sacct(start=start, allocations_only=True)
        if row["JobName"].startswith("orig_")
    ]
    steps = _sacct(start=start, allocations_only=False)
    batch_rss = {
        row["JobIDRaw"].removesuffix(".batch"): row["MaxRSS"]
        for row in steps
        if row["JobIDRaw"].endswith(".batch")
    }
    records: list[dict[str, Any]] = []
    for row in allocations:
        tres = row["AllocTRES"]
        match = re.search(r"(?:gres/gpu|gpu)=(\d+)", tres)
        gpu_count = int(match.group(1)) if match else 0
        elapsed = int(row["ElapsedRaw"] or 0)
        job_id = row["JobIDRaw"]
        records.append(
            {
                **row,
                "elapsed_seconds": elapsed,
                "allocated_gpu_count": gpu_count,
                "allocated_gpu_seconds": elapsed * gpu_count,
                "batch_max_rss": batch_rss.get(job_id),
            }
        )
    result = {
        "source": "Slurm sacct",
        "query_start": start,
        "user": getpass.getuser(),
        "records": records,
        "total_completed_gpu_seconds": sum(
            row["allocated_gpu_seconds"]
            for row in records
            if str(row["State"]).startswith("COMPLETED")
        ),
        "total_failed_or_cancelled_gpu_seconds": sum(
            row["allocated_gpu_seconds"]
            for row in records
            if not str(row["State"]).startswith(("COMPLETED", "RUNNING", "PENDING"))
        ),
    }
    atomic_json(RESULTS_DIR / "slurm_accounting.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2026-10-12")
    args = parser.parse_args()
    result = collect(args.start)
    print(json.dumps({"jobs": len(result["records"]), **{
        key: result[key] for key in (
            "total_completed_gpu_seconds",
            "total_failed_or_cancelled_gpu_seconds",
        )
    }}, sort_keys=True))


if __name__ == "__main__":
    main()
