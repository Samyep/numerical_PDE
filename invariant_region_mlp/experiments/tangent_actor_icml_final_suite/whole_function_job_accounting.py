"""Collect authoritative Slurm allocation cost for ``icml_whole*`` jobs."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import getpass
import json
from pathlib import Path
import re
import subprocess
from typing import Any

from .protocol import RESULTS_DIR


FIELDS = ("JobIDRaw", "JobName%64", "Partition", "State", "ElapsedRaw", "ExitCode", "AllocTRES%100", "Submit", "Start", "End")


def _gpu_count(value: str) -> int:
    match = re.search(r"(?:gres/gpu|gpu)=(\d+)", value or "")
    return int(match.group(1)) if match else 0


def collect(start: str = "2026-10-19") -> dict[str, Any]:
    command = [
        "sacct", "-X", "-S", start, "-u", getpass.getuser(), "-n", "-P",
        "--format=" + ",".join(FIELDS),
    ]
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    names = [item.split("%", 1)[0] for item in FIELDS]
    rows = []
    for line in output.splitlines():
        values = line.split("|")
        raw = dict(zip(names, values, strict=False))
        if not str(raw.get("JobName", "")).startswith("icml_whole"):
            continue
        elapsed = int(float(raw.get("ElapsedRaw") or 0))
        gpus = _gpu_count(str(raw.get("AllocTRES") or ""))
        rows.append({
            "job_id": raw.get("JobIDRaw"), "job_name": raw.get("JobName"),
            "partition": raw.get("Partition"), "state": raw.get("State"),
            "elapsed_seconds": elapsed, "allocated_gpu_count": gpus,
            "allocated_gpu_seconds": elapsed * gpus,
            "allocated_gpu_hours": elapsed * gpus / 3600.0,
            "exit_code": raw.get("ExitCode"), "alloc_tres": raw.get("AllocTRES"),
            "submit": raw.get("Submit"), "start": raw.get("Start"), "end": raw.get("End"),
        })
    rows.sort(key=lambda item: str(item["job_id"]))
    path = RESULTS_DIR / "WHOLE_FUNCTION_JOB_ACCOUNTING.csv"
    columns = list(rows[0]) if rows else [
        "job_id", "job_name", "partition", "state", "elapsed_seconds",
        "allocated_gpu_count", "allocated_gpu_seconds", "allocated_gpu_hours",
        "exit_code", "alloc_tres", "submit", "start", "end",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    allocated = sum(float(item["allocated_gpu_seconds"]) for item in rows)
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "query_start": start, "job_name_prefix": "icml_whole",
        "jobs": len(rows), "allocated_gpu_seconds": allocated,
        "allocated_gpu_hours": allocated / 3600.0,
        "cap_gpu_hours": 120.0,
        "within_cap": allocated <= 120 * 3600,
        "source": "Slurm sacct allocation records; includes failed and timed-out jobs",
    }
    (RESULTS_DIR / "WHOLE_FUNCTION_JOB_ACCOUNTING.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(collect(), indent=2, sort_keys=True))
