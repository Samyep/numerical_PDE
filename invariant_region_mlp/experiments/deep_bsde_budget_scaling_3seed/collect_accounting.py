"""Collect allocation-level Slurm accounting, including failures and retries."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import re
import subprocess
from typing import Any

from .protocol import MAX_H100_GPU_HOURS, RESULTS_DIR, write_json_atomic


def _gpu_count(tres: str) -> int:
    matches = re.findall(r"(?:gres/)?gpu(?::[^=,]+)?=(\d+)", tres or "")
    return max((int(value) for value in matches), default=0)


def collect() -> dict[str, Any]:
    submitted = json.loads((RESULTS_DIR / "submitted_jobs.json").read_text(encoding="utf-8"))
    job_ids = [str(value) for value in submitted.get("all_job_ids", [])]
    if not job_ids:
        job_ids = [
            str(value) for key, value in submitted.items()
            if key.endswith("job_id") and value
        ]
    command = [
        "sacct", "-X", "-j", ",".join(job_ids), "-P", "-n",
        "--format=JobIDRaw,JobName,State,ElapsedRaw,AllocTRES,ReqTRES,ExitCode,Start,End",
    ]
    output = subprocess.check_output(command, text=True)
    fields = ["job_id_raw", "job_name", "state", "elapsed_seconds", "alloc_tres", "req_tres", "exit_code", "start", "end"]
    rows = []
    for line in output.splitlines():
        values = line.split("|")
        if len(values) < 9:
            continue
        raw = dict(zip(fields, values[:9]))
        if "." in raw["job_id_raw"]:
            continue
        elapsed = int(raw["elapsed_seconds"] or 0)
        gpus = _gpu_count(raw["alloc_tres"]) or _gpu_count(raw["req_tres"])
        raw["gpu_count"] = gpus
        raw["allocated_gpu_seconds"] = elapsed * gpus
        rows.append(raw)
    csv_path = RESULTS_DIR / "slurm_resource_accounting.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else fields + ["gpu_count", "allocated_gpu_seconds"])
        writer.writeheader()
        writer.writerows(rows)
    total = sum(int(row["allocated_gpu_seconds"]) for row in rows)
    payload = {
        "status": "PASS" if total / 3600.0 <= MAX_H100_GPU_HOURS else "CAP_EXCEEDED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "submitted_job_ids": job_ids,
        "allocation_rows": len(rows),
        "failed_or_cancelled_rows": [
            row for row in rows if not str(row["state"]).startswith("COMPLETED")
        ],
        "total_allocated_h100_gpu_seconds": total,
        "total_allocated_h100_gpu_hours": total / 3600.0,
        "maximum_h100_gpu_hours": MAX_H100_GPU_HOURS,
        "rows": rows,
    }
    write_json_atomic(RESULTS_DIR / "resource_accounting.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(collect(), indent=2, sort_keys=True))

