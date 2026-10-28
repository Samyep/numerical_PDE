"""Record Slurm submissions and collect allocation-level H100 accounting."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import re
import subprocess

from .protocol import (
    ALLOCATED_GPU_HOUR_CAP,
    RESULTS_DIR,
    REVISION,
    write_json_atomic,
)


SUBMITTED = RESULTS_DIR / "submitted_jobs.json"


def record(stage: str, job_id: str) -> dict:
    payload = (
        json.loads(SUBMITTED.read_text(encoding="utf-8"))
        if SUBMITTED.is_file()
        else {"revision": REVISION, "jobs": []}
    )
    item = {
        "stage": stage,
        "job_id": str(job_id),
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
    }
    if not any(row["job_id"] == item["job_id"] for row in payload["jobs"]):
        payload["jobs"].append(item)
    write_json_atomic(SUBMITTED, payload)
    return payload


def _gpu_count(tres: str) -> int:
    match = re.search(r"gres/gpu(?:[:][^,=]+)?=(\d+)", tres)
    return int(match.group(1)) if match else 0


def _sacct(job_id: str) -> list[dict]:
    command = [
        "sacct", "-X", "-n", "-P", "-j", str(job_id),
        "--format=JobIDRaw,State,ElapsedRaw,AllocTRES,NodeList,ExitCode",
    ]
    output = subprocess.check_output(command, text=True)
    parsed = []
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split("|")
        if len(parts) < 6:
            continue
        raw, state, elapsed, tres, node, exit_code = parts[:6]
        if "." in raw:
            continue
        parsed.append({
            "job_id_raw": raw,
            "state": state,
            "elapsed_raw_seconds": int(elapsed or 0),
            "alloc_tres": tres,
            "node": node,
            "exit_code": exit_code,
            "gpu_count": _gpu_count(tres),
        })
    children = [row for row in parsed if row["job_id_raw"].startswith(f"{job_id}_")]
    if children:
        return children
    return [row for row in parsed if row["job_id_raw"] == str(job_id)]


def collect() -> dict:
    submitted = json.loads(SUBMITTED.read_text(encoding="utf-8"))
    rows = []
    for job in submitted["jobs"]:
        for allocation in _sacct(job["job_id"]):
            allocated = (
                allocation["elapsed_raw_seconds"] * allocation["gpu_count"] / 3600.0
            )
            rows.append({
                "stage": job["stage"],
                "submitted_job_id": job["job_id"],
                **allocation,
                "allocated_h100_gpu_hours": allocated,
            })
    fields = list(rows[0]) if rows else []
    csv_path = RESULTS_DIR / "resource_accounting.csv"
    temporary = csv_path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(csv_path)
    total = sum(float(row["allocated_h100_gpu_hours"]) for row in rows)
    payload = {
        "revision": REVISION,
        "status": "PASS" if total <= ALLOCATED_GPU_HOUR_CAP else "CAP_EXCEEDED",
        "total_allocated_h100_gpu_hours": total,
        "cap_h100_gpu_hours": ALLOCATED_GPU_HOUR_CAP,
        "allocation_count": len(rows),
        "rows": rows,
        "collected_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RESULTS_DIR / "resource_accounting.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    item = sub.add_parser("record")
    item.add_argument("--stage", required=True)
    item.add_argument("--job-id", required=True)
    sub.add_parser("collect")
    args = parser.parse_args()
    result = record(args.stage, args.job_id) if args.action == "record" else collect()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
