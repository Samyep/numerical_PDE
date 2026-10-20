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
CAP_SECONDS = 120 * 3600


def _gpu_count(value: str) -> int:
    match = re.search(r"(?:gres/gpu|gpu)=(\d+)", value or "")
    return int(match.group(1)) if match else 0


def _billing_units(value: str) -> int:
    match = re.search(r"(?:^|,)billing=(\d+)(?:,|$)", value or "")
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
        billing = _billing_units(str(raw.get("AllocTRES") or ""))
        rows.append({
            "job_id": raw.get("JobIDRaw"), "job_name": raw.get("JobName"),
            "partition": raw.get("Partition"), "state": raw.get("State"),
            "elapsed_seconds": elapsed, "allocated_gpu_count": gpus,
            "allocated_gpu_seconds": elapsed * gpus,
            "allocated_gpu_hours": elapsed * gpus / 3600.0,
            "allocated_billing_units": billing,
            "account_credit_hours": elapsed * billing / (3600.0 * 1000.0),
            "exit_code": raw.get("ExitCode"), "alloc_tres": raw.get("AllocTRES"),
            "submit": raw.get("Submit"), "start": raw.get("Start"), "end": raw.get("End"),
        })
    rows.sort(key=lambda item: str(item["job_id"]))
    path = RESULTS_DIR / "WHOLE_FUNCTION_JOB_ACCOUNTING.csv"
    columns = list(rows[0]) if rows else [
        "job_id", "job_name", "partition", "state", "elapsed_seconds",
        "allocated_gpu_count", "allocated_gpu_seconds", "allocated_gpu_hours",
        "allocated_billing_units", "account_credit_hours", "exit_code",
        "alloc_tres", "submit", "start", "end",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    allocated = sum(float(item["allocated_gpu_seconds"]) for item in rows)
    credits = sum(float(item["account_credit_hours"]) for item in rows)
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "query_start": start, "job_name_prefix": "icml_whole",
        "jobs": len(rows), "allocated_gpu_seconds": allocated,
        "allocated_gpu_hours": allocated / 3600.0,
        "account_credit_hours": credits,
        "cap_gpu_hours": 120.0,
        "within_cap": allocated <= 120 * 3600,
        "source": (
            "Slurm sacct allocation records; includes failed and timed-out jobs. "
            "Account credit hours are elapsed hours times AllocTRES billing units / 1000."
        ),
    }
    (RESULTS_DIR / "WHOLE_FUNCTION_JOB_ACCOUNTING.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    return payload


def _duration_seconds(value: str) -> int:
    days = 0
    if "-" in value:
        day, value = value.split("-", 1)
        days = int(day)
    fields = [int(item) for item in value.split(":")]
    if len(fields) == 3:
        hours, minutes, seconds = fields
    elif len(fields) == 2:
        hours, minutes, seconds = 0, fields[0], fields[1]
    else:
        raise ValueError(f"unsupported Slurm duration {value!r}")
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def verify_reservation(job_id: str) -> dict[str, Any]:
    """Fail closed unless this allocation's full GPU request fits the suite cap."""
    collect()
    accounting_path = RESULTS_DIR / "WHOLE_FUNCTION_JOB_ACCOUNTING.csv"
    with accounting_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    prior = sum(
        float(row.get("allocated_gpu_seconds") or 0.0)
        for row in rows if str(row.get("job_id")) != str(job_id)
    )
    info = subprocess.run(
        ["scontrol", "show", "job", "-o", str(job_id)],
        check=True, capture_output=True, text=True,
    ).stdout
    time_match = re.search(r"\bTimeLimit=(\S+)", info)
    tres_match = re.search(r"\bReqTRES=([^ ]+)", info)
    if not time_match or not tres_match:
        raise RuntimeError("could not read TimeLimit and ReqTRES for budget reservation")
    gpus = _gpu_count(tres_match.group(1))
    if gpus <= 0:
        raise RuntimeError("allocation has no GPU in ReqTRES")
    wall_seconds = _duration_seconds(time_match.group(1))
    reserved = wall_seconds * gpus
    passed = prior + reserved <= CAP_SECONDS
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "job_id": str(job_id), "gpu_count": gpus,
        "wall_seconds": wall_seconds, "reserved_gpu_seconds": reserved,
        "prior_allocated_gpu_seconds": prior,
        "cap_gpu_seconds": CAP_SECONDS,
        "status": "PASS" if passed else "FAIL",
    }
    directory = RESULTS_DIR / "whole_function_dispatch_reservations"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"job_{job_id}.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    if not passed:
        raise RuntimeError(
            f"job reservation exceeds suite cap: prior={prior}, reserved={reserved}, cap={CAP_SECONDS}"
        )
    return payload


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reserve-job-id")
    args = parser.parse_args()
    payload = verify_reservation(args.reserve_job_id) if args.reserve_job_id else collect()
    print(json.dumps(payload, indent=2, sort_keys=True))
