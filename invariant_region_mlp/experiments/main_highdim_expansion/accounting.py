"""Authoritative Slurm allocation accounting for the expansion suite."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import getpass
import json
import re
import subprocess
from typing import Any

from .protocol import GLOBAL_GPU_SECONDS_CAP, RESULTS_DIR, write_json_atomic


FIELDS = (
    "JobIDRaw", "JobName%64", "Partition", "State", "ElapsedRaw", "ExitCode",
    "AllocTRES%100", "Submit", "Start", "End",
)


def _gpu_count(value: str) -> int:
    match = re.search(r"(?:gres/gpu|gpu)=(\d+)", value or "")
    return int(match.group(1)) if match else 0


def _billing_units(value: str) -> int:
    match = re.search(r"(?:^|,)billing=(\d+)(?:,|$)", value or "")
    return int(match.group(1)) if match else 0


def collect(start: str = "2026-10-20") -> dict[str, Any]:
    command = [
        "sacct", "-X", "-S", start, "-u", getpass.getuser(), "-n", "-P",
        "--format=" + ",".join(FIELDS),
    ]
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    names = [item.split("%", 1)[0] for item in FIELDS]
    rows: list[dict[str, Any]] = []
    for line in output.splitlines():
        raw = dict(zip(names, line.split("|"), strict=False))
        if not str(raw.get("JobName", "")).startswith("icml_exp_"):
            continue
        elapsed = int(float(raw.get("ElapsedRaw") or 0))
        gpus = _gpu_count(str(raw.get("AllocTRES") or ""))
        billing = _billing_units(str(raw.get("AllocTRES") or ""))
        rows.append(
            {
                "job_id": raw.get("JobIDRaw"),
                "job_name": raw.get("JobName"),
                "partition": raw.get("Partition"),
                "state": raw.get("State"),
                "elapsed_seconds": elapsed,
                "allocated_gpu_count": gpus,
                "allocated_gpu_seconds": elapsed * gpus,
                "allocated_gpu_hours": elapsed * gpus / 3600.0,
                "allocated_billing_units": billing,
                "account_credit_hours": elapsed * billing / 3_600_000.0,
                "exit_code": raw.get("ExitCode"),
                "alloc_tres": raw.get("AllocTRES"),
                "submit": raw.get("Submit"),
                "start": raw.get("Start"),
                "end": raw.get("End"),
            }
        )
    rows.sort(key=lambda item: str(item["job_id"]))
    path = RESULTS_DIR / "resource_accounting.csv"
    columns = list(rows[0]) if rows else [
        "job_id", "job_name", "partition", "state", "elapsed_seconds",
        "allocated_gpu_count", "allocated_gpu_seconds", "allocated_gpu_hours",
        "allocated_billing_units", "account_credit_hours", "exit_code", "alloc_tres",
        "submit", "start", "end",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    allocated = sum(float(row["allocated_gpu_seconds"]) for row in rows)
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "query_start": start,
        "job_name_prefix": "icml_exp_",
        "jobs": len(rows),
        "allocated_gpu_seconds": allocated,
        "allocated_gpu_hours": allocated / 3600.0,
        "account_credit_hours": sum(float(row["account_credit_hours"]) for row in rows),
        "cap_gpu_hours": GLOBAL_GPU_SECONDS_CAP / 3600.0,
        "within_cap": allocated <= GLOBAL_GPU_SECONDS_CAP,
        "source": "Slurm sacct allocations, including failed and timed-out expansion jobs",
    }
    write_json_atomic(RESULTS_DIR / "resource_accounting.json", result)
    return result


def _duration_seconds(value: str) -> int:
    days = 0
    if "-" in value:
        prefix, value = value.split("-", 1)
        days = int(prefix)
    fields = [int(item) for item in value.split(":")]
    if len(fields) == 3:
        hours, minutes, seconds = fields
    elif len(fields) == 2:
        hours, minutes, seconds = 0, fields[0], fields[1]
    else:
        raise ValueError(value)
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def verify_reservation(job_id: str) -> dict[str, Any]:
    collect()
    with (RESULTS_DIR / "resource_accounting.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    prior = sum(
        float(row.get("allocated_gpu_seconds") or 0.0)
        for row in rows
        if str(row.get("job_id")) != str(job_id)
    )
    info = subprocess.run(
        ["scontrol", "show", "job", "-o", str(job_id)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    time_match = re.search(r"\bTimeLimit=(\S+)", info)
    tres_match = re.search(r"\bReqTRES=([^ ]+)", info)
    if not time_match or not tres_match:
        raise RuntimeError("cannot inspect Slurm reservation")
    gpus = _gpu_count(tres_match.group(1))
    wall_seconds = _duration_seconds(time_match.group(1))
    reserved = gpus * wall_seconds
    passed = gpus > 0 and prior + reserved <= GLOBAL_GPU_SECONDS_CAP
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "job_id": str(job_id),
        "gpu_count": gpus,
        "wall_seconds": wall_seconds,
        "reserved_gpu_seconds": reserved,
        "prior_allocated_gpu_seconds": prior,
        "cap_gpu_seconds": GLOBAL_GPU_SECONDS_CAP,
        "status": "PASS" if passed else "FAIL",
    }
    directory = RESULTS_DIR / "dispatch_reservations"
    directory.mkdir(parents=True, exist_ok=True)
    write_json_atomic(directory / f"job_{job_id}.json", result)
    if not passed:
        raise RuntimeError(f"expansion reservation exceeds cap: {result}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reserve-job-id")
    args = parser.parse_args()
    result = verify_reservation(args.reserve_job_id) if args.reserve_job_id else collect()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
