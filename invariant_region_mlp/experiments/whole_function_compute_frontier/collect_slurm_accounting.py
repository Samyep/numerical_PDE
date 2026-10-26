"""Collect scheduler accounting for every wcf_* job, including failures."""

from __future__ import annotations

import argparse
import csv
import getpass
import re
import subprocess

from .protocol import RESULTS_DIR


def collect(start: str) -> int:
    # Job names are experiment-specific, so widening an ISO timestamp to its
    # calendar day safely includes preflight failures/retries that preceded the
    # formal dependency chain without capturing unrelated work.
    start = start.split("T", 1)[0]
    command = [
        "sacct", "-u", getpass.getuser(), "-S", start, "-X", "-n", "-P",
        "-o", "JobIDRaw,JobName,Partition,Account,QOS,AllocTRES,Elapsed,State,ExitCode",
    ]
    output = subprocess.check_output(command, text=True)
    rows = []
    for line in output.splitlines():
        values = line.split("|")
        if len(values) < 9 or not values[1].startswith("wcf_"):
            continue
        match = re.search(r"gres/gpu(?::h100)?=(\d+)", values[5])
        rows.append({
            "job_id": values[0], "job_name": values[1], "partition": values[2],
            "account": values[3], "qos": values[4], "alloc_tres": values[5],
            "elapsed": values[6], "state": values[7], "exit_code": values[8],
            "gpu_count": int(match.group(1)) if match else 0,
        })
    path = RESULTS_DIR / "slurm_accounting_raw.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [
            "job_id", "job_name", "partition", "account", "qos", "alloc_tres",
            "elapsed", "state", "exit_code", "gpu_count",
        ])
        writer.writeheader(); writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", required=True, help="Slurm -S date/time")
    args = parser.parse_args()
    print(f"collected {collect(args.start)} scheduler rows", flush=True)


if __name__ == "__main__":
    main()
