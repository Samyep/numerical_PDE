"""Collect Slurm allocation, including failed attempts, for declared jobs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess

from invariant_region_mlp.experiments.whole_function_compute_frontier.common import write_json_atomic
from .protocol import GLOBAL_GPU_SECONDS_CAP, RESULTS_DIR


def _seconds(text: str) -> float:
    days = 0
    if "-" in text:
        day, text = text.split("-", 1)
        days = int(day)
    parts = [int(value) for value in text.split(":")]
    if len(parts) == 2:
        hours, minutes, seconds = 0, parts[0], parts[1]
    else:
        hours, minutes, seconds = parts
    return float(((days * 24 + hours) * 60 + minutes) * 60 + seconds)


def collect(job_ids: list[str]) -> dict:
    rows = []
    for job_id in job_ids:
        output = subprocess.check_output([
            "sacct", "-j", job_id, "--array", "-X", "-n", "-P",
            "--format=JobIDRaw,JobID,State,Elapsed,AllocTRES%100,NodeList,ExitCode",
        ], text=True)
        for line in output.splitlines():
            if not line.strip():
                continue
            raw, display, state, elapsed, tres, nodes, exit_code = line.split("|", 6)
            if "." in display:
                continue
            gpus = 1
            for item in tres.split(","):
                if item.startswith("gres/gpu="):
                    gpus = int(item.split("=", 1)[1])
            seconds = _seconds(elapsed) * gpus
            rows.append({
                "submitted_job": job_id, "job_id_raw": raw, "job_id": display,
                "state": state.strip(), "elapsed": elapsed, "allocated_gpu_count": gpus,
                "allocated_gpu_seconds": seconds, "node_list": nodes,
                "exit_code": exit_code,
            })
    path = RESULTS_DIR / "resource_accounting.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader(); writer.writerows(rows)
    total = sum(float(row["allocated_gpu_seconds"]) for row in rows)
    payload = {
        "job_ids": job_ids, "scheduler_rows": len(rows),
        "allocated_gpu_seconds": total, "allocated_gpu_hours": total / 3600.0,
        "within_90_h100_gpu_hour_cap": total <= GLOBAL_GPU_SECONDS_CAP,
        "failed_or_cancelled_rows": [row for row in rows if row["state"] not in {"COMPLETED"}],
    }
    write_json_atomic(RESULTS_DIR / "resource_accounting.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_ids", nargs="+")
    args = parser.parse_args()
    print(json.dumps(collect(args.job_ids), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
