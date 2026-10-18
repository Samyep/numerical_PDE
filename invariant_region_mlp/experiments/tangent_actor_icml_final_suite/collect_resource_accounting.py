"""Collect scheduler allocation records for the final ICML suite.

The output separates newly charged final-suite jobs from the historical source
experiments whose numerical artifacts are reused.  Failed and cancelled jobs
remain in the ledger; pending jobs are omitted until Slurm reports an
allocation record.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import getpass
import json
from pathlib import Path
import re
import subprocess
from typing import Any

import pandas as pd

from .protocol import RESULTS_DIR


SOURCE_ROOT = Path(__file__).resolve().parents[3]
ORIGINAL_ACCOUNTING = (
    SOURCE_ROOT
    / "invariant_region_mlp/experiments/original_pde_benchmark_comparison"
    / "runtime_resource_accounting.csv"
)
FRONTIER_SLURM = (
    SOURCE_ROOT / "invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm"
)

SACCT_FIELDS = [
    "JobIDRaw", "JobName%64", "Partition", "State", "ElapsedRaw",
    "ExitCode", "AllocTRES%100", "Submit", "Start", "End", "MaxRSS",
]


def _parse_sacct(command: list[str]) -> list[dict[str, str]]:
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    names = [field.split("%", 1)[0] for field in SACCT_FIELDS]
    rows: list[dict[str, str]] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        values = line.split("|")
        values += [""] * (len(names) - len(values))
        rows.append(dict(zip(names, values, strict=False)))
    return rows


def _sacct(*, start: str | None = None, job_ids: list[str] | None = None,
           allocations_only: bool = True) -> list[dict[str, str]]:
    command = ["sacct"]
    if allocations_only:
        command.append("-X")
    if start:
        command.extend(("-S", start, "-u", getpass.getuser()))
    if job_ids:
        command.extend(("-j", ",".join(job_ids)))
    command.extend(("-n", "-P", f"--format={','.join(SACCT_FIELDS)}"))
    return _parse_sacct(command)


def _gpu_count(alloc_tres: str) -> int:
    match = re.search(r"(?:gres/gpu|gpu)=(\d+)", alloc_tres or "")
    return int(match.group(1)) if match else 0


def _record(row: dict[str, Any], *, origin: str, role: str,
            incremental: bool, batch_max_rss: str | None = None) -> dict[str, Any]:
    elapsed = int(float(row.get("ElapsedRaw") or row.get("elapsed_seconds") or 0))
    count = int(row.get("allocated_gpu_count") or _gpu_count(str(row.get("AllocTRES", ""))))
    seconds = elapsed * count
    return {
        "accounting_origin": origin,
        "job_id": str(row.get("JobIDRaw", "")),
        "job_name": row.get("JobName"),
        "partition": row.get("Partition"),
        "state": row.get("State"),
        "elapsed_seconds": elapsed,
        "allocated_gpu_count": count,
        "allocated_gpu_seconds": seconds,
        "allocated_gpu_hours": seconds / 3600.0,
        "max_rss": batch_max_rss or row.get("batch_max_rss") or row.get("MaxRSS"),
        "submit": row.get("Submit"),
        "start": row.get("Start"),
        "end": row.get("End"),
        "exit_code": row.get("ExitCode"),
        "alloc_tres": row.get("AllocTRES"),
        "evidence_role": role,
        "charged_to_incremental_final_suite": incremental,
        "notes": (
            "All scheduler outcomes retained, including failed and cancelled attempts. "
            "Historical rows are supporting source-experiment accounting and are not "
            "charged again as incremental final-suite compute."
        ),
    }


def _frontier_job_ids() -> list[str]:
    ids: set[str] = set()
    if FRONTIER_SLURM.exists():
        for path in FRONTIER_SLURM.iterdir():
            match = re.search(r"_(\d+)\.(?:out|err)$", path.name)
            if match:
                ids.add(match.group(1))
    return sorted(ids, key=int)


def _batch_rss(rows: list[dict[str, str]]) -> dict[str, str]:
    return {
        row["JobIDRaw"].removesuffix(".batch"): row["MaxRSS"]
        for row in rows
        if row["JobIDRaw"].endswith(".batch") and row.get("MaxRSS")
    }


def collect(start: str) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    if ORIGINAL_ACCOUNTING.exists():
        for row in csv.DictReader(ORIGINAL_ACCOUNTING.open()):
            records.append(_record(
                row, origin="historical_original_pde_suite",
                role="source experiment for audited historical reuse",
                incremental=False,
            ))

    frontier_ids = _frontier_job_ids()
    if frontier_ids:
        allocations = _sacct(job_ids=frontier_ids, allocations_only=True)
        rss = _batch_rss(_sacct(job_ids=frontier_ids, allocations_only=False))
        for row in allocations:
            records.append(_record(
                row, origin="historical_accuracy_frontier",
                role="source experiment for audited historical reuse",
                incremental=False, batch_max_rss=rss.get(row["JobIDRaw"]),
            ))

    all_recent = _sacct(start=start, allocations_only=True)
    final_allocations = [row for row in all_recent if row["JobName"].startswith("icml_")]
    recent_steps = _sacct(start=start, allocations_only=False)
    recent_rss = _batch_rss(recent_steps)
    for row in final_allocations:
        records.append(_record(
            row, origin="icml_final_suite_incremental",
            role="new final-suite execution",
            incremental=True, batch_max_rss=recent_rss.get(row["JobIDRaw"]),
        ))

    records.sort(key=lambda row: (row["accounting_origin"], int(row["job_id"] or 0)))
    frame = pd.DataFrame(records)
    output = RESULTS_DIR / "RESOURCE_ACCOUNTING.csv"
    frame.to_csv(output, index=False)

    incremental = [row for row in records if row["charged_to_incremental_final_suite"]]
    historical = [row for row in records if not row["charged_to_incremental_final_suite"]]
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Slurm sacct plus immutable historical accounting CSV",
        "query_start": start,
        "record_count": len(records),
        "incremental_record_count": len(incremental),
        "incremental_allocated_gpu_hours": sum(row["allocated_gpu_hours"] for row in incremental),
        "incremental_completed_gpu_hours": sum(
            row["allocated_gpu_hours"] for row in incremental
            if str(row["state"]).startswith("COMPLETED")
        ),
        "incremental_failed_cancelled_gpu_hours": sum(
            row["allocated_gpu_hours"] for row in incremental
            if not str(row["state"]).startswith(("COMPLETED", "RUNNING", "PENDING"))
        ),
        "historical_source_experiment_gpu_hours": sum(
            row["allocated_gpu_hours"] for row in historical
        ),
        "caveat": (
            "Historical source-experiment hours include development and negative runs; "
            "they are reported for auditability and are not double-counted as incremental cost."
        ),
    }
    (RESULTS_DIR / "RESOURCE_ACCOUNTING_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2026-10-18")
    args = parser.parse_args()
    print(json.dumps(collect(args.start), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
