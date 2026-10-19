"""Generate budget-safe task waves for the whole-function Slurm array."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess

from .whole_function_artifacts import build_results_frame
from .whole_function_status import DIMENSIONS, METHODS, RESULTS_DIR
from .whole_function_job_accounting import collect as collect_job_accounting


def _active_suite_jobs() -> list[str]:
    """Return active suite jobs so two budget waves cannot overlap."""
    result = subprocess.run(
        ["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%i|%j|%T"],
        check=False, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"could not query active Slurm jobs: {result.stderr.strip()}")
    return [
        line for line in result.stdout.splitlines()
        if "|icml_whole" in line
    ]


def make_tasks(stage: str, requested: int) -> tuple[list[tuple[str, str, int]], dict]:
    active_jobs = _active_suite_jobs()
    if active_jobs:
        raise RuntimeError(
            "refusing to generate a new budget wave while suite jobs are active: "
            + "; ".join(active_jobs)
        )
    frame, attempt_rows = build_results_frame()
    accounting = collect_job_accounting()
    used_seconds = float(accounting["allocated_gpu_seconds"])
    remaining_seconds = max(0.0, 120 * 3600 - used_seconds)
    guaranteed_slots = int(math.floor(remaining_seconds / (6 * 3600)))
    dimensions = {
        "stage2": {20, 50, 200},
        "stage3": {400},
        "all": set(DIMENSIONS) - {100},
    }[stage]
    pending = frame[
        frame["dimension"].isin(dimensions)
        & (frame["status"] == "PENDING")
        & (~frame["active_slurm_attempt"])
    ].copy()
    method_rank = {method: index for index, method in enumerate(METHODS)}
    pde_rank = {pde: index for index, pde in enumerate(frame["pde"].drop_duplicates())}
    pending["_method"] = pending["method"].map(method_rank)
    pending["_pde"] = pending["pde"].map(pde_rank)
    pending = pending.sort_values(["dimension", "_pde", "_method"])
    limit = min(requested, guaranteed_slots, len(pending))
    selected = [
        (row.method, row.pde, int(row.dimension))
        for row in pending.head(limit).itertuples()
    ]
    metadata = {
        "stage": stage,
        "requested_cells": requested,
        "selected_cells": len(selected),
        "pending_eligible_cells": len(pending),
        "used_new_gpu_hours": used_seconds / 3600.0,
        "remaining_gpu_hours": remaining_seconds / 3600.0,
        "guaranteed_six_hour_slots": guaranteed_slots,
        "selection_rule": "dimension, PDE, frozen method order; no accuracy input",
    }
    return selected, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("stage2", "stage3", "all"), required=True)
    parser.add_argument("--max-cells", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tasks, metadata = make_tasks(args.stage, args.max_cells)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(f"{m},{p},{d}\n" for m, p, d in tasks), encoding="utf-8")
    metadata_path = args.output.with_suffix(args.output.suffix + ".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({**metadata, "task_file": str(args.output)}, indent=2, sort_keys=True))
    if not tasks:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
