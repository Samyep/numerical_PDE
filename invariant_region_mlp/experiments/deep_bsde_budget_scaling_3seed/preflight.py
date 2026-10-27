"""Truth-free H100 timing/memory preflight and allocation gate."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
import time

import torch

from .protocol import (
    FINAL_COUNT,
    MAX_H100_GPU_HOURS,
    PDES,
    PREFLIGHT_CONTINGENCY,
    PREFLIGHT_UPDATES,
    RESULTS_DIR,
    TRAINING_SEEDS,
    TRAJECTORIES,
    file_sha256,
    write_json_atomic,
)
from .solver import load_input, require_h100, train_query


def run(device: torch.device) -> dict:
    gpu_name = require_h100(device)
    rows = []
    started = time.perf_counter()
    for pde in PDES:
        original_index, t_start, x_start = load_input(pde, "calibration", 0)
        for trajectory in TRAJECTORIES:
            record = train_query(
                pde, trajectory, 991, 0, original_index, t_start, x_start,
                device, preflight_updates=PREFLIGHT_UPDATES,
            )
            if record["status"] != "complete":
                raise RuntimeError(f"preflight numerical failure: {pde} {trajectory.trajectory_id}")
            seconds_per_update = record["active_training_seconds"] / PREFLIGHT_UPDATES
            estimated_query_seconds = record["initialization_seconds"] + seconds_per_update * trajectory.updates
            rows.append({
                "pde": pde,
                "trajectory_id": trajectory.trajectory_id,
                "trajectory_signature": trajectory.signature,
                "preflight_updates": PREFLIGHT_UPDATES,
                "initialization_seconds": record["initialization_seconds"],
                "active_training_seconds": record["active_training_seconds"],
                "seconds_per_update": seconds_per_update,
                "official_updates": trajectory.updates,
                "estimated_seconds_per_query": estimated_query_seconds,
                "projected_final_seconds_all_seeds": estimated_query_seconds * FINAL_COUNT * len(TRAINING_SEEDS),
                "peak_gpu_bytes": record["peak_gpu_bytes"],
            })
    measured_preflight_seconds = time.perf_counter() - started
    projected_formal_seconds = sum(row["projected_final_seconds_all_seeds"] for row in rows)
    projected_with_contingency = projected_formal_seconds * PREFLIGHT_CONTINGENCY + measured_preflight_seconds
    projected_hours = projected_with_contingency / 3600.0
    passed = projected_hours <= MAX_H100_GPU_HOURS
    payload = {
        "status": "ALLOCATION_GATE_PASS" if passed else "ALLOCATION_GATE_FAIL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gpu_name": gpu_name,
        "hostname": os.uname().nodename,
        "preflight_truth_loaded": False,
        "preflight_updates": PREFLIGHT_UPDATES,
        "measured_preflight_seconds": measured_preflight_seconds,
        "projected_formal_active_gpu_seconds": projected_formal_seconds,
        "contingency_multiplier": PREFLIGHT_CONTINGENCY,
        "projected_total_gpu_hours_with_contingency": projected_hours,
        "maximum_h100_gpu_hours": MAX_H100_GPU_HOURS,
        "formal_submission_authorized": passed,
        "rows": rows,
        "query_manifest_sha256": file_sha256(RESULTS_DIR / "query_manifest.json"),
        "task_manifest_sha256": file_sha256(RESULTS_DIR / "task_manifest.json"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    write_json_atomic(RESULTS_DIR / "preflight_allocation_gate.json", payload)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = run(torch.device(args.device))
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["formal_submission_authorized"]:
        raise SystemExit(2)

