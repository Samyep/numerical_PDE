"""Run one resumable bundle: 64 queries for one PDE/seed/trajectory."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import traceback

import torch

from .protocol import (
    RESULTS_DIR,
    STATUS_DIR,
    TRAJECTORY_BY_ID,
    file_sha256,
    task_specs,
    write_json_atomic,
)
from .solver import load_input, require_h100, train_query


def _is_cuda_infrastructure_failure(error: BaseException) -> bool:
    """Identify node/device failures that make every later query invalid."""

    message = f"{type(error).__name__}: {error}".lower()
    return any(marker in message for marker in (
        "cuda-capable device(s) is/are busy or unavailable",
        "cuda initialization: cuda unknown error",
        "cuda unknown error",
        "no cuda gpus are available",
    ))


def run_task(index: int, device: torch.device) -> dict:
    gate_path = RESULTS_DIR / "preflight_allocation_gate.json"
    if not gate_path.is_file():
        raise RuntimeError("preflight allocation gate is missing")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if not gate.get("formal_submission_authorized", False):
        raise RuntimeError("preflight allocation gate did not pass")
    require_h100(device)
    tasks = task_specs()
    task = tasks[int(index)]
    trajectory = TRAJECTORY_BY_ID[task["trajectory_id"]]
    successes = 0
    failures = 0
    for position in range(task["query_count"]):
        original_index, t_start, x_start = load_input(task["pde"], "final", position)
        try:
            record = train_query(
                task["pde"], trajectory, task["seed"], position,
                original_index, t_start, x_start, device,
            )
            if record["status"] == "complete":
                successes += 1
            else:
                failures += 1
        except Exception as error:
            failures += 1
            fatal_cuda_infrastructure = _is_cuda_infrastructure_failure(error)
            failure_path = (
                RESULTS_DIR / "raw" / "formal" / task["pde"]
                / task["trajectory_id"] / f"s{task['seed']}"
                / f"q{position:02d}.infrastructure_failure.json"
            )
            write_json_atomic(failure_path, {
                "status": "infrastructure_failure",
                "task": task,
                "position": position,
                "error_type": type(error).__name__,
                "error": repr(error),
                "traceback": traceback.format_exc(),
                "fatal_cuda_infrastructure": fatal_cuda_infrastructure,
            })
            # A broken or unavailable CUDA device poisons the allocation, not
            # the query.  Fail the Slurm task immediately so it can be retried
            # on a healthy node; do not emit 63 redundant failure markers.
            if fatal_cuda_infrastructure:
                raise
    payload = {
        "status": "complete" if successes + failures == task["query_count"] else "incomplete",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "task": task,
        "successful_queries": successes,
        "failed_queries": failures,
        "allocation_gate_sha256": file_sha256(gate_path),
    }
    write_json_atomic(STATUS_DIR / f"task_{index:03d}.json", payload)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    print(json.dumps(run_task(args.index, torch.device(args.device)), indent=2, sort_keys=True))

