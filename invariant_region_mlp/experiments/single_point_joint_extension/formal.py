"""Restart-safe formal bundle runner for the five preregistered dimensions."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import time

import torch

from .deep_bsde import execute as deep_execute
from .fh_mlp import execute as fh_execute
from .protocol import (
    DIMENSIONS,
    FH_CONFIGURATIONS,
    FH_REPLICATES,
    QUERY_COUNT,
    STATUS_DIR,
    TRAINING_SEEDS,
    require_smoke,
    write_json_atomic,
)
from .tangent_actor import evaluate as ta_evaluate
from .tangent_actor import train as ta_train


def tasks() -> list[dict]:
    result: list[dict] = []
    for dimension in DIMENSIONS:
        result.append({"method": "tangent_actor", "dimension": dimension})
        result.append({"method": "fh_mlp_standard", "dimension": dimension})
        for seed in TRAINING_SEEDS:
            result.append({
                "method": "deep_bsde", "dimension": dimension, "seed": seed
            })
    return result


def run(index: int, device: torch.device) -> dict:
    require_smoke()
    task_list = tasks()
    if not 0 <= int(index) < len(task_list):
        raise IndexError(index)
    task = task_list[int(index)]
    method = task["method"]
    dimension = int(task["dimension"])
    task_id = f"{index:02d}_{method}_d{dimension}"
    if "seed" in task:
        task_id += f"_s{task['seed']}"
    output = STATUS_DIR / f"{task_id}.json"
    if output.exists():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "complete":
            return existing
    started = time.perf_counter()
    completed = 0
    if method == "tangent_actor":
        ta_train(dimension, 0, device)
        completed = len(
            ta_evaluate(
                dimension,
                0,
                device,
                actor_names=("selected",),
            )
        )
    elif method == "fh_mlp_standard":
        # Heaviest recursion first: a timeout then leaves only cheaper work.
        cells = sorted(
            (
                (n, M, replicate)
                for n, M in FH_CONFIGURATIONS
                for replicate in range(FH_REPLICATES)
            ),
            key=lambda item: (-(item[1] ** item[0]), -item[0], -item[1], item[2]),
        )
        for n, M, replicate in cells:
            fh_execute(dimension, n, M, replicate, device)
            completed += 1
    elif method == "deep_bsde":
        seed = int(task["seed"])
        for point in range(QUERY_COUNT):
            deep_execute(dimension, point, seed, device)
            completed += 1
    else:
        raise KeyError(method)
    result = {
        **task,
        "task_id": task_id,
        "status": "complete",
        "completed_runs": completed,
        "elapsed_seconds": time.perf_counter() - started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(output, result)
    return result


def pending() -> list[int]:
    result: list[int] = []
    for index, task in enumerate(tasks()):
        method, dimension = task["method"], task["dimension"]
        task_id = f"{index:02d}_{method}_d{dimension}"
        if "seed" in task:
            task_id += f"_s{task['seed']}"
        path = STATUS_DIR / f"{task_id}.json"
        if not path.exists() or json.loads(path.read_text(encoding="utf-8")).get("status") != "complete":
            result.append(index)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("list")
    subparsers.add_parser("pending")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--index", type=int, required=True)
    run_parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.action == "list":
        print(json.dumps(tasks(), indent=2))
        return
    if args.action == "pending":
        print(",".join(map(str, pending())))
        return
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    print(json.dumps(run(args.index, device), sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["pending", "run", "tasks"]
