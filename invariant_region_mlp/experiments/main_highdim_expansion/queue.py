"""Build immutable, accuracy-blind queues for four independent GPUs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
from typing import Any

from .protocol import DIMENSIONS, METHODS, PDES, RESULTS_DIR
from .status import cell_budget


ESTIMATE_PATH = RESULTS_DIR / "runtime_estimates.json"
STAGES = {
    "anchor": lambda method, pde, d: d == 100,
    "actor_scaling": lambda method, pde, d: method == "tangent_actor" and d != 100,
    "baseline_scaling": lambda method, pde, d: method != "tangent_actor" and d != 100,
}


def _active_jobs() -> list[str]:
    output = subprocess.run(
        ["squeue", "-h", "-u", __import__("getpass").getuser(), "-o", "%i|%j|%T"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return [line for line in output.splitlines() if "|icml_exp_" in line]


def _estimates() -> dict[tuple[str, str, int], tuple[float, float]]:
    result: dict[tuple[str, str, int], tuple[float, float]] = {}
    if ESTIMATE_PATH.is_file():
        payload = json.loads(ESTIMATE_PATH.read_text(encoding="utf-8"))
        for row in payload.get("cells", []):
            result[(row["method"], row["pde"], int(row["dimension"]))] = (
                min(6.0, float(row["estimated_gpu_hours"])),
                min(6.0, float(row["high_gpu_hours"])),
            )
    return result


def generate(stage: str, slots: int = 4) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if stage not in STAGES:
        raise KeyError(stage)
    active = _active_jobs()
    if active:
        raise RuntimeError("refusing to freeze queue while expansion jobs are active: " + "; ".join(active))
    estimates = _estimates()
    tasks: list[tuple[float, float, int, int, int, str, str]] = []
    for pde_index, pde in enumerate(PDES):
        for d in DIMENSIONS:
            for method_index, method in enumerate(METHODS):
                if not STAGES[stage](method, pde, d):
                    continue
                status = cell_budget(method, pde, d)["status"]
                if status not in {"PENDING", "CHECKPOINTED"}:
                    continue
                point, high = estimates.get(
                    (method, pde, d),
                    ((0.75 if method == "tangent_actor" else 2.0), (1.5 if method == "tangent_actor" else 6.0)),
                )
                tasks.append((point, high, pde_index, method_index, d, method, pde))
    tasks.sort(key=lambda row: (-row[0], row[2], row[3], row[4]))
    queues: list[list[dict[str, Any]]] = [[] for _ in range(slots)]
    loads = [0.0] * slots
    for point, high, _p, _m, d, method, pde in tasks:
        slot = min(range(slots), key=lambda index: (loads[index], index))
        queues[slot].append(
            {
                "slot": slot,
                "sequence": len(queues[slot]),
                "method": method,
                "pde": pde,
                "dimension": d,
                "estimated_gpu_hours": point,
                "high_gpu_hours": high,
            }
        )
        loads[slot] += point
    nonresumable = {"dpi": 0, "pinn_hte": 1}
    for queue in queues:
        queue.sort(
            key=lambda row: (
                0 if row["method"] in nonresumable else 1,
                nonresumable.get(row["method"], 0),
                row["sequence"],
            )
        )
        for sequence, row in enumerate(queue):
            row["sequence"] = sequence
    records = [row for queue in queues for row in queue]
    metadata = {
        "stage": stage,
        "slots": slots,
        "cells": len(records),
        "estimated_gpu_hours": sum(loads),
        "slot_estimated_gpu_hours": loads,
        "ordering": "frozen runtime LPT; no validation/test accuracy enters dispatch",
    }
    return records, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=tuple(STAGES), required=True)
    parser.add_argument("--slots", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records, metadata = generate(args.stage, args.slots)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "slot", "sequence", "method", "pde", "dimension",
                "estimated_gpu_hours", "high_gpu_hours",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)
    args.output.with_suffix(args.output.suffix + ".json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({**metadata, "task_file": str(args.output)}, indent=2, sort_keys=True))
    if not records:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
