"""Create deterministic, estimate-balanced task queues for multi-GPU waves."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd

from .whole_function_artifacts import build_results_frame
from .whole_function_status import METHODS, RESULTS_DIR
from .whole_function_tasks import _active_suite_jobs


ESTIMATES = RESULTS_DIR / "REMAINING_GPU_HOUR_ESTIMATE.csv"


def generate(stage: str, slots: int = 4) -> tuple[list[dict], dict]:
    active = _active_suite_jobs()
    if active:
        raise RuntimeError(
            "refusing to freeze a queue while suite jobs are active: " + "; ".join(active)
        )
    dimensions = {
        "stage2": (20, 50, 200),
        "stage3": (400,),
        "remaining": (20, 50, 200, 400),
    }[stage]
    frame, _ = build_results_frame()
    pending = frame[
        frame["dimension"].isin(dimensions)
        & (frame["status"] == "PENDING")
        & (~frame["active_slurm_attempt"])
    ].copy()
    estimate_frame = pd.read_csv(ESTIMATES)
    estimate = {}
    for row in estimate_frame.itertuples():
        if row.method not in METHODS:
            continue
        point = min(6.0, float(row.estimated_gpu_hours))
        high = min(6.0, float(getattr(row, "high_gpu_hours", point * 2.2)))
        estimate[(row.pde, int(row.dimension), row.method)] = (point, high)
    method_order = {name: index for index, name in enumerate(METHODS)}
    pde_order = {name: index for index, name in enumerate(frame["pde"].drop_duplicates())}
    queues: list[list[dict]] = [[] for _ in range(slots)]
    loads = [0.0] * slots
    dimension_groups = (dimensions,) if stage == "remaining" else tuple(
        (dimension,) for dimension in dimensions
    )
    for dimension_group in dimension_groups:
        rows = []
        for row in pending[pending["dimension"].isin(dimension_group)].itertuples():
            hours, high_hours = estimate.get(
                (row.pde, int(row.dimension), row.method), (6.0, 6.0),
            )
            rows.append((hours, high_hours, pde_order[row.pde], method_order[row.method], row))
        # Longest-processing-time placement is deterministic and does not use accuracy.
        rows.sort(key=lambda item: (-item[0], item[2], item[3]))
        for hours, high_hours, _, _, row in rows:
            slot = min(range(slots), key=lambda index: (loads[index], index))
            record = {
                "slot": slot,
                "sequence": len(queues[slot]),
                "method": row.method,
                "pde": row.pde,
                "dimension": int(row.dimension),
                "estimated_gpu_hours": hours,
                "high_gpu_hours": high_hours,
            }
            queues[slot].append(record)
            loads[slot] += hours
    # Keep external baselines at the front of each allocation. PINN-HTE is a
    # single nonresumable call; DPI resumes only at exact Picard-round
    # boundaries. This changes dispatch order only; LPT placement and every
    # scientific setting stay frozen.
    baseline_rank = {"dpi": 0, "pinn_hte": 1}
    for queue in queues:
        queue.sort(key=lambda record: (
            0 if record["method"] in baseline_rank else 1,
            baseline_rank.get(record["method"], 0),
            record["sequence"],
        ))
        for sequence, record in enumerate(queue):
            record["sequence"] = sequence
    records = [record for queue in queues for record in queue]
    metadata = {
        "stage": stage,
        "slots": slots,
        "cells": len(records),
        "dimensions": list(dimensions),
        "estimated_active_gpu_hours": sum(loads),
        "estimated_slot_gpu_hours": loads,
        "maximum_estimated_slot_gpu_hours": max(loads, default=0.0),
        "ordering": (
            "deterministic LPT placement using frozen pre-run runtime estimates; within "
            "each slot external DPI/PINN baseline first, followed by stable LPT order; "
            "no accuracy or test metric used"
        ),
    }
    return records, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("stage2", "stage3", "remaining"), required=True,
    )
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
        writer.writeheader(); writer.writerows(records)
    metadata_path = args.output.with_suffix(args.output.suffix + ".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({**metadata, "task_file": str(args.output)}, indent=2, sort_keys=True))
    if not records:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
