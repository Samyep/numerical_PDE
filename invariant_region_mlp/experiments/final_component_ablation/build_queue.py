"""Build the immutable, load-balanced four-GPU formal execution queue."""

from __future__ import annotations

import csv
from pathlib import Path

from .protocol import PACKAGE_DIR, PDES, TRAINING_SEEDS


ESTIMATE = {
    "hjb_log": 0.38,
    "joint_concave": 0.60,
    "joint_convex": 0.63,
    "nonconvex_coupled": 0.92,
}


def main() -> None:
    cells: list[dict[str, object]] = []
    for variant in ("a0_full", "a1_no_mcv", "a2_no_terminal_anchor"):
        for pde in PDES:
            for seed in TRAINING_SEEDS:
                reused = variant == "a0_full" and not (
                    pde in ("joint_concave", "joint_convex") and seed in (1, 2)
                )
                estimate = 0.10 if reused else ESTIMATE[pde]
                cells.append({
                    "pde": pde,
                    "variant": variant,
                    "seed": seed,
                    "estimated_gpu_hours": estimate,
                })
    for pde in PDES:
        cells.append({
            "pde": pde,
            "variant": "a3_no_updates",
            "seed": 0,
            "estimated_gpu_hours": 0.10,
        })
    loads = [0.0] * 4
    assigned: list[list[dict[str, object]]] = [[] for _ in range(4)]
    for cell in sorted(cells, key=lambda row: float(row["estimated_gpu_hours"]), reverse=True):
        slot = min(range(4), key=lambda index: loads[index])
        assigned[slot].append(cell)
        loads[slot] += float(cell["estimated_gpu_hours"])
    rows = []
    for slot, local in enumerate(assigned):
        for sequence, cell in enumerate(local):
            rows.append({"slot": slot, "sequence": sequence, **cell})
    output = PACKAGE_DIR / "formal_queue.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print({"cells": len(rows), "slot_estimated_gpu_hours": loads, "path": str(output)})


if __name__ == "__main__":
    main()
