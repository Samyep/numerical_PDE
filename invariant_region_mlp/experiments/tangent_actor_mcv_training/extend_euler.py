"""Extend the Euler diagnostic to the validation-selected Plain-trained Actors."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import torch

from .common import RAW_DIR, RESULTS_DIR, atomic_json, ensure_dirs, write_csv
from .evaluate import _evaluate_with_mcv_inference


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def run() -> None:
    ensure_dirs()
    output_path = RESULTS_DIR / "euler_sensitivity.csv"
    rows: list[dict[str, Any]] = _read_csv(output_path)
    summary = json.loads((RAW_DIR / "evaluation_summary.json").read_text())

    # The primary four-way comparison already contains Plain-trained Actors
    # with 32-step MCV inference. Include those rows in the symmetric Euler
    # table, then compute the two missing time grids.
    selected_plain = [
        row for row in summary["paired_final_rows"]
        if row["training_variant"] == "plain"
        and row["inference_variant"] == "plain"
        and row.get("b5_validation_selected")
    ]
    existing = {
        (
            row["training_variant"], row["benchmark"], int(row["seed"]),
            row["checkpoint"], int(row["inference_time_steps"]),
        )
        for row in rows
    }
    for source in selected_plain:
        matching_32 = next(
            row for row in summary["paired_final_rows"]
            if row["training_variant"] == "plain"
            and row["inference_variant"] == "mcv"
            and row["benchmark"] == source["benchmark"]
            and int(row["seed"]) == int(source["seed"])
            and row["checkpoint"] == source["checkpoint"]
        )
        key_32 = (
            "plain", source["benchmark"], int(source["seed"]),
            source["checkpoint"], 32,
        )
        if key_32 not in existing:
            rows.append(matching_32)
            existing.add(key_32)

    device = torch.device("cuda")
    for source in selected_plain:
        for steps in (16, 64):
            key = (
                "plain", source["benchmark"], int(source["seed"]),
                source["checkpoint"], steps,
            )
            if key in existing:
                continue
            row = _evaluate_with_mcv_inference(
                source["checkpoint_path"],
                source["benchmark"],
                int(source["seed"]),
                checkpoint=source["checkpoint"],
                training_variant="plain",
                training_seconds=float(source["training_seconds"]),
                steps=steps,
                device=device,
            )
            rows.append(row)
            existing.add(key)
            rows.sort(key=lambda item: (
                item["benchmark"], int(item["seed"]),
                item["training_variant"], int(item["inference_time_steps"]),
            ))
            write_csv(output_path, rows)
            atomic_json(RAW_DIR / "plain_euler_extension.json", {"rows": rows})
            print(
                f"Euler extension {source['benchmark']} s{source['seed']} "
                f"Plain training, steps={steps}: rel={row['relative_l2']:.6g} "
                f"sys={row['systematic_rmse_proxy']:.6g}",
                flush=True,
            )

    summary["euler_sensitivity"] = rows
    summary["euler_training_variants"] = ["plain", "mcv"]
    atomic_json(RAW_DIR / "evaluation_summary.json", summary)


if __name__ == "__main__":
    run()
