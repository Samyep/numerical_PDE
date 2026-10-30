#!/usr/bin/env python3
"""Consolidate per-PDE JSON artifacts without changing source results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ORDER = ("fourier_tanh", "joint_concave", "nonconvex_coupled", "ou_gmm")


def load_many(pattern: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pde in ORDER:
        path = ROOT / pattern.format(pde=pde)
        if path.exists():
            value = json.loads(path.read_text())
            rows.extend(value if isinstance(value, list) else [value])
    return rows


def flatten(row: dict[str, Any]) -> dict[str, Any]:
    output = {}
    for key, value in row.items():
        if isinstance(value, (dict, list)):
            output[key] = json.dumps(value, sort_keys=True, separators=(",", ":"))
        else:
            output[key] = value
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    flat = [flatten(row) for row in rows]
    keys = sorted({key for row in flat for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(flat)


def main() -> None:
    audit = load_many("audit/checkpoint_audit_{pde}.json")
    canonical = load_many("audit/canonical_metric_audit_{pde}.json")
    slices = load_many("metrics/slice_metrics_{pde}.json")
    (ROOT / "checkpoint_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    write_csv(ROOT / "checkpoint_audit.csv", audit)
    write_csv(ROOT / "canonical_metric_audit.csv", canonical)
    write_csv(ROOT / "slice_metrics.csv", slices)
    print(
        json.dumps(
            {
                "checkpoint_rows": len(audit),
                "canonical_rows": len(canonical),
                "slice_rows": len(slices),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
