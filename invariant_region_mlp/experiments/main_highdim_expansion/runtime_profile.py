"""Freeze accuracy-blind runtime estimates from smoke and prior GH200 runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median
from typing import Any

from .protocol import DIMENSIONS, METHODS, PDES, RESULTS_DIR, SOURCE_ROOT, write_json_atomic


OLD_RESULTS = SOURCE_ROOT.parent / "source" / "invariant_region_mlp" / "results" / "tangent_actor_icml_final_suite" / "WHOLE_FUNCTION_RESULTS.csv"


def _old_times() -> dict[tuple[str, int], list[float]]:
    result: dict[tuple[str, int], list[float]] = {}
    if not OLD_RESULTS.is_file():
        return result
    with OLD_RESULTS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "COMPLETE":
                continue
            method = row.get("method")
            if method not in METHODS:
                continue
            try:
                dimension = int(row["dimension"])
                seconds = float(row.get("total_seconds") or row.get("gpu_seconds") or 0.0)
            except (TypeError, ValueError):
                continue
            if seconds > 0:
                result.setdefault((method, dimension), []).append(seconds / 3600.0)
    return result


def profile() -> dict[str, Any]:
    old = _old_times()
    smoke_evidence: list[dict[str, Any]] = []
    for path in sorted((RESULTS_DIR / "raw").glob("*smoke*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        smoke_evidence.append({
            "artifact": str(path),
            "run_id": payload.get("run_id"),
            "status": payload.get("status"),
            "training_seconds": payload.get("training_seconds")
            or payload.get("training_seconds_this_invocation"),
            "role": "integration/runtime floor only; never extrapolated into accuracy or formal budget",
        })
    defaults = {
        "tangent_actor": (0.7, 1.4),
        "sdgd_pinn": (1.3, 3.0),
        "dpi": (2.0, 5.0),
        "pinn_hte": (1.5, 4.0),
    }
    cells = []
    for pde in PDES:
        for d in DIMENSIONS:
            for method in METHODS:
                samples = old.get((method, d), [])
                if samples:
                    point = min(6.0, median(samples))
                    high = min(6.0, max(max(samples), point * 1.75))
                    source = "median and conservative high bound from completed original-suite GH200 cells"
                else:
                    point, high = defaults[method]
                    dimension_factor = max(0.6, (d / 100.0) ** 0.35)
                    point, high = min(6.0, point * dimension_factor), min(6.0, high * dimension_factor)
                    source = "predeclared method default with sublinear dimension scaling"
                cells.append(
                    {
                        "method": method,
                        "pde": pde,
                        "dimension": d,
                        "estimated_gpu_hours": point,
                        "high_gpu_hours": max(point, high),
                        "source": source,
                        "test_accuracy_used": False,
                    }
                )
    payload = {
        "status": "FROZEN_BEFORE_FORMAL_RUNS",
        "cells": cells,
        "estimated_total_gpu_hours": sum(row["estimated_gpu_hours"] for row in cells),
        "sum_of_cellwise_high_bounds_gpu_hours": sum(row["high_gpu_hours"] for row in cells),
        "global_cap_gpu_hours": 120.0,
        "point_estimate_within_global_cap": sum(row["estimated_gpu_hours"] for row in cells) <= 120.0,
        "smoke_evidence": smoke_evidence,
        "old_results_path": str(OLD_RESULTS),
        "rule": "method/dimension runtime evidence only; PDE accuracy never enters estimates",
    }
    write_json_atomic(RESULTS_DIR / "runtime_estimates.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    result = profile()
    print(json.dumps({"status": result["status"], "cells": len(result["cells"])}, sort_keys=True))


if __name__ == "__main__":
    main()
