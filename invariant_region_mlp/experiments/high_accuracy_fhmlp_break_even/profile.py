"""Truth-free FH-MLP performance aggregation and frozen cost gate."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from .protocol import (
    ALLOCATED_GPU_HOUR_CAP,
    BANK_SIZES,
    CALIBRATION_TASK_CAP_SECONDS,
    COST_GATE,
    FH_ACCURACY_REPLICATES,
    FH_CONFIGS,
    PROFILE_REPLICATES,
    PROFILE_RESULTS,
    RAW_DIR,
    REVISION,
    TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE,
    TIMING_GPU_HOUR_RESERVE,
    require_preflight,
    write_json_atomic,
)


SAFETY_FACTOR = 1.35


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def aggregate() -> dict[str, Any]:
    require_preflight()
    rows: list[dict[str, Any]] = []
    profile_active_hours = 0.0
    for n, M in FH_CONFIGS:
        records = []
        missing = []
        for replicate in range(PROFILE_REPLICATES):
            path = RAW_DIR / "profile" / "fh_mlp" / f"fh__profile__n{n}_M{M}__r{replicate}.json"
            if path.is_file():
                records.append(_read(path))
            else:
                missing.append(str(path))
        key = f"n{n}_M{M}"
        if len(records) != PROFILE_REPLICATES or any(
            row.get("status") != "complete" or int(row.get("nonfinite_count", 1))
            for row in records
        ):
            rows.append({
                "candidate": key,
                "n": n,
                "M": M,
                "status": "NOT_EVALUATED_PROFILE_FAILURE",
                "missing_records": missing,
                "observed_records": len(records),
            })
            continue
        active = np.asarray([float(row["active_seconds"]) for row in records])
        seconds_per_query = np.asarray(
            [float(row["seconds_per_query"]) for row in records]
        )
        profile_active_hours += float(active.sum()) / 3600.0
        mean_seconds = float(seconds_per_query.mean())
        projected_bundle_seconds = (
            mean_seconds
            * BANK_SIZES["calibration"]
            * FH_ACCURACY_REPLICATES
            * SAFETY_FACTOR
        )
        generator_calls = []
        for row in records:
            stats = row["stats"]
            generator_calls.append(
                (
                    float(stats["terminal_samples"])
                    + float(stats["transition_samples"])
                )
                / float(row["queries"])
            )
        rows.append({
            "candidate": key,
            "n": n,
            "M": M,
            "status": "PROFILE_COMPLETE",
            "observed_records": len(records),
            "mean_seconds_per_query": mean_seconds,
            "std_seconds_per_query": float(seconds_per_query.std(ddof=1)),
            "peak_gpu_bytes": max(int(row["peak_gpu_bytes"]) for row in records),
            "mean_generator_calls_per_query": float(np.mean(generator_calls)),
            "projected_calibration_bundle_seconds": projected_bundle_seconds,
            "projected_calibration_gpu_hours": projected_bundle_seconds / 3600.0,
        })

    calibration_budget = (
        ALLOCATED_GPU_HOUR_CAP
        - TIMING_GPU_HOUR_RESERVE
        - TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE
        - SAFETY_FACTOR * profile_active_hours
    )
    if calibration_budget <= 0:
        raise RuntimeError("profile consumed the entire predeclared allocation")
    complete = sorted(
        (row for row in rows if row["status"] == "PROFILE_COMPLETE"),
        key=lambda row: (
            row["projected_calibration_gpu_hours"], row["n"], row["M"]
        ),
    )
    eligible: dict[str, dict[str, Any]] = {}
    projected_calibration_hours = 0.0
    for row in complete:
        projected = float(row["projected_calibration_gpu_hours"])
        single_task_ok = (
            float(row["projected_calibration_bundle_seconds"])
            <= CALIBRATION_TASK_CAP_SECONDS
        )
        cumulative_ok = projected_calibration_hours + projected <= calibration_budget
        if single_task_ok and cumulative_ok:
            row["cost_gate_status"] = "ELIGIBLE"
            projected_calibration_hours += projected
            eligible[row["candidate"]] = {
                key: value for key, value in row.items()
                if key not in {"cost_gate_status"}
            }
        else:
            row["cost_gate_status"] = "NOT_EVALUATED_COMPUTE_CAP"
            row["cost_gate_reason"] = (
                "projected_single_task_cap"
                if not single_task_ok
                else "projected_total_allocation_cap"
            )
    for row in rows:
        if row["status"] != "PROFILE_COMPLETE":
            row["cost_gate_status"] = row["status"]

    projected_total = (
        SAFETY_FACTOR * profile_active_hours
        + projected_calibration_hours
        + TIMING_GPU_HOUR_RESERVE
        + TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE
    )
    payload = {
        "revision": REVISION,
        "status": "PASS" if eligible and projected_total <= ALLOCATED_GPU_HOUR_CAP else "FAIL",
        "profile_rows": rows,
        "profile_active_gpu_hours": profile_active_hours,
        "safety_factor": SAFETY_FACTOR,
        "calibration_gpu_hour_budget": calibration_budget,
        "projected_calibration_gpu_hours": projected_calibration_hours,
        "timing_gpu_hour_reserve": TIMING_GPU_HOUR_RESERVE,
        "ta_and_final_accuracy_gpu_hour_reserve": TA_AND_FINAL_ACCURACY_GPU_HOUR_RESERVE,
        "projected_complete_study_gpu_hours": projected_total,
        "eligible_fh_configs": eligible,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "truth_used": False,
    }
    write_json_atomic(PROFILE_RESULTS, {
        "revision": REVISION,
        "rows": rows,
        "profile_active_gpu_hours": profile_active_hours,
        "truth_used": False,
    })
    write_json_atomic(COST_GATE, payload)
    if payload["status"] != "PASS":
        raise RuntimeError("projected cost gate failed")
    return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))
