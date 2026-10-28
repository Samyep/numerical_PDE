"""Final artifact and scientific-protocol audit."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from .analyze import _records
from .protocol import (
    ALLOCATED_GPU_HOUR_CAP,
    BANK_SIZES,
    FH_CONFIGS,
    FIGURE_DIR,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    TARGETS,
    array_sha256,
    file_sha256,
    load_bank,
    require_frozen_selections,
    require_frozen_timing_grid,
    write_json_atomic,
)


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def audit() -> dict[str, Any]:
    frozen = require_frozen_selections()
    timing_grid = require_frozen_timing_grid()
    final = _read(RESULTS_DIR / "final_accuracy_summary.json")
    analysis = _read(RESULTS_DIR / "analysis_summary.json")
    coverage = _read(RESULTS_DIR / "timing_coverage.json")
    resources = _read(RESULTS_DIR / "resource_accounting.json")
    preflight = _read(RESULTS_DIR / "preflight_audit.json")
    reference = _read(RESULTS_DIR / "reference_audit.json")
    profile = _read(RESULTS_DIR / "performance_profile.json")
    timing = load_bank("timing", require_truth=False)

    required = [
        RESULTS_DIR / name
        for name in (
            "FROZEN_PROTOCOL.json", "bank_manifest.json", "reference_audit.json",
            "preflight_audit.json", "performance_profile.json", "cost_gate.json",
            "calibration_rows.csv", "FROZEN_SELECTIONS.json",
            "final_accuracy_rows.csv", "final_accuracy_summary.json",
            "FROZEN_TIMING_GRID.json", "timing_coverage.json",
            "timing_replicates.csv", "timing_rows.csv", "break_even_summary.csv",
            "break_even_summary.json", "analysis_summary.json",
            "resource_accounting.csv", "resource_accounting.json",
            "ICML_FIGURE_CAPTION.tex", "RESULTS_PARAGRAPH.md",
        )
    ]
    primary_pair = (
        coverage["targets"]["primary"]["ta_direct"] == "RUN"
        and coverage["targets"]["primary"]["fh_mlp_standard"] == "RUN"
    )
    if primary_pair:
        required += [
            FIGURE_DIR / f"high_accuracy_break_even_primary.{suffix}"
            for suffix in ("pdf", "png", "svg")
        ]
    secondary_pair = (
        coverage["targets"]["secondary"]["ta_direct"] == "RUN"
        and coverage["targets"]["secondary"]["fh_mlp_standard"] == "RUN"
    )
    if secondary_pair:
        required += [
            FIGURE_DIR / f"high_accuracy_break_even_secondary.{suffix}"
            for suffix in ("pdf", "png", "svg")
        ]
    if coverage["targets"]["primary"]["ta_fk_mcv"] == "RUN" and primary_pair:
        required += [
            FIGURE_DIR / f"ta_inference_modes_timing.{suffix}"
            for suffix in ("pdf", "png", "svg")
        ]

    needed: set[tuple[str, str]] = set()
    for target in TARGETS:
        for method in ("ta_direct", "ta_fk_mcv", "fh_mlp_standard"):
            item = final["targets"][target][method]
            if coverage["targets"][target][method] == "RUN":
                needed.add((method, item["selected_candidate"]))
    timing_detail = []
    frozen_counts = set(int(value) for value in timing_grid["timing_counts"])
    for method, candidate in sorted(needed):
        records = _records(method, candidate)
        for replicate, record in enumerate(records):
            completed = int(record["completed_queries"])
            prediction_path = Path(record["prediction_file"])
            with np.load(prediction_path, allow_pickle=False) as data:
                prediction = np.asarray(data["prediction"], dtype=np.float64)
            counts = [int(row["queries"]) for row in record["prefix_rows"]]
            seconds = [float(row["total_seconds"]) for row in record["prefix_rows"]]
            checks = {
                "actual_prefix_execution": record.get("actual_prefix_execution") is True,
                "truth_not_loaded": record.get("truth_loaded") is False,
                "h100": "H100" in str(record.get("gpu", "")).upper(),
                "valid_status": record.get("status") in {"complete", "TIME_CAP"},
                "no_nonfinite": int(record.get("nonfinite_count", -1)) == 0,
                "input_prefix_hash": record.get("query_input_sha256")
                == array_sha256(timing["t"][:completed], timing["x"][:completed]),
                "prefixes_frozen_and_monotone": bool(counts)
                and set(counts).issubset(frozen_counts)
                and all(b > a for a, b in zip(counts, counts[1:])),
                "times_monotone": bool(seconds)
                and all(b >= a for a, b in zip(seconds, seconds[1:])),
                "prediction_shape": prediction.shape == (completed,),
                "prediction_finite": bool(np.isfinite(prediction).all()),
                "prediction_hash": record.get("prediction_sha256") == array_sha256(prediction),
            }
            timing_detail.append({
                "method": method,
                "candidate": candidate,
                "replicate": replicate,
                "checks": checks,
                "overall_pass": all(checks.values()),
            })

    calibration_candidates = {
        row["candidate"]
        for row in frozen["candidate_summaries"]
        if row["method"] == "fh_mlp_standard"
    }
    expected_candidates = {f"n{n}_M{M}" for n, M in FH_CONFIGS}
    result_names = [path.name.lower() for path in RESULTS_DIR.rglob("*") if path.is_file()]
    checks = {
        "all_required_artifacts_exist": all(path.is_file() for path in required),
        "timing_bank_has_no_truth": "reference" not in timing,
        "timing_bank_shape": timing["x"].shape == (BANK_SIZES["timing"], 100),
        "reference_validated": reference.get("status") == "VALIDATED",
        "preflight_validated": preflight.get("status") == "VALIDATED",
        "ta_source_artifacts_exact": preflight["checks"]["ta_source_artifacts_exact"],
        "original_fh_recursion_retained": preflight["fh_batched_execution_audit"]["original_recursion_retained"],
        "complete_candidate_status_table": calibration_candidates == expected_candidates,
        "profile_covers_candidate_grid": {
            row["candidate"] for row in profile["rows"]
        } == expected_candidates,
        "final_configurations_not_changed": final.get("configuration_changes_after_final_errors") is False,
        "three_timing_replicates": all(len(_records(*key)) == 3 for key in needed),
        "timing_records_integrity": bool(timing_detail)
        and all(row["overall_pass"] for row in timing_detail),
        "resource_accounting_pass": resources.get("status") == "PASS",
        "within_30_h100_hour_cap": float(resources["total_allocated_h100_gpu_hours"])
        <= ALLOCATED_GPU_HOUR_CAP,
        "analysis_complete": analysis.get("status") == "COMPLETE",
        "no_out_of_scope_baseline_artifacts": not any(
            token in name for name in result_names
            for token in ("deep_bsde", "pinn_hte", "sdgd", "dpi_")
        ),
    }
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(checks.values()) else "FAILED",
        "overall_pass": all(checks.values()),
        "checks": checks,
        "allocated_h100_gpu_hours": resources["total_allocated_h100_gpu_hours"],
        "timing_record_audit": timing_detail,
        "artifacts": {
            str(path.relative_to(SOURCE_ROOT)): file_sha256(path)
            for path in required if path.is_file()
        },
        "audited_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    if not payload["overall_pass"]:
        raise RuntimeError(f"artifact audit failed: {checks}")
    return payload


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
