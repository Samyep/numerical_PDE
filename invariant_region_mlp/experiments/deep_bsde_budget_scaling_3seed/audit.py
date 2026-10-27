"""Final compact-artifact, completeness, reuse, and resource audit."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from .protocol import (
    FINAL_COUNT,
    INPUT_DIR,
    MAX_H100_GPU_HOURS,
    PDES,
    RESULTS_DIR,
    TABLE_DIR,
    TRAINING_SEEDS,
    file_sha256,
    scaling_task_specs,
    write_json_atomic,
)


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def audit() -> dict:
    failures = []
    manifest = json.loads((RESULTS_DIR / "query_manifest.json").read_text(encoding="utf-8"))
    for pde in PDES:
        entry = manifest["pdes"][pde]
        cal = set(entry["calibration_indices"])
        final = set(entry["final_indices"])
        if len(cal) != 32 or len(final) != 64 or cal & final:
            failures.append(f"{pde}: invalid query subsets")
        with np.load(INPUT_DIR / f"{pde}_final_inputs.npz", allow_pickle=False) as payload:
            if payload["x"].shape != (64, 100) or payload["t"].shape != (64,):
                failures.append(f"{pde}: invalid final input shapes")
    math_audit = json.loads((RESULTS_DIR / "mathematical_adaptation_audit.json").read_text(encoding="utf-8"))
    if math_audit["status"] != "PASS":
        failures.append("mathematical adaptation audit failed")
    preflight = json.loads((RESULTS_DIR / "preflight_allocation_gate.json").read_text(encoding="utf-8"))
    if not preflight.get("formal_submission_authorized", False):
        failures.append("allocation gate did not pass")
    completeness = json.loads((RESULTS_DIR / "training_completeness.json").read_text(encoding="utf-8"))
    if completeness["status"] != "ALL_FORMAL_TRAINING_COMPLETE_TRUTH_GATE_OPEN":
        failures.append("formal training incomplete")
    active_tasks = scaling_task_specs()
    expected_records = len(active_tasks) * FINAL_COUNT
    if completeness.get("terminal_records") != expected_records:
        failures.append("terminal record count mismatch")
    seed_rows = _csv_rows(RESULTS_DIR / "deep_bsde_seed_metrics.csv")
    per_query_rows = _csv_rows(RESULTS_DIR / "deep_bsde_per_query.csv")
    aggregate_rows = _csv_rows(RESULTS_DIR / "deep_bsde_scaling_summary.csv")
    if len(seed_rows) != len(PDES) * 9 * len(TRAINING_SEEDS):
        failures.append("seed metric row count mismatch")
    if len(aggregate_rows) != len(PDES) * 9:
        failures.append("sweep aggregate row count mismatch")
    if len(per_query_rows) != len(PDES) * 9 * len(TRAINING_SEEDS) * FINAL_COUNT:
        failures.append("display per-query row count mismatch")
    matched_status = json.loads(
        (RESULTS_DIR / "matched_existing_baselines_status.json").read_text(encoding="utf-8")
    )
    if matched_status.get("status") != "COMPLETE" or matched_status.get("row_count") != 108:
        failures.append("frozen DPI/PINN matched-64 evaluation incomplete")
    matched_rows = _csv_rows(RESULTS_DIR / "matched_existing_baselines_64_seed_metrics.csv")
    if len(matched_rows) != 108:
        failures.append("matched-64 seed row count mismatch")
    reuse = json.loads((RESULTS_DIR / "checkpoint_reuse_audit.json").read_text(encoding="utf-8"))
    if reuse["status"] != "PASS" or not all(row["exact_single_artifact_reuse"] for row in reuse["rows"]):
        failures.append("checkpoint reuse audit failed")
    accounting = json.loads((RESULTS_DIR / "resource_accounting.json").read_text(encoding="utf-8"))
    if accounting["total_allocated_h100_gpu_hours"] > MAX_H100_GPU_HOURS:
        failures.append("H100 allocation cap exceeded")
    submitted = json.loads((RESULTS_DIR / "submitted_jobs.json").read_text(encoding="utf-8"))
    accounted_ids = [str(row["job_id_raw"]) for row in accounting.get("rows", [])]
    for job_id in map(str, submitted.get("all_job_ids", [])):
        if not any(value == job_id or value.startswith(job_id + "_") for value in accounted_ids):
            failures.append(f"submitted Slurm job missing from accounting: {job_id}")
    scheduling = json.loads(
        (RESULTS_DIR / "scheduling_amendment.json").read_text(encoding="utf-8")
    )
    if scheduling.get("numerical_protocol_changed") is not False or scheduling.get("truth_opened") is not False:
        failures.append("scheduling amendment changed numerical protocol or opened truth")
    scope = json.loads((RESULTS_DIR / "scope_amendment.json").read_text(encoding="utf-8"))
    active_indices = [int(task["task_index"]) for task in active_tasks]
    if (
        scope.get("status") != "APPLIED_BEFORE_TRUTH_GATE"
        or scope.get("active_scope") != "C0_C1_C2_ONLY"
        or scope.get("truth_gate_open_at_application") is not False
        or scope.get("numerical_protocol_changed") is not False
        or scope.get("active_array_task_indices") != active_indices
        or scope.get("expected_terminal_records") != expected_records
    ):
        failures.append("C0/C1/C2 scope amendment mismatch")
    incidents = json.loads(
        (RESULTS_DIR / "infrastructure_incidents.json").read_text(encoding="utf-8")
    )
    for incident in incidents.get("incidents", []):
        if (
            incident.get("accuracy_or_truth_observed") is not False
            or incident.get("allocation_included_in_final_accounting") is not True
            or incident.get("retry_array_job_id") not in submitted.get("all_job_ids", [])
        ):
            failures.append("infrastructure incident provenance mismatch")
    required_tables = [
        TABLE_DIR / "deep_bsde_budget_scaling_main.tex",
        TABLE_DIR / "deep_bsde_budget_scaling_main_preview.pdf",
        TABLE_DIR / "matched_existing_baselines_64.tex",
        TABLE_DIR / "matched_existing_baselines_64_preview.pdf",
    ]
    for path in required_tables:
        if not path.is_file():
            failures.append(f"missing table artifact: {path.name}")
    compact = [
        RESULTS_DIR / "query_manifest.json",
        RESULTS_DIR / "task_manifest.json",
        RESULTS_DIR / "mathematical_adaptation_audit.json",
        RESULTS_DIR / "preflight_allocation_gate.json",
        RESULTS_DIR / "training_completeness.json",
        RESULTS_DIR / "deep_bsde_seed_metrics.csv",
        RESULTS_DIR / "deep_bsde_per_query.csv",
        RESULTS_DIR / "deep_bsde_scaling_summary.csv",
        RESULTS_DIR / "matched_existing_baselines_status.json",
        RESULTS_DIR / "matched_existing_baselines_64_seed_metrics.csv",
        RESULTS_DIR / "submitted_jobs.json",
        RESULTS_DIR / "scheduling_amendment.json",
        RESULTS_DIR / "scope_amendment.json",
        RESULTS_DIR / "checkpoint_storage_archive.json",
        RESULTS_DIR / "infrastructure_incidents.json",
        RESULTS_DIR / "resource_accounting.json",
        *required_tables,
    ]
    artifacts = [
        {"path": str(path), "bytes": path.stat().st_size, "sha256": file_sha256(path)}
        for path in compact if path.is_file()
    ]
    payload = {
        "status": "PASS" if not failures else "FAIL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "failures": failures,
        "expected_formal_tasks": len(active_tasks),
        "expected_terminal_records": expected_records,
        "observed_terminal_records": completeness.get("terminal_records"),
        "numerical_failure_count": completeness.get("numerical_failure_count"),
        "total_allocated_h100_gpu_hours": accounting["total_allocated_h100_gpu_hours"],
        "scope": "C0_C1_C2_ONLY",
        "artifacts": artifacts,
    }
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    return payload


if __name__ == "__main__":
    result = audit()
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)

