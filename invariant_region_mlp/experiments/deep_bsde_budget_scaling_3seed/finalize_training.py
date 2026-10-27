"""Open the truth gate only after every formal query has a terminal record."""

from __future__ import annotations

from datetime import datetime, timezone
import json

from .protocol import FINAL_COUNT, RAW_DIR, RESULTS_DIR, STATUS_DIR, file_sha256, scaling_task_specs, write_json_atomic


def finalize() -> dict:
    incomplete = []
    numerical_failures = []
    completed = []
    active_tasks = scaling_task_specs()
    for task in active_tasks:
        status_path = STATUS_DIR / f"task_{task['task_index']:03d}.json"
        if not status_path.is_file():
            incomplete.append({"task_index": task["task_index"], "reason": "missing_task_status"})
        for position in range(FINAL_COUNT):
            path = (
                RAW_DIR / "formal" / task["pde"] / task["trajectory_id"]
                / f"s{task['seed']}" / f"q{position:02d}.json"
            )
            if not path.is_file():
                incomplete.append({
                    "task_index": task["task_index"], "position": position,
                    "reason": "missing_terminal_record",
                })
                continue
            record = json.loads(path.read_text(encoding="utf-8"))
            status = record.get("status")
            if status not in {"complete", "numerical_failure"}:
                incomplete.append({
                    "task_index": task["task_index"], "position": position,
                    "reason": f"nonterminal_status:{status}",
                })
            elif status == "numerical_failure":
                numerical_failures.append({
                    "task_index": task["task_index"], "position": position,
                    "failure_reason": record.get("failure_reason"),
                })
            completed.append({"path": str(path), "sha256": file_sha256(path)})
    if incomplete:
        payload = {
            "status": "INCOMPLETE_FORMAL_TRAINING",
            "checked_utc": datetime.now(timezone.utc).isoformat(),
            "incomplete": incomplete,
            "terminal_records": len(completed),
        }
        write_json_atomic(RESULTS_DIR / "training_completeness.json", payload)
        return payload
    payload = {
        "status": "ALL_FORMAL_TRAINING_COMPLETE_TRUTH_GATE_OPEN",
        "opened_utc": datetime.now(timezone.utc).isoformat(),
        "task_count": len(active_tasks),
        "expected_terminal_records": len(active_tasks) * FINAL_COUNT,
        "terminal_records": len(completed),
        "numerical_failure_count": len(numerical_failures),
        "numerical_failures": numerical_failures,
        "query_manifest_sha256": file_sha256(RESULTS_DIR / "query_manifest.json"),
        "task_manifest_sha256": file_sha256(RESULTS_DIR / "task_manifest.json"),
        "records": completed,
        "truth_used_before_gate": False,
        "scope": "C0_C1_C2_ONLY",
        "joint_scaling_excluded_by_user_amendment": True,
    }
    write_json_atomic(RESULTS_DIR / "training_completeness.json", payload)
    write_json_atomic(STATUS_DIR / "truth_gate_open.json", payload)
    return payload


if __name__ == "__main__":
    result = finalize()
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "ALL_FORMAL_TRAINING_COMPLETE_TRUTH_GATE_OPEN":
        raise SystemExit(1)

