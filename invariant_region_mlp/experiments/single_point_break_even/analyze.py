"""Aggregate measured prefixes, determine crossings, plot, and audit artifacts."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .protocol import (
    FIGURE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    TARGETS,
    array_sha256,
    file_sha256,
    load_bank,
    require_frozen_selections,
    write_json_atomic,
)


METHOD_LABELS = {
    "tangent_actor": "Tangent-Actor",
    "fh_mlp_standard": "Full-History MLP",
    "deep_bsde": "Deep BSDE",
}
METHOD_COLORS = {
    "tangent_actor": "#1665D8",
    "fh_mlp_standard": "#E76F00",
    "deep_bsde": "#2A9D62",
}


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _timing_records(method: str, candidate: str) -> list[dict[str, Any]]:
    root = RAW_DIR / "timing"
    if method == "tangent_actor":
        paths = [root / "tangent_actor" / f"ta__timing__s{seed}__{candidate}.json"
                 for seed in (10, 11, 12)]
    elif method == "fh_mlp_standard":
        n_text, M_text = candidate.split("_")
        n, M = int(n_text[1:]), int(M_text[1:])
        paths = [root / "fh_mlp" / f"fh__timing__n{n}_M{M}__r{rep}.json"
                 for rep in range(3)]
    else:
        updates = int(candidate.split("_")[1])
        paths = [root / "deep_bsde" / f"db__timing__u{updates}__r{rep}.json"
                 for rep in range(3)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing timing records: {missing}")
    return [_read(path) for path in paths]


def _flatten_timing(
    method: str, candidate: str, records: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows = []
    for replicate, record in enumerate(records):
        for prefix in record["prefix_rows"]:
            rows.append({
                "method": method,
                "candidate": candidate,
                "replicate": replicate,
                "queries": int(prefix["queries"]),
                "training_seconds": float(prefix.get("training_seconds", 0.0)),
                "checkpoint_load_seconds": float(prefix.get("checkpoint_load_seconds", 0.0)),
                "inference_seconds": float(prefix.get(
                    "inference_seconds", prefix.get("total_seconds", 0.0)
                )),
                "total_seconds": float(prefix["total_seconds"]),
                "run_status": record["status"],
                "gpu": record["gpu"],
            })
    return rows


def _mean_curve(rows: list[dict[str, Any]]) -> dict[int, dict[str, float]]:
    values: dict[int, list[float]] = {}
    for row in rows:
        values.setdefault(int(row["queries"]), []).append(float(row["total_seconds"]))
    result = {}
    for queries, local in values.items():
        if len(local) != 3:
            continue
        result[queries] = {
            "mean": float(np.mean(local)),
            "std": float(np.std(local, ddof=1)),
            "minimum": float(np.min(local)),
            "maximum": float(np.max(local)),
        }
    return result


def _crossing(
    ta_rows: list[dict[str, Any]], baseline_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    ta_curve, base_curve = _mean_curve(ta_rows), _mean_curve(baseline_rows)
    common = sorted(set(ta_curve) & set(base_curve))
    if not common:
        return {"status": "NOT_EVALUATED_NO_COMMON_MEASURED_PREFIX"}
    crossing = next(
        (queries for queries in common if base_curve[queries]["mean"] >= ta_curve[queries]["mean"]),
        None,
    )
    per_replicate = []
    for replicate in range(3):
        ta = {int(row["queries"]): float(row["total_seconds"]) for row in ta_rows
              if int(row["replicate"]) == replicate}
        base = {int(row["queries"]): float(row["total_seconds"]) for row in baseline_rows
                if int(row["replicate"]) == replicate}
        local_common = sorted(set(ta) & set(base))
        local_crossing = next(
            (queries for queries in local_common if base[queries] >= ta[queries]), None
        )
        per_replicate.append({
            "replicate": replicate,
            "crossing_queries": local_crossing,
            "previous_queries": (
                max((q for q in local_common if local_crossing is not None and q < local_crossing),
                    default=None)
            ),
            "largest_common_queries": local_common[-1] if local_common else None,
        })
    if crossing is None:
        largest = common[-1]
        return {
            "status": "NOT_REACHED_WITHIN_TESTED_BUDGET",
            "largest_common_queries": largest,
            "tangent_actor_mean_seconds": ta_curve[largest]["mean"],
            "baseline_mean_seconds": base_curve[largest]["mean"],
            "baseline_over_tangent_actor_time": (
                base_curve[largest]["mean"] / ta_curve[largest]["mean"]
            ),
            "per_replicate": per_replicate,
        }
    previous = max((q for q in common if q < crossing), default=None)
    return {
        "status": "OBSERVED_BREAK_EVEN",
        "previous_queries": previous,
        "crossing_queries": crossing,
        "interval": [previous, crossing],
        "tangent_actor_mean_seconds_at_crossing": ta_curve[crossing]["mean"],
        "baseline_mean_seconds_at_crossing": base_curve[crossing]["mean"],
        "per_replicate": per_replicate,
    }


def _plot(
    final: dict[str, Any],
    target_rows: dict[str, dict[str, list[dict[str, Any]]]],
    crossings: dict[str, Any],
    timing_coverage: dict[str, Any],
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.75), sharey=True)
    for axis, (target_name, threshold) in zip(axes, TARGETS.items()):
        details = final["targets"][target_name]
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            item = details[method]
            coverage_status = timing_coverage["targets"][target_name][method]
            if item["status"] != "FINAL_ACCURACY_PASSED" or coverage_status != "RUN":
                label = METHOD_LABELS[method]
                status = (
                    item["status"]
                    if item["status"] != "FINAL_ACCURACY_PASSED"
                    else coverage_status
                )
                reason = status.replace("_", " ").lower()
                axis.text(
                    0.98,
                    0.06 + 0.07 * list(METHOD_LABELS).index(method),
                    f"{label}: {reason}",
                    color=METHOD_COLORS[method],
                    transform=axis.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=7.5,
                )
                continue
            rows = target_rows[target_name][method]
            curve = _mean_curve(rows)
            q = np.asarray(sorted(curve), dtype=np.float64)
            mean = np.asarray([curve[int(value)]["mean"] for value in q])
            std = np.asarray([curve[int(value)]["std"] for value in q])
            axis.plot(
                q,
                mean,
                marker="o",
                markersize=3.6,
                linewidth=1.7,
                color=METHOD_COLORS[method],
                label=METHOD_LABELS[method],
            )
            axis.fill_between(
                q,
                np.maximum(mean - std, np.finfo(float).tiny),
                mean + std,
                color=METHOD_COLORS[method],
                alpha=0.14,
                linewidth=0,
            )
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.grid(True, which="major", linewidth=0.45, alpha=0.28)
        axis.set_xlabel("Number of distinct queries")
        axis.set_title(rf"RMSE $\leq {threshold:.1e}$", fontsize=10)
        axis.tick_params(labelsize=8)
        for baseline in ("fh_mlp_standard", "deep_bsde"):
            crossing = crossings[target_name].get(baseline, {})
            if crossing.get("status") == "OBSERVED_BREAK_EVEN":
                value = crossing["crossing_queries"]
                axis.axvline(
                    value,
                    color=METHOD_COLORS[baseline],
                    linestyle=":",
                    linewidth=1.0,
                    alpha=0.8,
                )
                axis.text(
                    value,
                    0.97,
                    f"{value} queries",
                    color=METHOD_COLORS[baseline],
                    transform=axis.get_xaxis_transform(),
                    rotation=90,
                    ha="right",
                    va="top",
                    fontsize=7,
                )
    axes[0].set_ylabel("Measured total H100 time (s)")
    handles, labels = axes[0].get_legend_handles_labels()
    if len(handles) < 3:
        other_h, other_l = axes[1].get_legend_handles_labels()
        for handle, label in zip(other_h, other_l):
            if label not in labels:
                handles.append(handle); labels.append(label)
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, fontsize=8.5)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.17, top=0.82, wspace=0.10)
    for suffix in ("png", "pdf"):
        fig.savefig(
            FIGURE_DIR / f"matched_accuracy_break_even.{suffix}",
            dpi=300 if suffix == "png" else None,
        )
    plt.close(fig)


def analyze() -> dict[str, Any]:
    frozen = require_frozen_selections()
    analysis_code_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()
    final = _read(RESULTS_DIR / "final_accuracy_summary.json")
    timing_coverage = _read(RESULTS_DIR / "timing_coverage.json")
    all_rows: list[dict[str, Any]] = []
    cache: dict[tuple[str, str], list[dict[str, Any]]] = {}
    target_rows: dict[str, dict[str, list[dict[str, Any]]]] = {}
    crossings: dict[str, Any] = {}
    break_even_rows: list[dict[str, Any]] = []
    for target_name in TARGETS:
        target_rows[target_name] = {}
        crossings[target_name] = {}
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            item = final["targets"][target_name][method]
            if item["status"] != "FINAL_ACCURACY_PASSED":
                continue
            if timing_coverage["targets"][target_name][method] != "RUN":
                continue
            candidate = item["selected_candidate"]
            key = (method, candidate)
            if key not in cache:
                cache[key] = _flatten_timing(method, candidate, _timing_records(method, candidate))
                all_rows.extend(cache[key])
            target_rows[target_name][method] = cache[key]
        if "tangent_actor" not in target_rows[target_name]:
            for baseline in ("fh_mlp_standard", "deep_bsde"):
                baseline_final_status = final["targets"][target_name][baseline]["status"]
                status = (
                    baseline_final_status
                    if baseline_final_status != "FINAL_ACCURACY_PASSED"
                    else timing_coverage["targets"][target_name]["tangent_actor"]
                )
                local = {"status": status}
                crossings[target_name][baseline] = local
                break_even_rows.append({
                    "target": target_name,
                    "threshold": TARGETS[target_name],
                    "baseline": baseline,
                    "status": local["status"],
                    "previous_queries": None,
                    "crossing_queries": None,
                    "largest_common_queries": None,
                    "tangent_actor_seconds": None,
                    "baseline_seconds": None,
                })
            continue
        for baseline in ("fh_mlp_standard", "deep_bsde"):
            if baseline not in target_rows[target_name]:
                status = final["targets"][target_name][baseline]["status"]
                crossings[target_name][baseline] = {"status": status}
            else:
                crossings[target_name][baseline] = _crossing(
                    target_rows[target_name]["tangent_actor"],
                    target_rows[target_name][baseline],
                )
            local = crossings[target_name][baseline]
            break_even_rows.append({
                "target": target_name,
                "threshold": TARGETS[target_name],
                "baseline": baseline,
                "status": local["status"],
                "previous_queries": local.get("previous_queries"),
                "crossing_queries": local.get("crossing_queries"),
                "largest_common_queries": local.get("largest_common_queries"),
                "tangent_actor_seconds": local.get(
                    "tangent_actor_mean_seconds_at_crossing",
                    local.get("tangent_actor_mean_seconds"),
                ),
                "baseline_seconds": local.get(
                    "baseline_mean_seconds_at_crossing",
                    local.get("baseline_mean_seconds"),
                ),
            })
    _write_csv(RESULTS_DIR / "timing_rows.csv", all_rows)
    _write_csv(RESULTS_DIR / "break_even_summary.csv", break_even_rows)
    _plot(final, target_rows, crossings, timing_coverage)
    payload = {
        "revision": REVISION,
        "status": "COMPLETE",
        "analyzed_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_selections_canonical_sha256": frozen["canonical_sha256"],
        "analysis_code_commit": analysis_code_commit,
        "targets": final["targets"],
        "break_even": crossings,
        "timing_row_count": len(all_rows),
        "timing_coverage": timing_coverage,
        "no_extrapolated_crossings": True,
        "hardware": "NVIDIA H100 80GB",
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload)
    return payload


def resource_accounting() -> list[dict[str, Any]]:
    rows = []
    for path in sorted((RAW_DIR / "task_status").glob("*.json")):
        task = _read(path)
        rows.append({
            "task_id": task["task_id"],
            "stage": task["stage"],
            "kind": task["task"]["kind"],
            "status": task["status"],
            "task_wall_seconds": task["task_wall_seconds"],
            "slurm_job_id": task.get("slurm_job_id"),
            "slurm_array_task_id": task.get("slurm_array_task_id"),
            "gpu_count": 1,
            "error": task.get("error"),
        })
    infrastructure_path = RESULTS_DIR / "infrastructure_attempts.json"
    if infrastructure_path.is_file():
        for attempt in _read(infrastructure_path)["attempts"]:
            rows.append({
                "task_id": attempt["task_id"],
                "stage": attempt["stage"],
                "kind": attempt["kind"],
                "status": attempt["status"],
                "task_wall_seconds": attempt["task_wall_seconds"],
                "slurm_job_id": attempt["slurm_job_id"],
                "slurm_array_task_id": None,
                "gpu_count": attempt["gpu_count"],
                "error": attempt.get("error"),
            })
    _write_csv(RESULTS_DIR / "resource_accounting.csv", rows)
    return rows


def _timing_integrity(
    timing: dict[str, np.ndarray],
) -> tuple[dict[str, bool], list[dict[str, Any]]]:
    """Independently audit every timing record selected for a reported curve."""
    final = _read(RESULTS_DIR / "final_accuracy_summary.json")
    coverage = _read(RESULTS_DIR / "timing_coverage.json")
    selected: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for target_name in TARGETS:
        for method in METHOD_LABELS:
            item = final["targets"][target_name][method]
            if (
                item["status"] == "FINAL_ACCURACY_PASSED"
                and coverage["targets"][target_name][method] == "RUN"
            ):
                key = (method, item["selected_candidate"])
                if key not in selected:
                    selected[key] = _timing_records(*key)

    details: list[dict[str, Any]] = []
    for (method, candidate), records in selected.items():
        for replicate, record in enumerate(records):
            completed = int(record["completed_queries"])
            expected_input_hash = array_sha256(
                timing["t"][:completed], timing["x"][:completed]
            )
            prefix_rows = record["prefix_rows"]
            query_counts = [int(row["queries"]) for row in prefix_rows]
            total_seconds = [float(row["total_seconds"]) for row in prefix_rows]
            prediction_path = Path(record["prediction_file"])
            prediction_exists = prediction_path.is_file()
            prediction_shape_ok = False
            prediction_finite = False
            prediction_hash_ok = False
            if prediction_exists:
                with np.load(prediction_path, allow_pickle=False) as payload:
                    prediction = np.asarray(payload["prediction"], dtype=np.float64)
                prediction_shape_ok = prediction.shape == (completed,)
                prediction_finite = bool(np.isfinite(prediction).all())
                prediction_hash_ok = (
                    array_sha256(prediction) == record["prediction_sha256"]
                )
            local_checks = {
                "actual_prefix_execution": record.get("actual_prefix_execution") is True,
                "truth_not_loaded": record.get("truth_loaded") is False,
                "h100_80gb": "H100 80GB" in str(record.get("gpu", "")),
                "valid_status": record.get("status") in {"complete", "TIME_CAP"},
                "no_nonfinite": int(record.get("nonfinite_count", -1)) == 0,
                "input_prefix_hash": record.get("query_input_sha256") == expected_input_hash,
                "prefixes_consistent_with_completed_queries": (
                    bool(query_counts)
                    and all(b > a for a, b in zip(query_counts, query_counts[1:]))
                    and query_counts[0] > 0
                    and query_counts[-1] <= completed
                ),
                "time_is_monotone": (
                    bool(total_seconds)
                    and all(b >= a for a, b in zip(total_seconds, total_seconds[1:]))
                ),
                "prediction_exists": prediction_exists,
                "prediction_shape": prediction_shape_ok,
                "prediction_finite": prediction_finite,
                "prediction_array_hash": prediction_hash_ok,
            }
            details.append({
                "method": method,
                "candidate": candidate,
                "replicate": replicate,
                "completed_queries": completed,
                "record_status": record["status"],
                "checks": local_checks,
                "overall_pass": all(local_checks.values()),
            })
    checks = {
        "reported_timing_configurations_present": bool(selected),
        "three_replicates_per_timing_configuration": all(
            len(records) == 3 for records in selected.values()
        ),
        "all_timing_records_integrity_checked": (
            bool(details) and all(item["overall_pass"] for item in details)
        ),
    }
    return checks, details


def artifact_audit() -> dict[str, Any]:
    required = [
        RESULTS_DIR / "bank_manifest.json",
        RESULTS_DIR / "reference_audit.json",
        RESULTS_DIR / "calibration_rows.csv",
        RESULTS_DIR / "FROZEN_SELECTIONS.json",
        RESULTS_DIR / "final_accuracy_rows.csv",
        RESULTS_DIR / "timing_rows.csv",
        RESULTS_DIR / "break_even_summary.csv",
        RESULTS_DIR / "resource_accounting.csv",
        RESULTS_DIR / "analysis_summary.json",
        RESULTS_DIR / "infrastructure_attempts.json",
        RESULTS_DIR / "timing_coverage.json",
        FIGURE_DIR / "matched_accuracy_break_even.png",
        FIGURE_DIR / "matched_accuracy_break_even.pdf",
    ]
    timing = load_bank("timing", require_truth=False)
    timing_checks, timing_details = _timing_integrity(timing)
    checks = {
        "all_required_artifacts_exist": all(path.is_file() for path in required),
        "timing_bank_has_no_truth": "reference" not in timing,
        "timing_bank_size": timing["x"].shape == (524288, 100),
        "reference_gate_validated": _read(RESULTS_DIR / "reference_audit.json").get("status") == "VALIDATED",
        "preflight_validated": _read(RESULTS_DIR / "preflight_audit.json").get("status") == "VALIDATED",
        "frozen_selection_hash_valid": bool(require_frozen_selections()),
        "analysis_complete": _read(RESULTS_DIR / "analysis_summary.json").get("status") == "COMPLETE",
        **timing_checks,
    }
    with (RESULTS_DIR / "resource_accounting.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        resource_rows = list(csv.DictReader(handle))
    allocated_gpu_hours = sum(
        float(row["task_wall_seconds"]) * float(row.get("gpu_count") or 0)
        for row in resource_rows
    ) / 3600.0
    checks["within_18_h100_hour_cap"] = allocated_gpu_hours <= 18.0
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(checks.values()) else "FAILED",
        "overall_pass": all(checks.values()),
        "checks": checks,
        "artifacts": {
            str(path.relative_to(SOURCE_ROOT)): file_sha256(path)
            for path in required if path.is_file()
        },
        "timing_input_sha256": array_sha256(timing["t"], timing["x"]),
        "timing_record_audit": timing_details,
        "allocated_h100_gpu_hours": allocated_gpu_hours,
        "audited_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    if not payload["overall_pass"]:
        raise RuntimeError("artifact audit failed")
    return payload


if __name__ == "__main__":
    result = analyze()
    resource_accounting()
    artifact_audit()
    print(json.dumps(result, indent=2, sort_keys=True))

