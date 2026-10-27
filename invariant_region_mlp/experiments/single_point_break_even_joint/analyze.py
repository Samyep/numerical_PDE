"""Measured break-even analysis, figure, resource ledger, and artifact audit."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .protocol import (
    ALLOCATED_GPU_HOUR_CAP, BANK_SIZES, FIGURE_DIR, FINAL_SEEDS,
    FROZEN_SELECTIONS, RESULTS_DIR, REVISION, SOURCE_ROOT, TARGETS,
    array_sha256, file_sha256, load_bank, require_frozen_selections,
    write_json_atomic,
)


METHODS = ("tangent_actor", "fh_mlp_standard", "deep_bsde")
LABELS = {"tangent_actor": "Tangent-Actor", "fh_mlp_standard": "Full-History MLP", "deep_bsde": "Deep BSDE"}
COLORS = {"tangent_actor": "#2166ac", "fh_mlp_standard": "#d95f02", "deep_bsde": "#1b9e77"}


def _read(path: Path) -> dict[str, Any]: return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    temporary.replace(path)


def _records(method: str, candidate: str) -> list[dict[str, Any]]:
    if method == "tangent_actor":
        paths = [RESULTS_DIR / "raw" / "timing" / method / f"ta__timing__s{seed}__{candidate}.json" for seed in FINAL_SEEDS]
    elif method == "fh_mlp_standard":
        n_text, m_text = candidate.split("_"); n, M = int(n_text[1:]), int(m_text[1:])
        paths = [RESULTS_DIR / "raw" / "timing" / "fh_mlp" / f"fh__timing__n{n}_M{M}__r{rep}.json" for rep in range(3)]
    else:
        updates = int(candidate.split("_")[1])
        paths = [RESULTS_DIR / "raw" / "timing" / "deep_bsde" / f"db__timing__u{updates}__r{rep}.json" for rep in range(3)]
    if not all(path.is_file() for path in paths): raise RuntimeError(f"missing timing records for {method}/{candidate}: {paths}")
    return [_read(path) for path in paths]


def _mean_curve(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    mappings = [{int(row["queries"]): float(row["total_seconds"]) for row in record["prefix_rows"]} for record in records]
    common = sorted(set.intersection(*(set(mapping) for mapping in mappings)))
    rows = []
    for queries in common:
        values = np.asarray([mapping[queries] for mapping in mappings])
        rows.append({"queries": queries, "mean_seconds": float(values.mean()), "std_seconds": float(values.std(ddof=1)), "min_seconds": float(values.min()), "max_seconds": float(values.max()), "replicate_seconds": values.tolist()})
    return rows


def _crossing(ta: list[dict[str, Any]], baseline: list[dict[str, Any]]) -> dict[str, Any]:
    left = {row["queries"]: row for row in ta}; right = {row["queries"]: row for row in baseline}
    common = sorted(set(left) & set(right))
    if not common: return {"status": "NO_COMMON_EXECUTED_PREFIX"}
    replicate_crossings = []
    for replicate in range(3):
        crossed = None; previous = None
        for query in common:
            if right[query]["replicate_seconds"][replicate] >= left[query]["replicate_seconds"][replicate]: crossed = query; break
            previous = query
        replicate_crossings.append({"replicate": replicate, "previous_queries": previous if crossed is not None else common[-1], "crossing_queries": crossed, "status": "OBSERVED_BREAK_EVEN" if crossed is not None else "NOT_REACHED_WITHIN_TESTED_BUDGET"})
    previous = None
    for query in common:
        if right[query]["mean_seconds"] >= left[query]["mean_seconds"]:
            return {"status": "OBSERVED_BREAK_EVEN", "previous_queries": previous, "crossing_queries": query, "tangent_actor_mean_seconds_at_crossing": left[query]["mean_seconds"], "tangent_actor_std_seconds_at_crossing": left[query]["std_seconds"], "baseline_mean_seconds_at_crossing": right[query]["mean_seconds"], "baseline_std_seconds_at_crossing": right[query]["std_seconds"], "per_replicate_crossings": replicate_crossings}
        previous = query
    query = common[-1]
    return {"status": "NOT_REACHED_WITHIN_TESTED_BUDGET", "largest_common_queries": query, "tangent_actor_mean_seconds": left[query]["mean_seconds"], "baseline_mean_seconds": right[query]["mean_seconds"], "baseline_over_tangent_actor_time": right[query]["mean_seconds"] / left[query]["mean_seconds"], "per_replicate_crossings": replicate_crossings}


def _plot(curves: dict[tuple[str, str], list[dict[str, Any]]], final: dict[str, Any], coverage: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    # The secondary target has no timed baseline under the frozen compute cap.
    # Plotting a TA-only panel creates a visually dominant but scientifically
    # empty axis, so the publication figure shows the primary matched-accuracy
    # comparison and leaves the secondary outcome in the report/table.
    target = "epsilon_1"
    threshold = TARGETS[target]
    fig, axis = plt.subplots(figsize=(7.25, 4.55), constrained_layout=True)
    for method in METHODS:
        item = final["targets"][target][method]
        if coverage["targets"][target][method] != "RUN" or item.get("selected_candidate") is None:
            continue
        curve = curves[(method, item["selected_candidate"])]
        x = np.asarray([row["queries"] for row in curve])
        mean = np.asarray([row["mean_seconds"] for row in curve])
        std = np.asarray([row["std_seconds"] for row in curve])
        positive = np.asarray([value for row in curve for value in row["replicate_seconds"] if value > 0])
        floor = float(positive.min()) if len(positive) else 1e-12
        axis.plot(x, mean, marker="o", ms=4.2, lw=2.0, color=COLORS[method], label=LABELS[method], zorder=3)
        axis.fill_between(x, np.maximum(mean - std, floor), mean + std, color=COLORS[method], alpha=0.14, linewidth=0, zorder=1)

    ta_item = final["targets"][target]["tangent_actor"]
    fh_item = final["targets"][target]["fh_mlp_standard"]
    crossing = _crossing(
        curves[("tangent_actor", ta_item["selected_candidate"])],
        curves[("fh_mlp_standard", fh_item["selected_candidate"])],
    )
    if crossing["status"] == "OBSERVED_BREAK_EVEN":
        previous = crossing["previous_queries"]
        observed = crossing["crossing_queries"]
        axis.axvspan(previous, observed, color="#6f6f6f", alpha=0.09, zorder=0)
        axis.axvline(observed, color="#555555", lw=1.0, ls=":", zorder=2)
        axis.text(
            np.sqrt(previous * observed), 0.035, "observed break-even interval",
            transform=axis.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=8.2, color="#444444",
        )

    axis.set_xscale("log", base=2)
    axis.set_yscale("log")
    ticks = [2**power for power in (0, 4, 8, 12, 16, 20)]
    axis.set_xticks(ticks, [rf"$2^{{{power}}}$" for power in (0, 4, 8, 12, 16, 20)])
    axis.set_xlim(0.8, 1.16 * 2**20)
    axis.set_xlabel("Distinct initial-state queries")
    axis.set_ylabel("Measured total H100 time (s)")
    axis.set_title(rf"Joint nonseparable, $d=100$  |  RMSE $\leq {threshold:.0e}$", pad=10)
    axis.grid(True, which="both", alpha=0.20, linewidth=0.6)
    axis.legend(loc="center left", bbox_to_anchor=(0.025, 0.52), frameon=False)
    axis.text(
        0.025, 0.965, "Deep BSDE: target not reached",
        transform=axis.transAxes, ha="left", va="top", fontsize=8.6,
        color="#444444",
    )
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURE_DIR / f"joint_matched_accuracy_break_even.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def analyze() -> dict[str, Any]:
    frozen = require_frozen_selections(); final = _read(RESULTS_DIR / "final_accuracy_summary.json"); coverage = _read(RESULTS_DIR / "timing_coverage.json")
    curves: dict[tuple[str, str], list[dict[str, Any]]] = {}; timing_rows = []
    for target in TARGETS:
        for method in METHODS:
            item = final["targets"][target][method]
            if coverage["targets"][target][method] != "RUN" or not item.get("selected_candidate"): continue
            key = (method, item["selected_candidate"])
            if key not in curves: curves[key] = _mean_curve(_records(*key))
    for (method, candidate), rows in curves.items():
        for row in rows:
            timing_rows.append({"method": method, "candidate": candidate, "queries": row["queries"], "mean_seconds": row["mean_seconds"], "std_seconds": row["std_seconds"], "min_seconds": row["min_seconds"], "max_seconds": row["max_seconds"]})
    _write_csv(RESULTS_DIR / "timing_rows.csv", timing_rows)
    break_even: dict[str, Any] = {}
    flat = []
    for target in TARGETS:
        break_even[target] = {}
        ta_item = final["targets"][target]["tangent_actor"]
        for method in ("fh_mlp_standard", "deep_bsde"):
            baseline_item = final["targets"][target][method]
            if coverage["targets"][target]["tangent_actor"] != "RUN" or coverage["targets"][target][method] != "RUN": result = {"status": coverage["targets"][target][method] if coverage["targets"][target][method] != "RUN" else coverage["targets"][target]["tangent_actor"]}
            else: result = _crossing(curves[("tangent_actor", ta_item["selected_candidate"])], curves[(method, baseline_item["selected_candidate"])])
            break_even[target][method] = result; flat.append({"target": target, "baseline": method, **result})
    _write_csv(RESULTS_DIR / "break_even_summary.csv", flat)
    _plot(curves, final, coverage)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True).strip()
    payload = {"revision": REVISION, "status": "COMPLETE", "analysis_code_commit": commit, "frozen_selections_canonical_sha256": frozen["canonical_sha256"], "targets": final["targets"], "break_even": break_even, "timing_coverage": coverage, "curve_keys": [list(key) for key in curves], "analyzed_utc": datetime.now(timezone.utc).isoformat()}
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", payload); return payload


def resource_accounting() -> list[dict[str, Any]]:
    rows = []
    preflight_path = RESULTS_DIR / "preflight_audit.json"
    if preflight_path.is_file():
        preflight = _read(preflight_path)
        rows.append({
            "task_id": "h100_preflight",
            "stage": "preflight",
            "kind": "mathematical_and_source_preflight",
            "status": "complete" if preflight.get("status") == "VALIDATED" else "failed",
            "task_wall_seconds": float(preflight.get("preflight_wall_seconds", 0.0)),
            "slurm_job_id": None,
            "slurm_array_task_id": None,
            "gpu_count": 1,
            "error": None,
        })
    for path in sorted((RESULTS_DIR / "raw" / "task_status").glob("*.json")):
        item = _read(path); rows.append({"task_id": item["task_id"], "stage": item["stage"], "kind": item["task"]["kind"], "status": item["status"], "task_wall_seconds": item["task_wall_seconds"], "slurm_job_id": item.get("slurm_job_id"), "slurm_array_task_id": item.get("slurm_array_task_id"), "gpu_count": item.get("gpu_count", 1), "error": item.get("error")})
    path = RESULTS_DIR / "infrastructure_attempts.json"
    if path.is_file(): rows.extend(_read(path).get("attempts", []))
    _write_csv(RESULTS_DIR / "resource_accounting.csv", rows); return rows


def artifact_audit() -> dict[str, Any]:
    required = [RESULTS_DIR / name for name in ("bank_manifest.json", "reference_audit.json", "preflight_audit.json", "calibration_rows.csv", "FROZEN_SELECTIONS.json", "final_accuracy_rows.csv", "final_accuracy_summary.json", "timing_coverage.json", "timing_rows.csv", "break_even_summary.csv", "resource_accounting.csv", "analysis_summary.json")]
    required += [FIGURE_DIR / "joint_matched_accuracy_break_even.png", FIGURE_DIR / "joint_matched_accuracy_break_even.pdf"]
    timing = load_bank("timing", require_truth=False); final = _read(RESULTS_DIR / "final_accuracy_summary.json"); coverage = _read(RESULTS_DIR / "timing_coverage.json")
    detail = []; selected = set()
    for target in TARGETS:
        for method in METHODS:
            item = final["targets"][target][method]
            if coverage["targets"][target][method] == "RUN": selected.add((method, item["selected_candidate"]))
    for method, candidate in sorted(selected):
        records = _records(method, candidate)
        for replicate, record in enumerate(records):
            completed = int(record["completed_queries"]); prediction_path = Path(record["prediction_file"])
            with np.load(prediction_path, allow_pickle=False) as data: prediction = np.asarray(data["prediction"], dtype=np.float64)
            counts = [int(row["queries"]) for row in record["prefix_rows"]]; seconds = [float(row["total_seconds"]) for row in record["prefix_rows"]]
            checks = {"actual_prefix_execution": record.get("actual_prefix_execution") is True, "truth_not_loaded": record.get("truth_loaded") is False, "h100": "H100" in str(record.get("gpu", "")).upper(), "valid_status": record.get("status") in {"complete", "TIME_CAP"}, "no_nonfinite": int(record.get("nonfinite_count", -1)) == 0, "input_prefix_hash": record.get("query_input_sha256") == array_sha256(timing["t"][:completed], timing["x"][:completed]), "prefixes_monotone": bool(counts) and counts[-1] <= completed and all(b > a for a, b in zip(counts, counts[1:])), "times_monotone": bool(seconds) and all(b >= a for a, b in zip(seconds, seconds[1:])), "prediction_shape": prediction.shape == (completed,), "prediction_finite": bool(np.isfinite(prediction).all()), "prediction_hash": record.get("prediction_sha256") == array_sha256(prediction)}
            detail.append({"method": method, "candidate": candidate, "replicate": replicate, "checks": checks, "overall_pass": all(checks.values())})
    resources = resource_accounting(); gpu_hours = sum(float(row["task_wall_seconds"]) * float(row.get("gpu_count") or 0) for row in resources) / 3600.0
    checks = {"all_required_artifacts_exist": all(path.is_file() for path in required), "timing_bank_has_no_truth": "reference" not in timing, "timing_bank_shape": timing["x"].shape == (BANK_SIZES["timing"], 100), "reference_validated": _read(RESULTS_DIR / "reference_audit.json").get("status") == "VALIDATED", "preflight_validated": _read(RESULTS_DIR / "preflight_audit.json").get("status") == "VALIDATED", "frozen_selection_hash_valid": bool(require_frozen_selections()), "analysis_complete": _read(RESULTS_DIR / "analysis_summary.json").get("status") == "COMPLETE", "three_replicates_per_timed_configuration": all(len(_records(*key)) == 3 for key in selected), "timing_records_integrity": bool(detail) and all(row["overall_pass"] for row in detail), "within_18_h100_hour_cap": gpu_hours <= ALLOCATED_GPU_HOUR_CAP}
    payload = {"revision": REVISION, "status": "VALIDATED" if all(checks.values()) else "FAILED", "overall_pass": all(checks.values()), "checks": checks, "allocated_h100_gpu_hours": gpu_hours, "timing_record_audit": detail, "artifacts": {str(path.relative_to(SOURCE_ROOT)): file_sha256(path) for path in required if path.is_file()}, "audited_utc": datetime.now(timezone.utc).isoformat()}
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", payload)
    if not payload["overall_pass"]: raise RuntimeError("artifact audit failed")
    return payload


if __name__ == "__main__":
    result = analyze(); resource_accounting(); audit = artifact_audit(); print(json.dumps({"analysis": result["status"], "audit": audit["status"]}, indent=2))
