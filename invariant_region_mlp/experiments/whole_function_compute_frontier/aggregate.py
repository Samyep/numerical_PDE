"""Aggregate measured checkpoints and perform the artifact/fidelity audit."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import array_sha256

from .common import file_sha256, write_json_atomic
from .prepare import historical_results_root
from .protocol import (
    BASE_REVISION,
    C0_CONSISTENCY_ABSOLUTE_RMSE,
    C0_CONSISTENCY_RELATIVE_RMSE,
    CONFIGURATIONS,
    GLOBAL_GPU_SECONDS_CAP,
    METHODS,
    PDES,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    STATUS_DIR,
    UPSTREAM_DPI_REVISION,
    protocol_payload,
)
from .run_frontier import training_inventory


STANDARD_FIELDS = (
    "revision", "pde", "dimension", "method", "configuration_id", "seed",
    "trajectory_id", "checkpoint_id", "algorithm_iteration",
    "training_seconds", "selection_seconds", "main_axis_seconds",
    "inference_seconds", "inference_seconds_mean", "total_compute_seconds",
    "validation_metric", "rmse", "rmse_sd", "relmae", "relative_l2",
    "mae", "bias", "n_inference_replicates", "causal_selected_checkpoint",
    "is_original_budget",
    "is_final", "checkpoint", "checkpoint_sha256", "prediction_files",
    "prediction_sha256s", "test_bank_sha256", "test_input_sha256",
    "truth_value_sha256", "method_source_revision", "nonfinite_count", "status",
)


def _write_csv(path: Path, rows: list[dict[str, Any]], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _metric_rows() -> list[dict[str, Any]]:
    result = []
    for path in sorted(PREDICTION_DIR.rglob("*.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        if row.get("status") == "COMPLETE" and "rmse" in row:
            row["_metric_path"] = str(path)
            result.append(row)
    return result


def _aggregate_replicates(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for row in rows:
        key = (
            row["pde"], row["method"], row["configuration_id"],
            row["checkpoint_id"], int(row["algorithm_iteration"]),
        )
        groups.setdefault(key, []).append(row)
    output: list[dict[str, Any]] = []
    for key, members in sorted(groups.items()):
        first = members[0]
        numeric_means = ("rmse", "relmae", "relative_l2", "mae", "bias")
        item = {name: first.get(name) for name in STANDARD_FIELDS}
        for name in numeric_means:
            values = [float(row[name]) for row in members]
            item[name] = float(np.mean(values))
        rmse_values = [float(row["rmse"]) for row in members]
        item.update({
            "revision": REVISION,
            "dimension": 100,
            "seed": 0,
            "training_seconds": float(first["training_seconds"]),
            "selection_seconds": float(first.get("selection_seconds", 0.0)),
            "main_axis_seconds": float(first["main_axis_seconds"]),
            "inference_seconds": float(sum(float(row["inference_seconds"]) for row in members)),
            "inference_seconds_mean": float(np.mean([float(row["inference_seconds"]) for row in members])),
            "rmse_sd": float(np.std(rmse_values, ddof=1)) if len(members) > 1 else 0.0,
            "n_inference_replicates": len(members),
            "nonfinite_count": int(sum(int(row.get("nonfinite_count", 0)) for row in members)),
            "prediction_files": ";".join(str(row["prediction_file"]) for row in members),
            "prediction_sha256s": ";".join(str(row["prediction_sha256"]) for row in members),
            "method_source_revision": (
                BASE_REVISION if first["method"] == "tangent_actor"
                else UPSTREAM_DPI_REVISION
            ),
            "status": "COMPLETE",
        })
        item["total_compute_seconds"] = item["main_axis_seconds"] + item["inference_seconds"]
        output.append(item)
    return output


def _trajectory_records() -> list[dict[str, Any]]:
    records = []
    for row in training_inventory():
        path = Path(row["path"])
        if not path.is_file():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        records.append(payload)
    return records


def _trajectory_summary(
    trajectories: list[dict[str, Any]], metrics: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    lookup: dict[str, list[dict[str, Any]]] = {}
    for row in metrics:
        lookup.setdefault(str(row["trajectory_id"]), []).append(row)
    output = []
    for trajectory in trajectories:
        rows = lookup.get(str(trajectory["trajectory_id"]), [])
        original = [row for row in rows if str(row.get("is_original_budget")).lower() == "true"]
        final = max(rows, key=lambda row: float(row["main_axis_seconds"])) if rows else None
        output.append({
            "trajectory_id": trajectory["trajectory_id"],
            "pde": trajectory["pde"], "method": trajectory["method"],
            "configuration_id": trajectory["configuration_id"],
            "seed": trajectory["seed"], "status": trajectory["status"],
            "stop_reason": trajectory.get("stop_reason"),
            "checkpoint_count": len(trajectory.get("checkpoints", [])),
            "reported_checkpoint_count": len(rows),
            "training_seconds": trajectory.get("training_seconds"),
            "selection_seconds": trajectory.get("selection_seconds", 0.0),
            "final_rmse": final.get("rmse") if final else None,
            "original_budget_rmse": original[0].get("rmse") if len(original) == 1 else None,
            "peak_gpu_bytes": trajectory.get("peak_gpu_bytes"),
            "failure": trajectory.get("failure"),
        })
    return output


def _elapsed_seconds(text: str) -> float:
    # Slurm format is [[days-]hours:]minutes:seconds.
    days = 0
    if "-" in text:
        day_text, text = text.split("-", 1)
        days = int(day_text)
    parts = [float(value) for value in text.split(":")]
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours, minutes, seconds = 0.0, *parts
    else:
        hours, minutes, seconds = 0.0, 0.0, parts[0]
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def _compute_accounting(trajectories: list[dict[str, Any]]) -> list[dict[str, Any]]:
    slurm = RESULTS_DIR / "slurm_accounting_raw.csv"
    rows: list[dict[str, Any]] = []
    if slurm.is_file():
        with slurm.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                elapsed = _elapsed_seconds(row["elapsed"])
                gpu_count = int(row.get("gpu_count") or 0)
                row["job_wall_seconds"] = elapsed
                row["allocated_gpu_seconds"] = elapsed * gpu_count
                rows.append(row)
    active_total = sum(float(row.get("training_seconds", 0.0)) for row in trajectories)
    allocated_total = sum(float(row.get("allocated_gpu_seconds", 0.0)) for row in rows)
    rows.append({
        "job_id": "TOTAL", "job_name": "all_wcf_jobs",
        "state": "SUMMARY", "gpu_count": "",
        "elapsed": "", "job_wall_seconds": "",
        "allocated_gpu_seconds": allocated_total,
        "active_training_seconds": active_total,
        "allocated_gpu_hours": allocated_total / 3600.0,
        "active_training_hours": active_total / 3600.0,
        "within_18_gpu_hour_cap": allocated_total <= GLOBAL_GPU_SECONDS_CAP,
    })
    return rows


def _overlap(metrics: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {"definition": "intersection of observed C0 main-axis time ranges"}
    for pde in PDES:
        ranges = {}
        for method in METHODS:
            times = [
                float(row["main_axis_seconds"]) / 60.0 for row in metrics
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == "C0"
            ]
            ranges[method] = [min(times), max(times)] if times else None
        valid = [value for value in ranges.values() if value is not None]
        low = max(value[0] for value in valid) if len(valid) == len(METHODS) else None
        high = min(value[1] for value in valid) if len(valid) == len(METHODS) else None
        result[pde] = {
            "ranges_minutes": ranges,
            "intersection_minutes": [low, high] if low is not None and low <= high else None,
            "meaningful_observed_overlap": bool(low is not None and low <= high),
            "no_interpolation_used": True,
        }
    return result


def _historical_rows() -> list[dict[str, str]]:
    path = historical_results_root() / "main_whole_function_results.csv"
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _historical_rmse(row: dict[str, str]) -> float | None:
    for key in ("value_rmse", "rmse", "Value RMSE"):
        value = row.get(key)
        if value not in (None, "", "NA", "nan"):
            return float(value)
    return None


def _fidelity(
    trajectories: list[dict[str, Any]], metrics: list[dict[str, Any]], accounting: list[dict[str, Any]]
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    monotone_failures = []
    checkpoint_failures = []
    for trajectory in trajectories:
        times = [float(row["main_axis_seconds"]) for row in trajectory.get("checkpoints", [])]
        if any(right < left for left, right in zip(times, times[1:])):
            monotone_failures.append(trajectory["trajectory_id"])
        for checkpoint in trajectory.get("checkpoints", []):
            path_text = checkpoint.get("checkpoint")
            digest = checkpoint.get("checkpoint_sha256")
            if path_text and digest:
                path = Path(path_text)
                if not path.is_file() or file_sha256(path) != digest:
                    checkpoint_failures.append(checkpoint["checkpoint_id"])
    checks.append({"id": 1, "name": "monotone checkpoint times", "pass": not monotone_failures,
                   "failures": monotone_failures})
    checks.append({"id": 2, "name": "checkpoint weights match recorded hashes", "pass": not checkpoint_failures,
                   "failures": checkpoint_failures})
    prediction_failures = []
    for row in metrics:
        paths = str(row["prediction_files"]).split(";")
        hashes = str(row["prediction_sha256s"]).split(";")
        if len(paths) != len(hashes) or any(
            not Path(path).is_file() or file_sha256(Path(path)) != digest
            for path, digest in zip(paths, hashes)
        ):
            prediction_failures.append(f"{row['trajectory_id']}:{row['checkpoint_id']}")
    checks.append({"id": 3, "name": "predictions bind to recorded checkpoints/artifacts",
                   "pass": not prediction_failures, "failures": prediction_failures})
    ta_missing = [
        row["trajectory_id"] for row in trajectories
        if row["method"] == "tangent_actor" and not row.get("phases")
    ]
    checks.append({"id": 4, "name": "TA timing includes label generation and fitting",
                   "pass": not ta_missing, "failures": ta_missing})
    truth_failures = [
        row["trajectory_id"] for row in trajectories
        if row.get("ground_truth_used_for_training_selection_or_stopping") is not False
    ]
    gate = STATUS_DIR / "all_training_complete.json"
    checks.append({"id": 5, "name": "truth-free training/selection and global truth gate",
                   "pass": not truth_failures and gate.is_file(), "failures": truth_failures})

    historical = _historical_rows()
    consistency = []
    aliases = {
        "tangent_actor": {"tangent_actor", "TA", "ta"},
        "dpi": {"dpi", "DPI"},
        "pinn_hte": {"pinn_hte", "PINN-HTE", "pinn"},
    }
    for pde in PDES:
        for method in METHODS:
            candidates = [
                row for row in metrics
                if row["pde"] == pde and row["method"] == method
                and row["configuration_id"] == "C0" and bool(row.get("is_original_budget"))
            ]
            if method == "tangent_actor" and len(candidates) == 1:
                selected_id = candidates[0].get("causal_selected_checkpoint")
                selected = [
                    row for row in metrics
                    if row["pde"] == pde and row["method"] == method
                    and row["configuration_id"] == "C0"
                    and row["checkpoint_id"] == selected_id
                ]
                if len(selected) == 1:
                    candidates = selected
            old_candidates = [
                row for row in historical
                if row.get("pde") == pde
                and row.get("method") in aliases[method]
                and int(float(row.get("dimension", row.get("d", 0)))) == 100
            ]
            new = float(candidates[0]["rmse"]) if len(candidates) == 1 else None
            old = _historical_rmse(old_candidates[0]) if len(old_candidates) == 1 else None
            if new is None or old is None:
                passed = False
                tolerance = None
            else:
                tolerance = max(C0_CONSISTENCY_ABSOLUTE_RMSE, C0_CONSISTENCY_RELATIVE_RMSE * old)
                passed = abs(new - old) <= tolerance
            consistency.append({
                "pde": pde, "method": method, "new_rmse": new,
                "historical_rmse": old, "tolerance": tolerance, "pass": passed,
                "absolute_difference": (
                    abs(new - old) if new is not None and old is not None else None
                ),
                "relative_difference_to_historical": (
                    abs(new - old) / abs(old)
                    if new is not None and old not in (None, 0.0) else None
                ),
                "new_rmse_sd_across_inference_replicates": (
                    float(candidates[0].get("rmse_sd", 0.0))
                    if len(candidates) == 1 else None
                ),
                "new_inference_replicates": (
                    int(candidates[0].get("n_inference_replicates", 0))
                    if len(candidates) == 1 else None
                ),
                "new_training_plus_selection_seconds": (
                    float(candidates[0].get("main_axis_seconds", 0.0))
                    if len(candidates) == 1 else None
                ),
                "historical_training_seconds": (
                    float(old_candidates[0]["training_seconds"])
                    if len(old_candidates) == 1
                    and old_candidates[0].get("training_seconds") not in (None, "")
                    else None
                ),
                "uncertainty_scope": (
                    "inference-only; one training seed in both studies"
                    if method == "tangent_actor"
                    else "no replicate uncertainty estimate; one training seed in both studies"
                ),
            })
    checks.append({"id": 6, "name": "original-budget C0 consistency",
                   "pass": all(row["pass"] for row in consistency), "details": consistency})
    recorded_failures = [
        row["trajectory_id"] for row in trajectories if row.get("failure")
    ]
    incomplete_trajectories = [
        row["trajectory_id"] for row in trajectories
        if row.get("status") != "COMPLETE"
    ]
    trajectory_ids = {str(row["trajectory_id"]) for row in trajectories}
    missing_metric_trajectories = sorted(
        trajectory_ids - {str(row["trajectory_id"]) for row in metrics}
    )
    nonfinite_metrics = []
    for row in metrics:
        label = f"{row['trajectory_id']}:{row['checkpoint_id']}"
        numeric = ("rmse", "relmae", "relative_l2", "mae", "bias")
        if int(row.get("nonfinite_count", 0)) > 0 or any(
            not math.isfinite(float(row[name])) for name in numeric
        ):
            nonfinite_metrics.append(label)
    checks.append({
        "id": 7,
        "name": "completion and numerical-failure records agree with predictions",
        "pass": not (
            incomplete_trajectories or missing_metric_trajectories or nonfinite_metrics
        ),
        "details": {
            "recorded_failure_payloads": recorded_failures,
            "incomplete_trajectories": incomplete_trajectories,
            "trajectories_without_reported_metrics": missing_metric_trajectories,
            "nonfinite_prediction_or_metric_rows": nonfinite_metrics,
        },
    })
    totals = [row for row in accounting if row.get("job_id") == "TOTAL"]
    within = bool(totals and str(totals[0].get("within_18_gpu_hour_cap")).lower() == "true")
    checkpoint_timing_failures = []
    metric_timing_failures = []
    checkpoints_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    for trajectory in trajectories:
        checkpoints = trajectory.get("checkpoints", [])
        if checkpoints:
            last = checkpoints[-1]
            if not math.isclose(
                float(trajectory.get("training_seconds", 0.0)),
                float(last.get("training_seconds", 0.0)), rel_tol=0.0, abs_tol=1e-6,
            ) or not math.isclose(
                float(trajectory.get("selection_seconds", 0.0)),
                float(last.get("selection_seconds", 0.0)), rel_tol=0.0, abs_tol=1e-6,
            ):
                checkpoint_timing_failures.append(f"{trajectory['trajectory_id']}:terminal")
        for checkpoint in checkpoints:
            expected_main = (
                float(checkpoint.get("training_seconds", 0.0))
                + float(checkpoint.get("selection_seconds", 0.0))
            )
            if not math.isclose(
                expected_main, float(checkpoint.get("main_axis_seconds", 0.0)),
                rel_tol=0.0, abs_tol=1e-6,
            ):
                checkpoint_timing_failures.append(
                    f"{trajectory['trajectory_id']}:{checkpoint['checkpoint_id']}"
                )
            checkpoints_by_key[(str(trajectory["trajectory_id"]), str(checkpoint["checkpoint_id"]))] = checkpoint
    for row in metrics:
        key = (str(row["trajectory_id"]), str(row["checkpoint_id"]))
        checkpoint = checkpoints_by_key.get(key)
        if checkpoint is None or not math.isclose(
            float(row["main_axis_seconds"]), float(checkpoint["main_axis_seconds"]),
            rel_tol=0.0, abs_tol=1e-6,
        ):
            metric_timing_failures.append(":".join(key))
    active_from_trajectories = sum(
        float(row.get("training_seconds", 0.0)) for row in trajectories
    )
    active_from_accounting = (
        float(totals[0].get("active_training_seconds", 0.0)) if totals else None
    )
    allocated = (
        float(totals[0].get("allocated_gpu_seconds", 0.0)) if totals else None
    )
    active_reconciled = bool(
        totals and math.isclose(
            active_from_trajectories, float(active_from_accounting),
            rel_tol=0.0, abs_tol=1e-6,
        )
    )
    allocated_covers_active = bool(
        totals and float(allocated) >= float(active_from_accounting)
    )
    live_gpu_jobs = [
        str(row.get("job_id")) for row in accounting
        if row.get("job_id") != "TOTAL"
        and int(row.get("gpu_count") or 0) > 0
        and str(row.get("state", "")).upper() in {"RUNNING", "PENDING", "CONFIGURING"}
    ]
    retry_gpu_seconds = sum(
        float(row.get("allocated_gpu_seconds") or 0.0) for row in accounting
        if row.get("job_id") != "TOTAL"
        and str(row.get("state", "")).upper().startswith(("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY"))
    )
    accounting_pass = bool(
        within and active_reconciled and allocated_covers_active
        and not checkpoint_timing_failures and not metric_timing_failures
        and not live_gpu_jobs
    )
    checks.append({
        "id": 8,
        "name": "checkpoint timing components and scheduler accounting reconcile",
        "pass": accounting_pass,
        "details": {
            "totals": totals,
            "checkpoint_timing_failures": checkpoint_timing_failures,
            "metric_timing_failures": metric_timing_failures,
            "active_seconds_from_trajectories": active_from_trajectories,
            "active_seconds_from_accounting": active_from_accounting,
            "active_seconds_reconciled": active_reconciled,
            "allocated_seconds_cover_active_training": allocated_covers_active,
            "failed_or_cancelled_retry_gpu_seconds_included": retry_gpu_seconds,
            "live_gpu_jobs_at_accounting_snapshot": live_gpu_jobs,
            "within_18_gpu_hour_cap": within,
        },
    })
    protocol = protocol_payload()
    expected_cells = {
        (pde, method, configuration)
        for pde in PDES for method in METHODS for configuration in CONFIGURATIONS
    }
    observed_cells = {
        (str(row["pde"]), str(row["method"]), str(row["configuration_id"]))
        for row in trajectories
    }
    protocol_failures = []
    for row in trajectories:
        expected_grid = protocol["configuration_grid"][row["method"]][row["configuration_id"]]
        snapshot = row.get("configuration", {})
        if row["method"] == "tangent_actor":
            # TA records the effective mature configuration directly, whereas
            # the imported baselines wrap their changed fields under
            # ``frontier_changes``.  Compare the frozen numerical knobs in the
            # representation actually emitted by each runner.
            mismatch = any(
                snapshot.get(key) != expected_grid[key]
                for key in ("paths_per_label", "optimizer_steps")
            )
        else:
            mismatch = snapshot.get("frontier_changes") != expected_grid
        if mismatch:
            protocol_failures.append(row["trajectory_id"])
    input_hash_failures = []
    direct_hash_rows = 0
    embedded_array_hash_rows = 0
    for row in metrics:
        expected = protocol["expected_inputs"][row["pde"]]
        direct = (
            row.get("test_input_sha256") == expected["test_input_sha256"]
            and row.get("truth_value_sha256") == expected["truth_value_sha256"]
        )
        if direct:
            direct_hash_rows += 1
            continue
        embedded_valid = (
            row.get("test_bank_sha256") == expected["test_bank_file_sha256"]
        )
        for path_text in str(row.get("prediction_files", "")).split(";"):
            path = Path(path_text)
            if not path.is_file():
                embedded_valid = False
                break
            with np.load(path, allow_pickle=False) as payload:
                if not {"t", "x", "truth"}.issubset(payload.files):
                    embedded_valid = False
                    break
                embedded_valid = embedded_valid and (
                    array_sha256(np.asarray(payload["t"]), np.asarray(payload["x"]))
                    == expected["test_input_sha256"]
                    and array_sha256(np.asarray(payload["truth"]))
                    == expected["truth_value_sha256"]
                )
        if embedded_valid:
            embedded_array_hash_rows += 1
        else:
            input_hash_failures.append(
                f"{row['trajectory_id']}:{row['checkpoint_id']}"
            )
    checks.append({
        "id": 9, "name": "common frozen protocols and input hashes on both PDEs",
        "pass": (
            observed_cells == expected_cells
            and not protocol_failures and not input_hash_failures
        ),
        "protocol_sha256": protocol["protocol_sha256"],
        "details": {
            "missing_cells": sorted(expected_cells - observed_cells),
            "unexpected_cells": sorted(observed_cells - expected_cells),
            "configuration_snapshot_failures": protocol_failures,
            "input_hash_failures": input_hash_failures,
            "direct_input_hash_rows": direct_hash_rows,
            "embedded_prediction_array_hash_rows": embedded_array_hash_rows,
        },
    })
    submitted = RESULTS_DIR / "submitted_jobs.json"
    submitted_payload: dict[str, Any] = {}
    if submitted.is_file():
        submitted_payload = json.loads(submitted.read_text(encoding="utf-8"))
    submitted_jobs = submitted_payload.get("jobs", {})
    manifest_says_no_sdgd = submitted_payload.get("no_sdgd_jobs") is True
    submitted_job_mentions_sdgd = "sdgd" in json.dumps(
        {key: value for key, value in submitted_payload.items() if key != "no_sdgd_jobs"},
        sort_keys=True,
    ).lower()
    accounted_job_mentions_sdgd = any(
        "sdgd" in str(row.get("job_name", "")).lower()
        for row in accounting
    )
    checks.append({"id": 10, "name": "no new SDGD-PINN job submitted",
                   "pass": (
                       manifest_says_no_sdgd
                       and not submitted_job_mentions_sdgd
                       and not accounted_job_mentions_sdgd
                   ),
                   "details": {
                       "manifest_no_sdgd_jobs": manifest_says_no_sdgd,
                       "submitted_job_mentions_sdgd": submitted_job_mentions_sdgd,
                       "accounted_job_mentions_sdgd": accounted_job_mentions_sdgd,
                   }})
    diagnostic_path = RESULTS_DIR / "diagnostics" / "dense_dpi_callback_replay.json"
    post_freeze_diagnostics: list[dict[str, Any]] = []
    if diagnostic_path.is_file():
        diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
        prediction_path = Path(diagnostic.get("prediction_file", ""))
        prediction_hash_valid = bool(
            prediction_path.is_file()
            and file_sha256(prediction_path) == diagnostic.get("prediction_sha256")
        )
        post_freeze_diagnostics.append({
            "name": "Dense Fourier DPI historical-callback RNG replay",
            "classification": diagnostic.get("classification"),
            "status": diagnostic.get("status"),
            "historical_rmse": diagnostic.get("historical_rmse"),
            "frontier_rmse": next(
                (
                    row["new_rmse"] for row in consistency
                    if row["pde"] == "dense_fourier" and row["method"] == "dpi"
                ),
                None,
            ),
            "rng_replay_rmse": diagnostic.get("replay_rmse"),
            "relative_rmse_difference_from_historical": diagnostic.get("relative_rmse_difference"),
            "rng_replay_callback_calls": diagnostic.get("rng_only_callback_calls"),
            "prediction_artifact_hash_valid": prediction_hash_valid,
            "excluded_from_frontier": True,
            "changes_frozen_check_6_verdict": False,
        })
    return {
        "revision": REVISION,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "overall_pass": all(bool(row["pass"]) for row in checks),
        "checks": checks,
        "post_freeze_diagnostics": post_freeze_diagnostics,
    }


def _write_fidelity_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Fidelity audit", "",
        f"Overall: **{'PASS' if payload['overall_pass'] else 'FAIL'}**", "",
        "The checks below are mechanical; a failed check is not silently waived.", "",
        "| # | Check | Result |", "|---:|---|---|",
    ]
    for row in payload["checks"]:
        lines.append(f"| {row['id']} | {row['name']} | {'PASS' if row['pass'] else 'FAIL'} |")
    if payload.get("post_freeze_diagnostics"):
        lines += [
            "", "## Post-freeze diagnostics", "",
            "These diagnostics explain failures but never change a frozen verdict or enter a frontier.", "",
        ]
        for diagnostic in payload["post_freeze_diagnostics"]:
            lines.append(
                f"- {diagnostic['name']}: historical RMSE `{diagnostic['historical_rmse']:.9g}`, "
                f"frontier RMSE `{diagnostic['frontier_rmse']:.9g}`, RNG-replay RMSE "
                f"`{diagnostic['rng_replay_rmse']:.9g}`; frozen Check 6 remains unchanged."
            )
    lines += ["", "```json", json.dumps(payload, indent=2, sort_keys=True), "```", ""]
    (RESULTS_DIR / "fidelity_audit.md").write_text("\n".join(lines), encoding="utf-8")


def aggregate() -> dict[str, Any]:
    inventory = training_inventory()
    if len(inventory) != 18 or any(row["status"] != "COMPLETE" for row in inventory):
        raise RuntimeError(f"all evaluated trajectories are required: {inventory}")
    raw_metrics = _metric_rows()
    metrics = _aggregate_replicates(raw_metrics)
    trajectories = _trajectory_records()
    summaries = _trajectory_summary(trajectories, metrics)
    accounting = _compute_accounting(trajectories)
    overlap = _overlap(metrics)
    fidelity = _fidelity(trajectories, metrics, accounting)
    _write_csv(RESULTS_DIR / "checkpoint_metrics.csv", metrics, STANDARD_FIELDS)
    _write_csv(RESULTS_DIR / "trajectory_summary.csv", summaries)
    _write_csv(RESULTS_DIR / "compute_accounting.csv", accounting)
    write_json_atomic(RESULTS_DIR / "overlap_audit.json", overlap)
    write_json_atomic(RESULTS_DIR / "artifact_audit.json", fidelity)
    _write_fidelity_markdown(fidelity)
    summary = {
        "revision": REVISION, "status": "COMPLETE",
        "checkpoint_rows": len(metrics), "trajectory_rows": len(summaries),
        "overlap": overlap, "fidelity_pass": fidelity["overall_pass"],
    }
    write_json_atomic(RESULTS_DIR / "analysis_summary.json", summary)
    return summary


def main() -> None:
    print(json.dumps(aggregate(), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
