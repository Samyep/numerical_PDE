"""Integrity audit for paired VB repetition artifacts."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
RESULT_ROOT = PROJECT_ROOT / "results" / "active_vb_high_budget"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def audit(stages: list[str]) -> dict[str, Any]:
    grouped: dict[tuple[int, int, int, int], list[Path]] = defaultdict(list)
    stage_counts = {}
    manifest_status = {}
    for stage in stages:
        paths = sorted((RESULT_ROOT / "repetitions" / stage).glob("*.npz"))
        stage_counts[stage] = len(paths)
        manifest = RESULT_ROOT / f"{stage}_manifest.json"
        if manifest.exists():
            manifest_status[stage] = json.loads(manifest.read_text(encoding="utf-8"))["status"]
        for path in paths:
            with np.load(path, allow_pickle=False) as data:
                metadata = json.loads(str(data["metadata_json"]))
            key = (
                int(metadata["dimension"]),
                int(metadata["n"]),
                int(metadata["M"]),
                int(metadata["repetition"]),
            )
            grouped[key].append(path)

    max_root_terminal_difference = 0.0
    max_truth_difference = 0.0
    max_z_zero_f_zero_difference = 0.0
    max_f_zero_nonlinear_correction = 0.0
    nonfinite_prediction_values = 0
    nonfinite_truth_values = 0
    missing_pairs = []
    fzero_nonzero_f_calls = []
    recursive_work_mismatches = []
    checked_files = 0
    for key, paths in sorted(grouped.items()):
        reference_root = None
        reference_truth = None
        method_payload: dict[str, tuple[np.ndarray, np.ndarray, dict[str, Any]]] = {}
        recursive_f_counts = set()
        for path in paths:
            with np.load(path, allow_pickle=False) as data:
                prediction = np.array(data["prediction"], copy=True)
                truth = np.array(data["truth"], copy=True)
                root = np.array(data["root_terminal"], copy=True)
                nonlinear = np.array(data["nonlinear_u_correction"], copy=True)
                metadata = json.loads(str(data["metadata_json"]))
            checked_files += 1
            method = metadata["method"]["name"]
            method_payload[method] = (prediction, nonlinear, metadata)
            nonfinite_prediction_values += int(np.size(prediction) - np.count_nonzero(np.isfinite(prediction)))
            nonfinite_truth_values += int(np.size(truth) - np.count_nonzero(np.isfinite(truth)))
            if reference_root is None:
                reference_root = root
                reference_truth = truth
            else:
                max_root_terminal_difference = max(
                    max_root_terminal_difference, float(np.max(np.abs(root - reference_root)))
                )
                max_truth_difference = max(
                    max_truth_difference, float(np.max(np.abs(truth - reference_truth)))
                )
            if method == "f_zero":
                if metadata["work"]["f_evals"] != 0:
                    fzero_nonzero_f_calls.append({"key": key, "value": metadata["work"]["f_evals"]})
                max_f_zero_nonlinear_correction = max(
                    max_f_zero_nonlinear_correction, float(np.max(np.abs(nonlinear)))
                )
            else:
                recursive_f_counts.add(int(metadata["work"]["f_evals"]))
        if "z_zero" not in method_payload or "f_zero" not in method_payload:
            missing_pairs.append({"key": key, "methods": sorted(method_payload)})
        else:
            difference = method_payload["z_zero"][0] - method_payload["f_zero"][0]
            max_z_zero_f_zero_difference = max(
                max_z_zero_f_zero_difference, float(np.max(np.abs(difference)))
            )
        if len(recursive_f_counts) > 1:
            recursive_work_mismatches.append({"key": key, "f_counts": sorted(recursive_f_counts)})

    checks = {
        "all_manifests_complete": all(value == "complete" for value in manifest_status.values()),
        "root_terminal_bitwise_paired": max_root_terminal_difference == 0.0,
        "truth_bitwise_identical_within_pair": max_truth_difference == 0.0,
        "z_zero_equals_f_zero_output": max_z_zero_f_zero_difference == 0.0,
        "f_zero_has_zero_nonlinear_correction": max_f_zero_nonlinear_correction == 0.0,
        "f_zero_has_zero_f_calls": not fzero_nonzero_f_calls,
        "recursive_methods_have_equal_f_work": not recursive_work_mismatches,
        "all_predictions_finite": nonfinite_prediction_values == 0,
        "all_truth_finite": nonfinite_truth_values == 0,
        "all_groups_have_z_zero_f_zero_pair": not missing_pairs,
    }
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "stages": stages,
        "stage_file_counts": stage_counts,
        "manifest_status": manifest_status,
        "paired_groups": len(grouped),
        "checked_files": checked_files,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "max_root_terminal_abs_difference": max_root_terminal_difference,
        "max_truth_abs_difference": max_truth_difference,
        "max_z_zero_f_zero_abs_difference": max_z_zero_f_zero_difference,
        "max_f_zero_nonlinear_correction": max_f_zero_nonlinear_correction,
        "nonfinite_prediction_values": nonfinite_prediction_values,
        "nonfinite_truth_values": nonfinite_truth_values,
        "missing_pairs": missing_pairs,
        "fzero_nonzero_f_calls": fzero_nonzero_f_calls,
        "recursive_work_mismatches": recursive_work_mismatches,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stages", default="main,high_budget_extension")
    args = parser.parse_args()
    payload = audit([stage for stage in args.stages.split(",") if stage])
    output = RESULT_ROOT / "validation_audit.json"
    atomic_json(output, payload)
    print(json.dumps(payload["checks"], indent=2))
    if not payload["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
