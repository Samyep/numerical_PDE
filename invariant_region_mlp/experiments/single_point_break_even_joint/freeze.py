"""Aggregate calibration and freeze target-specific configurations."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np

from .protocol import (
    DB_CANDIDATE_UPDATES, FH_CONFIGS, FH_REPLICATES, FROZEN_SELECTIONS,
    RESULTS_DIR, REVISION, SOURCE_ROOT, TA_INFERENCE_REPLICATES, TA_PATHS,
    TARGETS, array_sha256, canonical_sha256, file_sha256, load_bank,
    write_json_atomic,
)


RAW_DIR = RESULTS_DIR / "raw"


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    temporary.replace(path)


def _pooled(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    prediction, reference = [], []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            prediction.append(np.asarray(data["prediction"], dtype=np.float64).reshape(-1))
            local = np.asarray(data["reference"], dtype=np.float64)
            if local.ndim == 0: local = np.repeat(local.reshape(1), len(prediction[-1]))
            reference.append(local.reshape(-1))
    return np.concatenate(prediction), np.concatenate(reference)


def _summary(method: str, candidate: str, records: list[dict[str, Any]], prediction: np.ndarray, reference: np.ndarray, *, active_key: str, **extra: Any) -> dict[str, Any]:
    error = prediction - reference
    finite = error[np.isfinite(error)]
    active = np.asarray([float(row[active_key]) for row in records])
    points = np.asarray([float(row.get("points", 1)) for row in records])
    nonfinite = sum(int(row.get("nonfinite_count", 0)) for row in records) + int((~np.isfinite(prediction)).sum())
    result = {
        "method": method, "candidate": candidate, "run_count": len(records),
        "prediction_count": len(prediction),
        "pooled_rmse": float(np.sqrt(np.mean(finite**2))) if len(finite) else 1e300,
        "pooled_mae": float(np.mean(np.abs(finite))) if len(finite) else 1e300,
        "pooled_bias": float(np.mean(finite)) if len(finite) else 0.0,
        "mean_active_seconds": float(active.mean()),
        "std_active_seconds": float(active.std(ddof=1)) if len(active) > 1 else 0.0,
        "mean_seconds_per_query": float(np.mean(active / points)),
        "nonfinite_count": int(nonfinite),
        "gpu": ";".join(sorted({str(row.get("gpu")) for row in records})),
        **extra,
    }
    for target, threshold in TARGETS.items(): result[f"eligible_{target}"] = bool(nonfinite == 0 and result["pooled_rmse"] <= threshold)
    return result


def _ta_rows(bank_hash: str) -> list[dict[str, Any]]:
    rows = []
    for candidate in ("direct", *(f"mcv_{paths}" for paths in TA_PATHS)):
        records = []
        for seed in (0, 1, 2):
            for replicate in ((0,) if candidate == "direct" else range(TA_INFERENCE_REPLICATES)):
                path = RAW_DIR / "calibration" / "tangent_actor" / f"ta__calibration__s{seed}__{candidate}__r{replicate}.json"
                record = _read(path)
                if record["query_input_sha256"] != bank_hash: raise RuntimeError(f"TA bank hash mismatch: {path}")
                records.append(record)
        prediction, reference = _pooled([Path(row["prediction_file"]) for row in records])
        rows.append(_summary("tangent_actor", candidate, records, prediction, reference, active_key="inference_seconds", paths=0 if candidate == "direct" else int(candidate.split("_")[1]), time_steps=0 if candidate == "direct" else 64))
    return rows


def _fh_rows(bank_hash: str) -> list[dict[str, Any]]:
    rows = []
    for n, M in FH_CONFIGS:
        records = []
        for replicate in range(FH_REPLICATES):
            path = RAW_DIR / "calibration" / "fh_mlp" / f"fh__calibration__n{n}_M{M}__r{replicate}.json"
            record = _read(path)
            if record["query_input_sha256"] != bank_hash: raise RuntimeError(f"FH bank hash mismatch: {path}")
            records.append(record)
        prediction, reference = _pooled([Path(row["prediction_file"]) for row in records])
        generator_calls = [
            (float(row["stats"]["terminal_samples"]) + float(row["stats"]["transition_samples"])) / float(row["points"])
            for row in records
        ]
        rows.append(_summary("fh_mlp_standard", f"n{n}_M{M}", records, prediction, reference, active_key="inference_seconds", n=n, M=M, mean_generator_calls_per_query=float(np.mean(generator_calls))))
    return rows


def _db_rows(bank_hash: str) -> list[dict[str, Any]]:
    records = []
    bank = load_bank("calibration")
    for point in range(32):
        for seed in (0, 1, 2):
            path = RAW_DIR / "calibration" / "deep_bsde" / f"db__calibration__q{point:03d}__s{seed}__u4000__c500-1000-2000-4000.json"
            record = _read(path)
            expected = array_sha256(bank["t"][point:point + 1], bank["x"][point:point + 1])
            if record["query_input_sha256"] != expected: raise RuntimeError(f"Deep bank hash mismatch: {path}")
            records.append(record)
    rows = []
    for updates in DB_CANDIDATE_UPDATES:
        prediction, reference, active = [], [], []
        nonfinite = 0
        for record in records:
            item = {int(row["updates"]): row for row in record["candidates"]}.get(updates)
            if item is None: nonfinite += 1; continue
            prediction.append(float(item["y0"])); reference.append(float(item["reference"])); active.append(float(item["active_seconds"])); nonfinite += int(item.get("nonfinite_count_cumulative", 0))
        prediction_array, reference_array = np.asarray(prediction), np.asarray(reference)
        error = prediction_array - reference_array; finite = error[np.isfinite(error)]
        nonfinite += int((~np.isfinite(prediction_array)).sum())
        row = {
            "method": "deep_bsde", "candidate": f"updates_{updates}", "updates": updates,
            "run_count": len(records), "prediction_count": len(prediction_array),
            "pooled_rmse": float(np.sqrt(np.mean(finite**2))) if len(finite) else 1e300,
            "pooled_mae": float(np.mean(np.abs(finite))) if len(finite) else 1e300,
            "pooled_bias": float(np.mean(finite)) if len(finite) else 0.0,
            "mean_active_seconds": float(np.mean(active)) if active else 1e300,
            "std_active_seconds": float(np.std(active, ddof=1)) if len(active) > 1 else 0.0,
            "mean_seconds_per_query": float(np.mean(active)) if active else 1e300,
            "nonfinite_count": int(nonfinite),
            "gpu": ";".join(sorted({str(record.get("gpu")) for record in records})),
        }
        for target, threshold in TARGETS.items(): row[f"eligible_{target}"] = bool(len(prediction_array) == 96 and nonfinite == 0 and row["pooled_rmse"] <= threshold)
        rows.append(row)
    return rows


def _select(rows: list[dict[str, Any]], method: str, target: str) -> dict[str, Any]:
    local = [row for row in rows if row["method"] == method]
    eligible = [row for row in local if row[f"eligible_{target}"]]
    if not eligible:
        best = min(local, key=lambda row: row["pooled_rmse"])
        return {"status": "TARGET_NOT_REACHED", "selected_candidate": None, "best_candidate": best["candidate"], "best_calibration_rmse": best["pooled_rmse"], "best_actual_mean_seconds": best["mean_active_seconds"]}
    if method == "deep_bsde": selected = min(eligible, key=lambda row: row["updates"])
    else:
        fastest = min(row["mean_seconds_per_query"] for row in eligible)
        near = [row for row in eligible if row["mean_seconds_per_query"] <= 1.02 * fastest]
        if method == "tangent_actor":
            order = {"direct": 0, **{f"mcv_{paths}": index + 1 for index, paths in enumerate(TA_PATHS)}}
            selected = min(near, key=lambda row: order[row["candidate"]])
        else:
            selected = min(near, key=lambda row: (row["mean_generator_calls_per_query"], row["n"], row["M"]))
    return {"status": "SELECTED", "selected_candidate": selected["candidate"], "calibration_rmse": selected["pooled_rmse"], "calibration_mean_seconds_per_query": selected["mean_seconds_per_query"], "candidate": selected}


def freeze() -> dict[str, Any]:
    if FROZEN_SELECTIONS.exists(): raise RuntimeError("FROZEN_SELECTIONS already exists; refusing to overwrite")
    bank = load_bank("calibration", require_truth=True); bank_hash = array_sha256(bank["t"], bank["x"])
    rows = _ta_rows(bank_hash) + _fh_rows(bank_hash) + _db_rows(bank_hash)
    _write_csv(RESULTS_DIR / "calibration_rows.csv", rows)
    selections = {target: {"threshold": TARGETS[target], **{method: _select(rows, method, target) for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde")}} for target in TARGETS}
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True).strip()
    payload = {"revision": REVISION, "status": "FROZEN", "frozen_utc": datetime.now(timezone.utc).isoformat(), "code_commit": code_commit, "calibration_bank_input_sha256": bank_hash, "calibration_bank_file_sha256": file_sha256(RESULTS_DIR / "banks" / "calibration.npz"), "calibration_rows_file_sha256": file_sha256(RESULTS_DIR / "calibration_rows.csv"), "targets": selections, "candidate_summaries": rows, "truth_policy": "read only after every predeclared calibration run completed"}
    payload["canonical_sha256"] = canonical_sha256(payload); write_json_atomic(FROZEN_SELECTIONS, payload); return payload


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2, sort_keys=True))
