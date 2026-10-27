"""Mechanical calibration aggregation and immutable selection barrier."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np

from .protocol import (
    DB_CANDIDATE_UPDATES,
    FH_CONFIGS,
    FH_REPLICATES,
    FROZEN_SELECTIONS,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    TA_INFERENCE_REPLICATES,
    TA_PATHS,
    TARGETS,
    SOURCE_ROOT,
    array_sha256,
    canonical_sha256,
    file_sha256,
    load_bank,
    write_json_atomic,
)


def _read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _pooled_prediction(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    predictions, references = [], []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            predictions.append(np.asarray(data["prediction"], dtype=np.float64).reshape(-1))
            reference = np.asarray(data["reference"], dtype=np.float64)
            if reference.ndim == 0:
                reference = np.repeat(reference.reshape(1), len(predictions[-1]))
            references.append(reference.reshape(-1))
    return np.concatenate(predictions), np.concatenate(references)


def _summary(
    method: str,
    candidate: str,
    records: list[dict[str, Any]],
    predictions: np.ndarray,
    references: np.ndarray,
    *,
    active_key: str,
    **extra: Any,
) -> dict[str, Any]:
    error = predictions - references
    active = np.asarray([float(row[active_key]) for row in records], dtype=np.float64)
    point_count = np.asarray([float(row.get("points", 1)) for row in records])
    per_query = active / point_count
    nonfinite = int(sum(int(row.get("nonfinite_count", 0)) for row in records))
    nonfinite += int((~np.isfinite(predictions)).sum())
    row = {
        "method": method,
        "candidate": candidate,
        "run_count": len(records),
        "prediction_count": len(predictions),
        "pooled_rmse": float(np.sqrt(np.mean(error**2))),
        "pooled_mae": float(np.mean(np.abs(error))),
        "pooled_bias": float(np.mean(error)),
        "mean_active_seconds": float(np.mean(active)),
        "std_active_seconds": float(np.std(active, ddof=1)) if len(active) > 1 else 0.0,
        "mean_seconds_per_query": float(np.mean(per_query)),
        "nonfinite_count": nonfinite,
        "gpu": ";".join(sorted({str(value.get("gpu")) for value in records})),
        **extra,
    }
    for name, threshold in TARGETS.items():
        row[f"eligible_{name}"] = bool(nonfinite == 0 and row["pooled_rmse"] <= threshold)
    return row


def _ta_rows(bank_hash: str) -> list[dict[str, Any]]:
    rows = []
    for candidate in ("direct", *(f"mcv_{paths}" for paths in TA_PATHS)):
        expected_reps = 1 if candidate == "direct" else TA_INFERENCE_REPLICATES
        records = []
        for seed in (0, 1, 2):
            for replicate in range(expected_reps):
                path = RAW_DIR / "calibration" / "tangent_actor" / (
                    f"ta__calibration__s{seed}__{candidate}__r{replicate}.json"
                )
                if not path.is_file():
                    raise RuntimeError(f"missing Tangent-Actor calibration record: {path}")
                record = _read(path)
                if record["query_input_sha256"] != bank_hash:
                    raise RuntimeError(f"Tangent-Actor bank hash mismatch: {path}")
                records.append(record)
        predictions, references = _pooled_prediction(
            [Path(row["prediction_file"]) for row in records]
        )
        rows.append(_summary(
            "tangent_actor",
            candidate,
            records,
            predictions,
            references,
            active_key="inference_seconds",
            paths=0 if candidate == "direct" else int(candidate.split("_")[1]),
            time_steps=0 if candidate == "direct" else 64,
        ))
    return rows


def _fh_rows(bank_hash: str) -> list[dict[str, Any]]:
    rows = []
    for n, M in FH_CONFIGS:
        records = []
        for replicate in range(FH_REPLICATES):
            path = RAW_DIR / "calibration" / "fh_mlp" / (
                f"fh__calibration__n{n}_M{M}__r{replicate}.json"
            )
            if not path.is_file():
                raise RuntimeError(f"missing FH calibration record: {path}")
            record = _read(path)
            if record["query_input_sha256"] != bank_hash:
                raise RuntimeError(f"FH bank hash mismatch: {path}")
            records.append(record)
        predictions, references = _pooled_prediction(
            [Path(row["prediction_file"]) for row in records]
        )
        rows.append(_summary(
            "fh_mlp_standard",
            f"n{n}_M{M}",
            records,
            predictions,
            references,
            active_key="inference_seconds",
            n=n,
            M=M,
        ))
    return rows


def _db_rows(bank_hash: str) -> list[dict[str, Any]]:
    records = []
    for point in range(32):
        for seed in (0, 1, 2):
            path = RAW_DIR / "calibration" / "deep_bsde" / (
                f"db__calibration__q{point:03d}__s{seed}__u8000__c2000-4000-8000.json"
            )
            if not path.is_file():
                raise RuntimeError(f"missing Deep BSDE calibration record: {path}")
            record = _read(path)
            expected_point_hash = array_sha256(
                load_bank("calibration")["t"][point:point + 1],
                load_bank("calibration")["x"][point:point + 1],
            )
            if record["query_input_sha256"] != expected_point_hash:
                raise RuntimeError(f"Deep BSDE point hash mismatch: {path}")
            records.append(record)
    rows = []
    for updates in DB_CANDIDATE_UPDATES:
        predictions, references, active = [], [], []
        nonfinite = 0
        for record in records:
            candidates = {int(row["updates"]): row for row in record["candidates"]}
            if updates not in candidates:
                nonfinite += 1
                continue
            item = candidates[updates]
            predictions.append(float(item["y0"]))
            references.append(float(item["reference"]))
            active.append(float(item["active_seconds"]))
            nonfinite += int(item.get("nonfinite_count_cumulative", 0))
        prediction = np.asarray(predictions, dtype=np.float64)
        reference = np.asarray(references, dtype=np.float64)
        error = prediction - reference
        nonfinite += int((~np.isfinite(prediction)).sum())
        row = {
            "method": "deep_bsde",
            "candidate": f"updates_{updates}",
            "updates": updates,
            "run_count": len(records),
            "prediction_count": len(prediction),
            "pooled_rmse": float(np.sqrt(np.mean(error**2))) if len(error) else float("inf"),
            "pooled_mae": float(np.mean(np.abs(error))) if len(error) else float("inf"),
            "pooled_bias": float(np.mean(error)) if len(error) else float("nan"),
            "mean_active_seconds": float(np.mean(active)) if active else float("inf"),
            "std_active_seconds": float(np.std(active, ddof=1)) if len(active) > 1 else 0.0,
            "mean_seconds_per_query": float(np.mean(active)) if active else float("inf"),
            "nonfinite_count": nonfinite,
            "gpu": ";".join(sorted({str(value.get("gpu")) for value in records})),
        }
        for name, threshold in TARGETS.items():
            row[f"eligible_{name}"] = bool(
                len(prediction) == 96 and nonfinite == 0 and row["pooled_rmse"] <= threshold
            )
        rows.append(row)
    return rows


def _select(rows: list[dict[str, Any]], method: str, target: str) -> dict[str, Any]:
    local = [row for row in rows if row["method"] == method]
    eligible = [row for row in local if row[f"eligible_{target}"]]
    if not eligible:
        best = min(local, key=lambda row: row["pooled_rmse"])
        return {
            "status": "TARGET_NOT_REACHED",
            "selected_candidate": None,
            "best_candidate": best["candidate"],
            "best_calibration_rmse": best["pooled_rmse"],
            "best_actual_mean_seconds": best["mean_active_seconds"],
        }
    if method == "deep_bsde":
        selected = min(eligible, key=lambda row: row["updates"])
    elif method == "tangent_actor":
        minimum = min(row["mean_seconds_per_query"] for row in eligible)
        near = [
            row for row in eligible if row["mean_seconds_per_query"] <= 1.02 * minimum
        ]
        order = {"direct": 0, **{f"mcv_{p}": i + 1 for i, p in enumerate(TA_PATHS)}}
        selected = min(near, key=lambda row: order[row["candidate"]])
    else:
        order = {f"n{n}_M{M}": index for index, (n, M) in enumerate(FH_CONFIGS)}
        selected = min(
            eligible,
            key=lambda row: (row["mean_seconds_per_query"], order[row["candidate"]]),
        )
    return {
        "status": "SELECTED",
        "selected_candidate": selected["candidate"],
        "calibration_rmse": selected["pooled_rmse"],
        "calibration_mean_seconds_per_query": selected["mean_seconds_per_query"],
        "candidate": selected,
    }


def freeze() -> dict[str, Any]:
    if FROZEN_SELECTIONS.exists():
        raise RuntimeError("FROZEN_SELECTIONS already exists; refusing to overwrite")
    bank = load_bank("calibration", require_truth=True)
    bank_hash = array_sha256(bank["t"], bank["x"])
    rows = _ta_rows(bank_hash) + _fh_rows(bank_hash) + _db_rows(bank_hash)
    _write_csv(RESULTS_DIR / "calibration_rows.csv", rows)
    selections = {}
    for target in TARGETS:
        selections[target] = {
            "threshold": TARGETS[target],
            "tangent_actor": _select(rows, "tangent_actor", target),
            "fh_mlp_standard": _select(rows, "fh_mlp_standard", target),
            "deep_bsde": _select(rows, "deep_bsde", target),
        }
    code_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
    ).strip()
    payload = {
        "revision": REVISION,
        "status": "FROZEN",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "code_commit": code_commit,
        "calibration_bank_input_sha256": bank_hash,
        "calibration_bank_file_sha256": file_sha256(
            RESULTS_DIR / "banks" / "calibration.npz"
        ),
        "calibration_rows_file_sha256": file_sha256(
            RESULTS_DIR / "calibration_rows.csv"
        ),
        "targets": selections,
        "candidate_summaries": rows,
        "truth_policy": "read only after each frozen calibration run completed",
    }
    payload["canonical_sha256"] = canonical_sha256(payload)
    write_json_atomic(FROZEN_SELECTIONS, payload)
    return payload


if __name__ == "__main__":
    result = freeze()
    print(json.dumps({
        "status": result["status"],
        "canonical_sha256": result["canonical_sha256"],
        "targets": result["targets"],
    }, sort_keys=True))

