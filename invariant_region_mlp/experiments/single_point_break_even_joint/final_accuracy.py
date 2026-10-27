"""One-shot final-bank qualification of frozen joint-PDE configurations."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from .protocol import FINAL_SEEDS, RAW_DIR, RESULTS_DIR, REVISION, TARGETS, array_sha256, load_bank, require_frozen_selections, write_json_atomic


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


def _pooled(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    prediction, reference = [], []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            prediction.append(np.asarray(data["prediction"], dtype=np.float64).reshape(-1))
            local = np.asarray(data["reference"], dtype=np.float64)
            if local.ndim == 0: local = np.repeat(local.reshape(1), len(prediction[-1]))
            reference.append(local.reshape(-1))
    return np.concatenate(prediction), np.concatenate(reference)


def _metrics(method: str, candidate: str, records: list[dict[str, Any]], prediction: np.ndarray, reference: np.ndarray, *, extra_nonfinite: int = 0) -> dict[str, Any]:
    error = prediction - reference; finite = error[np.isfinite(error)]
    nonfinite = extra_nonfinite + int((~np.isfinite(prediction)).sum()) + sum(int(row.get("nonfinite_count", 0)) for row in records)
    return {"method": method, "candidate": candidate, "run_count": len(records), "prediction_count": len(prediction), "pooled_rmse": float(np.sqrt(np.mean(finite**2))) if len(finite) else 1e300, "pooled_mae": float(np.mean(np.abs(finite))) if len(finite) else 1e300, "pooled_bias": float(np.mean(finite)) if len(finite) else 0.0, "nonfinite_count": int(nonfinite), "gpu": ";".join(sorted({str(row.get("gpu")) for row in records}))}


def _ta(candidate: str, bank_hash: str) -> dict[str, Any]:
    records, paths = [], []
    for seed in FINAL_SEEDS:
        for replicate in ((0,) if candidate == "direct" else (0, 1)):
            path = RAW_DIR / "final" / "tangent_actor" / f"ta__final__s{seed}__{candidate}__r{replicate}.json"
            record = _read(path)
            if record["query_input_sha256"] != bank_hash: raise RuntimeError(f"TA final bank hash mismatch: {path}")
            records.append(record); paths.append(Path(record["prediction_file"]))
    prediction, reference = _pooled(paths); return _metrics("tangent_actor", candidate, records, prediction, reference)


def _fh(candidate: str, bank_hash: str) -> dict[str, Any]:
    n_text, m_text = candidate.split("_"); n, M = int(n_text[1:]), int(m_text[1:])
    records, paths = [], []
    for replicate in range(5):
        path = RAW_DIR / "final" / "fh_mlp" / f"fh__final__n{n}_M{M}__r{replicate}.json"
        record = _read(path)
        if record["query_input_sha256"] != bank_hash: raise RuntimeError(f"FH final bank hash mismatch: {path}")
        records.append(record); paths.append(Path(record["prediction_file"]))
    prediction, reference = _pooled(paths); return _metrics("fh_mlp_standard", candidate, records, prediction, reference)


def _db(candidate: str, bank: dict[str, np.ndarray]) -> dict[str, Any]:
    updates = int(candidate.split("_")[1]); prediction, reference, records = [], [], []; nonfinite = 0
    for point in range(len(bank["x"])):
        for seed in FINAL_SEEDS:
            matches = list((RAW_DIR / "final" / "deep_bsde").glob(f"db__final__q{point:03d}__s{seed}__u*__c*.json"))
            if len(matches) != 1: raise RuntimeError(f"expected one Deep BSDE final record, got {matches}")
            record = _read(matches[0]); item = {int(row["updates"]): row for row in record["candidates"]}.get(updates)
            if item is None: nonfinite += 1; continue
            prediction.append(float(item["y0"])); reference.append(float(item["reference"])); nonfinite += int(item.get("nonfinite_count_cumulative", 0)); records.append(record)
    return _metrics("deep_bsde", candidate, records, np.asarray(prediction), np.asarray(reference), extra_nonfinite=nonfinite)


def aggregate() -> dict[str, Any]:
    frozen = require_frozen_selections(); bank = load_bank("final", require_truth=True); bank_hash = array_sha256(bank["t"], bank["x"])
    cache: dict[tuple[str, str], dict[str, Any]] = {}; rows: list[dict[str, Any]] = []; targets: dict[str, Any] = {}
    for target, threshold in TARGETS.items():
        targets[target] = {"threshold": threshold}
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            selection = frozen["targets"][target][method]; candidate = selection.get("selected_candidate")
            if selection["status"] != "SELECTED" or candidate is None:
                targets[target][method] = {"status": "TARGET_NOT_REACHED", "selected_candidate": None, "calibration": selection}; continue
            key = (method, candidate)
            if key not in cache:
                cache[key] = _ta(candidate, bank_hash) if method == "tangent_actor" else (_fh(candidate, bank_hash) if method == "fh_mlp_standard" else _db(candidate, bank))
                rows.append(cache[key])
            metrics = cache[key]; passed = metrics["nonfinite_count"] == 0 and metrics["pooled_rmse"] <= threshold
            targets[target][method] = {"status": "FINAL_ACCURACY_PASSED" if passed else "FINAL_ACCURACY_FAILED", "selected_candidate": candidate, "threshold": threshold, "final_rmse": metrics["pooled_rmse"], "final_prediction_count": metrics["prediction_count"], "nonfinite_count": metrics["nonfinite_count"], "calibration_rmse": selection["calibration_rmse"]}
    _write_csv(RESULTS_DIR / "final_accuracy_rows.csv", rows)
    payload = {"revision": REVISION, "status": "COMPLETE", "aggregated_utc": datetime.now(timezone.utc).isoformat(), "frozen_selections_canonical_sha256": frozen["canonical_sha256"], "final_bank_input_sha256": bank_hash, "targets": targets, "rows": rows, "configuration_changes_after_final_errors": False}
    write_json_atomic(RESULTS_DIR / "final_accuracy_summary.json", payload); return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))
