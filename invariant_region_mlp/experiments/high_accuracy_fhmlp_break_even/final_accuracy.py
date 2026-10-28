"""One-shot independent final-bank qualification of frozen choices."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from .protocol import (
    FH_ACCURACY_REPLICATES,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    TA_MCV_ACCURACY_REPLICATES,
    TA_TRAINING_SEEDS,
    TARGETS,
    array_sha256,
    load_bank,
    require_frozen_selections,
    write_json_atomic,
)


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


def _pooled(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    predictions, references = [], []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            predictions.append(np.asarray(data["prediction"], dtype=np.float64).reshape(-1))
            references.append(np.asarray(data["reference"], dtype=np.float64).reshape(-1))
    return np.concatenate(predictions), np.concatenate(references)


def _metrics(
    method: str,
    candidate: str,
    records: list[dict[str, Any]],
    prediction: np.ndarray,
    reference: np.ndarray,
) -> dict[str, Any]:
    error = prediction - reference
    nonfinite = int((~np.isfinite(prediction)).sum()) + sum(
        int(row.get("nonfinite_count", 0)) for row in records
    )
    return {
        "method": method,
        "candidate": candidate,
        "run_count": len(records),
        "prediction_count": len(prediction),
        "pooled_rmse": float(np.sqrt(np.mean(error**2))),
        "pooled_mae": float(np.mean(np.abs(error))),
        "pooled_bias": float(np.mean(error)),
        "nonfinite_count": nonfinite,
        "gpu": ";".join(sorted({str(row.get("gpu")) for row in records})),
    }


def _ta(method: str, candidate: str, bank_hash: str) -> dict[str, Any]:
    records, paths = [], []
    for seed in TA_TRAINING_SEEDS:
        replicates = (0,) if candidate == "direct" else range(TA_MCV_ACCURACY_REPLICATES)
        for replicate in replicates:
            path = RAW_DIR / "final" / "tangent_actor" / f"ta__final__s{seed}__{candidate}__r{replicate}.json"
            record = _read(path)
            if record["query_input_sha256"] != bank_hash:
                raise RuntimeError(f"TA final bank hash mismatch: {path}")
            records.append(record)
            paths.append(Path(record["prediction_file"]))
    prediction, reference = _pooled(paths)
    return _metrics(method, candidate, records, prediction, reference)


def _fh(candidate: str, bank_hash: str) -> dict[str, Any]:
    n_text, m_text = candidate.split("_")
    n, M = int(n_text[1:]), int(m_text[1:])
    records, paths = [], []
    for replicate in range(FH_ACCURACY_REPLICATES):
        path = RAW_DIR / "final" / "fh_mlp" / f"fh__final__n{n}_M{M}__r{replicate}.json"
        record = _read(path)
        if record["query_input_sha256"] != bank_hash:
            raise RuntimeError(f"FH final bank hash mismatch: {path}")
        records.append(record)
        paths.append(Path(record["prediction_file"]))
    prediction, reference = _pooled(paths)
    return _metrics("fh_mlp_standard", candidate, records, prediction, reference)


def aggregate() -> dict[str, Any]:
    frozen = require_frozen_selections()
    bank = load_bank("final", require_truth=True)
    bank_hash = array_sha256(bank["t"], bank["x"])
    cache: dict[tuple[str, str], dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    targets: dict[str, Any] = {}
    for target, threshold in TARGETS.items():
        targets[target] = {"threshold": threshold}
        for method in ("ta_direct", "ta_fk_mcv", "fh_mlp_standard"):
            selection = frozen["targets"][target][method]
            candidate = selection.get("selected_candidate")
            if selection["status"] != "SELECTED" or candidate is None:
                targets[target][method] = {
                    "status": "TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET",
                    "selected_candidate": None,
                    "calibration": selection,
                }
                continue
            key = (method, candidate)
            if key not in cache:
                cache[key] = (
                    _ta(method, candidate, bank_hash)
                    if method.startswith("ta_")
                    else _fh(candidate, bank_hash)
                )
                rows.append(cache[key])
            metrics = cache[key]
            passed = (
                metrics["nonfinite_count"] == 0
                and metrics["pooled_rmse"] <= threshold
            )
            targets[target][method] = {
                "status": "FINAL_ACCURACY_PASSED" if passed else "FINAL_ACCURACY_FAILED",
                "selected_candidate": candidate,
                "threshold": threshold,
                "final_rmse": metrics["pooled_rmse"],
                "final_prediction_count": metrics["prediction_count"],
                "nonfinite_count": metrics["nonfinite_count"],
                "calibration_rmse": selection["calibration_rmse"],
            }
    _write_csv(RESULTS_DIR / "final_accuracy_rows.csv", rows)
    payload = {
        "revision": REVISION,
        "status": "COMPLETE",
        "aggregated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_selections_canonical_sha256": frozen["canonical_sha256"],
        "final_bank_input_sha256": bank_hash,
        "targets": targets,
        "rows": rows,
        "configuration_changes_after_final_errors": False,
    }
    write_json_atomic(RESULTS_DIR / "final_accuracy_summary.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(aggregate(), indent=2, sort_keys=True))
