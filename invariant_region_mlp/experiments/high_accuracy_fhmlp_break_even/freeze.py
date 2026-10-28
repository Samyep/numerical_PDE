"""Aggregate calibration and irreversibly freeze target-specific choices."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np

from .protocol import (
    COST_GATE,
    FH_ACCURACY_REPLICATES,
    FROZEN_SELECTIONS,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    TA_MCV_ACCURACY_REPLICATES,
    TA_PATHS,
    TA_TRAINING_SEEDS,
    TARGETS,
    array_sha256,
    canonical_sha256,
    file_sha256,
    load_bank,
    require_cost_gate,
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
    prediction, reference = [], []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            prediction.append(np.asarray(data["prediction"], dtype=np.float64).reshape(-1))
            reference.append(np.asarray(data["reference"], dtype=np.float64).reshape(-1))
    return np.concatenate(prediction), np.concatenate(reference)


def _summary(
    method: str,
    candidate: str,
    records: list[dict[str, Any]],
    prediction: np.ndarray,
    reference: np.ndarray,
    **extra: Any,
) -> dict[str, Any]:
    error = prediction - reference
    nonfinite = int((~np.isfinite(prediction)).sum()) + sum(
        int(row.get("nonfinite_count", 0)) for row in records
    )
    active = np.asarray([float(row["inference_seconds"]) for row in records])
    points = np.asarray([float(row["points"]) for row in records])
    result = {
        "method": method,
        "candidate": candidate,
        "status": "COMPLETE" if nonfinite == 0 else "NUMERICAL_FAILURE",
        "run_count": len(records),
        "prediction_count": len(prediction),
        "pooled_rmse": float(np.sqrt(np.mean(error**2))),
        "pooled_mae": float(np.mean(np.abs(error))),
        "pooled_bias": float(np.mean(error)),
        "mean_inference_seconds": float(active.mean()),
        "std_inference_seconds": float(active.std(ddof=1)) if len(active) > 1 else 0.0,
        "mean_seconds_per_query": float(np.mean(active / points)),
        "nonfinite_count": nonfinite,
        "gpu": ";".join(sorted({str(row.get("gpu")) for row in records})),
        **extra,
    }
    for target, threshold in TARGETS.items():
        result[f"eligible_{target}"] = bool(
            nonfinite == 0 and result["pooled_rmse"] <= threshold
        )
    return result


def _required_paths(gate: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    for seed in TA_TRAINING_SEEDS:
        paths.append(RAW_DIR / "calibration" / "tangent_actor" / f"ta__calibration__s{seed}__direct__r0.json")
        for budget in TA_PATHS:
            for replicate in range(TA_MCV_ACCURACY_REPLICATES):
                paths.append(RAW_DIR / "calibration" / "tangent_actor" / f"ta__calibration__s{seed}__mcv_{budget}__r{replicate}.json")
    for row in gate["eligible_fh_configs"].values():
        for replicate in range(FH_ACCURACY_REPLICATES):
            paths.append(RAW_DIR / "calibration" / "fh_mlp" / f"fh__calibration__n{row['n']}_M{row['M']}__r{replicate}.json")
    return paths


def _ta_rows(bank_hash: str) -> list[dict[str, Any]]:
    rows = []
    for candidate in ("direct", *(f"mcv_{paths}" for paths in TA_PATHS)):
        records = []
        for seed in TA_TRAINING_SEEDS:
            replicates = (0,) if candidate == "direct" else range(TA_MCV_ACCURACY_REPLICATES)
            for replicate in replicates:
                path = RAW_DIR / "calibration" / "tangent_actor" / f"ta__calibration__s{seed}__{candidate}__r{replicate}.json"
                record = _read(path)
                if record["query_input_sha256"] != bank_hash:
                    raise RuntimeError(f"TA bank hash mismatch: {path}")
                records.append(record)
        prediction, reference = _pooled([Path(row["prediction_file"]) for row in records])
        rows.append(_summary(
            "ta_direct" if candidate == "direct" else "ta_fk_mcv",
            candidate,
            records,
            prediction,
            reference,
            paths=0 if candidate == "direct" else int(candidate.split("_")[1]),
            time_steps=0 if candidate == "direct" else 64,
            charged_training_seconds_mean=float(np.mean([
                row["charged_training_seconds"] for row in records
            ])),
        ))
    return rows


def _fh_rows(bank_hash: str, gate: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    profile_lookup = {
        row["candidate"]: row for row in gate["profile_rows"]
    }
    for candidate, profile in profile_lookup.items():
        if candidate not in gate["eligible_fh_configs"]:
            rows.append({
                "method": "fh_mlp_standard",
                "candidate": candidate,
                "status": profile.get("cost_gate_status", profile["status"]),
                "n": profile["n"],
                "M": profile["M"],
                "pooled_rmse": None,
                "mean_seconds_per_query": profile.get("mean_seconds_per_query"),
                **{f"eligible_{target}": False for target in TARGETS},
            })
            continue
        n, M = int(profile["n"]), int(profile["M"])
        records = []
        for replicate in range(FH_ACCURACY_REPLICATES):
            path = RAW_DIR / "calibration" / "fh_mlp" / f"fh__calibration__n{n}_M{M}__r{replicate}.json"
            record = _read(path)
            if record["query_input_sha256"] != bank_hash:
                raise RuntimeError(f"FH bank hash mismatch: {path}")
            records.append(record)
        prediction, reference = _pooled([Path(row["prediction_file"]) for row in records])
        generator_calls = [
            (
                float(row["stats"]["terminal_samples"])
                + float(row["stats"]["transition_samples"])
            ) / float(row["points"])
            for row in records
        ]
        rows.append(_summary(
            "fh_mlp_standard",
            candidate,
            records,
            prediction,
            reference,
            n=n,
            M=M,
            mean_generator_calls_per_query=float(np.mean(generator_calls)),
            profile_mean_seconds_per_query=profile["mean_seconds_per_query"],
        ))
    return rows


def _select(rows: list[dict[str, Any]], method: str, target: str) -> dict[str, Any]:
    local = [row for row in rows if row["method"] == method]
    eligible = [row for row in local if row.get(f"eligible_{target}")]
    evaluated = [row for row in local if row.get("pooled_rmse") is not None]
    if not eligible:
        best = min(evaluated, key=lambda row: row["pooled_rmse"])
        return {
            "status": "TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET",
            "selected_candidate": None,
            "best_candidate": best["candidate"],
            "best_calibration_rmse": best["pooled_rmse"],
        }
    selected = min(
        eligible,
        key=lambda row: (row["mean_seconds_per_query"], row["candidate"]),
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
    gate = require_cost_gate()
    missing = [str(path) for path in _required_paths(gate) if not path.is_file()]
    if missing:
        raise RuntimeError(f"calibration incomplete; truth remains closed: {missing[:5]}")
    bank = load_bank("calibration", require_truth=True)
    bank_hash = array_sha256(bank["t"], bank["x"])
    rows = _ta_rows(bank_hash) + _fh_rows(bank_hash, gate)
    _write_csv(RESULTS_DIR / "calibration_rows.csv", rows)
    selections = {}
    for target, threshold in TARGETS.items():
        selections[target] = {
            "threshold": threshold,
            "ta_direct": _select(rows, "ta_direct", target),
            "ta_fk_mcv": _select(rows, "ta_fk_mcv", target),
            "fh_mlp_standard": _select(rows, "fh_mlp_standard", target),
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
        "calibration_bank_file_sha256": file_sha256(RESULTS_DIR / "banks" / "calibration.npz"),
        "calibration_rows_file_sha256": file_sha256(RESULTS_DIR / "calibration_rows.csv"),
        "cost_gate_sha256": file_sha256(COST_GATE),
        "targets": selections,
        "candidate_summaries": rows,
        "truth_policy": "calibration truth opened only after every cost-gate-eligible run completed",
    }
    payload["canonical_sha256"] = canonical_sha256(payload)
    write_json_atomic(FROZEN_SELECTIONS, payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2, sort_keys=True))
