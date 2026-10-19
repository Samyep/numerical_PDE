#!/usr/bin/env python3
"""Strictly verify frozen seed-0 whole-function evaluation arrays and metrics.

Unlike the original metadata-only guard, this audit opens every canonical
prediction artifact and hashes the numerical t, x, and truth arrays.  It also
recomputes all reported value metrics and verifies the frozen source artifact
hash, so equal strings in a CSV cannot make the audit pass by themselves.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .paired_bank_repair import (
    METHODS,
    PDES,
    RESULTS_DIR,
    canonical_bank_arrays,
    file_sha256,
    source_prediction,
    verify_canonical_prediction,
)


SOURCE_ROOT = Path(__file__).resolve().parents[2]
METRIC_NAMES = ("value_rmse", "value_relative_l2", "value_bias", "value_mae")


def _absolute(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else SOURCE_ROOT / path


def _metrics(path: Path) -> dict[str, float]:
    with np.load(path, allow_pickle=False) as data:
        truth = np.asarray(data["truth"], dtype=np.float64).reshape(-1)
        prediction = np.asarray(data["prediction"], dtype=np.float64).reshape(-1)
    error = prediction - truth
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(
            np.linalg.norm(error) / max(np.linalg.norm(truth), np.finfo(float).tiny)
        ),
        "value_bias": float(np.mean(error)),
        "value_mae": float(np.mean(np.abs(error))),
    }


def _canonical_arrays(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return tuple(np.asarray(data[name]) for name in ("t", "x", "truth", "prediction"))


def _close(left: Any, right: Any) -> bool:
    try:
        return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=1e-14)
    except (TypeError, ValueError):
        return False


def audit(main_results: Path, dimension: int = 100) -> dict[str, Any]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    with main_results.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                int(float(row["dimension"])) == dimension
                and int(float(row["seed"])) == 0
                and row["method"] in METHODS
                and row["pde"] in PDES
                and row.get("profile_role") == "primary"
            ):
                groups[row["pde"]].append(row)

    results: dict[str, Any] = {}
    paired_wins: list[str] = []
    for pde in PDES:
        rows = groups[pde]
        records: list[dict[str, Any]] = []
        for row in rows:
            prediction_path = _absolute(row.get("prediction_file", ""))
            source_path = _absolute(row.get("source_prediction_file", ""))
            errors: list[str] = []
            try:
                actual = verify_canonical_prediction(prediction_path)
                calculated = _metrics(prediction_path)
            except Exception as exc:  # audit records the exact failed artifact
                actual = {"valid": False}
                calculated = {}
                errors.append(f"canonical prediction unreadable: {exc}")

            expected_fields = {
                "states": actual.get("states"),
                "test_bank_sha256": actual.get("test_bank_array_sha256"),
                "test_bank_array_sha256": actual.get("test_bank_array_sha256"),
                "reference_array_sha256": actual.get("reference_array_sha256"),
                "prediction_file_sha256": actual.get("prediction_file_sha256"),
            }
            for name, observed in expected_fields.items():
                expected = row.get(name, "")
                if name == "states":
                    matches = _close(expected, observed)
                else:
                    matches = bool(expected) and expected == observed
                if not matches:
                    errors.append(f"{name}: CSV={expected!r}, actual={observed!r}")

            for name in METRIC_NAMES:
                if name not in calculated or not _close(row.get(name), calculated.get(name)):
                    errors.append(
                        f"{name}: CSV={row.get(name)!r}, recomputed={calculated.get(name)!r}"
                    )

            source_hash = None
            if source_path.is_file():
                source_hash = file_sha256(source_path)
                if source_hash != row.get("source_prediction_file_sha256"):
                    errors.append("frozen source prediction SHA-256 mismatch")
                try:
                    canonical_t, canonical_x, canonical_truth, canonical_prediction = _canonical_arrays(
                        prediction_path
                    )
                    extracted = source_prediction(row, source_path)
                    if not np.array_equal(np.asarray(extracted).reshape(-1), canonical_prediction):
                        errors.append("canonical prediction differs from frozen source prediction")
                    authoritative = canonical_bank_arrays(pde, dimension)
                    if not np.array_equal(canonical_t, authoritative["t"]):
                        errors.append("canonical times differ from authoritative frozen bank")
                    if not np.array_equal(canonical_x, authoritative["x"]):
                        errors.append("canonical states differ from authoritative frozen bank")
                    if not np.array_equal(canonical_truth, authoritative["truth"]):
                        errors.append("canonical truth differs from authoritative reference")
                except Exception as exc:
                    errors.append(f"source/authoritative reconstruction failed: {exc}")
            else:
                errors.append(f"frozen source prediction missing: {source_path}")

            records.append({
                "method": row["method"],
                "status": "PASS" if bool(actual.get("valid")) and not errors else "FAIL",
                "canonical_prediction_file": str(prediction_path),
                "source_prediction_file": str(source_path),
                "source_prediction_file_sha256": source_hash,
                "states": actual.get("states"),
                "test_bank_array_sha256": actual.get("test_bank_array_sha256"),
                "reference_array_sha256": actual.get("reference_array_sha256"),
                "recomputed_metrics": calculated,
                "errors": errors,
            })

        methods = sorted(record["method"] for record in records)
        bank_hashes = {record.get("test_bank_array_sha256") for record in records}
        reference_hashes = {record.get("reference_array_sha256") for record in records}
        state_counts = {record.get("states") for record in records}
        complete = methods == sorted(METHODS)
        same_arrays = (
            complete
            and len(bank_hashes) == 1
            and None not in bank_hashes
            and len(reference_hashes) == 1
            and None not in reference_hashes
            and len(state_counts) == 1
            and None not in state_counts
        )
        passed = complete and same_arrays and all(r["status"] == "PASS" for r in records)
        metric_by_method = {
            record["method"]: record["recomputed_metrics"].get("value_rmse")
            for record in records
        }
        finite_metrics = {
            method: value for method, value in metric_by_method.items()
            if value is not None and math.isfinite(value)
        }
        winner = min(finite_metrics, key=finite_metrics.get) if passed else None
        if winner == "tangent_actor":
            paired_wins.append(pde)
        results[pde] = {
            "status": "PASS" if passed else "FAIL",
            "methods_found": methods,
            "same_actual_t_x_arrays": same_arrays,
            "same_actual_reference_array": same_arrays,
            "test_bank_array_sha256": sorted(str(v) for v in bank_hashes),
            "reference_array_sha256": sorted(str(v) for v in reference_hashes),
            "states": sorted(int(v) for v in state_counts if v is not None),
            "winner_if_paired": winner,
            "records": records,
        }

    status = "PASS" if all(item["status"] == "PASS" for item in results.values()) else "FAIL"
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "dimension": dimension,
        "seed": 0,
        "pde_count": len(PDES),
        "method_count_per_pde": len(METHODS),
        "verification_level": "actual arrays, source/canonical files, and recomputed metrics",
        "tangent_actor_strict_paired_wins": paired_wins,
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--dimension", type=int, default=100)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.results_dir / "MAIN_RESULTS.csv", args.dimension)
    output = args.output or args.results_dir / f"PAIRWISE_BANK_VERIFICATION_{args.dimension}D.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "output": str(output),
        "paired_wins": result["tangent_actor_strict_paired_wins"],
    }, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
