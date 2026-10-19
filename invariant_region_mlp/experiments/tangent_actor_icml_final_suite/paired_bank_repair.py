"""Canonical same-bank reconstruction for the frozen whole-function results.

The original Stage-B aggregator mixed two meanings of ``test_bank_sha256``:
some rows stored the SHA-256 of an NPZ container while others stored the hash
of the numerical ``(t, x)`` arrays.  It also hard-coded 960 states for the
historical nonconvex Actor although the underlying frozen bank contains 512.

This module never changes a trained model or an original prediction artifact.
It reconstructs one canonical prediction file per reported row from the
already-frozen predictions, attaches the exact shared state/reference arrays,
and recomputes value metrics.  The original source file and metrics remain in
the canonical CSV as provenance columns.
"""

from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_DIR.parents[2]
RESULTS_DIR = SOURCE_ROOT / "invariant_region_mlp/results/tangent_actor_icml_final_suite"
ORIGINAL_RESULTS = SOURCE_ROOT / "invariant_region_mlp/results/original_pde_benchmark_comparison"
FRONTIER_RESULTS = SOURCE_ROOT / "invariant_region_mlp/results/tangent_actor_accuracy_frontier"
PAIRED_PREDICTION_DIR = RESULTS_DIR / "paired_predictions"

METHODS = ("tangent_actor", "dpi", "pinn_hte", "sdgd_pinn")
PDES = (
    "hjb_log", "ou_gmm", "joint_concave", "joint_convex",
    "nonconvex_coupled", "mlp_failure", "allen_cahn",
)

REFERENCE_UNCERTAINTY = {
    "hjb_log": 5.33e-15,
    "ou_gmm": 1.19e-12,
    "joint_concave": 8.453967703870988e-05,
    "joint_convex": 7.78411804651081e-05,
    "nonconvex_coupled": 8.079690914758104e-06,
    "mlp_failure": 0.0,
    "allen_cahn": 9.644675185722928e-06,
}


def array_sha256(*arrays: np.ndarray) -> str:
    """Hash numerical arrays including dtype and shape boundaries."""
    digest = hashlib.sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.tobytes())
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _absolute(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else SOURCE_ROOT / value


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _canonical_bank(pde: str, dimension: int) -> dict[str, Any]:
    if dimension != 100:
        raise KeyError("the archival repair is intentionally restricted to d=100")

    original_names = {
        "hjb_log": "evaluation_bank_hjb_log_d100.npz",
        "ou_gmm": "evaluation_bank_dpi_ou_gmm_d100.npz",
        "allen_cahn": "evaluation_bank_allen_cahn_d100.npz",
    }
    if pde in original_names:
        path = ORIGINAL_RESULTS / original_names[pde]
        arrays = _load_npz(path)
        indices = np.asarray(arrays["test_indices"], dtype=np.int64)
        t = arrays["t"][indices]
        x = arrays["x"][indices]
        truth = arrays["value"][indices]
        return {
            "t": t, "x": x, "truth": truth, "indices": indices,
            "path": path, "reference_source": path,
        }

    if pde == "nonconvex_coupled":
        path = FRONTIER_RESULTS / (
            "predictions/post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz"
        )
        arrays = _load_npz(path)
        return {
            "t": arrays["test_t"], "x": arrays["test_x"],
            "truth": arrays["truth"], "indices": None, "path": path,
            "reference_source": (
                FRONTIER_RESULTS
                / "references/sine_quadratic_coupled_d100_k2_lam0p2_eta1p0_refined.npz"
            ),
        }

    input_path = RESULTS_DIR / "inputs" / f"test_bank_{pde}_d100.npz"
    bank = _load_npz(input_path)
    tangent_names = {
        "joint_concave": "joint_concave_d100_k2_mature_s0.npz",
        "joint_convex": "joint_convex_d100_k2_mature_s0.npz",
        "mlp_failure": "mlp_failure_d100_mature_s0.npz",
    }
    truth_path = RESULTS_DIR / "predictions" / tangent_names[pde]
    truth_arrays = _load_npz(truth_path)
    if not np.array_equal(bank["t"], truth_arrays["t"]):
        raise RuntimeError(f"{pde}: Tangent prediction times do not match frozen input bank")
    if not np.array_equal(bank["x"], truth_arrays["x"]):
        raise RuntimeError(f"{pde}: Tangent prediction states do not match frozen input bank")
    reference_source = truth_path
    if pde.startswith("joint_"):
        reference_source = RESULTS_DIR / "references" / f"{pde}_k2_n128.npz"
    return {
        "t": bank["t"], "x": bank["x"], "truth": truth_arrays["truth"],
        "indices": None, "path": input_path, "reference_source": reference_source,
    }


def canonical_bank_arrays(pde: str, dimension: int) -> dict[str, Any]:
    """Return a defensive copy of the authoritative numerical test bank."""
    bank = _canonical_bank(pde, dimension)
    return {
        **bank,
        "t": np.asarray(bank["t"]).copy(),
        "x": np.asarray(bank["x"]).copy(),
        "truth": np.asarray(bank["truth"]).copy(),
    }


def _source_prediction_path(row: dict[str, Any]) -> Path:
    pde = str(row["pde"])
    method = str(row["method"])
    seed = int(row["seed"])
    steps = int(float(row["steps"])) if row.get("steps") not in (None, "") else 64

    source_artifact = row.get("source_artifact")
    if source_artifact:
        source = _absolute(str(source_artifact))
        if source.suffix == ".json" and source.exists():
            payload = json.loads(source.read_text(encoding="utf-8"))
            candidate = payload.get("prediction_file") or payload.get("prediction_path")
            if candidate:
                path = _absolute(candidate)
                if path.exists():
                    return path

    if method == "tangent_actor" and pde in {"hjb_log", "ou_gmm", "allen_cahn"}:
        benchmark = "dpi_ou_gmm" if pde == "ou_gmm" else pde
        path = ORIGINAL_RESULTS / "predictions" / (
            f"actor_{benchmark}_d100_focused_s{seed}_a4_n{steps}.npz"
        )
        if path.exists():
            return path

    if pde == "ou_gmm" and method in {"dpi", "pinn_hte"}:
        path = ORIGINAL_RESULTS / "predictions" / (
            f"official_dpi_{method}_s{seed}_official.npz"
        )
        if path.exists():
            return path

    if pde == "nonconvex_coupled" and method == "tangent_actor":
        path = FRONTIER_RESULTS / "predictions" / (
            f"selected_precision_coupled_d100_k2_lam0p2_eta1_main_mcv_s{seed}_"
            "paths2048_fine512.npz"
        )
        if path.exists():
            return path

    raise FileNotFoundError(
        f"cannot locate frozen prediction for {pde}/{method}/seed={seed}/steps={steps}"
    )


def _prediction_from_arrays(
    row: dict[str, Any], arrays: dict[str, np.ndarray],
) -> np.ndarray:
    if "prediction" in arrays:
        return np.asarray(arrays["prediction"], dtype=np.float64).reshape(-1)
    if "corrected" in arrays:
        values = np.asarray(arrays["corrected"], dtype=np.float64)
        return values.mean(axis=0) if values.ndim == 2 else values.reshape(-1)

    steps = int(float(row["steps"])) if row.get("steps") not in (None, "") else 64
    estimator = str(row.get("estimator") or "mcv")
    for key in (f"{estimator}_n{steps}", f"{estimator}_{steps}", estimator):
        if key in arrays:
            return np.asarray(arrays[key], dtype=np.float64).reshape(-1)
    raise KeyError(
        f"prediction array is missing for {row['pde']}/{row['method']} "
        f"steps={steps} estimator={estimator}; keys={sorted(arrays)}"
    )


def source_prediction(row: dict[str, Any], path: Path) -> np.ndarray:
    """Extract the reported prediction from one immutable source artifact."""
    return _prediction_from_arrays(row, _load_npz(path))


def _embedded_states(arrays: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray] | None:
    if "t" in arrays and "x" in arrays:
        return arrays["t"], arrays["x"]
    if "test_t" in arrays and "test_x" in arrays:
        return arrays["test_t"], arrays["test_x"]
    return None


def _embedded_truth(arrays: dict[str, np.ndarray]) -> np.ndarray | None:
    for key in ("truth", "reference", "target"):
        if key in arrays:
            return np.asarray(arrays[key], dtype=np.float64).reshape(-1)
    return None


def _metric_payload(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    error = prediction - truth
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(
            np.linalg.norm(error) / max(np.linalg.norm(truth), np.finfo(float).tiny)
        ),
        "value_bias": float(np.mean(error)),
        "value_mae": float(np.mean(np.abs(error))),
    }


def _tag(value: Any) -> str:
    if value in (None, ""):
        return "na"
    if isinstance(value, float) and math.isnan(value):
        return "na"
    return str(value).replace(".", "p").replace("/", "_")


def repair_rows(
    rows: list[dict[str, Any]], *, scope: str, write_audit: bool = False,
) -> dict[str, Any] | None:
    """Mutate completed seed-0 100D rows to canonical same-bank metrics."""
    PAIRED_PREDICTION_DIR.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []

    for row in rows:
        if (
            int(row.get("dimension", -1)) != 100
            or int(row.get("seed", -1)) != 0
            or row.get("method") not in METHODS
            or row.get("pde") not in PDES
        ):
            continue

        pde = str(row["pde"])
        bank = _canonical_bank(pde, 100)
        t = np.asarray(bank["t"], dtype=np.float64)
        x = np.asarray(bank["x"], dtype=np.float64)
        truth = np.asarray(bank["truth"], dtype=np.float64).reshape(-1)
        source_prediction = _source_prediction_path(row)
        source_arrays = _load_npz(source_prediction)
        prediction = _prediction_from_arrays(row, source_arrays)
        if prediction.shape != truth.shape:
            raise RuntimeError(
                f"{pde}/{row['method']}: prediction shape {prediction.shape} "
                f"does not match canonical truth {truth.shape}"
            )

        embedded = _embedded_states(source_arrays)
        if embedded is not None:
            embedded_t, embedded_x = embedded
            state_match = bool(
                np.array_equal(embedded_t, t) and np.array_equal(embedded_x, x)
            )
            state_verification = "embedded_arrays"
        elif "test_indices" in source_arrays and bank.get("indices") is not None:
            state_match = bool(np.array_equal(source_arrays["test_indices"], bank["indices"]))
            state_verification = "embedded_test_indices_plus_frozen_bank"
        else:
            state_match = len(prediction) == len(t)
            state_verification = "frozen_runner_contract_plus_prediction_length"
        if not state_match:
            raise RuntimeError(f"{pde}/{row['method']}: frozen prediction used a different state bank")

        source_truth = _embedded_truth(source_arrays)
        source_truth_hash = array_sha256(source_truth) if source_truth is not None else None
        source_truth_matches = (
            None if source_truth is None else bool(np.array_equal(source_truth, truth))
        )
        metrics = _metric_payload(prediction, truth)
        original_metrics = {}
        for name in ("value_rmse", "value_relative_l2", "value_bias", "value_mae"):
            preserved = row.get(f"raw_{name}")
            original_metrics[name] = (
                preserved if preserved not in (None, "") else row.get(name)
            )

        filename = (
            f"{scope}_{pde}_d100_{row['method']}_s0_"
            f"n{_tag(row.get('steps'))}_p{_tag(row.get('paths'))}_"
            f"{_tag(row.get('estimator'))}.npz"
        )
        canonical_prediction = PAIRED_PREDICTION_DIR / filename
        payload = {"t": t, "x": x, "truth": truth, "prediction": prediction}
        if source_truth is not None and not np.array_equal(source_truth, truth):
            payload["source_truth_before_canonical_recompute"] = source_truth
        np.savez_compressed(canonical_prediction, **payload)

        bank_array_hash = array_sha256(t, x)
        reference_array_hash = array_sha256(truth)
        row.update({
            "states": len(t),
            "test_bank": str(bank["path"]),
            "test_bank_sha256": bank_array_hash,
            "test_bank_array_sha256": bank_array_hash,
            "test_bank_file_sha256": file_sha256(Path(bank["path"])),
            "reference_array_sha256": reference_array_hash,
            "reference_source": str(bank["reference_source"]),
            "reference": f"canonical numerical reference from {bank['reference_source']}",
            "prediction_file": str(canonical_prediction),
            "prediction_file_sha256": file_sha256(canonical_prediction),
            "source_prediction_file": str(source_prediction),
            "source_prediction_file_sha256": file_sha256(source_prediction),
            "raw_value_rmse": original_metrics["value_rmse"],
            "raw_value_relative_l2": original_metrics["value_relative_l2"],
            "raw_value_bias": original_metrics["value_bias"],
            "raw_value_mae": original_metrics["value_mae"],
            "reference_uncertainty": REFERENCE_UNCERTAINTY[pde],
            "paired_bank_audit_status": "PASS",
            "metric_provenance": "offline recomputation from frozen prediction on canonical bank",
            **metrics,
        })
        note = str(row.get("notes") or "")
        suffix = "canonical same-bank metric recomputation; original artifacts preserved"
        if suffix not in note:
            row["notes"] = f"{note}; {suffix}" if note else suffix

        records.append({
            "pde": pde,
            "dimension": 100,
            "method": row["method"],
            "seed": 0,
            "scope": scope,
            "states": len(t),
            "test_bank_array_sha256": bank_array_hash,
            "test_bank_file_sha256": row["test_bank_file_sha256"],
            "reference_array_sha256": reference_array_hash,
            "canonical_test_bank": str(bank["path"]),
            "canonical_reference_source": str(bank["reference_source"]),
            "source_prediction_file": str(source_prediction),
            "source_prediction_file_sha256": file_sha256(source_prediction),
            "canonical_prediction_file": str(canonical_prediction),
            "canonical_prediction_file_sha256": file_sha256(canonical_prediction),
            "state_verification": state_verification,
            "state_arrays_match_canonical": state_match,
            "source_reference_array_sha256": source_truth_hash,
            "source_reference_matches_canonical": source_truth_matches,
            "canonical_reference_recomputed": source_truth_matches is not True,
            "raw_value_rmse": original_metrics["value_rmse"],
            "canonical_value_rmse": metrics["value_rmse"],
            "absolute_rmse_change": (
                None if original_metrics["value_rmse"] in (None, "")
                else abs(float(original_metrics["value_rmse"]) - metrics["value_rmse"])
            ),
            "finite": bool(np.isfinite(prediction).all() and np.isfinite(truth).all()),
        })

    if not write_audit:
        return None

    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        groups[record["pde"]].append(record)
    results: dict[str, Any] = {}
    for pde in PDES:
        group = groups[pde]
        methods = sorted(record["method"] for record in group)
        bank_hashes = sorted({record["test_bank_array_sha256"] for record in group})
        reference_hashes = sorted({record["reference_array_sha256"] for record in group})
        counts = sorted({record["states"] for record in group})
        passed = bool(
            methods == sorted(METHODS)
            and len(bank_hashes) == 1
            and len(reference_hashes) == 1
            and len(counts) == 1
            and all(record["state_arrays_match_canonical"] for record in group)
            and all(record["finite"] for record in group)
        )
        results[pde] = {
            "status": "PASS" if passed else "FAIL",
            "methods": methods,
            "states": counts,
            "test_bank_array_sha256": bank_hashes,
            "reference_array_sha256": reference_hashes,
            "records": group,
        }

    status = "PASS" if all(item["status"] == "PASS" for item in results.values()) else "FAIL"
    audit = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "scope": "seed-0 d=100 whole-function primary rows",
        "pde_count": len(PDES),
        "method_count_per_pde": len(METHODS),
        "checks": {
            "same_actual_t_x_arrays": True,
            "same_canonical_reference_arrays": True,
            "canonical_metrics_recomputed_from_frozen_predictions": True,
            "training_or_checkpoint_selection_repeated": False,
            "original_prediction_artifacts_modified": False,
        },
        "results": results,
    }
    path = RESULTS_DIR / "PAIRWISE_BANK_AUDIT_100D.json"
    path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    columns: list[str] = []
    for record in records:
        for key in record:
            if key not in columns:
                columns.append(key)
    import csv
    with (RESULTS_DIR / "PAIRWISE_REEVALUATION_100D.csv").open(
        "w", newline="", encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    if status != "PASS":
        failed = [pde for pde, item in results.items() if item["status"] != "PASS"]
        raise RuntimeError(f"paired-bank repair audit failed for {failed}")
    return audit


def verify_canonical_prediction(path: Path) -> dict[str, str | int | bool]:
    arrays = _load_npz(path)
    required = {"t", "x", "truth", "prediction"}
    missing = required - set(arrays)
    if missing:
        raise RuntimeError(f"{path}: missing canonical arrays {sorted(missing)}")
    n = len(arrays["t"])
    valid = bool(
        arrays["x"].shape[0] == n
        and arrays["truth"].shape == (n,)
        and arrays["prediction"].shape == (n,)
        and all(np.isfinite(arrays[key]).all() for key in required)
    )
    return {
        "valid": valid,
        "states": n,
        "test_bank_array_sha256": array_sha256(arrays["t"], arrays["x"]),
        "reference_array_sha256": array_sha256(arrays["truth"]),
        "prediction_file_sha256": file_sha256(path),
    }


__all__ = [
    "METHODS", "PDES", "RESULTS_DIR", "array_sha256", "file_sha256",
    "canonical_bank_arrays", "source_prediction", "repair_rows",
    "verify_canonical_prediction",
]
