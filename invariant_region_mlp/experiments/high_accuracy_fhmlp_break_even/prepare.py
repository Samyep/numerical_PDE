"""Create fresh banks and validate the independent periodic reference."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.single_point_joint_extension.prepare import (
    _query_inputs as previous_query_inputs,
    _reference_and_validation,
)

from .benchmark import frozen_basis, reduced_to_ambient_z
from .protocol import (
    BANK_SIZES,
    DIMENSION,
    REFERENCE_THRESHOLDS,
    RESULTS_DIR,
    REVISION,
    array_sha256,
    bank_path,
    ensure_directories,
    file_sha256,
    protocol_payload,
    source_results_dir,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


def _states(name: str, count: int) -> tuple[np.ndarray, int]:
    seed = stable_seed("query_bank", name)
    rng = np.random.default_rng(seed)
    return rng.standard_normal((count, DIMENSION), dtype=np.float64), seed


def _fingerprints(x: np.ndarray) -> np.ndarray:
    bits = np.ascontiguousarray(x, dtype=np.float64).view(np.uint64).reshape(len(x), -1)
    index = np.arange(1, bits.shape[1] + 1, dtype=np.uint64)
    first = np.bitwise_xor.reduce(
        bits * (index * np.uint64(0x9E3779B185EBCA87)), axis=1
    )
    second = np.bitwise_xor.reduce(
        (bits ^ (bits >> np.uint64(29)))
        * (index * np.uint64(0xC2B2AE3D27D4EB4F)),
        axis=1,
    )
    return np.rec.fromarrays((first, second), names=("a", "b"))


def _load_previous_small_banks() -> dict[str, np.ndarray]:
    previous: dict[str, np.ndarray] = {
        "previous_32": np.asarray(previous_query_inputs(DIMENSION)[1], dtype=np.float64),
    }
    old_bank_dir = source_results_dir() / "banks"
    for name in ("calibration", "final"):
        path = old_bank_dir / f"{name}.npz"
        if not path.is_file():
            raise RuntimeError(f"required audited prior bank is missing: {path}")
        with np.load(path, allow_pickle=False) as data:
            previous[f"previous_joint_{name}"] = np.asarray(data["x"], dtype=np.float64)
    return previous


def _exact_overlap(left: np.ndarray, right: np.ndarray) -> int:
    lookup = {row.tobytes() for row in np.ascontiguousarray(right)}
    return sum(row.tobytes() in lookup for row in np.ascontiguousarray(left))


def _separation_audit(states: dict[str, np.ndarray]) -> dict[str, Any]:
    old = _load_previous_small_banks()
    small = {"calibration": states["calibration"], "final": states["final"], **old}
    comparisons: list[dict[str, Any]] = []
    names = list(small)
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            overlap = _exact_overlap(small[left], small[right])
            comparisons.append({"left": left, "right": right, "exact_overlap_count": overlap})
            if overlap:
                raise RuntimeError(f"fresh-bank overlap detected: {left}/{right}={overlap}")

    timing = states["timing"]
    timing_fp = _fingerprints(timing)
    if len(np.unique(timing_fp)) != len(timing_fp):
        raise RuntimeError("timing stream contains duplicate row fingerprints")
    for name, values in small.items():
        candidate = np.intersect1d(timing_fp, _fingerprints(values))
        exact = _exact_overlap(timing, values) if len(candidate) else 0
        comparisons.append({
            "left": "timing",
            "right": name,
            "row_fingerprint_matches": int(len(candidate)),
            "exact_overlap_count": int(exact),
        })
        if exact:
            raise RuntimeError(f"timing stream overlaps {name}: {exact}")
    return {
        "status": "VALIDATED",
        "timing_unique_rows": int(len(np.unique(timing_fp))),
        "comparisons": comparisons,
    }


def prepare(*, force_verify: bool = False) -> dict[str, Any]:
    ensure_directories()
    manifest_path = RESULTS_DIR / "bank_manifest.json"
    audit_path = RESULTS_DIR / "reference_audit.json"
    protocol_path = RESULTS_DIR / "FROZEN_PROTOCOL.json"
    if manifest_path.is_file() and audit_path.is_file() and not force_verify:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") != "VALIDATED":
            raise RuntimeError("existing reference gate is not validated")
        return {
            "manifest": json.loads(manifest_path.read_text(encoding="utf-8")),
            "reference_audit": audit,
        }

    states: dict[str, np.ndarray] = {}
    seeds: dict[str, int] = {}
    for name, count in BANK_SIZES.items():
        states[name], seeds[name] = _states(name, count)
    separation = _separation_audit(states)

    primary, audit_solver, source_validation = _reference_and_validation()
    basis = frozen_basis(DIMENSION)
    combined_x = np.concatenate((states["calibration"], states["final"]), axis=0)
    combined_t = np.zeros(len(combined_x), dtype=np.float64)
    coordinates = combined_x @ basis
    reference = np.asarray(primary.value(combined_t, coordinates), dtype=np.float64)
    audit_reference = np.asarray(audit_solver.value(combined_t, coordinates), dtype=np.float64)
    uncertainty = np.abs(reference - audit_reference)
    primary_gradient = primary.gradient(combined_t, coordinates)
    audit_gradient = audit_solver.gradient(combined_t, coordinates)
    reference_z = reduced_to_ambient_z(DIMENSION, np.sqrt(2.0) * primary_gradient)
    audit_z = reduced_to_ambient_z(DIMENSION, np.sqrt(2.0) * audit_gradient)
    z_uncertainty = np.linalg.norm(reference_z - audit_z, axis=1)
    uncertainty_rmse = float(np.sqrt(np.mean(uncertainty**2)))
    uncertainty_max = float(np.max(uncertainty))
    finite = bool(
        np.isfinite(reference).all()
        and np.isfinite(reference_z).all()
        and np.isfinite(uncertainty).all()
    )
    gate = bool(
        finite
        and source_validation.get("overall_pass")
        and uncertainty_rmse <= REFERENCE_THRESHOLDS["value_refinement_rmse"]
        and uncertainty_max <= REFERENCE_THRESHOLDS["value_refinement_max_abs"]
    )

    cut = BANK_SIZES["calibration"]
    if gate:
        write_npz_atomic(
            bank_path("calibration"),
            t=np.zeros(cut, dtype=np.float64),
            x=states["calibration"],
            reference=reference[:cut],
            reference_z=reference_z[:cut],
            reference_uncertainty=uncertainty[:cut],
            reference_z_uncertainty=z_uncertainty[:cut],
        )
        write_npz_atomic(
            bank_path("final"),
            t=np.zeros(BANK_SIZES["final"], dtype=np.float64),
            x=states["final"],
            reference=reference[cut:],
            reference_z=reference_z[cut:],
            reference_uncertainty=uncertainty[cut:],
            reference_z_uncertainty=z_uncertainty[cut:],
        )
        write_npz_atomic(
            bank_path("timing"),
            t=np.zeros(BANK_SIZES["timing"], dtype=np.float64),
            x=states["timing"],
        )

    reference_audit = {
        "revision": REVISION,
        "status": "VALIDATED" if gate else "BLOCKED_REFERENCE_VALIDATION",
        "overall_pass": gate,
        "finite": finite,
        "fresh_bank_value_refinement_rmse": uncertainty_rmse,
        "fresh_bank_value_refinement_max_abs": uncertainty_max,
        "fresh_bank_gradient_refinement_vector_rmse": float(
            np.sqrt(np.mean(z_uncertainty**2))
        ),
        "thresholds": REFERENCE_THRESHOLDS,
        "source_reference_validation": source_validation,
        "combined_input_sha256": array_sha256(combined_t, combined_x),
        "combined_reference_sha256": array_sha256(reference),
        "combined_reference_z_sha256": array_sha256(reference_z),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(audit_path, reference_audit)
    if not gate:
        raise RuntimeError("BLOCKED_REFERENCE_VALIDATION")

    rows = {}
    for name in BANK_SIZES:
        path = bank_path(name)
        rows[name] = {
            "size": BANK_SIZES[name],
            "seed": seeds[name],
            "law": "t=0; x~N(0,I_100)",
            "path": str(path.relative_to(RESULTS_DIR.parents[2])),
            "file_sha256": file_sha256(path),
            "input_array_sha256": array_sha256(
                np.zeros(BANK_SIZES[name]), states[name]
            ),
            "contains_reference": name != "timing",
        }
    manifest = {
        "revision": REVISION,
        "protocol": protocol_payload(),
        "banks": rows,
        "basis_sha256": array_sha256(basis),
        "separation_audit": separation,
        "truth_policy": "calibration/final only; timing stream contains no truth",
    }
    write_json_atomic(manifest_path, manifest)
    write_json_atomic(protocol_path, protocol_payload())
    return {"manifest": manifest, "reference_audit": reference_audit}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-verify", action="store_true")
    args = parser.parse_args()
    result = prepare(force_verify=args.force_verify)
    print(json.dumps({
        "status": result["reference_audit"]["status"],
        "value_refinement_rmse": result["reference_audit"][
            "fresh_bank_value_refinement_rmse"
        ],
        "value_refinement_max_abs": result["reference_audit"][
            "fresh_bank_value_refinement_max_abs"
        ],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
