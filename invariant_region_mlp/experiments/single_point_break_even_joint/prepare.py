"""Create fresh banks and validate the independent two-dimensional reference."""

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
    RESULTS_DIR,
    REVISION,
    array_sha256,
    bank_path,
    ensure_directories,
    file_sha256,
    protocol_payload,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


def _states(name: str, count: int) -> tuple[np.ndarray, int]:
    seed = stable_seed("query_bank", name)
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((count, DIMENSION), dtype=np.float64)
    if name != "timing":
        x[0] = 0.0
    return x, seed


def _fingerprints(x: np.ndarray) -> np.ndarray:
    """Two vectorized uint64 row hashes, followed by exact checks on matches."""
    bits = np.ascontiguousarray(x, dtype=np.float64).view(np.uint64).reshape(len(x), -1)
    index = np.arange(1, bits.shape[1] + 1, dtype=np.uint64)
    first = np.bitwise_xor.reduce(bits * (index * np.uint64(0x9E3779B185EBCA87)), axis=1)
    second = np.bitwise_xor.reduce((bits ^ (bits >> np.uint64(29))) * (index * np.uint64(0xC2B2AE3D27D4EB4F)), axis=1)
    return np.rec.fromarrays((first, second), names=("a", "b"))


def _separation_audit(states: dict[str, np.ndarray]) -> dict[str, Any]:
    old_x = previous_query_inputs(DIMENSION)[1]
    origin = np.zeros(DIMENSION, dtype=np.float64)
    small = {
        "calibration": states["calibration"],
        "final": states["final"],
        "previous_32": old_x,
    }
    comparisons: list[dict[str, Any]] = []
    for left, right in (("calibration", "final"), ("calibration", "previous_32"), ("final", "previous_32")):
        matches = []
        for i, row in enumerate(small[left]):
            hit = np.flatnonzero(np.all(small[right] == row, axis=1))
            matches.extend((i, int(j)) for j in hit)
        unexpected = [pair for pair in matches if not (
            np.array_equal(small[left][pair[0]], origin)
            and np.array_equal(small[right][pair[1]], origin)
        )]
        comparisons.append({
            "left": left,
            "right": right,
            "exact_overlap_count": len(matches),
            "mandated_origin_overlap": any(
                np.array_equal(small[left][i], origin) and np.array_equal(small[right][j], origin)
                for i, j in matches
            ),
            "unexpected_overlap_count": len(unexpected),
        })
        if unexpected:
            raise RuntimeError(f"unexpected exact query overlap: {left}/{right}")

    timing_fp = _fingerprints(states["timing"])
    unique_timing = len(np.unique(timing_fp))
    if unique_timing != len(timing_fp):
        raise RuntimeError("timing bank contains a duplicate row fingerprint")
    timing_pairs = []
    for name, values in small.items():
        fingerprints = _fingerprints(values)
        candidate = np.intersect1d(timing_fp, fingerprints)
        exact = 0
        if len(candidate):
            # Hash matches are rare; confirm every candidate using exact bytes.
            local_lookup = {row.tobytes() for row in np.ascontiguousarray(values)}
            exact = sum(row.tobytes() in local_lookup for row in np.ascontiguousarray(states["timing"]))
        if exact:
            raise RuntimeError(f"timing bank overlaps {name}: {exact} rows")
        timing_pairs.append({
            "left": "timing",
            "right": name,
            "row_fingerprint_matches": int(len(candidate)),
            "exact_overlap_count": int(exact),
        })
    return {
        "status": "VALIDATED_WITH_MANDATED_ORIGIN_EXCEPTION",
        "hash_method": "two independent vectorized uint64 fingerprints; all matches exactly confirmed",
        "timing_unique_rows": unique_timing,
        "comparisons": comparisons + timing_pairs,
    }


def prepare(*, force_verify: bool = False) -> dict[str, Any]:
    ensure_directories()
    manifest_path = RESULTS_DIR / "bank_manifest.json"
    audit_path = RESULTS_DIR / "reference_audit.json"
    if manifest_path.is_file() and audit_path.is_file() and not force_verify:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") != "VALIDATED":
            raise RuntimeError("existing reference gate is not validated")
        return {"manifest": json.loads(manifest_path.read_text(encoding="utf-8")), "reference_audit": audit}

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
    finite = bool(np.isfinite(reference).all() and np.isfinite(uncertainty).all())
    source_pass = bool(source_validation.get("overall_pass"))
    gate = bool(finite and source_pass and uncertainty_rmse <= 5e-5 and uncertainty_max <= 2e-4)

    cut = BANK_SIZES["calibration"]
    if gate:
        write_npz_atomic(
            bank_path("calibration"),
            t=np.zeros(cut, dtype=np.float64), x=states["calibration"],
            reference=reference[:cut], reference_z=reference_z[:cut],
            reference_uncertainty=uncertainty[:cut],
            reference_z_uncertainty=z_uncertainty[:cut],
        )
        write_npz_atomic(
            bank_path("final"),
            t=np.zeros(BANK_SIZES["final"], dtype=np.float64), x=states["final"],
            reference=reference[cut:], reference_z=reference_z[cut:],
            reference_uncertainty=uncertainty[cut:],
            reference_z_uncertainty=z_uncertainty[cut:],
        )
        write_npz_atomic(
            bank_path("timing"),
            t=np.zeros(BANK_SIZES["timing"], dtype=np.float64), x=states["timing"],
        )

    reference_audit = {
        "revision": REVISION,
        "status": "VALIDATED" if gate else "BLOCKED_REFERENCE_VALIDATION",
        "overall_pass": gate,
        "finite": finite,
        "fresh_bank_value_refinement_rmse": uncertainty_rmse,
        "fresh_bank_value_refinement_max_abs": uncertainty_max,
        "fresh_bank_gradient_refinement_vector_rmse": float(np.sqrt(np.mean(z_uncertainty**2))),
        "thresholds": {"value_refinement_rmse": 5e-5, "value_refinement_max_abs": 2e-4},
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
            "law": "t=0; x~N(0,I_100)" + ("; row zero replaced by exact origin" if name != "timing" else ""),
            "path": str(path.relative_to(RESULTS_DIR.parents[2])),
            "file_sha256": file_sha256(path),
            "input_array_sha256": array_sha256(np.zeros(BANK_SIZES[name]), states[name]),
            "contains_reference": name != "timing",
        }
    manifest = {
        "revision": REVISION,
        "protocol": protocol_payload(),
        "banks": rows,
        "basis_sha256": array_sha256(basis),
        "separation_audit": separation,
        "truth_policy": "calibration/final only; timing bank contains no truth",
    }
    write_json_atomic(manifest_path, manifest)
    return {"manifest": manifest, "reference_audit": reference_audit}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-verify", action="store_true")
    args = parser.parse_args()
    result = prepare(force_verify=args.force_verify)
    print(json.dumps({
        "status": result["reference_audit"]["status"],
        "value_refinement_rmse": result["reference_audit"]["fresh_bank_value_refinement_rmse"],
        "value_refinement_max_abs": result["reference_audit"]["fresh_bank_value_refinement_max_abs"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
