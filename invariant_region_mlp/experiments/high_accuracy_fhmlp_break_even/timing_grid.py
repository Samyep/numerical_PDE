"""Freeze deterministic timing refinements from truth-free profile slopes."""

from __future__ import annotations

from datetime import datetime, timezone
import json

import numpy as np

from .protocol import (
    BANK_SIZES,
    EXPECTED_TA_ARTIFACTS,
    FROZEN_TIMING_GRID,
    GEOMETRIC_TIMING_COUNTS,
    RESULTS_DIR,
    REVISION,
    TA_TRAINING_SEEDS,
    TARGETS,
    canonical_sha256,
    refined_timing_counts,
    require_frozen_selections,
    write_json_atomic,
)


def freeze() -> dict:
    if FROZEN_TIMING_GRID.is_file():
        raise RuntimeError("FROZEN_TIMING_GRID already exists; refusing to overwrite")
    frozen = require_frozen_selections()
    final = json.loads(
        (RESULTS_DIR / "final_accuracy_summary.json").read_text(encoding="utf-8")
    )
    profile = json.loads(
        (RESULTS_DIR / "performance_profile.json").read_text(encoding="utf-8")
    )
    profile_lookup = {row["candidate"]: row for row in profile["rows"]}
    charged = float(np.mean([
        EXPECTED_TA_ARTIFACTS[int(seed)]["charged_training_seconds"]
        for seed in TA_TRAINING_SEEDS
    ]))
    # Conservatively ignore the Actor's tiny per-query term in the pilot
    # estimate.  Formal timing executes and records it; no accuracy-bearing
    # artifact is used to place the refinement points.
    direct_per_query = 0.0
    estimates = []
    counts = set(GEOMETRIC_TIMING_COUNTS)
    for target in TARGETS:
        local = final["targets"][target]
        if not (
            local["ta_direct"]["status"] == "FINAL_ACCURACY_PASSED"
            and local["fh_mlp_standard"]["status"] == "FINAL_ACCURACY_PASSED"
        ):
            estimates.append({"target": target, "status": "PAIR_NOT_ACCURACY_QUALIFIED"})
            continue
        candidate = local["fh_mlp_standard"]["selected_candidate"]
        fh_per_query = float(profile_lookup[candidate]["mean_seconds_per_query"])
        denominator = fh_per_query - direct_per_query
        estimate = (
            charged / denominator
            if denominator > 0
            else float(BANK_SIZES["timing"])
        )
        estimate = min(max(estimate, 1.0), float(BANK_SIZES["timing"]))
        refined = refined_timing_counts(estimate)
        counts.update(refined)
        estimates.append({
            "target": target,
            "status": "ESTIMATED_FROM_TRUTH_FREE_PROFILE",
            "fh_candidate": candidate,
            "charged_ta_fixed_seconds": charged,
            "ta_direct_profile_seconds_per_query": direct_per_query,
            "fh_profile_seconds_per_query": fh_per_query,
            "estimated_crossing_queries": estimate,
            "refined_count": len(refined),
        })
    payload = {
        "revision": REVISION,
        "status": "FROZEN",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "pilot_truth_used": False,
        "pilot_source": "actual H100 four-query profile, two replicates per FH candidate; zero-order TA direct slope",
        "refinement_rule": "power-of-two grid plus octiles of estimated bracket and its immediate lower/upper brackets",
        "estimates": estimates,
        "timing_counts": sorted(
            value for value in counts if 1 <= value <= BANK_SIZES["timing"]
        ),
    }
    payload["canonical_sha256"] = canonical_sha256(payload)
    write_json_atomic(FROZEN_TIMING_GRID, payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(freeze(), indent=2, sort_keys=True))
