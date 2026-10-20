"""Freeze the B6/Dense pilot compatibility and reuse decision.

The old B6 scalar records are valuable evidence, but neither its trained
checkpoints nor paired 960-state prediction arrays are tracked in the source
branch.  A dimension match therefore cannot be silently promoted to a formal
common-bank result.  This audit records that decision before formal training.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
from typing import Any

from .benchmarks import array_sha256, dense_parameters, fixed_b6_points
from .protocol import (
    DIMENSIONS,
    RESULTS_DIR,
    SOURCE_COMMITS,
    SOURCE_ROOT,
    write_json_atomic,
)


OLD_B6_FULL_BANK_HASHES = {
    100: "84b2132607bc6b17e44adbea838915f8ce657813f6a0987e34deca308660f1c1",
    400: "0df59abc155e90aa674532760578dc2563cbdb4aade97a849961c8a8edd9f2f9",
}
OLD_DENSE_D100_PARAMETER_HASH = (
    "b96804e2f7c20cbd4e72c0f961a4ebfc21bf46cae120d198b6afe6ce2e2684b6"
)


def _tracked_checkpoint_count(revision: str, prefix: str) -> int:
    output = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", revision, "--", prefix],
        cwd=SOURCE_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return sum(path.endswith(".pt") for path in output.splitlines())


def audit() -> dict[str, Any]:
    b6 = []
    for dimension in (100, 400):
        t, x, split = fixed_b6_points(dimension)
        actual = array_sha256(t, x, split)
        b6.append({
            "dimension": dimension,
            "old_full_bank_sha256": OLD_B6_FULL_BANK_HASHES[dimension],
            "reproduced_full_bank_sha256": actual,
            "bank_exactly_compatible": actual == OLD_B6_FULL_BANK_HASHES[dimension],
            "mature_training_configuration_compatible": True,
            "old_evaluation_paths": 256,
            "formal_evaluation_paths": 1024,
            "evaluation_protocol_exactly_compatible": False,
            "tracked_checkpoint_available": False,
            "tracked_common_bank_predictions_available": False,
            "reuse_decision": (
                "retrain and reevaluate under frozen formal protocol; retain old scalar result "
                "as historical evidence only"
            ),
        })
    omega, signs, phases = dense_parameters(100)
    dense_hash = array_sha256(omega, signs, phases)
    tracked_b6 = _tracked_checkpoint_count(
        SOURCE_COMMITS["b6_validated"], "invariant_region_mlp/results/highdim_nosym"
    )
    payload = {
        "status": "PASS",
        "b6": b6,
        "b6_source_commit": SOURCE_COMMITS["b6_validated"],
        "b6_tracked_checkpoint_count_in_source_commit": tracked_b6,
        "dense_fourier": {
            "dimension": 100,
            "old_parameter_sha256": OLD_DENSE_D100_PARAMETER_HASH,
            "reproduced_parameter_sha256": dense_hash,
            "parameter_exactly_compatible": dense_hash == OLD_DENSE_D100_PARAMETER_HASH,
            "pilot_budget_is_mature_protocol": False,
            "reuse_decision": "correctness evidence only; run new mature formal training",
        },
        "formal_dimensions": list(DIMENSIONS),
        "no_historical_accuracy_copied_into_formal_table": True,
    }
    if tracked_b6 != 0 or not all(row["bank_exactly_compatible"] for row in b6):
        payload["status"] = "FAIL"
    if dense_hash != OLD_DENSE_D100_PARAMETER_HASH:
        payload["status"] = "FAIL"
    write_json_atomic(RESULTS_DIR / "reuse_audit.json", payload)
    if payload["status"] != "PASS":
        raise RuntimeError(json.dumps(payload, indent=2))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
