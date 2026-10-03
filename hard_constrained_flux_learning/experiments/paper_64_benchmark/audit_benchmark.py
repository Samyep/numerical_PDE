"""Fail-closed integrity audit for the final 64-cell benchmark artifacts."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent


def finite_at_most(value: Any, tolerance: float, label: str) -> None:
    if value is None or not math.isfinite(float(value)):
        raise RuntimeError(f"{label} is missing or nonfinite: {value}")
    if float(value) > tolerance:
        raise RuntimeError(f"{label}={value} exceeds {tolerance}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    output = args.results_dir.resolve()
    summary_path = output / "benchmark64_summary.json"
    replicate_path = output / "hcfl64_replicate_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    replicates = json.loads(replicate_path.read_text(encoding="utf-8"))

    protocol = summary["protocol"]
    if protocol["candidate_cells"] != 64 or protocol["reference_cells"] != 2048:
        raise RuntimeError("Final benchmark is not the required 64-vs-2048 protocol")
    for system in ("Euler", "SWE"):
        if not replicates["aggregate"][system]["all_converged"]:
            raise RuntimeError(f"{system} replicate set contains a nonconverged run")
        if len(replicates["aggregate"][system]["seeds"]) != 3:
            raise RuntimeError(f"{system} does not contain exactly three HCFL seeds")

    required_methods = {
        "native_fvm64",
        "muscl_fvm64",
        "hcfl64",
        "hcfl64_no_low",
        "learned_flux64",
        "fno64",
    }
    present = {row["method"] for row in summary["aggregate"]}
    if present != required_methods:
        raise RuntimeError(f"Benchmark method set mismatch: {sorted(present)}")

    checked_rows = 0
    for row in summary["aggregate"]:
        if row["method"] not in ("hcfl64", "hcfl64_no_low"):
            continue
        checked_rows += 1
        if row["mean_completion_rate"] != 1.0:
            raise RuntimeError(
                f"{row['system']} {row['method']} seed {row['seed']} did not complete"
            )
        if row["minimum_density_or_depth"] <= 0.0:
            raise RuntimeError("HCFL lost positive density/depth")
        if row["system"] == "euler" and row["minimum_pressure"] <= 0.0:
            raise RuntimeError("Euler HCFL lost positive pressure")
        finite_at_most(
            row["maximum_interface_entropy_residual"],
            1.0e-5,
            "maximum interface entropy residual",
        )
        finite_at_most(
            row["maximum_internal_total_entropy_change"],
            1.0e-8,
            "maximum fully-discrete entropy change",
        )
        finite_at_most(
            row["maximum_relative_conservation_drift"],
            1.0e-5,
            "maximum conservation drift",
        )

    expected_hcfl_rows = 2 * 2 * 2 * 3
    if checked_rows != expected_hcfl_rows:
        raise RuntimeError(
            f"Expected {expected_hcfl_rows} HCFL aggregate rows, found {checked_rows}"
        )

    audit = {
        "status": "pass",
        "protocol": "64 candidate cells; 2048-cell scoring reference",
        "hcfl_aggregate_rows_checked": checked_rows,
        "checks": {
            "three_converged_seeds_per_system": "pass",
            "all_requested_same-task_baselines_present": "pass",
            "hcfl_physical_completion": "pass",
            "periodic_conservation_tolerance_1e-5": "pass",
            "interface_entropy_tolerance_1e-5": "pass",
            "fully_discrete_entropy_tolerance_1e-8": "pass",
        },
    }
    (output / "benchmark64_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
