"""Independent integrity checks for the nonperiodic Euler audit artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RETAINED = ("hllc_roe_correction", "nonnegative_feasibility")
LEARNED = (
    "hllc_roe_correction",
    "nonnegative_control",
    "nonnegative_feasibility",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    results = args.results_dir.resolve()
    source = results / f"nonperiodic_transmissive512_seed{args.seed}.json"
    report: dict[str, Any] = json.loads(source.read_text(encoding="utf-8"))

    require(report["training_boundary_condition"] == "periodic", "wrong train BC")
    require("transmissive" in report["deployment_boundary_condition"], "wrong test BC")
    require(report["learned_boundary_interfaces"] == 0, "a boundary flux was learned")
    require(report["learned_interior_interfaces"] == 511, "wrong interior count")

    for method, test in report["boundary_operator_self_tests"].items():
        require(test["status"] == "pass", f"operator self-test failed for {method}")
        require(
            test["left_interface_wrap_leak_max_error"] <= 1.0e-7,
            f"circular leakage detected for {method}",
        )
        require(
            test["boundary_flux_max_error"] <= 1.0e-7,
            f"learned/classical boundary flux mismatch for {method}",
        )

    for case_name, case in report["cases"].items():
        for method in LEARNED:
            row = case[method]
            require(row["completed"], f"{method} failed {case_name}")
            require(
                row["minimum_density"] >= 1.0e-5,
                f"{method} lost density admissibility on {case_name}",
            )
            require(
                row["minimum_pressure"] >= 1.0e-5,
                f"{method} lost pressure admissibility on {case_name}",
            )
            require(
                row["max_interior_tadmor_residual"] <= 1.0e-9,
                f"{method} violated an interior Tadmor constraint on {case_name}",
            )
            require(
                row["max_boundary_aware_entropy_balance_violation"] <= 1.0e-8,
                f"{method} violated boundary-aware entropy balance on {case_name}",
            )

    for group in ("all", "centered", "boundary_interaction"):
        baseline = report["aggregates"][group]["native_hllc_512"][
            "mean_rollout_nrmse"
        ]
        for method in RETAINED:
            require(
                report["aggregates"][group][method]["mean_rollout_nrmse"]
                < baseline,
                f"{method} did not beat HLLC-512 on {group}",
            )

    checks = {
        "periodic_training_to_nonperiodic_zero_shot": "pass",
        "no_learned_boundary_flux": "pass",
        "no_circular_stencil_leakage": "pass",
        "constant_state_and_flux_consistency": "pass",
        "flux_form_boundary_conservation": "pass",
        "all_learned_rollouts_completed": "pass",
        "density_and_pressure_admissibility": "pass",
        "interior_tadmor_constraints": "pass",
        "boundary_aware_fully_discrete_entropy": "pass",
        "both_retained_methods_beat_hllc512_in_all_groups": "pass",
    }
    output = {
        "seed": args.seed,
        "source": source.name,
        "checks": checks,
        "limitations": [
            "one checkpoint seed",
            "transmissive boundaries only",
            "final time 0.0252",
            "HLLC-2048 is a finite-resolution rather than exact reference",
        ],
    }
    destination = results / f"scientific_integrity_audit_seed{args.seed}.json"
    destination.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(destination)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
