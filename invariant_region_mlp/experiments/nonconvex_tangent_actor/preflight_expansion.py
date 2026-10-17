"""Freeze independent references for controlled lambda/dimension expansion."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json

from .common import REFERENCE_DIR, RESULTS_DIR, atomic_json, ensure_dirs
from .configuration import ExperimentSpec
from .drivers import SineQuadraticDriver
from .reference2d import save_reference, validate_reference


def prefix(spec: ExperimentSpec) -> str:
    return spec.run_id.split(f"_{spec.budget}")[0]


def run(spec: ExperimentSpec) -> dict:
    ensure_dirs()
    # References are reduced and ambient-dimension independent. Always build
    # them with d=100, while retaining the requested ambient d in metadata.
    reference_driver = SineQuadraticDriver(
        d=100, lam=spec.lam, eta=spec.eta if spec.coupled else 0.0
    )
    reference, audit = validate_reference(reference_driver)
    reference_path = REFERENCE_DIR / f"{reference_driver.name}.npz"
    save_reference(reference_path, reference)
    curvature_required = spec.lam > 0.0
    passed = bool(
        audit["inside_invariant_interval"]
        and audit["terminal_value_rmse"] < 2e-3
        and audit["fine_vs_production"]["value_rmse"] < 5e-4
        and (
            audit["curvature_coverage"]["meaningful_sign_change_gate"]
            if curvature_required
            else True
        )
    )
    payload = {
        "spec": asdict(spec),
        "reference_driver": reference_driver.metadata(),
        "reference": audit,
        "reference_path": str(reference_path),
        "curvature_sign_change_required": curvature_required,
        "lambda_zero_control": spec.lam == 0.0,
        "passed": passed,
        "ground_truth_visible_to_training": False,
    }
    output = RESULTS_DIR / f"preflight_{prefix(spec)}.json"
    atomic_json(output, payload)
    print(
        json.dumps(
            {
                "output": str(output),
                "passed": passed,
                "lambda": spec.lam,
                "eta": spec.eta,
                "reference_rmse_proxy": audit["fine_vs_production"]["value_rmse"],
                "curvature_gate": audit["curvature_coverage"][
                    "meaningful_sign_change_gate"
                ],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--lam", type=float, required=True)
    parser.add_argument("--eta", type=float, default=1.0)
    args = parser.parse_args()
    run(
        ExperimentSpec(
            family=args.family,
            d=args.d,
            lam=args.lam,
            eta=args.eta if args.family == "coupled" else 0.0,
            seed=0,
            budget="main",
            estimator="mcv",
        )
    )


if __name__ == "__main__":
    main()
