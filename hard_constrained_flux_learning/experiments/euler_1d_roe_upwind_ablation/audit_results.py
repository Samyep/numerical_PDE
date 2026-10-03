"""Recompute selection and structural checks for the Roe/upwind ablation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
AUDIT_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_convergence_audit"
CENTRAL_EXPERIMENT = (
    HCFL_ROOT / "experiments" / "euler_1d_consistent_central_flux"
)
for module_path in (AUDIT_EXPERIMENT, CENTRAL_EXPERIMENT, HERE):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import run_consistent_central_flux as central_baseline  # noqa: E402
import run_convergence_audit as convergence_audit  # noqa: E402
import run_roe_upwind_ablation as experiment  # noqa: E402


base = convergence_audit.base
shared = convergence_audit.shared


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_arm(
    arm: str,
    model_name: str,
    results: Path,
    seed: int,
    validation_data: torch.Tensor,
    state_std: torch.Tensor,
) -> dict[str, Any]:
    curve = read_csv(results / f"training_curve_{arm}_seed{seed}.csv")
    convergence = read_csv(
        results / f"convergence_{arm}_seed{seed}.csv"
    )[0]
    best_row = min(
        curve, key=lambda row: float(row["validation_rollout_nrmse"])
    )
    expected_update = int(best_row["update"])
    recorded_update = int(convergence["best_update"])
    expected_metric = float(best_row["validation_rollout_nrmse"])
    recorded_metric = float(convergence["best_validation_rollout_nrmse"])
    if expected_update != recorded_update or expected_metric != recorded_metric:
        raise AssertionError(f"{arm} was not selected at the curve minimum")
    if convergence["converged"].lower() != "true":
        raise AssertionError(f"{arm} was retained without convergence")

    checkpoint = results / f"{arm}_converged_best_seed{seed}.pt"
    state = torch.load(checkpoint, map_location="cpu", weights_only=True)
    model = base.Solver(model_name, np.zeros(3), np.ones(3), width=72)
    model.load_state_dict(state)
    model.eval()
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    if parameter_count != 6627:
        raise AssertionError(f"{arm} parameter count changed: {parameter_count}")
    recomputed = convergence_audit.validation_rollout_nrmse(
        model, validation_data, state_std
    )
    if abs(recomputed - recorded_metric) > 1.0e-12:
        raise AssertionError(
            f"{arm} checkpoint metric mismatch: {recomputed} vs "
            f"{recorded_metric}"
        )

    consistency = central_baseline.equal_interface_consistency(
        model, validation_data
    )
    multipliers = experiment.multiplier_diagnostics(
        model, model_name, validation_data
    )
    feasibility = experiment.proposal_feasibility_diagnostics(
        model, validation_data
    )
    if model_name.startswith("central_roe_"):
        if consistency["maximum_raw_absolute_error"] != 0.0:
            raise AssertionError(f"{arm} lost raw consistency")
        if consistency["maximum_projected_absolute_error"] != 0.0:
            raise AssertionError(f"{arm} lost projected consistency")
    if model_name == "central_roe_upwind":
        if multipliers is None or multipliers[
            "negative_multiplier_fraction"
        ] != 0.0:
            raise AssertionError("Automatic-upwind arm has a negative multiplier")

    return {
        "model": model_name,
        "parameter_count": parameter_count,
        "selected_update": recorded_update,
        "stop_update": int(convergence["stop_update"]),
        "recorded_validation_rollout_nrmse": recorded_metric,
        "recomputed_validation_rollout_nrmse": recomputed,
        "absolute_metric_difference": abs(recomputed - recorded_metric),
        "checkpoint_sha256": sha256(checkpoint),
        "equal_interface_consistency": consistency,
        "wave_multiplier_diagnostics": multipliers,
        "raw_proposal_feasibility": feasibility,
    }


def audit_512(results: Path, seed: int) -> dict[str, Any]:
    report = json.loads(
        (
            results / f"roe_upwind512_with_hllc2048_seed{seed}.json"
        ).read_text(encoding="utf-8")
    )
    new_methods = (
        "roe_complete_hcfl_512",
        "central_roe_signed_hcfl_512",
        "central_roe_upwind_hcfl_512",
    )
    maximum_conservation_drift = 0.0
    minimum_density = float("inf")
    minimum_pressure = float("inf")
    for case in report["cases"].values():
        for method in new_methods:
            metrics = case[method]
            if not metrics["completed"]:
                raise AssertionError(f"{method} has an incomplete 512 case")
            minimum_density = min(
                minimum_density, float(metrics["minimum_density"])
            )
            minimum_pressure = min(
                minimum_pressure, float(metrics["minimum_pressure"])
            )
            maximum_conservation_drift = max(
                maximum_conservation_drift,
                *(float(value) for value in metrics[
                    "max_conservation_drift"
                ].values()),
            )
    if minimum_density <= 0.0 or minimum_pressure <= 0.0:
        raise AssertionError("A completed 512-cell rollout is inadmissible")
    if maximum_conservation_drift > 2.0e-5:
        raise AssertionError(
            "512-cell conservation drift exceeded the float32 audit tolerance"
        )
    original = {
        "all_new_methods_completed_all_five_cases": True,
        "minimum_density": minimum_density,
        "minimum_pressure": minimum_pressure,
        "maximum_conservation_drift": maximum_conservation_drift,
        "mean_metrics": {
            method: report["mean_metrics"][method] for method in new_methods
        },
    }
    focused_report = json.loads(
        (results / f"upwind_feasibility512_seed{seed}.json").read_text(
            encoding="utf-8"
        )
    )
    focused_methods = ("control", "feas_1e4", "feas_1e3")
    focused_maximum_conservation_drift = 0.0
    focused_minimum_density = float("inf")
    focused_minimum_pressure = float("inf")
    for case in focused_report["cases"].values():
        for method in focused_methods:
            metrics = case[method]
            if not metrics["completed"]:
                raise AssertionError(
                    f"Focused feasibility arm {method} has an incomplete case"
                )
            focused_minimum_density = min(
                focused_minimum_density, float(metrics["minimum_density"])
            )
            focused_minimum_pressure = min(
                focused_minimum_pressure, float(metrics["minimum_pressure"])
            )
            focused_maximum_conservation_drift = max(
                focused_maximum_conservation_drift,
                *(float(value) for value in metrics[
                    "max_conservation_drift"
                ].values()),
            )
    if focused_minimum_density <= 0.0 or focused_minimum_pressure <= 0.0:
        raise AssertionError("A focused feasibility rollout is inadmissible")
    if focused_maximum_conservation_drift > 2.0e-5:
        raise AssertionError(
            "Focused feasibility conservation drift exceeded tolerance"
        )
    original_control = report["mean_metrics"][
        "central_roe_upwind_hcfl_512"
    ]["mean_rollout_nrmse"]
    focused_control = focused_report["mean_metrics"]["control"][
        "mean_rollout_nrmse"
    ]
    if original_control != focused_control:
        raise AssertionError(
            "Focused control does not reproduce the original 512 result"
        )
    return {
        "original_ablation": original,
        "focused_nonnegative_roe_feasibility": {
            "all_three_methods_completed_all_five_cases": True,
            "minimum_density": focused_minimum_density,
            "minimum_pressure": focused_minimum_pressure,
            "maximum_conservation_drift": (
                focused_maximum_conservation_drift
            ),
            "mean_metrics": focused_report["mean_metrics"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--results-dir", type=Path, default=HERE / "results"
    )
    args = parser.parse_args()
    results = args.results_dir.resolve()

    train_data = shared.make_baseline_training_data(args.seed)
    validation_data = convergence_audit.make_validation_data(args.seed)
    state_std = train_data.std(dim=(0, 1, 2))
    arms = {
        arm: audit_arm(
            arm,
            model_name,
            results,
            args.seed,
            validation_data,
            state_std,
        )
        for arm, model_name in experiment.ARMS.items()
    }
    if arms["roe_complete_broad"]["equal_interface_consistency"][
        "maximum_raw_absolute_error"
    ] == 0.0:
        raise AssertionError(
            "Direct Roe-complete unexpectedly passed the non-anchored "
            "consistency diagnostic"
        )
    signed_negative = arms["central_roe_signed_broad"][
        "wave_multiplier_diagnostics"
    ]["negative_multiplier_fraction"]
    if signed_negative <= 0.0:
        raise AssertionError("Signed arm did not exercise its negative branch")

    report = {
        "seed": args.seed,
        "verdict": (
            "PASS for checkpoint selection, metric reproduction, exact "
            "central-Roe consistency, automatic-upwind nonnegativity, and "
            "completed admissible conservative 512-cell rollouts, including "
            "both proposal-feasibility weights"
        ),
        "checks": {
            "validation_only_checkpoint_selection": "pass",
            "all_retained_new_checkpoints_converged": "pass",
            "checkpoint_metrics_exactly_reproduced": "pass",
            "central_roe_equal_interface_consistency": "pass",
            "automatic_upwind_multiplier_nonnegative": "pass",
            "signed_negative_branch_exercised": "pass",
            "all_new_512_rollouts_completed": "pass",
            "focused_feasibility_512_rollouts_completed": "pass",
            "512_admissibility_and_conservation": "pass",
        },
        "training_tensor_shape": list(train_data.shape),
        "validation_tensor_shape": list(validation_data.shape),
        "canonical_cases_used_for_selection": False,
        "arms": arms,
        "cell512": audit_512(results, args.seed),
        "limitations": [
            "one training seed",
            "five named canonical Riemann problems",
            "HLLC-2048 is a finite-resolution rather than exact reference",
            "512-cell evaluation is zero-shot transfer from 64-cell training",
        ],
    }
    output = results / f"scientific_integrity_audit_seed{args.seed}.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(output)
    print(json.dumps(report["checks"], indent=2))


if __name__ == "__main__":
    main()
