"""Independent mechanical integrity audit for the stencil-ablation outputs."""

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
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import run_stencil_ablation as experiment  # noqa: E402


base = experiment.base
convergence = experiment.convergence
nonperiodic = experiment.nonperiodic


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


@torch.no_grad()
def proposal_activity(
    model: base.Solver,
    validation_data: torch.Tensor,
    batch_size: int = 128,
) -> dict[str, float | int]:
    states = validation_data.reshape(-1, validation_data.shape[-2], 3)
    proposal_squared = 0.0
    baseline_squared = 0.0
    coefficient_squared = 0.0
    coefficient_absolute = 0.0
    coefficient_count = 0
    raw_violations = 0
    interfaces = 0
    projection_active = 0
    projection_squared = 0.0
    maximum_projected_residual = -float("inf")

    for start in range(0, states.shape[0], batch_size):
        state = states[start : start + batch_size]
        flux_net = model.flux_net
        primitive = base.primitive(state)
        features = torch.cat(
            [
                (torch.roll(primitive, shift, dims=-2) - flux_net.mean)
                / flux_net.std
                for shift in flux_net.stencil_shifts
            ],
            dim=-1,
        )
        coefficients = torch.tanh(flux_net.net(features))
        raw = flux_net(state)
        projected = model.flux(state)

        if flux_net.__class__.__name__ == "DissipationFlux":
            baseline = base.t_hllc(state)
        elif flux_net.__class__.__name__ == "CentralRoeUpwindFlux":
            right_state = torch.roll(state, -1, dims=-2)
            central_flux = 0.5 * (
                base.t_flux(state) + base.t_flux(right_state)
            )
            roe_vectors, wave_strengths, speeds = base.entropy_fixed_roe_waves(
                state
            )
            dissipation = torch.einsum(
                "...ij,...j->...i",
                roe_vectors,
                speeds * wave_strengths,
            )
            baseline = central_flux - 0.5 * dissipation
        else:  # pragma: no cover - fixed preregistered arm list
            raise TypeError(type(flux_net))

        delta = raw - baseline
        projection_delta = projected - raw
        tolerance = 1.0e-7 * (1.0 + raw.abs().amax(dim=-1))
        raw_residual = base.entropy_residual(raw.double(), state.double())
        projected_residual = base.entropy_residual(
            projected.double(), state.double()
        )

        proposal_squared += float(delta.double().square().sum())
        baseline_squared += float(baseline.double().square().sum())
        coefficient_squared += float(coefficients.double().square().sum())
        coefficient_absolute += float(coefficients.double().abs().sum())
        coefficient_count += coefficients.numel()
        raw_violations += int((raw_residual > 0.0).sum())
        interfaces += raw_residual.numel()
        projection_active += int(
            (
                torch.linalg.vector_norm(projection_delta, dim=-1)
                > tolerance
            ).sum()
        )
        projection_squared += float(
            projection_delta.double().square().sum()
        )
        maximum_projected_residual = max(
            maximum_projected_residual, float(projected_residual.max())
        )

    return {
        "evaluated_states": int(states.shape[0]),
        "evaluated_interfaces": interfaces,
        "proposal_change_relative_rms": float(
            np.sqrt(proposal_squared / max(baseline_squared, 1.0e-30))
        ),
        "learned_coefficient_rms": float(
            np.sqrt(coefficient_squared / max(coefficient_count, 1))
        ),
        "learned_coefficient_mean_absolute": (
            coefficient_absolute / max(coefficient_count, 1)
        ),
        "raw_entropy_violation_rate": raw_violations / max(interfaces, 1),
        "hard_projection_intervention_rate": (
            projection_active / max(interfaces, 1)
        ),
        "hard_projection_relative_rms": float(
            np.sqrt(projection_squared / max(baseline_squared, 1.0e-30))
        ),
        "maximum_projected_tadmor_residual_float64_check": (
            maximum_projected_residual
        ),
    }


@torch.no_grad()
def nonperiodic_index_alignment(
    model: base.Solver,
    state: torch.Tensor,
) -> dict[str, float | int]:
    """Match padded nonperiodic features to periodic indexing away from edges.

    The feature comparison is bitwise exact.  Raw network fluxes are also
    reported, but their matrix multiplications can differ by a few float32
    ulps when the surrounding tensor shape/stride changes.
    """
    shifts = tuple(model.flux_net.stencil_shifts)
    halo = max(abs(shift) for shift in shifts)
    padded = torch.cat(
        [
            state[:, :1].expand(-1, halo, -1),
            state,
            state[:, -1:].expand(-1, halo, -1),
        ],
        dim=1,
    )
    padded_raw = model.flux_net(padded)
    extracted = padded_raw[:, halo : halo + state.shape[1] - 1]
    periodic_raw = model.flux_net(state)[:, :-1]

    padded_primitive = base.primitive(padded)
    periodic_primitive = base.primitive(state)
    padded_features = torch.cat(
        [
            (torch.roll(padded_primitive, shift, dims=-2)
             - model.flux_net.mean)
            / model.flux_net.std
            for shift in shifts
        ],
        dim=-1,
    )[:, halo : halo + state.shape[1] - 1]
    periodic_features = torch.cat(
        [
            (torch.roll(periodic_primitive, shift, dims=-2)
             - model.flux_net.mean)
            / model.flux_net.std
            for shift in shifts
        ],
        dim=-1,
    )[:, :-1]

    left_extent = max(max(shifts), 0)
    right_extent = max(-min(shifts), 1)
    stop = state.shape[1] - right_extent
    feature_error = float(
        (
            padded_features[:, left_extent:stop]
            - periodic_features[:, left_extent:stop]
        ).abs().max()
    )
    raw_flux_error = float(
        (
            extracted[:, left_extent:stop]
            - periodic_raw[:, left_extent:stop]
        ).abs().max()
    )
    return {
        "compared_interface_count": int(stop - left_extent),
        "maximum_feature_error": feature_error,
        "maximum_raw_flux_error": raw_flux_error,
    }


def audit(args: argparse.Namespace) -> dict[str, Any]:
    results = args.results_dir.resolve()
    convergence_rows = {
        row["arm"]: row
        for row in read_csv(results / f"convergence_seed{args.seed}.csv")
    }
    periodic = json.loads(
        (results / f"periodic512_seed{args.seed}.json").read_text(
            encoding="utf-8"
        )
    )
    nonperiodic_result = json.loads(
        (results / f"nonperiodic512_seed{args.seed}.json").read_text(
            encoding="utf-8"
        )
    )

    train_data, mean, std, _ = experiment.prepare_training_statistics(args.seed)
    validation_data = convergence.make_validation_data(args.seed)
    probe = validation_data[:2, 0]
    arm_reports: dict[str, Any] = {}
    failures: list[str] = []

    for arm in experiment.ARMS:
        checkpoint = experiment.checkpoint_path(results, arm, args.seed)
        report_path = results / f"report_{arm.name}_seed{args.seed}.json"
        curve_path = results / f"training_curve_{arm.name}_seed{args.seed}.csv"
        model = experiment.load_model(
            arm, checkpoint, mean, std, args.width
        )
        curve = read_csv(curve_path)
        metrics = [float(row["validation_rollout_nrmse"]) for row in curve]
        minimum_index = int(np.argmin(metrics))
        curve_best_metric = metrics[minimum_index]
        curve_best_update = int(curve[minimum_index]["update"])
        recorded = convergence_rows[arm.name]
        recorded_metric = float(recorded["best_validation_rollout_nrmse"])
        recorded_update = int(recorded["best_update"])
        expected_parameters = (
            args.width * args.width
            + (3 * arm.stencil_size + 5) * args.width
            + 3
        )
        actual_parameters = sum(
            parameter.numel() for parameter in model.parameters()
        )
        activity = proposal_activity(model, validation_data)
        alignment = nonperiodic_index_alignment(model, probe)
        boundary = nonperiodic.boundary_operator_self_test(model)
        output_norm = float(
            torch.sqrt(
                sum(
                    parameter.detach().double().square().sum()
                    for parameter in model.flux_net.net[-1].parameters()
                )
            )
        )

        checks = {
            "converged": recorded["converged"] == "True",
            "best_checkpoint_matches_curve": (
                abs(curve_best_metric - recorded_metric) <= 1.0e-14
                and curve_best_update == recorded_update
            ),
            "parameter_count_matches_design": (
                actual_parameters == expected_parameters
            ),
            "learned_output_is_nonzero": (
                output_norm > 1.0e-8
                and activity["proposal_change_relative_rms"] > 1.0e-6
            ),
            "nonperiodic_feature_index_alignment_exact": (
                alignment["maximum_feature_error"] == 0.0
            ),
            "nonperiodic_raw_flux_alignment_within_float32_tolerance": (
                alignment["maximum_raw_flux_error"] <= 5.0e-7
            ),
            "boundary_operator_passed": boundary["status"] == "pass",
            "projected_ground_truth_states_satisfy_tadmor": (
                activity[
                    "maximum_projected_tadmor_residual_float64_check"
                ]
                <= 1.0e-4
            ),
        }
        for check_name, passed in checks.items():
            if not passed:
                failures.append(f"{arm.name}: {check_name}")
        arm_reports[arm.name] = {
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": sha256(checkpoint),
            "stencil_size": arm.stencil_size,
            "stencil_shifts": list(model.flux_net.stencil_shifts),
            "expected_parameter_count": expected_parameters,
            "actual_parameter_count": actual_parameters,
            "best_validation_metric_from_curve": curve_best_metric,
            "best_update_from_curve": curve_best_update,
            "final_layer_parameter_l2": output_norm,
            "proposal_activity_on_independent_validation_states": activity,
            "nonperiodic_index_alignment": alignment,
            "boundary_operator": boundary,
            "checks": checks,
        }

    periodic_safety: dict[str, Any] = {}
    for arm in experiment.ARMS:
        rows = [case[arm.name] for case in periodic["cases"].values()]
        check = {
            "minimum_density": min(row["minimum_density"] for row in rows),
            "minimum_pressure": min(row["minimum_pressure"] for row in rows),
            "maximum_entropy_violation_rate": max(
                row["entropy_violation_rate"] for row in rows
            ),
            "maximum_total_entropy_change": max(
                row["max_total_entropy_change"] for row in rows
            ),
        }
        check["passed"] = (
            check["minimum_density"] >= experiment.shared.RHO_FLOOR
            and check["minimum_pressure"] >= experiment.shared.PRESSURE_FLOOR
            and check["maximum_entropy_violation_rate"] == 0.0
            and check["maximum_total_entropy_change"]
            <= experiment.shared.FD_ENTROPY_TOLERANCE
        )
        if not check["passed"]:
            failures.append(f"{arm.name}: periodic rollout safety")
        periodic_safety[arm.name] = check

    nonperiodic_safety: dict[str, Any] = {}
    for arm in experiment.ARMS:
        rows = [
            case[arm.name] for case in nonperiodic_result["cases"].values()
        ]
        check = {
            "minimum_density": min(row["minimum_density"] for row in rows),
            "minimum_pressure": min(row["minimum_pressure"] for row in rows),
            "maximum_boundary_aware_entropy_balance_violation": max(
                row["max_boundary_aware_entropy_balance_violation"]
                for row in rows
            ),
            "maximum_interior_tadmor_residual": max(
                row["max_interior_tadmor_residual"] for row in rows
            ),
        }
        check["passed"] = (
            check["minimum_density"] >= experiment.shared.RHO_FLOOR
            and check["minimum_pressure"] >= experiment.shared.PRESSURE_FLOOR
            and check[
                "maximum_boundary_aware_entropy_balance_violation"
            ]
            <= experiment.shared.FD_ENTROPY_TOLERANCE
            and check["maximum_interior_tadmor_residual"] <= 1.0e-8
        )
        if not check["passed"]:
            failures.append(f"{arm.name}: nonperiodic rollout safety")
        nonperiodic_safety[arm.name] = check

    result = {
        "status": "pass" if not failures else "fail",
        "seed": args.seed,
        "failures": failures,
        "data_shapes": {
            "training": list(train_data.shape),
            "validation": list(validation_data.shape),
        },
        "arm_audits": arm_reports,
        "periodic_rollout_safety": periodic_safety,
        "nonperiodic_rollout_safety": nonperiodic_safety,
    }
    output = results / f"scientific_integrity_audit_seed{args.seed}.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    return parser.parse_args()


if __name__ == "__main__":
    audited = audit(parse_args())
    print(json.dumps({
        "status": audited["status"],
        "failures": audited["failures"],
    }, indent=2))
    if audited["status"] != "pass":
        raise SystemExit(1)
