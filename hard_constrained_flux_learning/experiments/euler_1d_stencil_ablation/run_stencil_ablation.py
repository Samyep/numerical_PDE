"""Controlled 2/4/6-cell stencil ablation for retained Euler HCFL models.

Every even stencil is centred on interface i+1/2:

* 2 cells: (i | i+1)
* 4 cells: (i-1, i | i+1, i+2)
* 6 cells: (i-2, i-1, i | i+1, i+2, i+3)

The models are trained on exactly the same periodic 64-cell trajectories and
selected by exactly the same independent validation-rollout rule.  Deployment
is audited both periodically and with unlearned transmissive boundaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
CONVERGENCE_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_convergence_audit"
NONPERIODIC_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_nonperiodic_audit"
for module_path in (CONVERGENCE_EXPERIMENT, NONPERIODIC_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import plot_hllc512_vs_hcfl512 as periodic_eval  # noqa: E402
import run_convergence_audit as convergence  # noqa: E402
import run_nonperiodic_audit as nonperiodic  # noqa: E402


base = convergence.base
shared = convergence.shared
STENCIL_SIZES = (2, 4, 6)
TARGET_CELLS = 512
REFERENCE_CELLS = 2048
FAMILY_COLORS = {
    "hllc_roe": "#0072B2",
    "central_nonnegative_feas": "#D55E00",
}
STENCIL_COLORS = {2: "#0072B2", 4: "#D55E00", 6: "#009E73"}
STENCIL_LINESTYLES = {2: "--", 4: "-", 6: "-."}


@dataclass(frozen=True)
class Arm:
    name: str
    family: str
    stencil_size: int
    model_name: str
    feasibility_weight: float


def make_arms() -> tuple[Arm, ...]:
    arms: list[Arm] = []
    for stencil_size in STENCIL_SIZES:
        arms.append(
            Arm(
                name=f"hllc_roe_s{stencil_size}",
                family="hllc_roe",
                stencil_size=stencil_size,
                model_name=f"dissipation_{stencil_size}",
                feasibility_weight=0.0,
            )
        )
        arms.append(
            Arm(
                name=f"central_nonnegative_feas_s{stencil_size}",
                family="central_nonnegative_feas",
                stencil_size=stencil_size,
                model_name=f"central_roe_upwind_{stencil_size}",
                feasibility_weight=1.0e-3,
            )
        )
    return tuple(arms)


ARMS = make_arms()
ARM_BY_NAME = {arm.name: arm for arm in ARMS}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def checkpoint_path(output: Path, arm: Arm, seed: int) -> Path:
    return output / f"{arm.name}_converged_best_seed{seed}.pt"


def load_model(
    arm: Arm,
    checkpoint: Path,
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
) -> base.Solver:
    model = base.Solver(arm.model_name, mean, std, width=width)
    state = torch.load(checkpoint, map_location="cpu", weights_only=True)
    model.load_state_dict(state)
    model.eval()
    return model


def prepare_training_statistics(
    seed: int,
) -> tuple[torch.Tensor, np.ndarray, np.ndarray, torch.Tensor]:
    train_data = shared.make_baseline_training_data(seed)
    primitive = base.primitive(train_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))
    return train_data, mean, std, state_std


def train_all(
    args: argparse.Namespace,
    output: Path,
    train_data: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    state_std: torch.Tensor,
) -> tuple[dict[str, base.Solver], list[dict[str, Any]]]:
    validation_data = convergence.make_validation_data(args.seed)
    evaluation_suite = shared.make_evaluation_suite(args.seed)
    models: dict[str, base.Solver] = {}
    convergence_rows: list[dict[str, Any]] = []
    periodic64_rows: list[dict[str, Any]] = []

    for arm in ARMS:
        checkpoint = checkpoint_path(output, arm, args.seed)
        report_path = output / f"report_{arm.name}_seed{args.seed}.json"
        if args.resume and checkpoint.exists() and report_path.exists():
            print(json.dumps({"stage": "resume_arm", "arm": arm.name}), flush=True)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            convergence_rows.append(report["convergence"])
            periodic64_rows.extend(report["periodic_64_cell_test_metrics"])
            models[arm.name] = load_model(
                arm, checkpoint, mean, std, args.width
            )
            continue

        print(json.dumps({"stage": "train_arm", "arm": arm.name}), flush=True)
        model, result, curve = convergence.train_to_convergence(
            arm=arm.name,
            model_name=arm.model_name,
            train_data=train_data,
            validation_data=validation_data,
            mean=mean,
            std=std,
            state_std=state_std,
            seed=args.seed,
            width=args.width,
            batch_size=args.batch_size,
            learning_rate=args.lr,
            max_updates=args.max_updates,
            validation_interval=args.validation_interval,
            fixed_budget_updates=args.fixed_budget_updates,
            plateau_patience=args.plateau_patience,
            minimum_relative_improvement=args.minimum_relative_improvement,
            learning_rate_factor=args.learning_rate_factor,
            minimum_learning_rate=args.minimum_learning_rate,
            output=output,
            proposal_feasibility_weight=arm.feasibility_weight,
        )
        summary = {
            key: value
            for key, value in result.items()
            if key not in ("fixed", "best")
        }
        if not bool(summary["converged"]):
            checkpoint.unlink(missing_ok=True)
            raise RuntimeError(
                f"{arm.name} reached the {args.max_updates}-update cap "
                "without satisfying the preregistered validation plateau."
            )

        metrics = convergence.evaluate_state(
            arm.model_name,
            result["best"],
            mean,
            std,
            args.width,
            state_std,
            evaluation_suite,
            args.seed,
            arm.name,
            "converged_best",
        )
        for row in metrics:
            row.update(
                {
                    "family": arm.family,
                    "stencil_size": arm.stencil_size,
                }
            )
        summary.update(
            {
                "family": arm.family,
                "stencil_size": arm.stencil_size,
                "stencil_shifts": list(model.flux_net.stencil_shifts),
            }
        )
        report = {
            "design": {
                "interface": "i+1/2",
                "stencil_size": arm.stencil_size,
                "stencil_shifts": list(model.flux_net.stencil_shifts),
                "hidden_width": args.width,
                "proposal_feasibility_weight": arm.feasibility_weight,
                "training_boundary_condition": "periodic",
                "training_trajectories": int(train_data.shape[0]),
                "training_cells": int(train_data.shape[-2]),
            },
            "convergence": summary,
            "periodic_64_cell_test_metrics": metrics,
        }
        write_csv(
            output / f"training_curve_{arm.name}_seed{args.seed}.csv",
            curve,
        )
        write_csv(
            output / f"periodic64_metrics_{arm.name}_seed{args.seed}.csv",
            metrics,
        )
        report_path.write_text(
            json.dumps(json_ready(report), indent=2), encoding="utf-8"
        )
        models[arm.name] = model
        convergence_rows.append(summary)
        periodic64_rows.extend(metrics)

    write_csv(output / f"convergence_seed{args.seed}.csv", convergence_rows)
    write_csv(output / f"periodic64_metrics_seed{args.seed}.csv", periodic64_rows)
    return models, convergence_rows


def load_all_models(
    args: argparse.Namespace,
    output: Path,
    mean: np.ndarray,
    std: np.ndarray,
) -> tuple[dict[str, base.Solver], list[dict[str, Any]]]:
    models: dict[str, base.Solver] = {}
    rows: list[dict[str, Any]] = []
    for arm in ARMS:
        checkpoint = checkpoint_path(output, arm, args.seed)
        report_path = output / f"report_{arm.name}_seed{args.seed}.json"
        if not checkpoint.exists() or not report_path.exists():
            raise FileNotFoundError(
                f"Missing converged artifacts for {arm.name}; run the "
                "training phase first."
            )
        models[arm.name] = load_model(
            arm, checkpoint, mean, std, args.width
        )
        rows.append(
            json.loads(report_path.read_text(encoding="utf-8"))["convergence"]
        )
    return models, rows


def total_variation(values: np.ndarray, periodic: bool) -> float:
    differences = np.diff(values)
    if periodic:
        differences = np.concatenate([differences, values[:1] - values[-1:]])
    return float(np.abs(differences).sum())


def significant_extrema(
    values: np.ndarray,
    scale: float,
    periodic: bool,
) -> int:
    threshold = 1.0e-3 * max(scale, 1.0e-8)
    if periodic:
        left = values - np.roll(values, 1)
        right = np.roll(values, -1) - values
    else:
        left = values[1:-1] - values[:-2]
        right = values[2:] - values[1:-1]
    return int(
        (
            (left * right < 0.0)
            & (np.minimum(np.abs(left), np.abs(right)) > threshold)
        ).sum()
    )


def oscillation_metrics(
    reference: torch.Tensor,
    candidate: torch.Tensor,
    periodic: bool,
) -> dict[str, float | int]:
    truth = base.primitive(reference)[0, -1].numpy()
    prediction = base.primitive(candidate)[0, -1].numpy()
    tv_excesses: list[float] = []
    range_violations: list[float] = []
    excess_extrema = 0
    for variable in range(3):
        exact = truth[:, variable]
        estimated = prediction[:, variable]
        scale = max(float(np.ptp(exact)), float(np.max(np.abs(exact))), 1.0)
        tv_excesses.append(
            max(
                total_variation(estimated, periodic)
                - total_variation(exact, periodic),
                0.0,
            )
            / scale
        )
        range_violations.append(
            (
                max(float(exact.min() - estimated.min()), 0.0)
                + max(float(estimated.max() - exact.max()), 0.0)
            )
            / scale
        )
        excess_extrema += max(
            significant_extrema(estimated, scale, periodic)
            - significant_extrema(exact, scale, periodic),
            0,
        )
    return {
        "mean_normalized_final_tv_excess": float(np.mean(tv_excesses)),
        "mean_normalized_global_range_violation": float(
            np.mean(range_violations)
        ),
        "total_excess_significant_extrema": excess_extrema,
    }


def aggregate_case_metrics(
    cases: dict[str, dict[str, Any]],
    methods: tuple[str, ...],
    selected_cases: tuple[str, ...] | None = None,
) -> dict[str, dict[str, Any]]:
    names = selected_cases or tuple(cases)
    aggregate: dict[str, dict[str, Any]] = {}
    for method in methods:
        rows = [cases[name][method] for name in names]
        aggregate[method] = {
            "case_count": len(rows),
            "mean_rollout_nrmse": float(
                np.mean([row["rollout_nrmse"] for row in rows])
            ),
            "mean_final_snapshot_nrmse": float(
                np.mean([row["final_snapshot_nrmse"] for row in rows])
            ),
            "mean_normalized_final_tv_excess": float(
                np.mean([
                    row["mean_normalized_final_tv_excess"] for row in rows
                ])
            ),
            "mean_normalized_global_range_violation": float(
                np.mean([
                    row["mean_normalized_global_range_violation"]
                    for row in rows
                ])
            ),
            "total_excess_significant_extrema": int(
                sum(row["total_excess_significant_extrema"] for row in rows)
            ),
            "minimum_density": float(
                min(row["minimum_density"] for row in rows)
            ),
            "minimum_pressure": float(
                min(row["minimum_pressure"] for row in rows)
            ),
        }
        for key in (
            "hard_projection_intervention_rate",
            "local_limiter_intervention_rate",
            "fd_entropy_intervention_rate",
        ):
            if all(key in row for row in rows):
                aggregate[method][f"mean_{key}"] = float(
                    np.mean([row[key] for row in rows])
                )
    return aggregate


def evaluate_periodic512(
    args: argparse.Namespace,
    models: dict[str, base.Solver],
    state_std: torch.Tensor,
) -> tuple[
    dict[str, Any],
    dict[str, torch.Tensor],
    dict[str, torch.Tensor],
    dict[str, dict[str, torch.Tensor]],
]:
    reference_native: dict[str, torch.Tensor] = {}
    native_512: dict[str, torch.Tensor] = {}
    predictions: dict[str, dict[str, torch.Tensor]] = {
        arm.name: {} for arm in ARMS
    }
    cases: dict[str, dict[str, Any]] = {}
    methods = ("native_hllc_512", *(arm.name for arm in ARMS))

    for name in periodic_eval.comparison.CASES:
        print(
            json.dumps({"stage": "periodic512", "case": name}),
            flush=True,
        )
        reference_native[name], reference_stats = (
            periodic_eval.strict_native_hllc_rollout(name, REFERENCE_CELLS)
        )
        restricted = periodic_eval.comparison.precision.conservative_restrict(
            reference_native[name].numpy(), TARGET_CELLS
        )
        scoring_reference = torch.from_numpy(restricted.astype(np.float32))
        native_512[name], native_stats = (
            periodic_eval.strict_native_hllc_rollout(name, TARGET_CELLS)
        )
        case: dict[str, Any] = {
            "reference_hllc_2048": reference_stats,
            "native_hllc_512": periodic_eval.enrich_metrics(
                scoring_reference, native_512[name], state_std, native_stats
            ),
        }
        case["native_hllc_512"].update(
            oscillation_metrics(scoring_reference, native_512[name], True)
        )
        for arm in ARMS:
            trajectory, run_stats = periodic_eval.matched_safe_rollout(
                models[arm.name], name, TARGET_CELLS
            )
            predictions[arm.name][name] = trajectory
            case[arm.name] = periodic_eval.enrich_metrics(
                scoring_reference, trajectory, state_std, run_stats
            )
            case[arm.name].update(
                oscillation_metrics(scoring_reference, trajectory, True)
            )
        cases[name] = case

    return (
        {
            "boundary_condition": "periodic",
            "reference": (
                "native HLLC + SSP-RK2 on 2048 cells, conservatively "
                "restricted to 512 cells for scoring"
            ),
            "cases": cases,
            "aggregate": aggregate_case_metrics(cases, methods),
        },
        reference_native,
        native_512,
        predictions,
    )


def evaluate_nonperiodic512(
    args: argparse.Namespace,
    models: dict[str, base.Solver],
    state_std: torch.Tensor,
) -> tuple[
    dict[str, Any],
    dict[str, torch.Tensor],
    dict[str, torch.Tensor],
    dict[str, dict[str, torch.Tensor]],
]:
    reference_native: dict[str, torch.Tensor] = {}
    native_512: dict[str, torch.Tensor] = {}
    predictions: dict[str, dict[str, torch.Tensor]] = {
        arm.name: {} for arm in ARMS
    }
    cases: dict[str, dict[str, Any]] = {}
    methods = ("native_hllc_512", *(arm.name for arm in ARMS))
    operator_tests = {
        arm.name: nonperiodic.boundary_operator_self_test(models[arm.name])
        for arm in ARMS
    }

    for name in nonperiodic.CASES:
        print(
            json.dumps({"stage": "nonperiodic512", "case": name}),
            flush=True,
        )
        reference_native[name], reference_stats = nonperiodic.strict_hllc_rollout(
            name, REFERENCE_CELLS
        )
        scoring_reference = nonperiodic.restrict_reference(
            reference_native[name], TARGET_CELLS
        )
        native_512[name], native_stats = nonperiodic.strict_hllc_rollout(
            name, TARGET_CELLS
        )
        case: dict[str, Any] = {
            "reference_hllc_2048": reference_stats,
            "native_hllc_512": dict(native_stats),
        }
        case["native_hllc_512"].update(
            periodic_eval.comparison.diagnostics(
                scoring_reference, native_512[name], state_std
            )
        )
        case["native_hllc_512"].update(
            nonperiodic.integrity_metrics(native_512[name])
        )
        case["native_hllc_512"].update(
            oscillation_metrics(scoring_reference, native_512[name], False)
        )

        for arm in ARMS:
            trajectory, run_stats = nonperiodic.learned_rollout(
                models[arm.name], name, TARGET_CELLS
            )
            predictions[arm.name][name] = trajectory
            case[arm.name] = dict(run_stats)
            case[arm.name].update(
                periodic_eval.comparison.diagnostics(
                    scoring_reference, trajectory, state_std
                )
            )
            case[arm.name].update(nonperiodic.integrity_metrics(trajectory))
            case[arm.name].update(
                oscillation_metrics(scoring_reference, trajectory, False)
            )
        cases[name] = case

    centered = tuple(
        name
        for name, spec in nonperiodic.CASES.items()
        if spec["group"] == "centered"
    )
    boundary = tuple(
        name
        for name, spec in nonperiodic.CASES.items()
        if spec["group"] == "boundary_interaction"
    )
    return (
        {
            "training_boundary_condition": "periodic",
            "deployment_boundary_condition": (
                "transmissive constant extrapolation; physical Euler flux "
                "at both domain boundaries; no learned boundary flux"
            ),
            "boundary_operator_self_tests": operator_tests,
            "reference": (
                "strict nonperiodic HLLC + SSP-RK2 on 2048 cells, "
                "conservatively restricted to 512 cells for scoring"
            ),
            "cases": cases,
            "aggregates": {
                "all": aggregate_case_metrics(cases, methods),
                "centered": aggregate_case_metrics(cases, methods, centered),
                "boundary_interaction": aggregate_case_metrics(
                    cases, methods, boundary
                ),
            },
        },
        reference_native,
        native_512,
        predictions,
    )


def flatten_case_rows(
    boundary: str,
    cases: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    methods = ("native_hllc_512", *(arm.name for arm in ARMS))
    for case_name, case in cases.items():
        for method in methods:
            metrics = case[method]
            arm = ARM_BY_NAME.get(method)
            rows.append(
                {
                    "boundary_condition": boundary,
                    "case": case_name,
                    "method": method,
                    "family": arm.family if arm else "native_hllc",
                    "stencil_size": arm.stencil_size if arm else None,
                    "rollout_nrmse": metrics["rollout_nrmse"],
                    "final_snapshot_nrmse": metrics["final_snapshot_nrmse"],
                    "mean_normalized_final_tv_excess": metrics[
                        "mean_normalized_final_tv_excess"
                    ],
                    "mean_normalized_global_range_violation": metrics[
                        "mean_normalized_global_range_violation"
                    ],
                    "total_excess_significant_extrema": metrics[
                        "total_excess_significant_extrema"
                    ],
                    "minimum_density": metrics["minimum_density"],
                    "minimum_pressure": metrics["minimum_pressure"],
                    "hard_projection_intervention_rate": metrics.get(
                        "hard_projection_intervention_rate"
                    ),
                    "local_limiter_intervention_rate": metrics.get(
                        "local_limiter_intervention_rate"
                    ),
                    "fd_entropy_intervention_rate": metrics.get(
                        "fd_entropy_intervention_rate"
                    ),
                }
            )
    return rows


def plot_family_profiles(
    boundary: str,
    family: str,
    case_names: tuple[str, ...],
    display_names: dict[str, str],
    references: dict[str, torch.Tensor],
    native_512: dict[str, torch.Tensor],
    predictions: dict[str, dict[str, torch.Tensor]],
    output: Path,
) -> None:
    columns = len(case_names)
    width = 18.0 if columns <= 5 else 24.0
    figure, axes = plt.subplots(
        3, columns, figsize=(width, 8.7), sharex=True, squeeze=False
    )
    figure.subplots_adjust(
        left=0.05, right=0.995, bottom=0.09, top=0.82,
        wspace=0.22, hspace=0.14,
    )
    x_reference = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    x_target = (np.arange(TARGET_CELLS) + 0.5) / TARGET_CELLS
    variables = (r"Density $\rho$", r"Velocity $u$", r"Pressure $p$")

    for column, name in enumerate(case_names):
        reference_primitive = base.primitive(references[name])[0, -1].numpy()
        native_primitive = base.primitive(native_512[name])[0, -1].numpy()
        learned = {
            arm.stencil_size: base.primitive(predictions[arm.name][name])[
                0, -1
            ].numpy()
            for arm in ARMS
            if arm.family == family
        }
        axes[0, column].set_title(
            display_names[name], fontsize=9.5, fontweight="semibold"
        )
        for row in range(3):
            axis = axes[row, column]
            axis.plot(
                x_reference,
                reference_primitive[:, row],
                color="#CBD5E1",
                linewidth=2.0,
                alpha=0.9,
                zorder=1,
            )
            axis.plot(
                x_target,
                native_primitive[:, row],
                color="#374151",
                linewidth=1.0,
                linestyle=":",
                zorder=2,
            )
            for stencil_size in STENCIL_SIZES:
                axis.plot(
                    x_target,
                    learned[stencil_size][:, row],
                    color=STENCIL_COLORS[stencil_size],
                    linewidth=1.25 if stencil_size != 4 else 1.55,
                    linestyle=STENCIL_LINESTYLES[stencil_size],
                    zorder=2 + stencil_size,
                )
            periodic_eval.comparison.style_axis(axis)
            axis.set_xlim(0.0, 1.0)
            if column == 0:
                axis.set_ylabel(variables[row], fontsize=9.5)
            if row == 2:
                axis.set_xlabel("x", fontsize=8.5)

    family_label = (
        "HLLC + Roe correction"
        if family == "hllc_roe"
        else r"Central + nonnegative Roe + feasibility ($\lambda=10^{-3}$)"
    )
    handles = [
        Line2D([0], [0], color="#CBD5E1", linewidth=2.0, label="HLLC-2048"),
        Line2D(
            [0], [0], color="#374151", linewidth=1.0,
            linestyle=":", label="native HLLC-512",
        ),
    ]
    handles.extend(
        Line2D(
            [0], [0], color=STENCIL_COLORS[size],
            linestyle=STENCIL_LINESTYLES[size], linewidth=1.5,
            label=f"{size}-cell NN stencil",
        )
        for size in STENCIL_SIZES
    )
    figure.suptitle(
        f"{family_label}: {boundary} deployment, t = 0.0252",
        y=0.975,
        fontsize=15.5,
        fontweight="bold",
    )
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.52, 0.91),
        ncol=5,
        frameon=False,
        fontsize=8.8,
    )
    figure.savefig(output, dpi=210, facecolor="white")
    plt.close(figure)


def build_summary_rows(
    convergence_rows: list[dict[str, Any]],
    periodic: dict[str, Any],
    nonperiodic_result: dict[str, Any],
) -> list[dict[str, Any]]:
    convergence_by_arm = {row["arm"]: row for row in convergence_rows}
    rows: list[dict[str, Any]] = []
    for arm in ARMS:
        conv = convergence_by_arm[arm.name]
        periodic_metrics = periodic["aggregate"][arm.name]
        nonperiodic_all = nonperiodic_result["aggregates"]["all"][arm.name]
        nonperiodic_centered = nonperiodic_result["aggregates"]["centered"][
            arm.name
        ]
        boundary = nonperiodic_result["aggregates"]["boundary_interaction"][
            arm.name
        ]
        periodic_native = periodic["aggregate"]["native_hllc_512"]
        nonperiodic_native = nonperiodic_result["aggregates"]["all"][
            "native_hllc_512"
        ]
        rows.append(
            {
                "arm": arm.name,
                "family": arm.family,
                "stencil_size": arm.stencil_size,
                "stencil_shifts": " ".join(
                    str(value)
                    for value in base.interface_stencil_shifts(
                        arm.stencil_size
                    )
                ),
                "parameter_count": conv["parameter_count"],
                "best_update": conv["best_update"],
                "stop_update": conv["stop_update"],
                "validation_rollout_nrmse": conv[
                    "best_validation_rollout_nrmse"
                ],
                "periodic512_mean_rollout_nrmse": periodic_metrics[
                    "mean_rollout_nrmse"
                ],
                "periodic512_mean_final_nrmse": periodic_metrics[
                    "mean_final_snapshot_nrmse"
                ],
                "periodic512_reduction_vs_native_hllc_percent": 100.0
                * (
                    1.0
                    - periodic_metrics["mean_rollout_nrmse"]
                    / periodic_native["mean_rollout_nrmse"]
                ),
                "periodic512_tv_excess": periodic_metrics[
                    "mean_normalized_final_tv_excess"
                ],
                "periodic512_excess_extrema": periodic_metrics[
                    "total_excess_significant_extrema"
                ],
                "nonperiodic512_mean_rollout_nrmse": nonperiodic_all[
                    "mean_rollout_nrmse"
                ],
                "nonperiodic512_mean_final_nrmse": nonperiodic_all[
                    "mean_final_snapshot_nrmse"
                ],
                "nonperiodic512_reduction_vs_native_hllc_percent": 100.0
                * (
                    1.0
                    - nonperiodic_all["mean_rollout_nrmse"]
                    / nonperiodic_native["mean_rollout_nrmse"]
                ),
                "nonperiodic512_centered_rollout_nrmse": (
                    nonperiodic_centered["mean_rollout_nrmse"]
                ),
                "nonperiodic512_boundary_rollout_nrmse": boundary[
                    "mean_rollout_nrmse"
                ],
                "nonperiodic512_tv_excess": nonperiodic_all[
                    "mean_normalized_final_tv_excess"
                ],
                "nonperiodic512_excess_extrema": nonperiodic_all[
                    "total_excess_significant_extrema"
                ],
                "nonperiodic512_centered_tv_excess": nonperiodic_centered[
                    "mean_normalized_final_tv_excess"
                ],
                "nonperiodic512_centered_excess_extrema": (
                    nonperiodic_centered["total_excess_significant_extrema"]
                ),
                "nonperiodic512_minimum_pressure": nonperiodic_all[
                    "minimum_pressure"
                ],
            }
        )
    return rows


def plot_summary(rows: list[dict[str, Any]], output: Path) -> None:
    by_key = {
        (row["family"], int(row["stencil_size"])): row for row in rows
    }
    panels = (
        ("validation_rollout_nrmse", "Periodic validation (64 cells)"),
        ("periodic512_mean_rollout_nrmse", "Periodic canonical (512 cells)"),
        ("nonperiodic512_mean_rollout_nrmse", "Nonperiodic all (512 cells)"),
        (
            "nonperiodic512_boundary_rollout_nrmse",
            "Nonperiodic boundary interaction (512 cells)",
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    for axis, (key, title) in zip(axes.flat, panels):
        for family in FAMILY_COLORS:
            values = [by_key[(family, size)][key] for size in STENCIL_SIZES]
            label = (
                "HLLC + Roe correction"
                if family == "hllc_roe"
                else "Central + nonnegative Roe + feasibility"
            )
            axis.plot(
                STENCIL_SIZES,
                values,
                marker="o",
                markersize=6,
                linewidth=1.8,
                color=FAMILY_COLORS[family],
                label=label,
            )
            for size, value in zip(STENCIL_SIZES, values):
                axis.annotate(
                    f"{value:.4f}",
                    (size, value),
                    xytext=(0, 7),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7.8,
                    color=FAMILY_COLORS[family],
                )
        axis.set_title(title, fontsize=10.5, fontweight="semibold")
        axis.set_xlabel("NN stencil cells")
        axis.set_ylabel("rollout NRMSE (lower is better)")
        axis.set_xticks(STENCIL_SIZES)
        axis.grid(True, color="#D7DEE7", linewidth=0.7, alpha=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    figure.suptitle(
        "Symmetric 2/4/6-cell interface-stencil ablation (seed 0)",
        y=1.025,
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(output, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_stability_summary(rows: list[dict[str, Any]], output: Path) -> None:
    by_key = {
        (row["family"], int(row["stencil_size"])): row for row in rows
    }
    panels = (
        ("periodic512_tv_excess", "Periodic final TV excess", False),
        ("periodic512_excess_extrema", "Periodic excess extrema", False),
        (
            "nonperiodic512_centered_tv_excess",
            "Nonperiodic centered final TV excess",
            False,
        ),
        (
            "nonperiodic512_minimum_pressure",
            "Nonperiodic minimum pressure",
            True,
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    for axis, (key, title, logarithmic) in zip(axes.flat, panels):
        for family in FAMILY_COLORS:
            values = [by_key[(family, size)][key] for size in STENCIL_SIZES]
            label = (
                "HLLC + Roe correction"
                if family == "hllc_roe"
                else "Central + nonnegative Roe + feasibility"
            )
            axis.plot(
                STENCIL_SIZES,
                values,
                marker="o",
                markersize=6,
                linewidth=1.8,
                color=FAMILY_COLORS[family],
                label=label,
            )
        if logarithmic:
            axis.set_yscale("log")
            axis.set_ylabel("pressure (higher margin is safer)")
        elif "extrema" in key:
            axis.set_ylabel("count (lower is better)")
        else:
            axis.set_ylabel("normalized excess (lower is better)")
        axis.set_title(title, fontsize=10.5, fontweight="semibold")
        axis.set_xlabel("NN stencil cells")
        axis.set_xticks(STENCIL_SIZES)
        axis.grid(True, color="#D7DEE7", linewidth=0.7, alpha=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    figure.suptitle(
        "Stencil stability and admissibility diagnostics (seed 0)",
        y=1.025,
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(output, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def evaluate_all(
    args: argparse.Namespace,
    output: Path,
    models: dict[str, base.Solver],
    state_std: torch.Tensor,
    convergence_rows: list[dict[str, Any]],
) -> None:
    periodic, p_reference, p_native, p_predictions = evaluate_periodic512(
        args, models, state_std
    )
    nonperiodic_result, n_reference, n_native, n_predictions = (
        evaluate_nonperiodic512(args, models, state_std)
    )
    periodic_path = output / f"periodic512_seed{args.seed}.json"
    nonperiodic_path = output / f"nonperiodic512_seed{args.seed}.json"
    periodic_path.write_text(
        json.dumps(json_ready(periodic), indent=2), encoding="utf-8"
    )
    nonperiodic_path.write_text(
        json.dumps(json_ready(nonperiodic_result), indent=2),
        encoding="utf-8",
    )
    case_rows = flatten_case_rows("periodic", periodic["cases"])
    case_rows.extend(
        flatten_case_rows("transmissive_nonperiodic", nonperiodic_result["cases"])
    )
    write_csv(output / f"deployment_metrics_seed{args.seed}.csv", case_rows)

    periodic_names = tuple(periodic_eval.comparison.CASES)
    periodic_display = dict(periodic_eval.comparison.DISPLAY)
    nonperiodic_names = tuple(nonperiodic.CASES)
    nonperiodic_display = {
        name: spec["display"] for name, spec in nonperiodic.CASES.items()
    }
    for family in FAMILY_COLORS:
        plot_family_profiles(
            "periodic",
            family,
            periodic_names,
            periodic_display,
            p_reference,
            p_native,
            p_predictions,
            output / f"periodic512_{family}_profiles_seed{args.seed}.png",
        )
        plot_family_profiles(
            "transmissive nonperiodic",
            family,
            nonperiodic_names,
            nonperiodic_display,
            n_reference,
            n_native,
            n_predictions,
            output / f"nonperiodic512_{family}_profiles_seed{args.seed}.png",
        )

    summary_rows = build_summary_rows(
        convergence_rows, periodic, nonperiodic_result
    )
    write_csv(output / f"summary_seed{args.seed}.csv", summary_rows)
    summary = {
        "scope": "single-seed controlled screening; not a multi-seed claim",
        "seed": args.seed,
        "stencil_definition": {
            str(size): list(base.interface_stencil_shifts(size))
            for size in STENCIL_SIZES
        },
        "training": (
            "identical periodic 580-trajectory, 64-cell data; independent "
            "136-trajectory validation; validation-plateau checkpointing"
        ),
        "native_hllc_512_baselines": {
            "periodic_mean_rollout_nrmse": periodic["aggregate"][
                "native_hllc_512"
            ]["mean_rollout_nrmse"],
            "nonperiodic_mean_rollout_nrmse": nonperiodic_result[
                "aggregates"
            ]["all"]["native_hllc_512"]["mean_rollout_nrmse"],
        },
        "screening_interpretation": {
            "recommended_general_stencil": 4,
            "reason": (
                "best validation and periodic-512 accuracy in both model "
                "families; best centered nonperiodic accuracy with far less "
                "oscillation than the two-cell models"
            ),
            "two_cell_boundary_finding": (
                "lowest boundary-interaction NRMSE, but severe shock-region "
                "TV/extrema elsewhere; motivates a boundary-local fallback, "
                "not a global two-cell replacement"
            ),
            "six_cell_finding": (
                "lower TV than four cells in several tests but consistently "
                "higher rollout error"
            ),
        },
        "rows": summary_rows,
    }
    (output / f"summary_seed{args.seed}.json").write_text(
        json.dumps(json_ready(summary), indent=2), encoding="utf-8"
    )
    plot_summary(
        summary_rows, output / f"stencil_ablation_summary_seed{args.seed}.png"
    )
    plot_stability_summary(
        summary_rows, output / f"stencil_stability_summary_seed{args.seed}.png"
    )
    print(json.dumps(json_ready(summary), indent=2), flush=True)


def self_test() -> None:
    expected = {
        2: (0, -1),
        4: (1, 0, -1, -2),
        6: (2, 1, 0, -1, -2, -3),
    }
    mean = np.array([1.0, 0.0, 1.0], dtype=np.float32)
    std = np.ones(3, dtype=np.float32)
    primitive = np.zeros((2, 24, 3), dtype=np.float32)
    primitive[..., 0] = 1.0
    primitive[..., 1] = np.linspace(-0.2, 0.2, 24)[None, :]
    primitive[..., 2] = 1.0
    state = torch.from_numpy(
        base.prim_to_cons(
            primitive[..., 0], primitive[..., 1], primitive[..., 2]
        ).astype(np.float32)
    )
    for arm in ARMS:
        model = base.Solver(arm.model_name, mean, std, width=16)
        if tuple(model.flux_net.stencil_shifts) != expected[arm.stencil_size]:
            raise RuntimeError(f"Incorrect stencil for {arm.name}")
        if model.flux_net.net[0].in_features != 3 * arm.stencil_size:
            raise RuntimeError(f"Incorrect feature count for {arm.name}")
        flux = model.flux(state)
        if flux.shape != state.shape or not bool(torch.isfinite(flux).all()):
            raise RuntimeError(f"Invalid periodic flux for {arm.name}")
        translated_flux = model.flux(torch.roll(state, 3, dims=-2))
        translation_error = float(
            (
                translated_flux - torch.roll(flux, 3, dims=-2)
            ).detach().abs().max()
        )
        if translation_error > 2.0e-6:
            raise RuntimeError(
                f"Periodic translation equivariance failed for {arm.name}: "
                f"{translation_error}"
            )
        boundary_test = nonperiodic.boundary_operator_self_test(model)
        if boundary_test["status"] != "pass":
            raise RuntimeError(f"Boundary test failed for {arm.name}")
    print("self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase", choices=("train", "evaluate", "all"), default="all"
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=56)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=1100)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--validation-interval", type=int, default=100)
    parser.add_argument("--plateau-patience", type=int, default=10)
    parser.add_argument("--minimum-relative-improvement", type=float, default=1.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    print("Generating shared strict training data...", flush=True)
    train_data, mean, std, state_std = prepare_training_statistics(args.seed)

    if args.phase in ("train", "all"):
        models, convergence_rows = train_all(
            args, output, train_data, mean, std, state_std
        )
    else:
        models, convergence_rows = load_all_models(
            args, output, mean, std
        )

    if args.phase in ("evaluate", "all"):
        evaluate_all(
            args, output, models, state_std, convergence_rows
        )
    print(
        json.dumps(
            {
                "stage": "complete",
                "phase": args.phase,
                "wall_seconds": time.perf_counter() - started,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
