"""Validation-converged direct full-flux baseline for 1D Euler HCFL."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
AUDIT_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_convergence_audit"
if str(AUDIT_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(AUDIT_EXPERIMENT))

import run_convergence_audit as audit  # noqa: E402


shared = audit.shared
base = audit.base
ARM = "full_flux_broad"
MODEL_NAME = "full"
COMPARATORS = ("direct_broad", "dissipation_broad")
COMPONENTS = ("density_flux", "momentum_flux", "energy_flux")


def read_comparator_metrics(
    seed: int,
) -> dict[tuple[str, str], dict[str, str]]:
    path = AUDIT_EXPERIMENT / "results" / f"metrics_seed{seed}.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {
        (row["arm"], row["split"]): row
        for row in rows
        if row["arm"] in COMPARATORS
        and row["stage"] == "converged_best"
    }
    expected = {
        (arm, split)
        for arm in COMPARATORS
        for split in audit.TEST_SPLITS
    }
    missing = sorted(expected - selected.keys())
    if missing:
        raise RuntimeError(f"Missing converged comparator rows: {missing}")
    return selected


def comparison_rows(
    candidate_metrics: list[dict[str, Any]],
    seed: int,
) -> list[dict[str, Any]]:
    candidate = {row["split"]: row for row in candidate_metrics}
    comparator = read_comparator_metrics(seed)
    rows: list[dict[str, Any]] = []
    for arm in COMPARATORS:
        for split in audit.TEST_SPLITS:
            old = float(comparator[(arm, split)]["rollout_nrmse"])
            new = float(candidate[split]["rollout_nrmse"])
            rows.append(
                {
                    "seed": seed,
                    "split": split,
                    "comparator": f"{arm}_converged_best",
                    "candidate": f"{ARM}_converged_best",
                    "comparator_rollout_nrmse": old,
                    "candidate_rollout_nrmse": new,
                    "relative_change_percent": 100.0 * (new / old - 1.0),
                    "candidate_local_limiter_intervention_rate": candidate[
                        split
                    ]["local_limiter_intervention_rate"],
                    "candidate_fd_entropy_intervention_rate": candidate[
                        split
                    ]["fd_entropy_intervention_rate"],
                }
            )
    return rows


@torch.no_grad()
def proposal_diagnostics(
    model: base.Solver,
    data: torch.Tensor,
    batch_size: int = 256,
) -> dict[str, Any]:
    """Measure proposal activity and projection reliance on fixed GT states."""
    states = data.reshape(-1, data.shape[-2], data.shape[-1])
    totals = {
        "values": 0,
        "interfaces": 0,
        "raw_squared": 0.0,
        "projected_squared": 0.0,
        "hllc_squared": 0.0,
        "projection_squared": 0.0,
        "projected_hllc_error_squared": 0.0,
        "raw_divergence_squared": 0.0,
        "projected_divergence_squared": 0.0,
        "hllc_divergence_squared": 0.0,
        "projection_active": 0,
        "raw_absolute_max": 0.0,
    }
    for start in range(0, states.shape[0], batch_size):
        state = states[start : start + batch_size]
        raw = model.flux_net(state)
        projected = model.flux(state)
        hllc = base.t_hllc(state)
        delta = projected - raw
        tolerance = 1.0e-7 * (1.0 + raw.abs().amax(dim=-1))
        raw_div = raw - torch.roll(raw, 1, dims=-2)
        projected_div = projected - torch.roll(projected, 1, dims=-2)
        hllc_div = hllc - torch.roll(hllc, 1, dims=-2)

        totals["values"] += raw.numel()
        totals["interfaces"] += raw.shape[0] * raw.shape[1]
        totals["raw_squared"] += float((raw.double() ** 2).sum())
        totals["projected_squared"] += float(
            (projected.double() ** 2).sum()
        )
        totals["hllc_squared"] += float((hllc.double() ** 2).sum())
        totals["projection_squared"] += float((delta.double() ** 2).sum())
        totals["projected_hllc_error_squared"] += float(
            ((projected - hllc).double() ** 2).sum()
        )
        totals["raw_divergence_squared"] += float(
            (raw_div.double() ** 2).sum()
        )
        totals["projected_divergence_squared"] += float(
            (projected_div.double() ** 2).sum()
        )
        totals["hllc_divergence_squared"] += float(
            (hllc_div.double() ** 2).sum()
        )
        totals["projection_active"] += int(
            (torch.linalg.vector_norm(delta, dim=-1) > tolerance).sum()
        )
        totals["raw_absolute_max"] = max(
            totals["raw_absolute_max"], float(raw.abs().max())
        )

    count = max(int(totals["values"]), 1)
    hllc_squared = max(float(totals["hllc_squared"]), 1.0e-30)
    hllc_div_squared = max(
        float(totals["hllc_divergence_squared"]), 1.0e-30
    )
    raw_squared = max(float(totals["raw_squared"]), 1.0e-30)
    return {
        "evaluated_ground_truth_states": int(states.shape[0]),
        "raw_flux_rms": np.sqrt(float(totals["raw_squared"]) / count),
        "projected_flux_rms": np.sqrt(
            float(totals["projected_squared"]) / count
        ),
        "hllc_flux_rms": np.sqrt(hllc_squared / count),
        "raw_flux_absolute_max": totals["raw_absolute_max"],
        "hard_projection_intervention_rate": (
            int(totals["projection_active"])
            / max(int(totals["interfaces"]), 1)
        ),
        "hard_projection_relative_flux_rms": np.sqrt(
            float(totals["projection_squared"]) / raw_squared
        ),
        "projected_vs_hllc_relative_flux_rmse": np.sqrt(
            float(totals["projected_hllc_error_squared"]) / hllc_squared
        ),
        "raw_flux_divergence_rms": np.sqrt(
            float(totals["raw_divergence_squared"]) / count
        ),
        "projected_flux_divergence_rms": np.sqrt(
            float(totals["projected_divergence_squared"]) / count
        ),
        "hllc_flux_divergence_rms": np.sqrt(hllc_div_squared / count),
        "projected_divergence_rms_ratio_vs_hllc": np.sqrt(
            float(totals["projected_divergence_squared"])
            / hllc_div_squared
        ),
    }


@torch.no_grad()
def constant_state_consistency(
    model: base.Solver,
    train_data: torch.Tensor,
    sample_count: int = 256,
) -> dict[str, Any]:
    """Check F_hat(U,U)=F(U), which trajectory loss does not identify."""
    states = train_data.reshape(-1, 3)
    count = min(sample_count, states.shape[0])
    indices = torch.linspace(
        0, states.shape[0] - 1, count, dtype=torch.float64
    ).round().long()
    sampled = states[indices]
    constant = sampled[:, None, :].repeat(1, 8, 1)
    raw = model.flux_net(constant)
    projected = model.flux(constant)
    exact = base.t_flux(sampled)
    raw_error = raw[:, 0] - exact
    projected_error = projected[:, 0] - exact
    component_scale = torch.sqrt((exact.double() ** 2).mean(dim=0)).clamp_min(
        1.0e-12
    )
    raw_component_rmse = torch.sqrt(
        (raw_error.double() ** 2).mean(dim=0)
    )
    projected_component_rmse = torch.sqrt(
        (projected_error.double() ** 2).mean(dim=0)
    )
    projected_divergence = projected - torch.roll(projected, 1, dims=-2)
    return {
        "sample_count": count,
        "raw_relative_consistency_nrmse": float(
            torch.sqrt(
                ((raw_error.double() / component_scale) ** 2).mean()
            )
        ),
        "projected_relative_consistency_nrmse": float(
            torch.sqrt(
                ((projected_error.double() / component_scale) ** 2).mean()
            )
        ),
        "raw_component_rmse": {
            name: float(value)
            for name, value in zip(COMPONENTS, raw_component_rmse)
        },
        "projected_component_rmse": {
            name: float(value)
            for name, value in zip(COMPONENTS, projected_component_rmse)
        },
        "maximum_constant_state_flux_divergence": float(
            projected_divergence.abs().max()
        ),
        "note": (
            "The periodic trajectory objective identifies flux divergence, "
            "not the additive flux gauge; this diagnostic is audit-only."
        ),
    }


def self_test() -> None:
    initial = base.generate_ic(4, base.NREF, 2991, ood=False)
    data = torch.from_numpy(shared.strict_rollout_reference(initial, nsnap=4))
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    model = base.Solver(MODEL_NAME, mean, std, width=16)
    loss = audit.deterministic_one_step_loss(
        model, data, state_std, batch_size=8
    )
    diagnostics = proposal_diagnostics(model, data)
    if not np.isfinite(loss):
        raise RuntimeError("Full-flux self-test produced a nonfinite loss")
    if diagnostics["raw_flux_absolute_max"] != 0.0:
        raise RuntimeError("Full-flux final layer is not zero initialized")
    print("self-test passed")


def run(args: argparse.Namespace) -> None:
    output = Path(args.outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    total_started = time.perf_counter()

    print("Generating strict train, validation, and test trajectories...")
    data_started = time.perf_counter()
    train_data = shared.make_baseline_training_data(args.seed)
    validation_data = audit.make_validation_data(args.seed)
    evaluation_suite = shared.make_evaluation_suite(args.seed)
    data_seconds = time.perf_counter() - data_started

    primitive = base.primitive(train_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))

    model, result, curve = audit.train_to_convergence(
        arm=ARM,
        model_name=MODEL_NAME,
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
    )

    convergence = {
        key: value
        for key, value in result.items()
        if key not in ("fixed", "best")
    }
    shared.write_csv(output / f"training_curve_seed{args.seed}.csv", curve)
    shared.write_csv(
        output / f"convergence_seed{args.seed}.csv", [convergence]
    )
    if not bool(result["converged"]):
        checkpoint = output / f"{ARM}_converged_best_seed{args.seed}.pt"
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            "Full-flux arm did not meet the validation convergence rule; "
            "the nonconverged checkpoint was removed."
        )

    metrics = audit.evaluate_state(
        MODEL_NAME,
        result["best"],
        mean,
        std,
        args.width,
        state_std,
        evaluation_suite,
        args.seed,
        ARM,
        "converged_best",
    )
    comparisons = comparison_rows(metrics, args.seed)
    proposal = proposal_diagnostics(model, validation_data)
    consistency = constant_state_consistency(model, train_data)

    shared.write_csv(output / f"metrics_seed{args.seed}.csv", metrics)
    shared.write_csv(
        output / f"comparison_seed{args.seed}.csv", comparisons
    )

    metadata = {
        "seed": args.seed,
        "pre_run_question": (
            "Can the same five-cell MLP learn the complete shared Euler "
            "interface flux without an HLLC base or correction parameterization?"
        ),
        "single_scientific_factor": "complete_flux_vs_hllc_plus_correction",
        "candidate": (
            "five-point normalized primitive stencil mapped directly to the "
            "complete three-component interface flux"
        ),
        "removed_from_candidate": [
            "HLLC base flux",
            "jump gate",
            "Roe characteristic basis",
            "analytic output scale",
        ],
        "held_fixed": [
            "MLP dimensions and trainable parameter count",
            "broad training and independent validation tensors",
            "initialization and minibatch seeds",
            "optimizer, learning-rate schedule, and convergence rule",
            "hard Tadmor projection",
            "local positivity limiter",
            "fully-discrete entropy limiter",
        ],
        "controlled_comparators": list(COMPARATORS),
        "parameter_count": convergence["parameter_count"],
        "training_trajectories": int(train_data.shape[0]),
        "validation_trajectories": int(validation_data.shape[0]),
        "fixed_budget_updates": args.fixed_budget_updates,
        "max_updates": args.max_updates,
        "validation_interval": args.validation_interval,
        "plateau_patience": args.plateau_patience,
        "minimum_relative_improvement": args.minimum_relative_improvement,
        "initial_learning_rate": args.lr,
        "learning_rate_factor": args.learning_rate_factor,
        "minimum_learning_rate": args.minimum_learning_rate,
        "width": args.width,
        "batch_size": args.batch_size,
        "reference": "512-cell Rusanov + SSP-RK2, strict/no repair",
        "training_and_deployment_cells_for_reported_metrics": 64,
        "boundary_condition": "periodic",
        "checkpoint_selection": "minimum independent validation rollout NRMSE",
        "constant_state_consistency_used_for_selection": False,
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
    }
    report = {
        "metadata": metadata,
        "convergence": convergence,
        "proposal_diagnostics_on_validation_ground_truth_states": proposal,
        "constant_state_consistency": consistency,
        "metrics": metrics,
        "comparisons": comparisons,
    }
    (output / f"report_seed{args.seed}.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    (output / f"metadata_seed{args.seed}.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
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
    parser.add_argument("--outdir", default=str(HERE / "results"))
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    parsed = parse_args()
    if parsed.self_test:
        self_test()
    else:
        run(parsed)
