"""Train an exactly consistent central-flux plus learned-correction arm."""

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
FULL_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_full_flux_baseline"
for module_path in (AUDIT_EXPERIMENT, FULL_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import run_convergence_audit as audit  # noqa: E402
import run_full_flux_baseline as full_baseline  # noqa: E402


shared = audit.shared
base = audit.base
ARM = "central_consistent_broad"
MODEL_NAME = "central_consistent"
COMPARATORS = (
    "full_flux_broad",
    "direct_broad",
    "dissipation_broad",
)


def read_comparator_metrics(
    seed: int,
) -> dict[tuple[str, str], dict[str, str]]:
    sources = (
        AUDIT_EXPERIMENT / "results" / f"metrics_seed{seed}.csv",
        FULL_EXPERIMENT / "results" / f"metrics_seed{seed}.csv",
    )
    selected: dict[tuple[str, str], dict[str, str]] = {}
    for path in sources:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if (
                    row["arm"] in COMPARATORS
                    and row["stage"] == "converged_best"
                ):
                    selected[(row["arm"], row["split"])] = row
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
def correction_diagnostics(
    model: base.Solver,
    data: torch.Tensor,
    batch_size: int = 256,
) -> dict[str, float]:
    states = data.reshape(-1, data.shape[-2], 3)
    correction_squared = 0.0
    central_squared = 0.0
    correction_divergence_squared = 0.0
    central_divergence_squared = 0.0
    maximum_correction = 0.0
    values = 0
    for start in range(0, states.shape[0], batch_size):
        state = states[start : start + batch_size]
        right = torch.roll(state, -1, dims=-2)
        central = 0.5 * (base.t_flux(state) + base.t_flux(right))
        correction = model.flux_net(state) - central
        correction_divergence = correction - torch.roll(
            correction, 1, dims=-2
        )
        central_divergence = central - torch.roll(central, 1, dims=-2)
        correction_squared += float((correction.double() ** 2).sum())
        central_squared += float((central.double() ** 2).sum())
        correction_divergence_squared += float(
            (correction_divergence.double() ** 2).sum()
        )
        central_divergence_squared += float(
            (central_divergence.double() ** 2).sum()
        )
        maximum_correction = max(
            maximum_correction, float(correction.abs().max())
        )
        values += correction.numel()
    return {
        "learned_correction_rms": np.sqrt(correction_squared / values),
        "central_flux_rms": np.sqrt(central_squared / values),
        "correction_relative_flux_rms": np.sqrt(
            correction_squared / max(central_squared, 1.0e-30)
        ),
        "correction_divergence_relative_rms": np.sqrt(
            correction_divergence_squared
            / max(central_divergence_squared, 1.0e-30)
        ),
        "maximum_absolute_correction": maximum_correction,
    }


@torch.no_grad()
def equal_interface_consistency(
    model: base.Solver,
    data: torch.Tensor,
    sample_count: int = 256,
) -> dict[str, Any]:
    """Test equal interface states while leaving the outer stencil varied."""
    states = data.reshape(-1, data.shape[-2], 3)
    count = min(sample_count, states.shape[0])
    indices = torch.linspace(
        0, states.shape[0] - 1, count, dtype=torch.float64
    ).round().long()
    samples = states[indices].clone()
    interface = 2
    samples[:, interface + 1] = samples[:, interface]
    exact = base.t_flux(samples[:, interface])
    raw = model.flux_net(samples)[:, interface]
    projected = model.flux(samples)[:, interface]
    raw_error = raw - exact
    projected_error = projected - exact
    scale = torch.sqrt((exact.double() ** 2).mean(dim=0)).clamp_min(
        1.0e-12
    )
    return {
        "sample_count": count,
        "outer_stencil_was_held_constant": False,
        "maximum_raw_absolute_error": float(raw_error.abs().max()),
        "maximum_projected_absolute_error": float(
            projected_error.abs().max()
        ),
        "raw_relative_consistency_nrmse": float(torch.sqrt(
            ((raw_error.double() / scale) ** 2).mean()
        )),
        "projected_relative_consistency_nrmse": float(torch.sqrt(
            ((projected_error.double() / scale) ** 2).mean()
        )),
    }


def self_test() -> None:
    initial = base.generate_ic(4, base.NREF, 3991, ood=False)
    data = torch.from_numpy(shared.strict_rollout_reference(initial, nsnap=4))
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    model = base.Solver(MODEL_NAME, mean, std, width=16)
    loss = audit.deterministic_one_step_loss(
        model, data, state_std, batch_size=8
    )
    consistency = equal_interface_consistency(model, data)
    if not np.isfinite(loss):
        raise RuntimeError("Central-consistent self-test loss is nonfinite")
    if consistency["maximum_raw_absolute_error"] != 0.0:
        raise RuntimeError("Raw central proposal is not exactly consistent")
    if consistency["maximum_projected_absolute_error"] != 0.0:
        raise RuntimeError("Projected central proposal lost consistency")
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
            "Central-consistent arm did not meet the convergence rule; "
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
    proposal = full_baseline.proposal_diagnostics(model, validation_data)
    constant = full_baseline.constant_state_consistency(model, train_data)
    equal_interface = equal_interface_consistency(model, validation_data)
    correction = correction_diagnostics(model, validation_data)

    shared.write_csv(output / f"metrics_seed{args.seed}.csv", metrics)
    shared.write_csv(
        output / f"comparison_seed{args.seed}.csv", comparisons
    )

    metadata = {
        "seed": args.seed,
        "pre_run_question": (
            "Does restoring the analytic central physical flux and exact "
            "F_hat(U,U)=F(U) consistency repair direct full-flux learning?"
        ),
        "formula": (
            "0.5*(F(U_L)+F(U_R)) + 0.18*normalized_jump*scale*"
            "tanh(MLP(five_point_primitive_stencil))"
        ),
        "exact_consistency_mechanism": (
            "epsilon-free interface jump norm makes the correction exactly "
            "zero whenever U_L equals U_R"
        ),
        "roe_decomposition_used": False,
        "hllc_used": False,
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
        "consistency_diagnostics_used_for_selection": False,
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
    }
    report = {
        "metadata": metadata,
        "convergence": convergence,
        "proposal_diagnostics_on_validation_ground_truth_states": proposal,
        "learned_correction_diagnostics": correction,
        "constant_state_consistency": constant,
        "equal_interface_with_varied_outer_stencil_consistency": (
            equal_interface
        ),
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
