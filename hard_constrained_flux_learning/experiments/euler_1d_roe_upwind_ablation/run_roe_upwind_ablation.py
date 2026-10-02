"""Train one preregistered Roe-coordinate/upwind Euler ablation arm."""

from __future__ import annotations

import argparse
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
CENTRAL_EXPERIMENT = (
    HCFL_ROOT / "experiments" / "euler_1d_consistent_central_flux"
)
for module_path in (AUDIT_EXPERIMENT, FULL_EXPERIMENT, CENTRAL_EXPERIMENT):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

import run_consistent_central_flux as central_baseline  # noqa: E402
import run_convergence_audit as audit  # noqa: E402
import run_full_flux_baseline as full_baseline  # noqa: E402


shared = audit.shared
base = audit.base
ARMS = {
    "roe_complete_broad": "roe_complete",
    "central_roe_signed_broad": "central_roe_signed",
    "central_roe_upwind_broad": "central_roe_upwind",
}


def stencil_logits(model: base.Solver, states: torch.Tensor) -> torch.Tensor:
    primitive = base.primitive(states)
    flux_net = model.flux_net
    features = torch.cat([
        (torch.roll(primitive, shift, dims=-2) - flux_net.mean)
        / flux_net.std
        for shift in [2, 1, 0, -1, -2]
    ], dim=-1)
    return flux_net.net(features)


@torch.no_grad()
def multiplier_diagnostics(
    model: base.Solver,
    model_name: str,
    data: torch.Tensor,
    batch_size: int = 256,
) -> dict[str, Any] | None:
    if model_name not in ("central_roe_signed", "central_roe_upwind"):
        return None
    states = data.reshape(-1, data.shape[-2], 3)
    chunks: list[torch.Tensor] = []
    for start in range(0, states.shape[0], batch_size):
        logits = stencil_logits(model, states[start : start + batch_size])
        if model_name == "central_roe_signed":
            multiplier = 1.0 + 2.0 * torch.tanh(logits)
        else:
            multiplier = 1.0 + torch.tanh(logits)
        chunks.append(multiplier.double().cpu())
    values = torch.cat(chunks, dim=0)
    flat = values.reshape(-1)
    return {
        "minimum_multiplier": float(flat.min()),
        "maximum_multiplier": float(flat.max()),
        "mean_multiplier": float(flat.mean()),
        "standard_deviation": float(flat.std()),
        "negative_multiplier_fraction": float((flat < 0.0).double().mean()),
        "near_zero_fraction": float((flat < 0.05).double().mean()),
        "near_upper_nominal_bound_fraction": float(
            (flat > 1.95).double().mean()
        ),
        "mean_by_wave": {
            "left_acoustic": float(values[..., 0].mean()),
            "contact": float(values[..., 1].mean()),
            "right_acoustic": float(values[..., 2].mean()),
        },
    }


def self_test() -> None:
    initial = base.generate_ic(4, base.NREF, 4491, ood=False)
    data = torch.from_numpy(shared.strict_rollout_reference(initial, nsnap=4))
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    for model_name in ARMS.values():
        model = base.Solver(model_name, mean, std, width=16)
        loss = audit.deterministic_one_step_loss(
            model, data, state_std, batch_size=8
        )
        if not np.isfinite(loss):
            raise RuntimeError(f"{model_name} self-test loss is nonfinite")
        if model_name.startswith("central_roe_"):
            consistency = central_baseline.equal_interface_consistency(
                model, data
            )
            if consistency["maximum_raw_absolute_error"] != 0.0:
                raise RuntimeError(
                    f"{model_name} raw proposal lost exact consistency"
                )
            if consistency["maximum_projected_absolute_error"] != 0.0:
                raise RuntimeError(
                    f"{model_name} projection lost exact consistency"
                )
    upwind = base.Solver("central_roe_upwind", mean, std, width=16)
    diagnostics = multiplier_diagnostics(
        upwind, "central_roe_upwind", data
    )
    if diagnostics is None or diagnostics["minimum_multiplier"] < 0.0:
        raise RuntimeError("Automatic-upwind multiplier became negative")
    print("self-test passed")


def run(args: argparse.Namespace) -> None:
    arm = args.arm
    model_name = ARMS[arm]
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
        arm=arm,
        model_name=model_name,
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
    checkpoint = output / f"{arm}_converged_best_seed{args.seed}.pt"
    if not bool(result["converged"]):
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"{arm} did not meet the convergence rule; its checkpoint and "
            "nonconverged result artifacts were not retained."
        )

    metrics = audit.evaluate_state(
        model_name,
        result["best"],
        mean,
        std,
        args.width,
        state_std,
        evaluation_suite,
        args.seed,
        arm,
        "converged_best",
    )
    proposal = full_baseline.proposal_diagnostics(model, validation_data)
    constant = full_baseline.constant_state_consistency(model, train_data)
    equal_interface = central_baseline.equal_interface_consistency(
        model, validation_data
    )
    correction = central_baseline.correction_diagnostics(
        model, validation_data
    )
    multipliers = multiplier_diagnostics(
        model, model_name, validation_data
    )

    curve_path = output / f"training_curve_{arm}_seed{args.seed}.csv"
    convergence_path = output / f"convergence_{arm}_seed{args.seed}.csv"
    metrics_path = output / f"metrics_{arm}_seed{args.seed}.csv"
    shared.write_csv(curve_path, curve)
    shared.write_csv(convergence_path, [convergence])
    shared.write_csv(metrics_path, metrics)

    metadata = {
        "seed": args.seed,
        "arm": arm,
        "model": model_name,
        "parameter_count": convergence["parameter_count"],
        "training_trajectories": int(train_data.shape[0]),
        "validation_trajectories": int(validation_data.shape[0]),
        "training_cells": int(train_data.shape[-2]),
        "fine_reference_cells": base.NREF,
        "reference": "512-cell Rusanov + SSP-RK2, strict/no repair",
        "boundary_condition": "periodic",
        "checkpoint_selection": "minimum independent validation rollout NRMSE",
        "canonical_cases_used_for_selection": False,
        "max_updates": args.max_updates,
        "validation_interval": args.validation_interval,
        "plateau_patience": args.plateau_patience,
        "minimum_relative_improvement": args.minimum_relative_improvement,
        "initial_learning_rate": args.lr,
        "learning_rate_factor": args.learning_rate_factor,
        "minimum_learning_rate": args.minimum_learning_rate,
        "width": args.width,
        "batch_size": args.batch_size,
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
    }
    report = {
        "metadata": metadata,
        "convergence": convergence,
        "proposal_diagnostics_on_validation_ground_truth_states": proposal,
        "raw_minus_central_flux_diagnostics": correction,
        "constant_state_consistency": constant,
        "equal_interface_with_varied_outer_stencil_consistency": (
            equal_interface
        ),
        "learned_wave_multiplier_diagnostics": multipliers,
        "metrics": metrics,
    }
    (output / f"report_{arm}_seed{args.seed}.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    (output / f"metadata_{arm}_seed{args.seed}.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=tuple(ARMS))
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
    args = parser.parse_args()
    if not args.self_test and args.arm is None:
        parser.error("--arm is required unless --self-test is used")
    return args


if __name__ == "__main__":
    parsed = parse_args()
    if parsed.self_test:
        self_test()
    else:
        run(parsed)
