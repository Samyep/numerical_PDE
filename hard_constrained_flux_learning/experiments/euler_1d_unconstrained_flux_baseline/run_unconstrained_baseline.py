"""Validation-converged unconstrained DNN flux baseline for 1D Euler.

The direct HLLC + MLP proposal is held fixed. Unlike HCFL, its proposal is used
without Tadmor projection, local admissibility limiting, or a fully-discrete
entropy line search. The finite-volume update remains conservative.
"""

from __future__ import annotations

import argparse
import copy
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
ARM = "unconstrained_direct_broad"
MODEL_NAME = "unconstrained_direct"
RHO_FLOOR = shared.RHO_FLOOR
PRESSURE_FLOOR = shared.PRESSURE_FLOOR
ENTROPY_RESIDUAL_TOLERANCE = shared.ENTROPY_RESIDUAL_TOLERANCE


class UnconstrainedDirectSolver(torch.nn.Module):
    """Existing direct proposal with no hard flux post-processing."""

    def __init__(
        self,
        mean: np.ndarray,
        std: np.ndarray,
        width: int = 72,
    ) -> None:
        super().__init__()
        self.flux_net = base.DirectFlux(mean, std, width=width)

    def flux(self, state: torch.Tensor) -> torch.Tensor:
        return self.flux_net(state)

    def one_step(self, state: torch.Tensor) -> torch.Tensor:
        return base.fv_step(state, self.flux(state))


def clone_state_dict(
    model: torch.nn.Module,
) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


def finite_minimum(values: torch.Tensor, default: float) -> float:
    finite = values[torch.isfinite(values)]
    if finite.numel() == 0:
        return default
    return float(finite.min())


@torch.no_grad()
def evaluate_raw_dataset(
    model: UnconstrainedDirectSolver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    split: str,
) -> dict[str, Any]:
    """Roll out the raw learned flux and record failures without repair."""
    started = time.perf_counter()
    was_training = model.training
    model.eval()

    state = data[:, 0].clone()
    trajectory_count = state.shape[0]
    active = torch.ones(trajectory_count, dtype=torch.bool)
    first_failure = torch.full(
        (trajectory_count,),
        -1,
        dtype=torch.int64,
    )
    squared_errors = torch.full(
        (trajectory_count, data.shape[1] - 1),
        float("nan"),
        dtype=state.dtype,
    )
    min_rho = float(state[..., 0].min())
    min_pressure = float(base.t_pressure(state).min())
    entropy_violations = 0
    entropy_total = 0
    max_entropy_residual = -float("inf")

    for snapshot in range(1, data.shape[1]):
        flux = model.flux(state)
        if bool(active.any()):
            residual = shared.entropy_residual64(
                flux[active],
                state[active],
            )
            entropy_violations += int(
                (residual > ENTROPY_RESIDUAL_TOLERANCE).sum()
            )
            entropy_total += residual.numel()
            max_entropy_residual = max(
                max_entropy_residual,
                float(residual.max()),
            )

        candidate = base.fv_step(state, flux)
        candidate_rho = candidate[..., 0]
        candidate_pressure = base.t_pressure(candidate)
        if bool(active.any()):
            min_rho = min(
                min_rho,
                finite_minimum(candidate_rho[active], min_rho),
            )
            min_pressure = min(
                min_pressure,
                finite_minimum(candidate_pressure[active], min_pressure),
            )

        finite = torch.isfinite(candidate).all(dim=(1, 2))
        admissible = (
            (candidate_rho >= RHO_FLOOR).all(dim=1)
            & (candidate_pressure >= PRESSURE_FLOOR).all(dim=1)
        )
        survives = active & finite & admissible
        newly_failed = active & ~survives
        first_failure[newly_failed] = snapshot

        if bool(survives.any()):
            normalized = (
                candidate[survives] - data[survives, snapshot]
            ) / state_std
            squared_errors[survives, snapshot - 1] = (
                normalized * normalized
            ).mean(dim=(1, 2))

        state = torch.where(
            survives[:, None, None],
            candidate,
            state,
        )
        active = survives

    completed = int(active.sum())
    failed = trajectory_count - completed
    partial_values = squared_errors[torch.isfinite(squared_errors)]
    partial_nrmse = (
        float(torch.sqrt(partial_values.mean()))
        if partial_values.numel()
        else None
    )
    completed_values = squared_errors[active]
    completed_nrmse = (
        float(torch.sqrt(completed_values.mean()))
        if completed_values.numel()
        else None
    )
    rollout_nrmse = completed_nrmse if failed == 0 else None
    failure_snapshot_values = first_failure[first_failure >= 0]
    first_failure_snapshot = (
        int(failure_snapshot_values.min())
        if failure_snapshot_values.numel()
        else None
    )
    selection_score = (
        float(rollout_nrmse)
        if rollout_nrmse is not None
        else 1.0
        + failed / trajectory_count
        + min(float(partial_nrmse) if partial_nrmse is not None else 1.0, 0.999)
    )

    model.train(was_training)
    return {
        "split": split,
        "ntraj": int(trajectory_count),
        "completed_trajectories": completed,
        "failed_trajectories": failed,
        "completion_rate": completed / trajectory_count,
        "first_failure_snapshot": first_failure_snapshot,
        "rollout_nrmse": rollout_nrmse,
        "completed_rollout_nrmse": completed_nrmse,
        "partial_rollout_nrmse": partial_nrmse,
        "selection_score": selection_score,
        "min_rho": min_rho,
        "min_pressure": min_pressure,
        "entropy_violation_rate": (
            entropy_violations / max(entropy_total, 1)
        ),
        "max_entropy_residual": max_entropy_residual,
        "local_limiter_intervention_rate": None,
        "fd_entropy_intervention_rate": None,
        "eval_seconds": time.perf_counter() - started,
    }


def train_to_convergence(
    train_data: torch.Tensor,
    validation_data: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    state_std: torch.Tensor,
    seed: int,
    width: int,
    batch_size: int,
    learning_rate: float,
    max_updates: int,
    validation_interval: int,
    fixed_budget_updates: int,
    plateau_patience: int,
    minimum_relative_improvement: float,
    learning_rate_factor: float,
    minimum_learning_rate: float,
    output: Path,
) -> tuple[
    UnconstrainedDirectSolver,
    dict[str, Any],
    list[dict[str, Any]],
]:
    if max_updates < fixed_budget_updates:
        raise ValueError("max_updates must include the fixed-budget checkpoint")

    torch.manual_seed(12000 + seed)
    model = UnconstrainedDirectSolver(mean, std, width=width)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    generator = torch.Generator().manual_seed(13000 + seed)

    best_score = float("inf")
    best_update = -1
    best_state: dict[str, torch.Tensor] | None = None
    fixed_state: dict[str, torch.Tensor] | None = None
    plateau_anchor = float("inf")
    checks_without_improvement = 0
    converged = False
    stop_reason = "max_updates"
    curve: list[dict[str, Any]] = []
    started = time.perf_counter()

    def validate(update: int, batch_loss: float | None) -> float:
        nonlocal best_score, best_update, best_state
        train_loss = audit.deterministic_one_step_loss(
            model,
            train_data,
            state_std,
        )
        validation_loss = audit.deterministic_one_step_loss(
            model,
            validation_data,
            state_std,
        )
        rollout = evaluate_raw_dataset(
            model,
            validation_data,
            state_std,
            "validation",
        )
        score = float(rollout["selection_score"])
        improved = score < best_score
        if improved:
            best_score = score
            best_update = update
            best_state = clone_state_dict(model)
            torch.save(
                best_state,
                output / f"{ARM}_converged_best_seed{seed}.pt",
            )
        row = {
            "seed": seed,
            "arm": ARM,
            "model": MODEL_NAME,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "full_train_one_step_loss": train_loss,
            "validation_one_step_loss": validation_loss,
            "validation_selection_score": score,
            "validation_rollout_nrmse": rollout["rollout_nrmse"],
            "validation_partial_rollout_nrmse": rollout[
                "partial_rollout_nrmse"
            ],
            "validation_completion_rate": rollout["completion_rate"],
            "validation_entropy_violation_rate": rollout[
                "entropy_violation_rate"
            ],
            "new_absolute_best": improved,
        }
        curve.append(row)
        print(json.dumps(row))
        return score

    metric = validate(0, None)
    last_batch_loss: float | None = None
    for update in range(1, max_updates + 1):
        model.train()
        indices = torch.randint(
            0,
            train_data.shape[0],
            (batch_size,),
            generator=generator,
        )
        times = torch.randint(
            0,
            train_data.shape[1] - 1,
            (batch_size,),
            generator=generator,
        )
        inputs = train_data[indices, times]
        targets = train_data[indices, times + 1]
        prediction = model.one_step(inputs)
        loss = (((prediction - targets) / state_std) ** 2).mean()
        if not bool(torch.isfinite(loss)):
            raise RuntimeError(f"Nonfinite training loss at update {update}")

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_batch_loss = float(loss.detach())

        if update == fixed_budget_updates:
            fixed_state = clone_state_dict(model)
            torch.save(
                fixed_state,
                output / f"{ARM}_fixed1100_seed{seed}.pt",
            )

        if update % validation_interval != 0 and update != max_updates:
            continue

        metric = validate(update, last_batch_loss)
        if update < fixed_budget_updates:
            continue
        if update == fixed_budget_updates:
            plateau_anchor = best_score
            checks_without_improvement = 0
            continue

        if metric <= plateau_anchor * (1.0 - minimum_relative_improvement):
            plateau_anchor = metric
            checks_without_improvement = 0
        else:
            checks_without_improvement += 1

        if checks_without_improvement < plateau_patience:
            continue

        current_lr = float(optimizer.param_groups[0]["lr"])
        if current_lr <= minimum_learning_rate * (1.0 + 1.0e-12):
            converged = True
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break

        new_lr = max(
            current_lr * learning_rate_factor,
            minimum_learning_rate,
        )
        for group in optimizer.param_groups:
            group["lr"] = new_lr
        checks_without_improvement = 0
        plateau_anchor = best_score
        print(
            json.dumps(
                {
                    "seed": seed,
                    "arm": ARM,
                    "event": "reduce_learning_rate",
                    "update": update,
                    "old_learning_rate": current_lr,
                    "new_learning_rate": new_lr,
                    "best_validation_selection_score": best_score,
                }
            )
        )

    if fixed_state is None:
        raise RuntimeError("The fixed 1,100-update checkpoint was not created")
    if best_state is None:
        raise RuntimeError("No validation checkpoint was created")

    summary = {
        "seed": seed,
        "arm": ARM,
        "model": MODEL_NAME,
        "parameter_count": sum(
            parameter.numel() for parameter in model.parameters()
        ),
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "best_validation_selection_score": best_score,
        "best_validation_rollout_nrmse": next(
            row["validation_rollout_nrmse"]
            for row in curve
            if row["update"] == best_update
        ),
        "best_validation_completion_rate": next(
            row["validation_completion_rate"]
            for row in curve
            if row["update"] == best_update
        ),
        "final_learning_rate": float(optimizer.param_groups[0]["lr"]),
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
    }
    model.load_state_dict(best_state)
    model.eval()
    return model, {"fixed": fixed_state, "best": best_state, **summary}, curve


def evaluate_state(
    state: dict[str, torch.Tensor],
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
    state_std: torch.Tensor,
    suite: dict[str, torch.Tensor],
    seed: int,
    stage: str,
) -> list[dict[str, Any]]:
    model = UnconstrainedDirectSolver(mean, std, width=width)
    model.load_state_dict(state)
    model.eval()
    parameter_count = sum(
        parameter.numel() for parameter in model.parameters()
    )
    rows: list[dict[str, Any]] = []
    for split, data in suite.items():
        row = evaluate_raw_dataset(model, data, state_std, split)
        row.update(
            {
                "seed": seed,
                "arm": ARM,
                "model": MODEL_NAME,
                "stage": stage,
                "parameter_count": parameter_count,
            }
        )
        rows.append(row)
        print(
            json.dumps(
                {
                    "seed": seed,
                    "arm": ARM,
                    "stage": stage,
                    "split": split,
                    "completion_rate": row["completion_rate"],
                    "rollout_nrmse": row["rollout_nrmse"],
                    "partial_rollout_nrmse": row["partial_rollout_nrmse"],
                }
            )
        )
    return rows


def read_constrained_direct_metrics(seed: int) -> dict[str, dict[str, str]]:
    path = AUDIT_EXPERIMENT / "results" / f"metrics_seed{seed}.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        row["split"]: row
        for row in rows
        if row["arm"] == "direct_broad"
        and row["stage"] == "converged_best"
    }


def make_comparison(
    metrics: list[dict[str, Any]],
    seed: int,
) -> list[dict[str, Any]]:
    plain = {
        row["split"]: row
        for row in metrics
        if row["stage"] == "converged_best"
    }
    constrained = read_constrained_direct_metrics(seed)
    rows: list[dict[str, Any]] = []
    for split in audit.TEST_SPLITS:
        candidate = plain[split]
        baseline = constrained[split]
        candidate_nrmse = candidate["rollout_nrmse"]
        baseline_nrmse = float(baseline["rollout_nrmse"])
        rows.append(
            {
                "seed": seed,
                "split": split,
                "comparator": "constrained_direct_broad_converged_best",
                "candidate": f"{ARM}_converged_best",
                "comparator_rollout_nrmse": baseline_nrmse,
                "candidate_rollout_nrmse": candidate_nrmse,
                "candidate_partial_rollout_nrmse": candidate[
                    "partial_rollout_nrmse"
                ],
                "candidate_completion_rate": candidate["completion_rate"],
                "candidate_entropy_violation_rate": candidate[
                    "entropy_violation_rate"
                ],
                "relative_change_percent": (
                    100.0 * (float(candidate_nrmse) / baseline_nrmse - 1.0)
                    if candidate_nrmse is not None
                    else None
                ),
            }
        )
    return rows


def self_test() -> None:
    initial = base.generate_ic(4, base.NREF, 1991, ood=False)
    data = torch.from_numpy(shared.strict_rollout_reference(initial, nsnap=4))
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    model = UnconstrainedDirectSolver(mean, std, width=16)
    loss = audit.deterministic_one_step_loss(
        model,
        data,
        state_std,
        batch_size=8,
    )
    metrics = evaluate_raw_dataset(model, data, state_std, "self_test")
    if not np.isfinite(loss) or metrics["completion_rate"] != 1.0:
        raise RuntimeError("Unconstrained baseline self-test failed")
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

    _model, result, curve = train_to_convergence(
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
    metrics = evaluate_state(
        result["fixed"],
        mean,
        std,
        args.width,
        state_std,
        evaluation_suite,
        args.seed,
        "fixed_1100",
    )
    metrics.extend(
        evaluate_state(
            result["best"],
            mean,
            std,
            args.width,
            state_std,
            evaluation_suite,
            args.seed,
            "converged_best",
        )
    )
    comparisons = make_comparison(metrics, args.seed)

    shared.write_csv(
        output / f"training_curve_seed{args.seed}.csv",
        curve,
    )
    shared.write_csv(
        output / f"convergence_seed{args.seed}.csv",
        [
            {
                key: value
                for key, value in result.items()
                if key not in ("fixed", "best")
            }
        ],
    )
    shared.write_csv(
        output / f"metrics_seed{args.seed}.csv",
        metrics,
    )
    shared.write_csv(
        output / f"comparison_vs_constrained_direct_seed{args.seed}.csv",
        comparisons,
    )

    metadata = {
        "seed": args.seed,
        "hypothesis": (
            "Removing hard output post-processing may slightly improve "
            "ordinary validation error, but will increase entropy violations "
            "and can lose admissibility on severe shock tubes."
        ),
        "single_scientific_factor": "hard_flux_output_safety_stack",
        "candidate": (
            "same DirectFlux proposal emitted without Tadmor projection, "
            "local admissibility limiting, or fully-discrete entropy limiting"
        ),
        "controlled_comparator": "direct_broad_converged_best",
        "conservation": "retained through one shared interface flux",
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
        "boundary_condition": "periodic",
        "raw_rollout_failure_definition": {
            "density_floor": RHO_FLOOR,
            "pressure_floor": PRESSURE_FLOOR,
            "nonfinite_is_failure": True,
            "repair_or_fallback": False,
        },
        "selection_rule": (
            "fully completed validation rollouts outrank any checkpoint with "
            "a failed validation trajectory; then minimize rollout NRMSE"
        ),
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
    }
    (output / f"metadata_seed{args.seed}.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=56)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=1100)
    parser.add_argument("--max-updates", type=int, default=20000)
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
