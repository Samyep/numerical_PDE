"""Validation-controlled convergence audit for 1D Euler HCFL variants."""

from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
COVERAGE_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_wave_coverage_ablation"
if str(COVERAGE_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(COVERAGE_EXPERIMENT))

import run_coverage_ablation as shared  # noqa: E402


base = shared.base
ARMS: dict[str, tuple[str, str]] = {
    "direct_broad": ("direct", "broad"),
    "invariant_broad": ("invariant", "broad"),
    "characteristic_broad": ("characteristic", "broad"),
    "dissipation_broad": ("dissipation", "broad"),
    "conv_broad": ("conv", "broad"),
    "direct_wave": ("direct", "wave"),
}
TEST_SPLITS = (
    "ordinary_id",
    "broad_random_in_support",
    "moderate_ood_high_frequency",
    *shared.CANONICAL_CASES,
)


def clone_state_dict(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


def make_validation_data(seed: int) -> torch.Tensor:
    """Independent common validation distribution for every audit arm."""
    ordinary = shared.strict_rollout_reference(
        base.generate_ic(44, base.NREF, 14000 + seed, ood=False)
    )
    broad = shared.strict_rollout_reference(
        base.generate_ic(52, base.NREF, 15000 + seed, ood=True)
    )
    extreme = shared.strict_rollout_reference(
        base.generate_extreme_ic(20, base.NREF, 16000 + seed)
    )
    wave = shared.strict_rollout_reference(
        shared.generate_wave_coverage_ic(20, base.NREF, 17000 + seed)
    )
    return torch.from_numpy(np.concatenate([ordinary, broad, extreme, wave], axis=0))


@torch.no_grad()
def deterministic_one_step_loss(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    batch_size: int = 256,
) -> float:
    inputs = data[:, :-1].reshape(-1, base.NCOARSE, 3)
    targets = data[:, 1:].reshape(-1, base.NCOARSE, 3)
    total_squared = 0.0
    total_values = 0
    was_training = model.training
    model.eval()
    for start in range(0, inputs.shape[0], batch_size):
        stop = min(start + batch_size, inputs.shape[0])
        prediction = model.one_step(inputs[start:stop])
        normalized = (prediction - targets[start:stop]) / state_std
        total_squared += float((normalized * normalized).sum())
        total_values += normalized.numel()
    model.train(was_training)
    return total_squared / total_values


@torch.no_grad()
def validation_rollout_nrmse(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
) -> float:
    was_training = model.training
    model.eval()
    value = float(
        shared.evaluate_dataset(model, data, state_std, "validation")[
            "rollout_nrmse"
        ]
    )
    model.train(was_training)
    return value


def train_to_convergence(
    arm: str,
    model_name: str,
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
) -> tuple[base.Solver, dict[str, Any], list[dict[str, Any]]]:
    """Train one arm until safe validation rollout NRMSE reaches a plateau."""
    if max_updates < fixed_budget_updates:
        raise ValueError("max_updates must include the fixed-budget checkpoint")

    torch.manual_seed(12000 + seed)
    model = base.Solver(model_name, mean, std, width=width)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    generator = torch.Generator().manual_seed(13000 + seed)

    best_metric = float("inf")
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
        nonlocal best_metric, best_update, best_state
        train_loss = deterministic_one_step_loss(model, train_data, state_std)
        validation_loss = deterministic_one_step_loss(
            model, validation_data, state_std
        )
        rollout_nrmse = validation_rollout_nrmse(
            model, validation_data, state_std
        )
        improved = rollout_nrmse < best_metric
        if improved:
            best_metric = rollout_nrmse
            best_update = update
            best_state = clone_state_dict(model)
            torch.save(
                best_state, output / f"{arm}_converged_best_seed{seed}.pt"
            )
        row = {
            "seed": seed,
            "arm": arm,
            "model": model_name,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "full_train_one_step_loss": train_loss,
            "validation_one_step_loss": validation_loss,
            "validation_rollout_nrmse": rollout_nrmse,
            "new_absolute_best": improved,
        }
        curve.append(row)
        print(json.dumps(row))
        return rollout_nrmse

    metric = validate(0, None)
    last_batch_loss: float | None = None
    for update in range(1, max_updates + 1):
        model.train()
        indices = torch.randint(
            0, train_data.shape[0], (batch_size,), generator=generator
        )
        times = torch.randint(
            0, train_data.shape[1] - 1, (batch_size,), generator=generator
        )
        inputs = train_data[indices, times]
        targets = train_data[indices, times + 1]
        prediction = model.one_step(inputs)
        loss = (((prediction - targets) / state_std) ** 2).mean()

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_batch_loss = float(loss.detach())

        if update == fixed_budget_updates:
            fixed_state = clone_state_dict(model)

        if update % validation_interval != 0 and update != max_updates:
            continue

        metric = validate(update, last_batch_loss)
        if update < fixed_budget_updates:
            continue
        if update == fixed_budget_updates:
            plateau_anchor = best_metric
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
            current_lr * learning_rate_factor, minimum_learning_rate
        )
        for group in optimizer.param_groups:
            group["lr"] = new_lr
        checks_without_improvement = 0
        plateau_anchor = best_metric
        print(
            json.dumps(
                {
                    "seed": seed,
                    "arm": arm,
                    "event": "reduce_learning_rate",
                    "update": update,
                    "old_learning_rate": current_lr,
                    "new_learning_rate": new_lr,
                    "best_validation_rollout_nrmse": best_metric,
                }
            )
        )

    if fixed_state is None:
        raise RuntimeError("The fixed 1,100-update state was not captured")
    if best_state is None:
        raise RuntimeError("No validation checkpoint was created")

    elapsed = time.perf_counter() - started
    stop_update = int(curve[-1]["update"])
    summary = {
        "seed": seed,
        "arm": arm,
        "model": model_name,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "best_update": best_update,
        "stop_update": stop_update,
        "best_validation_rollout_nrmse": best_metric,
        "final_learning_rate": float(optimizer.param_groups[0]["lr"]),
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": elapsed,
    }
    model.load_state_dict(best_state)
    model.eval()
    return model, {"fixed": fixed_state, "best": best_state, **summary}, curve


def evaluate_state(
    model_name: str,
    state: dict[str, torch.Tensor],
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
    state_std: torch.Tensor,
    suite: dict[str, torch.Tensor],
    seed: int,
    arm: str,
    stage: str,
) -> list[dict[str, Any]]:
    model = base.Solver(model_name, mean, std, width=width)
    model.load_state_dict(state)
    model.eval()
    rows: list[dict[str, Any]] = []
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    for split, data in suite.items():
        row = shared.evaluate_dataset(model, data, state_std, split)
        row.update(
            {
                "seed": seed,
                "arm": arm,
                "model": model_name,
                "stage": stage,
                "parameter_count": parameter_count,
            }
        )
        rows.append(row)
        print(
            json.dumps(
                {
                    "seed": seed,
                    "arm": arm,
                    "stage": stage,
                    "split": split,
                    "nrmse": row["rollout_nrmse"],
                }
            )
        )
    return rows


def comparison_rows(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {
        (row["arm"], row["stage"], row["split"]): row for row in metrics
    }
    rows: list[dict[str, Any]] = []

    for arm in ARMS:
        for split in TEST_SPLITS:
            fixed = indexed[(arm, "fixed_1100", split)]
            converged = indexed[(arm, "converged_best", split)]
            old = float(fixed["rollout_nrmse"])
            new = float(converged["rollout_nrmse"])
            rows.append(
                {
                    "comparison": "within_arm_convergence",
                    "arm": arm,
                    "comparator": "fixed_1100",
                    "candidate": "converged_best",
                    "split": split,
                    "comparator_nrmse": old,
                    "candidate_nrmse": new,
                    "relative_change_percent": 100.0 * (new / old - 1.0),
                }
            )

    for arm in ARMS:
        if arm == "direct_broad":
            continue
        comparison = (
            "data_distribution"
            if arm == "direct_wave"
            else "architecture"
        )
        for split in TEST_SPLITS:
            baseline = indexed[("direct_broad", "converged_best", split)]
            candidate = indexed[(arm, "converged_best", split)]
            old = float(baseline["rollout_nrmse"])
            new = float(candidate["rollout_nrmse"])
            rows.append(
                {
                    "comparison": comparison,
                    "arm": arm,
                    "comparator": "direct_broad_converged",
                    "candidate": f"{arm}_converged",
                    "split": split,
                    "comparator_nrmse": old,
                    "candidate_nrmse": new,
                    "relative_change_percent": 100.0 * (new / old - 1.0),
                }
            )
    return rows


def self_test() -> None:
    initial = base.generate_ic(4, base.NREF, 991, ood=False)
    data = torch.from_numpy(shared.strict_rollout_reference(initial, nsnap=3))
    primitive = base.primitive(data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = data.std(dim=(0, 1, 2))
    for model_name in ("full", "central_consistent", "direct", "invariant", "characteristic", "dissipation", "conv"):
        model = base.Solver(model_name, mean, std, width=16)
        value = deterministic_one_step_loss(model, data, state_std, batch_size=8)
        if not np.isfinite(value):
            raise RuntimeError(f"Nonfinite self-test loss for {model_name}")
    print("self-test passed")


def run(args: argparse.Namespace) -> None:
    output = Path(args.outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    total_started = time.perf_counter()

    print("Generating strict train, validation, and test trajectories...")
    data_started = time.perf_counter()
    broad_data = shared.make_baseline_training_data(args.seed)
    wave_data = shared.make_wave_coverage_training_data(broad_data, args.seed)
    validation_data = make_validation_data(args.seed)
    evaluation_suite = shared.make_evaluation_suite(args.seed)
    data_seconds = time.perf_counter() - data_started

    primitive = base.primitive(broad_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = broad_data.std(dim=(0, 1, 2))
    datasets = {"broad": broad_data, "wave": wave_data}

    curves: list[dict[str, Any]] = []
    convergence_rows: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []

    for arm, (model_name, data_name) in ARMS.items():
        print(json.dumps({"stage": "train_arm", "arm": arm}))
        _best_model, result, arm_curve = train_to_convergence(
            arm=arm,
            model_name=model_name,
            train_data=datasets[data_name],
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
        curves.extend(arm_curve)
        convergence_rows.append(
            {key: value for key, value in result.items() if key not in ("fixed", "best")}
        )
        metrics.extend(
            evaluate_state(
                model_name,
                result["fixed"],
                mean,
                std,
                args.width,
                state_std,
                evaluation_suite,
                args.seed,
                arm,
                "fixed_1100",
            )
        )
        metrics.extend(
            evaluate_state(
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
        )
        shared.write_csv(output / f"training_curve_seed{args.seed}.csv", curves)
        shared.write_csv(output / f"convergence_seed{args.seed}.csv", convergence_rows)
        shared.write_csv(output / f"metrics_seed{args.seed}.csv", metrics)

    comparisons = comparison_rows(metrics)
    shared.write_csv(output / f"comparison_seed{args.seed}.csv", comparisons)

    metadata = {
        "seed": args.seed,
        "purpose": "re-audit fixed-budget architecture and wave-coverage conclusions after validation convergence",
        "arms": ARMS,
        "training_trajectories_per_arm": 580,
        "validation_composition": {
            "ordinary": 44,
            "broad_random": 52,
            "random_extreme": 20,
            "structured_wave": 20,
        },
        "validation_used_canonical_cases": False,
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
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
    }
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
