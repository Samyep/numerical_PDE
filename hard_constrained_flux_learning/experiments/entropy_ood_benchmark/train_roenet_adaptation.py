"""Train the exploratory same-data RoeNet architecture adaptation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
PAPER_DIR = HERE.parent / "paper_64_benchmark"
if str(PAPER_DIR) not in sys.path:
    sys.path.insert(0, str(PAPER_DIR))

import run_64_benchmark as paper  # noqa: E402
from roenet_adaptation import RoeNetEuler1d  # noqa: E402


def clone_state(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def admissible(state: torch.Tensor) -> torch.Tensor:
    finite = torch.isfinite(state).all(dim=(-1, -2))
    safe = torch.where(torch.isfinite(state), state, torch.zeros_like(state))
    density = safe[..., 0]
    pressure = paper.euler.base.t_pressure(safe)
    return (
        finite
        & (density >= paper.euler.shared.RHO_FLOOR).all(dim=-1)
        & (pressure >= paper.euler.shared.PRESSURE_FLOOR).all(dim=-1)
    )


@torch.no_grad()
def validation_metrics(
    model: RoeNetEuler1d,
    data: torch.Tensor,
    state_std: torch.Tensor,
    batch_size: int,
) -> dict[str, float]:
    model.eval()
    squared = 0.0
    count = 0
    completed = 0
    minimum_density = float("inf")
    minimum_pressure = float("inf")
    for start in range(0, data.shape[0], batch_size):
        target = data[start : start + batch_size].to(model.state_mean.device)
        state = target[:, 0]
        alive = admissible(state)
        sample_squared = torch.zeros(state.shape[0], device=state.device)
        sample_count = 0
        for snapshot in range(1, target.shape[1]):
            state = model(state, paper.euler.base.DT_SNAPSHOT)
            alive = alive & admissible(state)
            finite_state = torch.where(
                torch.isfinite(state), state, torch.zeros_like(state)
            )
            error = (finite_state - target[:, snapshot]) / state_std
            sample_squared += error.square().sum(dim=(1, 2))
            sample_count += error.shape[1] * error.shape[2]
            minimum_density = min(minimum_density, float(finite_state[..., 0].min()))
            minimum_pressure = min(
                minimum_pressure,
                float(paper.euler.base.t_pressure(finite_state).min()),
            )
        squared += float(sample_squared[alive].sum())
        count += int(alive.sum()) * sample_count
        completed += int(alive.sum())
    return {
        "rollout_nrmse_completed_only": (
            math.sqrt(squared / count) if count else float("inf")
        ),
        "completion_rate": completed / data.shape[0],
        "completed_trajectories": float(completed),
        "minimum_density": minimum_density,
        "minimum_pressure": minimum_pressure,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--validation-batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=1.0e-3)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--fixed-budget-updates", type=int, default=2000)
    parser.add_argument("--validation-interval", type=int, default=250)
    parser.add_argument("--plateau-patience", type=int, default=8)
    parser.add_argument("--minimum-relative-improvement", type=float, default=1.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--hidden-waves", type=int, default=64)
    parser.add_argument("--internal-steps", type=int, default=1)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--results-dir", type=Path, default=HERE / "results" / "roenet")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    device = torch.device(args.device)
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checkpoint = output / f"roenet_adaptation64_best_seed{args.seed}.pt"
    report_path = output / f"roenet_adaptation64_report_seed{args.seed}.json"
    if args.resume and checkpoint.exists() and report_path.exists():
        print(report_path.read_text(encoding="utf-8"))
        return

    paper_results = PAPER_DIR / "results"
    train, validation = paper.replicates.euler_cache(paper_results, args.seed)
    state_mean = train.mean(dim=(0, 1, 2))
    state_std_cpu = train.std(dim=(0, 1, 2)).clamp_min(1.0e-6)
    state_std = state_std_cpu.reshape(1, 1, 3).to(device)
    train_inputs = train[:, :-1].reshape(-1, 64, 3)
    train_targets = train[:, 1:].reshape(-1, 64, 3)
    model = RoeNetEuler1d(
        state_mean,
        state_std_cpu,
        cells=64,
        hidden_waves=args.hidden_waves,
        internal_steps=args.internal_steps,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(81173 + args.seed)

    best_error = float("inf")
    best_completion = -1.0
    best_update = 0
    best_state: dict[str, torch.Tensor] | None = None
    stale = 0
    reductions = 0
    curve: list[dict[str, Any]] = []
    converged = False
    stop_reason = "maximum update cap"
    plateau_completion = -1.0
    plateau_error = float("inf")
    started = time.perf_counter()

    for update in range(1, args.max_updates + 1):
        indices = torch.randint(
            0, train_inputs.shape[0], (args.batch_size,), generator=generator
        )
        inputs = train_inputs[indices].to(device)
        targets = train_targets[indices].to(device)
        model.train()
        prediction = model(inputs, paper.euler.base.DT_SNAPSHOT)
        loss = (((prediction - targets) / state_std) ** 2).mean()
        if not bool(torch.isfinite(loss)):
            raise RuntimeError(f"nonfinite RoeNet training loss at update {update}")
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        if update % args.validation_interval != 0:
            continue
        validation_result = validation_metrics(
            model,
            validation,
            state_std,
            args.validation_batch_size,
        )
        # Checkpoint selection is an exact lexicographic comparison:
        # completion first, then completed-only error.  Do not fold completion
        # into a scalar score because doing so changes the effective error
        # tolerance and can skip the true validation winner.
        error = validation_result["rollout_nrmse_completed_only"]
        completion = validation_result["completion_rate"]
        improved = (
            best_state is None
            or completion > best_completion + 1.0e-12
            or (
                abs(completion - best_completion) <= 1.0e-12
                and error < best_error
            )
        )
        if improved:
            best_error = error
            best_completion = completion
            best_update = update
            best_state = clone_state(model)

        # The relative-improvement threshold is used only for plateau/LR
        # scheduling; it never changes which checkpoint is retained.
        plateau_relative_improvement = (
            (plateau_error - error) / max(abs(plateau_error), 1.0e-12)
            if math.isfinite(plateau_error)
            else float("inf")
        )
        significant_progress = (
            plateau_completion < 0.0
            or completion > plateau_completion + 1.0e-12
            or (
                abs(completion - plateau_completion) <= 1.0e-12
                and plateau_relative_improvement >= args.minimum_relative_improvement
            )
        )
        if significant_progress:
            plateau_completion = completion
            plateau_error = error
            stale = 0
        else:
            stale += 1
        row = {
            "update": update,
            "training_loss": float(loss.detach()),
            "validation_score": error,
            **validation_result,
            "learning_rate": optimizer.param_groups[0]["lr"],
            "best_update": best_update,
            "new_best": improved,
            "significant_plateau_progress": significant_progress,
            "stale_checks": stale,
        }
        curve.append(row)
        print(json.dumps(row), flush=True)

        if update < args.fixed_budget_updates or stale < args.plateau_patience:
            continue
        current_lr = optimizer.param_groups[0]["lr"]
        next_lr = current_lr * args.learning_rate_factor
        if next_lr >= args.minimum_learning_rate * (1.0 + 1.0e-12):
            for group in optimizer.param_groups:
                group["lr"] = next_lr
            reductions += 1
            stale = 0
        else:
            converged = True
            stop_reason = "validation plateau at minimum learning rate"
            break

    if best_state is None:
        raise RuntimeError("RoeNet adaptation produced no checkpoint")
    model.load_state_dict(best_state)
    final_validation = validation_metrics(
        model, validation, state_std, args.validation_batch_size
    )
    torch.save(best_state, checkpoint)
    report = {
        "baseline": "same-data RoeNet architecture adaptation",
        "is_exact_paper_reproduction": False,
        "official_source_commit": "ef877957c1c0ddb16eac17006d75bdf7bd786d45",
        "controlled_changes": [
            "64-cell periodic benchmark instead of the paper's 200-cell Sod task",
            "same 580/136 trajectory split used by HCFL",
            "state normalization inside the learned maps",
            "regularized least-squares inverse instead of an unregularized matrix inverse",
            "validation-plateau stopping with a 50000-update cap",
        ],
        "seed": args.seed,
        "training_data_shape": list(train.shape),
        "validation_data_shape": list(validation.shape),
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "hidden_waves": args.hidden_waves,
        "internal_steps": args.internal_steps,
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "converged": converged,
        "stop_reason": stop_reason,
        "learning_rate_reductions": reductions,
        "best_validation_score": best_error,
        "checkpoint_selection_rule": "lexicographic(completed_trajectories, -completed_only_nrmse)",
        "best_validation_completion_rate": best_completion,
        "validation_metrics": final_validation,
        "training_seconds": time.perf_counter() - started,
        "checkpoint": checkpoint.name,
    }
    write_csv(output / f"roenet_adaptation64_curve_seed{args.seed}.csv", curve)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
