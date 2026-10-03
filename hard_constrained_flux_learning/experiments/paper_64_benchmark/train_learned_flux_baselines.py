"""Train matched unconstrained learned-flux finite-volume baselines.

The network proposes one periodic interface flux from the same symmetric
four-cell stencil used by the retained HCFL model.  The state update is a
telescoping flux difference, so conservation is architectural, but there is no
Tadmor projection, positivity limiter, or fully-discrete entropy limiter.
This is a same-data learned-neural-FVM family baseline, not an exact
reproduction of any one published implementation.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn


HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
EULER_DIR = EXPERIMENTS / "euler_1d_stencil_ablation"
SWE_DIR = EXPERIMENTS / "swe_1d_consistent_hcfl"
for module_dir in (HERE, EULER_DIR, SWE_DIR):
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import run_stencil_ablation as euler  # noqa: E402
import run_swe_consistent as swe  # noqa: E402
import train_hcfl_replicates as replicates  # noqa: E402
import train_operator_baselines as common  # noqa: E402


STENCIL_SHIFTS = (1, 0, -1, -2)


class UnconstrainedLearnedFluxFV(nn.Module):
    """HLL(C) plus an unconstrained, consistency-scaled local correction."""

    def __init__(
        self,
        system: str,
        primitive_mean: torch.Tensor,
        primitive_std: torch.Tensor,
        width: int = 72,
    ) -> None:
        super().__init__()
        if system not in ("euler", "swe"):
            raise ValueError(system)
        channels = 3 if system == "euler" else 2
        self.system = system
        self.channels = channels
        self.width = width
        self.hidden1 = nn.Linear(channels * len(STENCIL_SHIFTS), width)
        self.hidden2 = nn.Linear(width, width)
        self.output = nn.Linear(width, channels)
        nn.init.zeros_(self.output.weight)
        nn.init.zeros_(self.output.bias)
        self.register_buffer(
            "primitive_mean", primitive_mean.reshape(1, 1, channels).float()
        )
        self.register_buffer(
            "primitive_std", primitive_std.reshape(1, 1, channels).float()
        )
        flux_scale = (
            torch.tensor([0.6, 1.2, 2.5])
            if system == "euler"
            else torch.tensor([1.0, 5.0])
        )
        self.register_buffer("flux_scale", flux_scale.reshape(1, 1, channels))

    def primitive(self, state: torch.Tensor) -> torch.Tensor:
        if self.system == "euler":
            return euler.base.primitive(state)
        return swe.primitive(state)

    def base_flux(self, state: torch.Tensor) -> torch.Tensor:
        if self.system == "euler":
            return euler.base.t_hllc(state)
        return swe.t_hll(state)

    def flux(self, state: torch.Tensor) -> torch.Tensor:
        primitive = self.primitive(state)
        normalized = (primitive - self.primitive_mean) / self.primitive_std
        stencil = torch.cat(
            [torch.roll(normalized, shift, dims=-2) for shift in STENCIL_SHIFTS],
            dim=-1,
        )
        hidden = torch.tanh(self.hidden1(stencil))
        hidden = torch.tanh(self.hidden2(hidden))
        correction = torch.tanh(self.output(hidden))
        jump = torch.linalg.vector_norm(
            (
                torch.roll(primitive, -1, dims=-2) - primitive
            )
            / self.primitive_std,
            dim=-1,
        )
        # The exact zero jump factor preserves F(U,U)=f(U), while otherwise
        # leaving the learned vector correction completely unconstrained.
        return (
            self.base_flux(state)
            + 0.18 * jump[..., None] * correction * self.flux_scale
        )

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        dt = (
            euler.base.DT_SNAPSHOT
            if self.system == "euler"
            else swe.DT_SNAPSHOT
        )
        lam = dt * state.shape[-2]
        flux = self.flux(state)
        return state - lam * (flux - torch.roll(flux, 1, dims=-2))


def invalid_sentinel(system: str, state: torch.Tensor) -> torch.Tensor:
    sentinel = torch.zeros_like(state)
    sentinel[..., 0] = -1.0
    return sentinel


@torch.no_grad()
def raw_rollout(
    system: str,
    model: UnconstrainedLearnedFluxFV,
    initial: torch.Tensor,
    snapshots: int,
) -> torch.Tensor:
    """Freeze a failed trajectory at a finite invalid marker, never repair it."""
    state = initial.clone()
    active = common.admissible(system, state).all(dim=-1)
    saved = [state.clone()]
    for _ in range(1, snapshots):
        next_state = state.clone()
        if bool(active.any()):
            proposal = model(state[active])
            finite = torch.isfinite(proposal).all(dim=(-1, -2))
            physical = common.admissible(system, proposal).all(dim=-1)
            survives = finite & physical
            replacement = torch.where(
                finite[:, None, None],
                proposal,
                invalid_sentinel(system, proposal),
            )
            next_state[active] = replacement
            active_indices = torch.nonzero(active, as_tuple=False).squeeze(-1)
            active[active_indices] = survives
        state = next_state
        saved.append(state.clone())
    return torch.stack(saved, dim=1)


@torch.no_grad()
def rollout_metrics(
    system: str,
    model: UnconstrainedLearnedFluxFV,
    data: torch.Tensor,
    state_std: torch.Tensor,
) -> dict[str, Any]:
    model.eval()
    trajectory = raw_rollout(system, model, data[:, 0], data.shape[1])
    alive = common.admissible(system, trajectory).all(dim=(-1, -2))
    normalized = (trajectory[:, 1:] - data[:, 1:]) / state_std
    per_trajectory = normalized.double().square().mean(dim=(1, 2, 3))
    completed = int(alive.sum())
    primary = float(trajectory[..., 0].min())
    pressure: float | None = (
        float(euler.base.t_pressure(trajectory).min())
        if system == "euler"
        else None
    )
    return {
        "rollout_nrmse_completed_only": (
            float(torch.sqrt(per_trajectory[alive].mean()))
            if completed
            else float("inf")
        ),
        "rollout_nrmse_all_finite_outputs": float(
            torch.sqrt(per_trajectory.mean())
        ),
        "completion_rate": completed / data.shape[0],
        "completed_trajectories": float(completed),
        "minimum_density_or_depth": primary,
        "minimum_pressure": pressure,
    }


def train_one(
    system: str,
    train_data: torch.Tensor,
    validation_data: torch.Tensor,
    args: argparse.Namespace,
    output: Path,
) -> dict[str, Any]:
    checkpoint = output / (
        f"{system}_learned_flux64_converged_best_seed{args.seed}.pt"
    )
    report_path = output / f"{system}_learned_flux64_report_seed{args.seed}.json"
    if args.resume and checkpoint.exists() and report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))

    primitives = (
        euler.base.primitive(train_data)
        if system == "euler"
        else swe.primitive(train_data)
    )
    primitive_mean = primitives.mean(dim=(0, 1, 2))
    primitive_std = primitives.std(dim=(0, 1, 2))
    state_std = train_data.std(dim=(0, 1, 2))
    torch.manual_seed(61000 + args.seed + (0 if system == "euler" else 1000))
    model = UnconstrainedLearnedFluxFV(
        system, primitive_mean, primitive_std, args.width
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(
        62000 + args.seed + (0 if system == "euler" else 1000)
    )
    best_metric = float("inf")
    best_update = -1
    best_state: dict[str, torch.Tensor] | None = None
    plateau_anchor = float("inf")
    stale = 0
    converged = False
    stop_reason = "max_updates"
    curve: list[dict[str, Any]] = []
    started = time.perf_counter()

    def validate(update: int, batch_loss: float | None) -> float:
        nonlocal best_metric, best_update, best_state
        metrics = rollout_metrics(
            system, model, validation_data, state_std
        )
        selection_metric = float(
            metrics["rollout_nrmse_all_finite_outputs"]
        )
        improved = selection_metric < best_metric
        if improved:
            best_metric = selection_metric
            best_update = update
            best_state = common.clone_state_dict(model)
            torch.save(best_state, checkpoint)
        row = {
            "system": system,
            "seed": args.seed,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "validation_rollout_nrmse_completed_only": metrics[
                "rollout_nrmse_completed_only"
            ],
            "validation_rollout_nrmse_all_finite_outputs": selection_metric,
            "validation_completion_rate": metrics["completion_rate"],
            "selection_metric": selection_metric,
            "new_absolute_best": improved,
        }
        curve.append(row)
        print(json.dumps(row), flush=True)
        return selection_metric

    validate(0, None)
    last_loss: float | None = None
    for update in range(1, args.max_updates + 1):
        model.train()
        indices = torch.randint(
            0, train_data.shape[0], (args.batch_size,), generator=generator
        )
        times = torch.randint(
            0,
            train_data.shape[1] - 1,
            (args.batch_size,),
            generator=generator,
        )
        inputs = train_data[indices, times]
        targets = train_data[indices, times + 1]
        prediction = model(inputs)
        loss = (((prediction - targets) / state_std) ** 2).mean()
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach())

        if update % args.validation_interval and update != args.max_updates:
            continue
        metric = validate(update, last_loss)
        if update < args.fixed_budget_updates:
            continue
        if plateau_anchor == float("inf"):
            plateau_anchor = best_metric
            stale = 0
            continue
        if metric <= plateau_anchor * (1.0 - args.minimum_relative_improvement):
            plateau_anchor = metric
            stale = 0
        else:
            stale += 1
        if stale < args.plateau_patience:
            continue
        current_lr = float(optimizer.param_groups[0]["lr"])
        if current_lr <= args.minimum_learning_rate * (1.0 + 1.0e-12):
            converged = True
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_lr = max(
            current_lr * args.learning_rate_factor,
            args.minimum_learning_rate,
        )
        for group in optimizer.param_groups:
            group["lr"] = new_lr
        plateau_anchor = best_metric
        stale = 0
        print(
            json.dumps(
                {
                    "system": system,
                    "event": "reduce_learning_rate",
                    "update": update,
                    "old_learning_rate": current_lr,
                    "new_learning_rate": new_lr,
                }
            ),
            flush=True,
        )

    if best_state is None:
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(f"{system} learned-flux baseline made no checkpoint")
    if not converged:
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"{system} learned-flux baseline did not converge by cap"
        )
    model.load_state_dict(best_state)
    final_metrics = rollout_metrics(system, model, validation_data, state_std)
    report = {
        "system": system,
        "baseline": "unconstrained conservative local learned flux",
        "scope": (
            "same-data learned-neural-FVM family baseline; not an exact "
            "reproduction of a published architecture"
        ),
        "training_cells": int(train_data.shape[-2]),
        "training_data_shape": list(train_data.shape),
        "validation_data_shape": list(validation_data.shape),
        "parameter_count": sum(p.numel() for p in model.parameters()),
        "stencil_cells": len(STENCIL_SHIFTS),
        "width": args.width,
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "best_validation_rollout_nrmse": best_metric,
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "validation_metrics": final_metrics,
    }
    common.write_csv(
        output / f"{system}_learned_flux64_curve_seed{args.seed}.csv",
        curve,
    )
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--systems",
        nargs="+",
        choices=("euler", "swe"),
        default=["euler", "swe"],
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=1200)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--validation-interval", type=int, default=200)
    parser.add_argument("--plateau-patience", type=int, default=8)
    parser.add_argument("--minimum-relative-improvement", type=float, default=1.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def self_test() -> None:
    specifications = (
        (
            "euler",
            torch.tensor([1.0, 0.0, 1.0]),
            torch.tensor([1.0, 1.0, 1.0]),
            torch.tensor([1.2, 0.24, 2.52]),
        ),
        (
            "swe",
            torch.tensor([1.0, 0.0]),
            torch.tensor([1.0, 1.0]),
            torch.tensor([1.2, 0.24]),
        ),
    )
    for system, mean, std, value in specifications:
        model = UnconstrainedLearnedFluxFV(system, mean, std, width=16)
        nn.init.normal_(model.output.weight)
        nn.init.normal_(model.output.bias)
        constant = value.reshape(1, 1, -1).expand(2, 64, -1).clone()
        flux = model.flux(constant)
        expected = (
            euler.base.t_flux(constant)
            if system == "euler"
            else swe.t_flux(constant)
        )
        if float((flux - expected).detach().abs().max()) > 2.0e-6:
            raise RuntimeError(f"{system} learned flux lost equal-state consistency")
        updated = model(constant)
        if float(
            (updated.sum(dim=-2) - constant.sum(dim=-2))
            .detach()
            .abs()
            .max()
        ) > 2.0e-5:
            raise RuntimeError(f"{system} learned FV update lost conservation")
    print("self-test passed")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reports: dict[str, Any] = {}
    if "euler" in args.systems:
        train, validation = replicates.euler_cache(output, args.seed)
        reports["euler"] = train_one(
            "euler", train, validation, args, output
        )
    if "swe" in args.systems:
        train, validation = swe.load_or_make_data(SWE_DIR / "results", args.seed)
        reports["swe"] = train_one(
            "swe", train, validation, args, output
        )
    print(json.dumps(reports, indent=2), flush=True)


if __name__ == "__main__":
    main()
