"""Train matched periodic Fourier-neural-operator baselines on 64-cell data.

This is a same-data architecture baseline, not a claim that this compact
implementation reproduces every detail of the original FNO paper.  It learns
the one-snapshot solution operator directly and is rolled out
autoregressively.  It has no finite-volume flux form, positivity layer, or
entropy projection.
"""

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
from torch import nn
import torch.nn.functional as functional


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


class SpectralConv1d(nn.Module):
    def __init__(self, channels: int, modes: int) -> None:
        super().__init__()
        self.channels = channels
        self.modes = modes
        scale = 1.0 / channels
        self.weight = nn.Parameter(
            scale
            * torch.randn(
                channels,
                channels,
                modes,
                dtype=torch.cfloat,
            )
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        cells = values.shape[-1]
        transformed = torch.fft.rfft(values, dim=-1)
        kept = min(self.modes, transformed.shape[-1])
        output = torch.zeros(
            values.shape[0],
            self.channels,
            transformed.shape[-1],
            dtype=torch.cfloat,
            device=values.device,
        )
        output[..., :kept] = torch.einsum(
            "bim,iom->bom",
            transformed[..., :kept],
            self.weight[..., :kept],
        )
        return torch.fft.irfft(output, n=cells, dim=-1)


class ResidualFNO1d(nn.Module):
    """A compact periodic one-step FNO with an identity initialization."""

    def __init__(
        self,
        channels: int,
        state_mean: torch.Tensor,
        state_std: torch.Tensor,
        width: int = 24,
        modes: int = 12,
        layers: int = 4,
    ) -> None:
        super().__init__()
        self.channels = channels
        self.width = width
        self.modes = modes
        self.layers = layers
        self.lift = nn.Conv1d(channels, width, 1)
        self.spectral = nn.ModuleList(
            [SpectralConv1d(width, modes) for _ in range(layers)]
        )
        self.local = nn.ModuleList(
            [nn.Conv1d(width, width, 1) for _ in range(layers)]
        )
        self.project_hidden = nn.Conv1d(width, width, 1)
        self.project = nn.Conv1d(width, channels, 1)
        nn.init.zeros_(self.project.weight)
        nn.init.zeros_(self.project.bias)
        self.register_buffer(
            "state_mean", state_mean.reshape(1, 1, channels).float()
        )
        self.register_buffer(
            "state_std", state_std.reshape(1, 1, channels).float()
        )

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        normalized = (state - self.state_mean) / self.state_std
        values = self.lift(normalized.transpose(1, 2))
        for index, (spectral, local) in enumerate(
            zip(self.spectral, self.local)
        ):
            values = spectral(values) + local(values)
            if index + 1 != self.layers:
                values = functional.gelu(values)
        delta = self.project(
            functional.gelu(self.project_hidden(values))
        ).transpose(1, 2)
        return state + delta * self.state_std


def clone_state_dict(model: nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


def admissible(system: str, state: torch.Tensor) -> torch.Tensor:
    finite = torch.isfinite(state).all(dim=-1)
    if system == "euler":
        return (
            finite
            & (state[..., 0] >= 1.0e-5)
            & (euler.base.t_pressure(state) >= 1.0e-5)
        )
    if system == "swe":
        return finite & (state[..., 0] >= swe.H_FLOOR)
    raise ValueError(system)


@torch.no_grad()
def rollout_metrics(
    system: str,
    model: ResidualFNO1d,
    data: torch.Tensor,
    state_std: torch.Tensor,
    batch_size: int = 64,
) -> dict[str, Any]:
    model.eval()
    total_squared = 0.0
    total_values = 0
    all_squared = 0.0
    all_values = 0
    all_finite = True
    completed = 0
    minimum_primary = float("inf")
    minimum_pressure = float("inf")
    for start in range(0, data.shape[0], batch_size):
        batch = data[start : start + batch_size]
        state = batch[:, 0].clone()
        alive = admissible(system, state).all(dim=-1)
        squared = torch.zeros(
            state.shape[0], dtype=torch.float64, device=state.device
        )
        values_per_trajectory = 0
        for snapshot in range(1, batch.shape[1]):
            state = model(state)
            all_finite = all_finite and bool(torch.isfinite(state).all())
            good = admissible(system, state).all(dim=-1)
            alive = alive & good
            finite_state = torch.where(
                torch.isfinite(state), state, torch.zeros_like(state)
            )
            error = (finite_state - batch[:, snapshot]) / state_std
            squared += error.double().square().sum(dim=(1, 2))
            values_per_trajectory += error.shape[1] * error.shape[2]
            minimum_primary = min(
                minimum_primary, float(finite_state[..., 0].min())
            )
            if system == "euler":
                minimum_pressure = min(
                    minimum_pressure,
                    float(euler.base.t_pressure(finite_state).min()),
                )
        if bool(alive.any()):
            total_squared += float(squared[alive].sum())
            total_values += int(alive.sum()) * values_per_trajectory
        all_squared += float(squared.sum())
        all_values += state.shape[0] * values_per_trajectory
        completed += int(alive.sum())
    return {
        "rollout_nrmse_completed_only": (
            float(np.sqrt(total_squared / total_values))
            if total_values
            else float("inf")
        ),
        "rollout_nrmse_all_finite_outputs": (
            float(np.sqrt(all_squared / all_values))
            if all_finite and all_values
            else float("inf")
        ),
        "completion_rate": completed / data.shape[0],
        "completed_trajectories": float(completed),
        "minimum_density_or_depth": minimum_primary,
        "minimum_pressure": minimum_pressure if system == "euler" else None,
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


def train_one(
    system: str,
    train_data: torch.Tensor,
    validation_data: torch.Tensor,
    args: argparse.Namespace,
    output: Path,
) -> dict[str, Any]:
    checkpoint = output / f"{system}_fno64_converged_best_seed{args.seed}.pt"
    report_path = output / f"{system}_fno64_report_seed{args.seed}.json"
    if args.resume and checkpoint.exists() and report_path.exists():
        return json.loads(report_path.read_text(encoding="utf-8"))

    state_mean = train_data.mean(dim=(0, 1, 2))
    state_std = train_data.std(dim=(0, 1, 2))
    torch.manual_seed(51000 + args.seed + (0 if system == "euler" else 1000))
    model = ResidualFNO1d(
        train_data.shape[-1],
        state_mean,
        state_std,
        width=args.width,
        modes=args.modes,
        layers=args.layers,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(
        52000 + args.seed + (0 if system == "euler" else 1000)
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
        metrics = rollout_metrics(system, model, validation_data, state_std)
        raw_metric = float(metrics["rollout_nrmse_completed_only"])
        selection_metric = float(metrics["rollout_nrmse_all_finite_outputs"])
        improved = selection_metric < best_metric
        if improved:
            best_metric = selection_metric
            best_update = update
            best_state = clone_state_dict(model)
            torch.save(best_state, checkpoint)
        row = {
            "system": system,
            "seed": args.seed,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "validation_rollout_nrmse_completed_only": raw_metric,
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
        raise RuntimeError(f"{system} FNO never produced an admissible checkpoint")
    if not converged:
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"{system} FNO reached {args.max_updates} updates without convergence"
        )
    model.load_state_dict(best_state)
    final_metrics = rollout_metrics(system, model, validation_data, state_std)
    report = {
        "system": system,
        "baseline": "compact periodic residual FNO",
        "scope": (
            "same-data architecture baseline; not an exact reproduction of "
            "the original FNO paper"
        ),
        "training_cells": int(train_data.shape[-2]),
        "training_data_shape": list(train_data.shape),
        "validation_data_shape": list(validation_data.shape),
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "width": args.width,
        "modes": args.modes,
        "layers": args.layers,
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "best_validation_rollout_nrmse": best_metric,
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "validation_metrics": final_metrics,
    }
    write_csv(output / f"{system}_fno64_curve_seed{args.seed}.csv", curve)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--systems", nargs="+", choices=("euler", "swe"), default=["euler", "swe"])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=24)
    parser.add_argument("--modes", type=int, default=12)
    parser.add_argument("--layers", type=int, default=4)
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reports: dict[str, Any] = {}
    if "euler" in args.systems:
        train, validation = replicates.euler_cache(output, args.seed)
        reports["euler"] = train_one(
            "euler", train, validation, args, output
        )
    if "swe" in args.systems:
        train, validation = swe.load_or_make_data(
            SWE_DIR / "results", args.seed
        )
        reports["swe"] = train_one("swe", train, validation, args, output)
    print(json.dumps(reports, indent=2), flush=True)


if __name__ == "__main__":
    main()
