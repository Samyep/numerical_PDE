"""Train and audit the retained Euler HCFL method directly on 512 cells.

The learned proposal is the symmetric four-cell Central + nonnegative Roe
architecture with proposal-feasibility loss.  Training and validation targets
are native periodic HLLC-2048 trajectories conservatively restricted to 512
cell averages.  Checkpoint selection uses only independent validation rollout
NRMSE through the same full hard-safety stack used at deployment.  Held-out
periodic and transmissive cases are evaluated against native HLLC-2048 and a
zero-network Roe-512 control using that same stack.  A separate projection-only
rollout removes F_low to measure what the low-order safety endpoint contributes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
STENCIL_EXPERIMENT = HCFL_ROOT / "experiments" / "euler_1d_stencil_ablation"
if str(STENCIL_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(STENCIL_EXPERIMENT))

import run_stencil_ablation as stencil  # noqa: E402


base = stencil.base
shared = stencil.shared
periodic_eval = stencil.periodic_eval
nonperiodic = stencil.nonperiodic

REFERENCE_CELLS = 2048
TRAIN_CELLS = 512
RESTRICTION_FACTOR = REFERENCE_CELLS // TRAIN_CELLS
INTERVAL_SUBSTEPS = TRAIN_CELLS // base.NCOARSE
TRAIN_LAMBDA = base.DT_SNAPSHOT * base.NCOARSE
ARM = stencil.Arm(
    name="central_nonnegative_feas_s4_train512",
    family="central_nonnegative_feas",
    stencil_size=4,
    model_name="central_roe_upwind_4",
    feasibility_weight=1.0e-3,
)

PERIODIC_CASES = {
    "sod": "Sod",
    "lax": "Lax",
    "collision": "Collision",
    "strong_pressure": "Strong pressure",
}
NONPERIODIC_CASES = {
    "near_vacuum_expansion": "Near-vacuum",
    "contact_left_exit": "Contact exits left",
    "contact_right_exit": "Contact exits right",
    "pressure_right_exit": "Pressure wave exits right",
}
VARIABLES = ((0, r"Density $\rho$"), (1, r"Velocity $u$"), (2, r"Pressure $p$"))


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
        raise ValueError(f"No rows for {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def device_from_argument(name: str) -> torch.device:
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    return device


def admissible(state: torch.Tensor) -> torch.Tensor:
    return (state[..., 0] > 0.0) & (base.t_pressure(state) > 0.0)


def assert_admissible(state: torch.Tensor, context: str) -> None:
    if not bool(torch.isfinite(state).all()) or not bool(admissible(state).all()):
        raise RuntimeError(f"Strict HLLC-2048 reference failed during {context}")


def restrict_reference(state: torch.Tensor) -> torch.Tensor:
    batch = state.shape[0]
    return state.reshape(batch, TRAIN_CELLS, RESTRICTION_FACTOR, 3).mean(dim=2)


@torch.no_grad()
def generate_reference_batch(
    initial: np.ndarray,
    device: torch.device,
    dtype: torch.dtype = torch.float32,
) -> torch.Tensor:
    """Generate strict periodic HLLC-2048 + SSP-RK2 trajectories."""
    state = torch.as_tensor(initial, dtype=dtype, device=device)
    dx = 1.0 / REFERENCE_CELLS
    saved = [restrict_reference(state).cpu()]
    for snapshot in range(1, base.NSNAP):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            density = state[..., 0]
            pressure = base.t_pressure(state)
            velocity = state[..., 1] / density
            sound = torch.sqrt(base.GAMMA * pressure / density)
            maximum_speed = float((torch.abs(velocity) + sound).max().item())
            dt = min(remaining, 0.20 * dx / max(maximum_speed, 1.0e-12))
            lam = dt / dx
            flux = base.t_hllc(state)
            stage_one = state - lam * (
                flux - torch.roll(flux, 1, dims=-2)
            )
            assert_admissible(
                stage_one,
                f"snapshot {snapshot}, SSP-RK2 stage one",
            )
            stage_flux = base.t_hllc(stage_one)
            stage_two = stage_one - lam * (
                stage_flux - torch.roll(stage_flux, 1, dims=-2)
            )
            state = 0.5 * state + 0.5 * stage_two
            assert_admissible(
                state,
                f"snapshot {snapshot}, SSP-RK2 completion",
            )
            remaining -= dt
        saved.append(restrict_reference(state).cpu())
    return torch.stack(saved, dim=1)


def initial_condition_groups(split: str, seed: int) -> list[tuple[str, np.ndarray]]:
    if split == "train":
        return [
            (
                "ordinary",
                base.generate_ic(220, REFERENCE_CELLS, 6000 + seed, ood=False),
            ),
            (
                "broad",
                base.generate_ic(260, REFERENCE_CELLS, 7000 + seed, ood=True),
            ),
            (
                "extreme",
                base.generate_extreme_ic(100, REFERENCE_CELLS, 8000 + seed),
            ),
        ]
    if split == "validation":
        return [
            (
                "ordinary",
                base.generate_ic(44, REFERENCE_CELLS, 14000 + seed, ood=False),
            ),
            (
                "broad",
                base.generate_ic(52, REFERENCE_CELLS, 15000 + seed, ood=True),
            ),
            (
                "extreme",
                base.generate_extreme_ic(20, REFERENCE_CELLS, 16000 + seed),
            ),
            (
                "wave",
                shared.generate_wave_coverage_ic(
                    20, REFERENCE_CELLS, 17000 + seed
                ),
            ),
        ]
    raise ValueError(split)


def generate_split(
    split: str,
    seed: int,
    device: torch.device,
    batch_size: int,
) -> torch.Tensor:
    pieces: list[torch.Tensor] = []
    for group_name, initial in initial_condition_groups(split, seed):
        for start in range(0, initial.shape[0], batch_size):
            stop = min(start + batch_size, initial.shape[0])
            started = time.perf_counter()
            trajectory = generate_reference_batch(
                initial[start:stop], device
            ).float()
            pieces.append(trajectory)
            print(
                json.dumps(
                    {
                        "stage": "reference_generation",
                        "split": split,
                        "group": group_name,
                        "start": start,
                        "stop": stop,
                        "seconds": time.perf_counter() - started,
                    }
                ),
                flush=True,
            )
    return torch.cat(pieces, dim=0)


def reference_precision_audit(
    seed: int,
    device: torch.device,
) -> dict[str, float]:
    sample = base.generate_ic(
        2, REFERENCE_CELLS, 14000 + seed, ood=False
    )
    gpu = generate_reference_batch(sample, device, torch.float32).double()
    cpu = generate_reference_batch(
        sample, torch.device("cpu"), torch.float64
    )
    error = gpu - cpu
    scale = cpu.std(dim=(0, 1, 2)).clamp_min(1.0e-12)
    normalized = error / scale
    return {
        "trajectories": 2,
        "maximum_absolute_state_difference": float(error.abs().max()),
        "rms_state_difference": float(torch.sqrt((error * error).mean())),
        "normalized_rms_state_difference": float(
            torch.sqrt((normalized * normalized).mean())
        ),
    }


def load_or_generate_data(
    output: Path,
    seed: int,
    device: torch.device,
    reference_batch_size: int,
) -> tuple[torch.Tensor, torch.Tensor, dict[str, float]]:
    cache = output / f"data_cache_train512_seed{seed}.pt"
    expected_train = (580, base.NSNAP, TRAIN_CELLS, 3)
    expected_validation = (136, base.NSNAP, TRAIN_CELLS, 3)
    if cache.exists():
        payload = torch.load(cache, map_location="cpu", weights_only=True)
        if (
            tuple(payload["train"].shape) == expected_train
            and tuple(payload["validation"].shape) == expected_validation
            and int(payload["reference_cells"]) == REFERENCE_CELLS
        ):
            print(f"Loading deterministic 512-cell cache: {cache}", flush=True)
            return (
                payload["train"],
                payload["validation"],
                payload["precision_audit"],
            )
        print(f"Ignoring stale cache: {cache}", flush=True)

    audit = reference_precision_audit(seed, device)
    print(json.dumps({"reference_precision_audit": audit}), flush=True)
    train = generate_split("train", seed, device, reference_batch_size)
    validation = generate_split(
        "validation", seed, device, reference_batch_size
    )
    torch.save(
        {
            "train": train,
            "validation": validation,
            "reference_cells": REFERENCE_CELLS,
            "target_cells": TRAIN_CELLS,
            "precision_audit": audit,
        },
        cache,
    )
    return train, validation, audit


def prepare_statistics(
    train_data: torch.Tensor,
) -> tuple[np.ndarray, np.ndarray, torch.Tensor]:
    primitive = base.primitive(train_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))
    return mean, std, state_std


def interval_forward(
    model: base.Solver,
    state: torch.Tensor,
    collect_feasibility: bool,
) -> tuple[torch.Tensor, torch.Tensor]:
    feasibility = torch.zeros((), dtype=state.dtype, device=state.device)
    for _ in range(INTERVAL_SUBSTEPS):
        raw_flux = model.flux_net(state)
        if collect_feasibility:
            raw_residual = base.entropy_residual(raw_flux, state)
            feasibility = feasibility + torch.relu(raw_residual).square().mean()
        projected = shared.strict_entropy_projection(raw_flux, state)
        state = state - TRAIN_LAMBDA * (
            projected - torch.roll(projected, 1, dims=-2)
        )
        if not collect_feasibility and (
            not bool(torch.isfinite(state).all())
            or not bool(admissible(state).all())
        ):
            break
    return state, feasibility / INTERVAL_SUBSTEPS


@torch.no_grad()
def validation_rollout(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    """Autoregressive validation through the actual deployed safe solver.

    ``advance_safe_snapshot`` uses the original 64-cell virtual spacing.  On a
    512-cell state one call therefore represents 1/8 of a saved physical
    interval; eight calls recover the 512-cell time scale while preserving the
    exact production limiter and CFL logic.
    """
    scale = state_std.to(device)
    squared_error_sum = 0.0
    squared_error_count = 0
    minimum_density = float("inf")
    minimum_pressure = float("inf")
    local_active = 0.0
    local_total = 0.0
    fd_active = 0.0
    fd_total = 0.0
    fd_beta_sum = 0.0
    fd_beta_min = 1.0
    # Keep the deterministic validation strata separate.  A batched adaptive
    # solver uses the largest wave speed in the batch; mixing a near-vacuum
    # extreme with every ordinary trajectory therefore changes all their time
    # grids and needlessly multiplies CPU work.  These boundaries exactly match
    # ``initial_condition_groups("validation", ...)`` above.
    if data.shape[0] == 136:
        batches = ((0, 44), (44, 96), (96, 116), (116, 136))
    else:
        batches = tuple(
            (start, min(start + 32, data.shape[0]))
            for start in range(0, data.shape[0], 32)
        )

    for start, stop in batches:
        state = data[start:stop, 0].to(device)
        minimum_density = min(minimum_density, float(state[..., 0].min()))
        minimum_pressure = min(
            minimum_pressure, float(base.t_pressure(state).min())
        )
        for snapshot in range(1, data.shape[1]):
            for _ in range(INTERVAL_SUBSTEPS):
                state, step = shared.advance_safe_snapshot(model, state)
                local_active += step["local_active"]
                local_total += step["local_total"]
                fd_active += step["fd_active"]
                fd_total += step["fd_total"]
                fd_beta_sum += step["fd_beta_sum"]
                fd_beta_min = min(fd_beta_min, step["fd_beta_min"])
            target = data[start:stop, snapshot].to(device)
            error = (state - target) / scale
            squared_error_sum += float((error * error).sum())
            squared_error_count += error.numel()
            minimum_density = min(
                minimum_density, float(state[..., 0].min())
            )
            minimum_pressure = min(
                minimum_pressure, float(base.t_pressure(state).min())
            )
    return {
        "rollout_nrmse": float(
            np.sqrt(squared_error_sum / max(squared_error_count, 1))
        ),
        "minimum_density": minimum_density,
        "minimum_pressure": minimum_pressure,
        "local_limiter_intervention_rate": local_active / max(local_total, 1.0),
        "fd_entropy_intervention_rate": fd_active / max(fd_total, 1.0),
        "mean_fd_beta": fd_beta_sum / max(fd_total, 1.0),
        "minimum_fd_beta": fd_beta_min,
    }


@torch.no_grad()
def teacher_forced_validation_loss(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    device: torch.device,
    pairs: int = 256,
    batch_size: int = 32,
) -> float:
    inputs = data[:, :-1].reshape(-1, TRAIN_CELLS, 3)
    targets = data[:, 1:].reshape(-1, TRAIN_CELLS, 3)
    indices = torch.linspace(0, inputs.shape[0] - 1, pairs).long()
    scale = state_std.to(device)
    total = 0.0
    count = 0
    for start in range(0, pairs, batch_size):
        chosen = indices[start : start + batch_size]
        prediction, _ = interval_forward(
            model, inputs[chosen].to(device), False
        )
        error = (prediction - targets[chosen].to(device)) / scale
        total += float((error * error).sum())
        count += error.numel()
    return total / count


def clone_state_dict(model: torch.nn.Module) -> dict[str, torch.Tensor]:
    return {
        key: value.detach().cpu().clone()
        for key, value in model.state_dict().items()
    }


def train_to_convergence(
    train_data: torch.Tensor,
    validation_data: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    state_std: torch.Tensor,
    args: argparse.Namespace,
    output: Path,
    device: torch.device,
) -> tuple[base.Solver, dict[str, Any], list[dict[str, Any]]]:
    torch.manual_seed(51000 + args.seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(51000 + args.seed)
    model = base.Solver(
        ARM.model_name, mean, std, width=args.width
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(52000 + args.seed)
    scale = state_std.to(device)
    best_metric = float("inf")
    best_update = -1
    best_state: dict[str, torch.Tensor] | None = None
    plateau_anchor = float("inf")
    checks_without_improvement = 0
    converged = False
    stop_reason = "max_updates"
    curve: list[dict[str, Any]] = []
    checkpoint = output / f"{ARM.name}_converged_best_seed{args.seed}.pt"
    started = time.perf_counter()

    def validate(update: int, batch_loss: float | None) -> float:
        nonlocal best_metric, best_update, best_state
        model.eval()
        teacher_loss = teacher_forced_validation_loss(
            model, validation_data, state_std, device
        )
        rollout = validation_rollout(
            model, validation_data, state_std, device
        )
        metric = rollout["rollout_nrmse"]
        improved = math.isfinite(metric) and metric < best_metric
        if improved:
            best_metric = metric
            best_update = update
            best_state = clone_state_dict(model)
            torch.save(best_state, checkpoint)
        row = {
            "seed": args.seed,
            "arm": ARM.name,
            "training_cells": TRAIN_CELLS,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": batch_loss,
            "validation_teacher_forced_loss": teacher_loss,
            "validation_rollout_nrmse": metric,
            "validation_minimum_density": rollout["minimum_density"],
            "validation_minimum_pressure": rollout["minimum_pressure"],
            "validation_local_limiter_intervention_rate": rollout[
                "local_limiter_intervention_rate"
            ],
            "validation_fd_entropy_intervention_rate": rollout[
                "fd_entropy_intervention_rate"
            ],
            "validation_mean_fd_beta": rollout["mean_fd_beta"],
            "validation_minimum_fd_beta": rollout["minimum_fd_beta"],
            "new_absolute_best": improved,
        }
        curve.append(row)
        print(json.dumps(row), flush=True)
        return metric

    validate(0, None)
    last_loss: float | None = None
    for update in range(1, args.max_updates + 1):
        model.train()
        indices = torch.randint(
            0, train_data.shape[0], (args.batch_size,), generator=generator
        )
        times = torch.randint(
            0, train_data.shape[1] - 1, (args.batch_size,), generator=generator
        )
        inputs = train_data[indices, times].to(device)
        targets = train_data[indices, times + 1].to(device)
        prediction, feasibility = interval_forward(model, inputs, True)
        trajectory_loss = (((prediction - targets) / scale) ** 2).mean()
        loss = trajectory_loss + ARM.feasibility_weight * feasibility
        if not bool(torch.isfinite(loss)):
            raise RuntimeError(f"Nonfinite training loss at update {update}")

        optimizer.zero_grad(set_to_none=True)
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
            checks_without_improvement = 0
            continue
        if metric <= plateau_anchor * (1.0 - args.minimum_relative_improvement):
            plateau_anchor = metric
            checks_without_improvement = 0
        else:
            checks_without_improvement += 1
        if checks_without_improvement < args.plateau_patience:
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
        checks_without_improvement = 0
        plateau_anchor = best_metric
        print(
            json.dumps(
                {
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
        raise RuntimeError("No finite direct-512 validation checkpoint was produced")
    if not converged:
        checkpoint.unlink(missing_ok=True)
        raise RuntimeError(
            f"Direct-512 Euler training reached the {args.max_updates}-update "
            "cap without satisfying the validation plateau rule"
        )
    model.load_state_dict(best_state)
    model.eval()
    final_validation = validation_rollout(
        model, validation_data, state_std, device
    )
    summary = {
        "seed": args.seed,
        "arm": ARM.name,
        "model": ARM.model_name,
        "stencil_cells": ARM.stencil_size,
        "training_cells": TRAIN_CELLS,
        "reference_cells": REFERENCE_CELLS,
        "interval_substeps": INTERVAL_SUBSTEPS,
        "substep_dt_over_dx": TRAIN_LAMBDA,
        "physical_supervision_interval": base.DT_SNAPSHOT,
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "best_update": best_update,
        "stop_update": int(curve[-1]["update"]),
        "best_validation_rollout_nrmse": best_metric,
        "best_checkpoint_validation": final_validation,
        "final_learning_rate": float(optimizer.param_groups[0]["lr"]),
        "converged": converged,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "proposal_feasibility_weight": ARM.feasibility_weight,
        "device": str(device),
    }
    model.cpu()
    return model, summary, curve


def load_model(
    output: Path,
    seed: int,
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
) -> base.Solver:
    checkpoint = output / f"{ARM.name}_converged_best_seed{seed}.pt"
    model = base.Solver(ARM.model_name, mean, std, width=width)
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


@torch.no_grad()
def multiplier_activity(
    model: base.Solver,
    data: torch.Tensor,
    pairs: int = 256,
) -> dict[str, float]:
    states = data[:, :-1].reshape(-1, TRAIN_CELLS, 3)
    indices = torch.linspace(0, states.shape[0] - 1, pairs).long()
    states = states[indices]
    primitive = base.primitive(states)
    flux_net = model.flux_net
    features = torch.cat(
        [
            (torch.roll(primitive, shift, dims=-2) - flux_net.mean)
            / flux_net.std
            for shift in flux_net.stencil_shifts
        ],
        dim=-1,
    )
    multipliers = 1.0 + torch.tanh(flux_net.net(features))
    return {
        "mean_absolute_multiplier_change_from_one": float(
            torch.mean(torch.abs(multipliers - 1.0))
        ),
        "minimum_multiplier": float(multipliers.min()),
        "maximum_multiplier": float(multipliers.max()),
        "mean_multiplier": float(multipliers.mean()),
    }


def method_metrics(
    boundary: str,
    case: str,
    method: str,
    reference: torch.Tensor,
    prediction: torch.Tensor,
    state_std: torch.Tensor,
    periodic: bool,
    run_stats: dict[str, Any],
) -> dict[str, Any]:
    metrics: dict[str, Any] = {
        "boundary": boundary,
        "case": case,
        "method": method,
        **periodic_eval.comparison.diagnostics(
            reference, prediction, state_std
        ),
        **stencil.oscillation_metrics(reference, prediction, periodic),
    }
    integrity = (
        periodic_eval.integrity_metrics(prediction)
        if periodic
        else nonperiodic.integrity_metrics(prediction)
    )
    metrics.update(integrity)
    metrics.update(
        {
            f"run_{key}": value
            for key, value in run_stats.items()
            if key != "seconds"
        }
    )
    return metrics


@torch.no_grad()
def projection_only_rollout(
    model: base.Solver,
    name: str,
    periodic: bool,
    cells: int = TRAIN_CELLS,
    cfl: float = 0.42,
) -> tuple[torch.Tensor | None, dict[str, Any]]:
    """Deploy hard interface projection without either F_low blend.

    This intentionally removes both the local positivity blend and the global
    fully-discrete entropy blend.  It is an empirical simplification ablation,
    not a hard-positive method.
    """
    if periodic:
        state = torch.from_numpy(
            periodic_eval.comparison.initial_condition(name, cells)
        ).float()
    else:
        state = nonperiodic.initial_condition(name, cells).float()
    dx = 1.0 / cells
    nominal_dt = base.DT_SNAPSHOT * base.NCOARSE / cells
    snapshots = [state.clone()]
    maximum_internal_substeps = 50_000
    maximum_substeps_per_snapshot = 4_096
    minimum_practical_dt = (
        base.DT_SNAPSHOT / maximum_substeps_per_snapshot
    )
    stats: dict[str, Any] = {
        "completed": True,
        "substeps": 0,
        "minimum_density": float(state[..., 0].min()),
        "minimum_pressure": float(base.t_pressure(state).min()),
        "fully_discrete_entropy_violations": 0,
        "fully_discrete_entropy_checks": 0,
        "max_fully_discrete_entropy_increase": -float("inf"),
    }

    for snapshot in range(1, shared.CANONICAL_NSNAP):
        remaining = base.DT_SNAPSHOT
        snapshot_substeps = 0
        while remaining > 1.0e-14:
            density = state[..., 0].clamp_min(1.0e-10)
            pressure = base.t_pressure(state).clamp_min(1.0e-10)
            velocity = state[..., 1] / density
            sound = torch.sqrt(base.GAMMA * pressure / density)
            maximum_speed = float((torch.abs(velocity) + sound).max())
            dt = min(
                remaining,
                nominal_dt,
                cfl * dx / max(maximum_speed, 1.0e-12),
            )
            if (
                stats["substeps"] >= maximum_internal_substeps
                or snapshot_substeps >= maximum_substeps_per_snapshot
                or dt < minimum_practical_dt
            ):
                stats.update(
                    {
                        "completed": False,
                        "failure_reason": "cfl_timestep_collapse",
                        "failure_snapshot": snapshot,
                        "failure_substep": stats["substeps"] + 1,
                        "required_dt": dt,
                        "minimum_practical_dt": minimum_practical_dt,
                        "maximum_internal_substeps": maximum_internal_substeps,
                        "maximum_substeps_per_snapshot": (
                            maximum_substeps_per_snapshot
                        ),
                    }
                )
                checks = max(stats["fully_discrete_entropy_checks"], 1)
                stats["fully_discrete_entropy_violation_rate"] = (
                    stats["fully_discrete_entropy_violations"] / checks
                )
                return None, stats
            lam = dt / dx
            if periodic:
                raw_flux = model.flux_net(state)
                flux = shared.strict_entropy_projection(raw_flux, state)
                next_state = state - lam * (
                    flux - torch.roll(flux, 1, dims=-2)
                )
                entropy_target = base.entropy(state.double()).sum(dim=-1)
            else:
                flux, _ = nonperiodic.model_flux(model, state)
                next_state = nonperiodic.update_from_flux(state, flux, lam)
                entropy_target = nonperiodic.entropy_target(state, lam)
            entropy_after = base.entropy(next_state.double()).sum(dim=-1)
            entropy_increase = entropy_after - entropy_target
            stats["fully_discrete_entropy_violations"] += int(
                (entropy_increase > shared.FD_ENTROPY_TOLERANCE).sum()
            )
            stats["fully_discrete_entropy_checks"] += entropy_increase.numel()
            stats["max_fully_discrete_entropy_increase"] = max(
                stats["max_fully_discrete_entropy_increase"],
                float(entropy_increase.max()),
            )
            stats["minimum_density"] = min(
                stats["minimum_density"], float(next_state[..., 0].min())
            )
            stats["minimum_pressure"] = min(
                stats["minimum_pressure"],
                float(base.t_pressure(next_state).min()),
            )
            if (
                not bool(torch.isfinite(next_state).all())
                or not bool(admissible(next_state).all())
            ):
                stats.update(
                    {
                        "completed": False,
                        "failure_reason": "non_admissible_state",
                        "failure_snapshot": snapshot,
                        "failure_substep": stats["substeps"] + 1,
                    }
                )
                checks = max(stats["fully_discrete_entropy_checks"], 1)
                stats["fully_discrete_entropy_violation_rate"] = (
                    stats["fully_discrete_entropy_violations"] / checks
                )
                return None, stats
            state = next_state
            remaining -= dt
            stats["substeps"] += 1
            snapshot_substeps += 1
        snapshots.append(state.clone())
    checks = max(stats["fully_discrete_entropy_checks"], 1)
    stats["fully_discrete_entropy_violation_rate"] = (
        stats["fully_discrete_entropy_violations"] / checks
    )
    return torch.stack(snapshots, dim=1), stats


def evaluate_deployment(
    model: base.Solver,
    baseline: base.Solver,
    state_std: torch.Tensor,
) -> tuple[
    list[dict[str, Any]],
    dict[tuple[str, str], dict[str, torch.Tensor]],
    list[dict[str, Any]],
]:
    rows: list[dict[str, Any]] = []
    trajectories: dict[tuple[str, str], dict[str, torch.Tensor]] = {}
    no_flow_rows: list[dict[str, Any]] = []

    for name, display in PERIODIC_CASES.items():
        print(json.dumps({"stage": "periodic512", "case": name}), flush=True)
        reference, reference_stats = periodic_eval.strict_native_hllc_rollout(
            name, REFERENCE_CELLS
        )
        restricted = periodic_eval.comparison.precision.conservative_restrict(
            reference.numpy(), TRAIN_CELLS
        )
        scoring_reference = torch.from_numpy(restricted.astype(np.float32))
        roe, roe_stats = periodic_eval.matched_safe_rollout(
            baseline, name, TRAIN_CELLS
        )
        hcfl, hcfl_stats = periodic_eval.matched_safe_rollout(
            model, name, TRAIN_CELLS
        )
        no_flow, no_flow_stats = projection_only_rollout(
            model, name, True
        )
        rows.append(
            method_metrics(
                "periodic",
                name,
                "roe_512",
                scoring_reference,
                roe,
                state_std,
                True,
                roe_stats,
            )
        )
        rows.append(
            method_metrics(
                "periodic",
                name,
                "hcfl_trained_512",
                scoring_reference,
                hcfl,
                state_std,
                True,
                hcfl_stats,
            )
        )
        no_flow_row: dict[str, Any] = {
            "boundary": "periodic",
            "case": name,
            "method": "hcfl_projection_only_no_f_low",
            **no_flow_stats,
        }
        if no_flow is not None:
            no_flow_row.update(
                periodic_eval.comparison.diagnostics(
                    scoring_reference, no_flow, state_std
                )
            )
            no_flow_row.update(
                stencil.oscillation_metrics(
                    scoring_reference, no_flow, True
                )
            )
            no_flow_row["max_absolute_state_difference_vs_full_safe"] = float(
                torch.max(torch.abs(no_flow - hcfl))
            )
        no_flow_rows.append(no_flow_row)
        trajectories[("periodic", name)] = {
            "display": display,
            "reference": reference,
            "roe": roe,
            "hcfl": hcfl,
            "reference_stats": reference_stats,
        }

    for name, display in NONPERIODIC_CASES.items():
        print(json.dumps({"stage": "nonperiodic512", "case": name}), flush=True)
        reference, reference_stats = nonperiodic.strict_hllc_rollout(
            name, REFERENCE_CELLS
        )
        scoring_reference = nonperiodic.restrict_reference(
            reference, TRAIN_CELLS
        )
        roe, roe_stats = nonperiodic.learned_rollout(
            baseline, name, TRAIN_CELLS
        )
        hcfl, hcfl_stats = nonperiodic.learned_rollout(
            model, name, TRAIN_CELLS
        )
        no_flow, no_flow_stats = projection_only_rollout(
            model, name, False
        )
        rows.append(
            method_metrics(
                "transmissive",
                name,
                "roe_512",
                scoring_reference,
                roe,
                state_std,
                False,
                roe_stats,
            )
        )
        rows.append(
            method_metrics(
                "transmissive",
                name,
                "hcfl_trained_512",
                scoring_reference,
                hcfl,
                state_std,
                False,
                hcfl_stats,
            )
        )
        no_flow_row = {
            "boundary": "transmissive",
            "case": name,
            "method": "hcfl_projection_only_no_f_low",
            **no_flow_stats,
        }
        if no_flow is not None:
            no_flow_row.update(
                periodic_eval.comparison.diagnostics(
                    scoring_reference, no_flow, state_std
                )
            )
            no_flow_row.update(
                stencil.oscillation_metrics(
                    scoring_reference, no_flow, False
                )
            )
            no_flow_row["max_absolute_state_difference_vs_full_safe"] = float(
                torch.max(torch.abs(no_flow - hcfl))
            )
        no_flow_rows.append(no_flow_row)
        trajectories[("transmissive", name)] = {
            "display": display,
            "reference": reference,
            "roe": roe,
            "hcfl": hcfl,
            "reference_stats": reference_stats,
        }
    return rows, trajectories, no_flow_rows


def aggregate(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for boundary in ("periodic", "transmissive"):
        output[boundary] = {}
        for method in ("roe_512", "hcfl_trained_512"):
            selected = [
                row
                for row in rows
                if row["boundary"] == boundary and row["method"] == method
            ]
            output[boundary][method] = {
                "cases": len(selected),
                "mean_rollout_nrmse": float(
                    np.mean([row["rollout_nrmse"] for row in selected])
                ),
                "mean_final_snapshot_nrmse": float(
                    np.mean([row["final_snapshot_nrmse"] for row in selected])
                ),
                "mean_normalized_final_tv_excess": float(
                    np.mean(
                        [
                            row["mean_normalized_final_tv_excess"]
                            for row in selected
                        ]
                    )
                ),
                "total_excess_significant_extrema": int(
                    np.sum(
                        [
                            row["total_excess_significant_extrema"]
                            for row in selected
                        ]
                    )
                ),
                "minimum_density": float(
                    min(row["minimum_density"] for row in selected)
                ),
                "minimum_pressure": float(
                    min(row["minimum_pressure"] for row in selected)
                ),
            }
    return output


def add_zoom_inset(
    axis: plt.Axes,
    fine_x: np.ndarray,
    coarse_x: np.ndarray,
    reference: np.ndarray,
    roe: np.ndarray,
    hcfl: np.ndarray,
) -> None:
    # Do not magnify roundoff-level motion in a physically constant component
    # (notably u and p for a pure contact).  Insets are reserved for actual
    # solution variation.
    magnitude = max(
        float(np.max(np.abs(reference))),
        float(np.max(np.abs(roe))),
        float(np.max(np.abs(hcfl))),
        1.0,
    )
    if float(np.ptp(reference)) <= 1.0e-5 * magnitude:
        return
    reference_on_coarse = np.interp(coarse_x, fine_x, reference)
    gradient = np.abs(np.gradient(reference_on_coarse, coarse_x))
    advantage = (
        np.abs(roe - reference_on_coarse)
        - np.abs(hcfl - reference_on_coarse)
    )
    window_cells = 24
    kernel = np.ones(window_cells, dtype=np.float64)
    window_gradient = np.convolve(gradient, kernel, mode="valid")
    window_advantage = np.convolve(advantage, kernel, mode="valid")
    threshold = np.quantile(window_gradient, 0.85)
    candidates = np.flatnonzero(window_gradient >= threshold)
    improved = candidates[window_advantage[candidates] > 0.0]
    start = int(
        improved[np.argmax(window_advantage[improved])]
        if improved.size
        else np.argmax(window_gradient)
    )
    stop = start + window_cells
    dx = coarse_x[1] - coarse_x[0]
    x_left = max(0.0, float(coarse_x[start] - 0.5 * dx))
    x_right = min(1.0, float(coarse_x[stop - 1] + 0.5 * dx))
    inset_left = 0.54 if 0.5 * (x_left + x_right) < 0.5 else 0.04
    inset = axis.inset_axes([inset_left, 0.08, 0.42, 0.42])
    inset.plot(fine_x, reference, color="#111111", linewidth=0.85)
    inset.plot(
        coarse_x,
        roe,
        color="#D55E00",
        linestyle=(0, (8, 3)),
        linewidth=0.95,
    )
    inset.plot(
        coarse_x,
        hcfl,
        color="#0072B2",
        linestyle="--",
        linewidth=1.0,
    )
    inset.set_xlim(x_left, x_right)
    fine_mask = (fine_x >= x_left) & (fine_x <= x_right)
    coarse_mask = (coarse_x >= x_left) & (coarse_x <= x_right)
    zoom_values = np.concatenate(
        (reference[fine_mask], roe[coarse_mask], hcfl[coarse_mask])
    )
    y_low = float(np.min(zoom_values))
    y_high = float(np.max(zoom_values))
    global_span = max(float(np.ptp(reference)), 1.0e-8)
    padding = max(0.08 * (y_high - y_low), 0.005 * global_span)
    inset.set_ylim(y_low - padding, y_high + padding)
    inset.set_xticks([])
    inset.set_yticks([])
    axis.indicate_inset_zoom(
        inset,
        edgecolor="#555555",
        alpha=0.50,
        linewidth=0.55,
    )


def plot_combined_profiles(
    trajectories: dict[tuple[str, str], dict[str, torch.Tensor]],
    output: Path,
) -> None:
    figure, axes = plt.subplots(6, 4, figsize=(14.8, 15.0), squeeze=False)
    fine_x = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    coarse_x = (np.arange(TRAIN_CELLS) + 0.5) / TRAIN_CELLS
    blocks = (
        (0, "periodic", PERIODIC_CASES),
        (3, "transmissive", NONPERIODIC_CASES),
    )
    for row_offset, boundary, cases in blocks:
        for column, (case, display) in enumerate(cases.items()):
            item = trajectories[(boundary, case)]
            reference = base.primitive(item["reference"])[0, -1].numpy()
            roe = base.primitive(item["roe"])[0, -1].numpy()
            hcfl = base.primitive(item["hcfl"])[0, -1].numpy()
            for local_row, (component, variable) in enumerate(VARIABLES):
                row = row_offset + local_row
                axis = axes[row, column]
                first = row == 0 and column == 0
                axis.plot(
                    fine_x,
                    reference[:, component],
                    color="#111111",
                    linewidth=1.05,
                    label="HLLC-2048" if first else None,
                )
                axis.plot(
                    coarse_x,
                    roe[:, component],
                    color="#D55E00",
                    linestyle=(0, (8, 3)),
                    linewidth=1.0,
                    label="Roe-512" if first else None,
                )
                axis.plot(
                    coarse_x,
                    hcfl[:, component],
                    color="#0072B2",
                    linestyle="--",
                    linewidth=1.1,
                    label="HCFL trained at 512" if first else None,
                )
                component_reference = reference[:, component]
                component_roe = roe[:, component]
                component_hcfl = hcfl[:, component]
                component_magnitude = max(
                    float(np.max(np.abs(component_reference))),
                    float(np.max(np.abs(component_roe))),
                    float(np.max(np.abs(component_hcfl))),
                    1.0,
                )
                if (
                    float(np.ptp(component_reference))
                    <= 1.0e-5 * component_magnitude
                ):
                    center = float(np.mean(component_reference))
                    visible_deviation = float(
                        np.max(
                            np.abs(
                                np.concatenate(
                                    (component_roe, component_hcfl)
                                )
                                - center
                            )
                        )
                    )
                    half_span = max(
                        5.0 * visible_deviation,
                        5.0e-3 * component_magnitude,
                    )
                    axis.set_ylim(center - half_span, center + half_span)
                    axis.ticklabel_format(
                        axis="y", style="plain", useOffset=False
                    )
                add_zoom_inset(
                    axis,
                    fine_x,
                    coarse_x,
                    reference[:, component],
                    roe[:, component],
                    hcfl[:, component],
                )
                axis.grid(alpha=0.18)
                axis.set_xlim(0.0, 1.0)
                if local_row == 0:
                    axis.set_title(display)
                if column == 0:
                    axis.set_ylabel(variable)
                if local_row == 2:
                    axis.set_xlabel("x")
    final_time = base.DT_SNAPSHOT * (shared.CANONICAL_NSNAP - 1)
    figure.suptitle(
        f"Euler, 512-grid deployment — method comparison, t={final_time:.4f}",
        fontsize=13,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=3,
        frameon=False,
    )
    figure.text(
        0.5,
        0.010,
        "Insets: high-gradient windows selected by local HCFL-versus-Roe error reduction.",
        ha="center",
        fontsize=7.5,
        color="#555555",
    )
    figure.tight_layout(rect=(0.0, 0.075, 1.0, 0.975))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_training_curve(curve: list[dict[str, Any]], output: Path) -> None:
    updates = np.array([int(row["update"]) for row in curve])
    rollout = np.array(
        [float(row["validation_rollout_nrmse"]) for row in curve]
    )
    teacher = np.sqrt(
        np.array(
            [float(row["validation_teacher_forced_loss"]) for row in curve]
        )
    )
    figure, axis = plt.subplots(figsize=(7.0, 4.2))
    axis.plot(updates, rollout, label="validation rollout NRMSE", linewidth=1.4)
    axis.plot(
        updates,
        teacher,
        label="teacher-forced validation RMSE",
        linewidth=1.2,
    )
    best = int(np.argmin(rollout))
    axis.scatter([updates[best]], [rollout[best]], color="#D55E00", zorder=3)
    axis.set_xlabel("optimizer updates")
    axis.set_ylabel("normalized RMSE")
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def self_test(device: torch.device) -> None:
    if not math.isclose(
        INTERVAL_SUBSTEPS * TRAIN_LAMBDA / TRAIN_CELLS,
        base.DT_SNAPSHOT,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        raise RuntimeError("512-cell temporal subdivision is inconsistent")
    mean = np.array([1.0, 0.0, 1.0], dtype=np.float32)
    std = np.ones(3, dtype=np.float32)
    model = base.Solver(ARM.model_name, mean, std, width=16).to(device)
    x = (torch.arange(TRAIN_CELLS, device=device) + 0.5) / TRAIN_CELLS
    density = 1.0 + 0.1 * torch.sin(2.0 * math.pi * x)
    velocity = 0.2 * torch.cos(2.0 * math.pi * x)
    pressure = 1.0 + 0.1 * torch.sin(4.0 * math.pi * x)
    state = torch.stack(
        (
            density,
            density * velocity,
            pressure / (base.GAMMA - 1.0)
            + 0.5 * density * velocity.square(),
        ),
        dim=-1,
    )[None]
    prediction, feasibility = interval_forward(model, state, True)
    loss = prediction.square().mean() + ARM.feasibility_weight * feasibility
    loss.backward()
    if not bool(torch.isfinite(prediction).all()) or not bool(
        admissible(prediction).all()
    ):
        raise RuntimeError("Euler direct-512 differentiable interval failed")
    if not any(parameter.grad is not None for parameter in model.parameters()):
        raise RuntimeError("No gradient crossed the Euler direct-512 map")
    print("Euler train512 self-test passed", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase", choices=("data", "train", "evaluate", "all"), default="all"
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--reference-batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--fixed-budget-updates", type=int, default=2000)
    parser.add_argument("--max-updates", type=int, default=50000)
    parser.add_argument("--validation-interval", type=int, default=200)
    parser.add_argument("--plateau-patience", type=int, default=10)
    parser.add_argument("--minimum-relative-improvement", type=float, default=2.0e-3)
    parser.add_argument("--learning-rate-factor", type=float, default=0.3)
    parser.add_argument("--minimum-learning-rate", type=float, default=3.0e-6)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = device_from_argument(args.device)
    if args.self_test:
        self_test(device)
        return
    output = args.results_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    print(
        json.dumps(
            {
                "device": str(device),
                "gpu": (
                    torch.cuda.get_device_name(device)
                    if device.type == "cuda"
                    else None
                ),
                "phase": args.phase,
            }
        ),
        flush=True,
    )

    train_data, validation_data, precision_audit = load_or_generate_data(
        output,
        args.seed,
        device,
        args.reference_batch_size,
    )
    mean, std, state_std = prepare_statistics(train_data)
    checkpoint = output / f"{ARM.name}_converged_best_seed{args.seed}.pt"
    report_path = output / f"report_{ARM.name}_seed{args.seed}.json"
    curve_path = output / f"training_curve_{ARM.name}_seed{args.seed}.csv"

    if args.phase == "data":
        return

    model: base.Solver
    training_report: dict[str, Any]
    curve: list[dict[str, Any]]
    can_resume = args.resume and checkpoint.exists() and report_path.exists()
    if args.phase == "evaluate" or can_resume:
        model = load_model(output, args.seed, mean, std, args.width)
        training_report = json.loads(report_path.read_text(encoding="utf-8"))
        if curve_path.exists():
            with curve_path.open(newline="", encoding="utf-8") as handle:
                curve = list(csv.DictReader(handle))
        else:
            curve = []
        print(f"Loading converged checkpoint: {checkpoint}", flush=True)
    else:
        model, training_report, curve = train_to_convergence(
            train_data,
            validation_data,
            mean,
            std,
            state_std,
            args,
            output,
            device,
        )
        training_report["reference_precision_audit"] = precision_audit
        training_report["training_trajectory_count"] = int(train_data.shape[0])
        training_report["validation_trajectory_count"] = int(
            validation_data.shape[0]
        )
        training_report["multiplier_activity"] = multiplier_activity(
            model, validation_data
        )
        report_path.write_text(
            json.dumps(json_ready(training_report), indent=2),
            encoding="utf-8",
        )
        write_csv(curve_path, curve)
        plot_training_curve(
            curve,
            output / f"validation_convergence_train512_seed{args.seed}.png",
        )

    if args.phase == "train":
        return

    # The deployment cases are single trajectories with many small limiter
    # kernels.  A single CPU thread avoids OpenMP launch overhead and makes the
    # reduction order deterministic; training batches retain the default.
    if device.type == "cpu":
        torch.set_num_threads(1)

    torch.manual_seed(53000 + args.seed)
    baseline = base.Solver(ARM.model_name, mean, std, width=args.width)
    baseline.eval()
    rows, trajectories, no_flow_rows = evaluate_deployment(
        model, baseline, state_std
    )
    metrics_path = output / f"deployment_roe_vs_hcfl_train512_seed{args.seed}.csv"
    write_csv(metrics_path, rows)
    write_csv(
        output / f"no_f_low_ablation_seed{args.seed}.csv",
        no_flow_rows,
    )
    figure_path = output / f"euler_roe512_vs_hcfl512_with_hllc2048_seed{args.seed}.png"
    plot_combined_profiles(trajectories, figure_path)
    summary = {
        "scope": "single-seed direct 512-cell Euler training validation",
        "seed": args.seed,
        "architecture": (
            "central + nonnegative Roe + feasibility, symmetric 4-cell stencil"
        ),
        "training_report": training_report,
        "reference_precision_audit": precision_audit,
        "held_out_final_time": base.DT_SNAPSHOT
        * (shared.CANONICAL_NSNAP - 1),
        "figure_cases": {
            "first_block": list(PERIODIC_CASES),
            "second_block": list(NONPERIODIC_CASES),
        },
        "aggregate": aggregate(rows),
        "no_f_low_ablation": no_flow_rows,
    }
    (output / f"summary_roe_vs_hcfl_train512_seed{args.seed}.json").write_text(
        json.dumps(json_ready(summary), indent=2),
        encoding="utf-8",
    )
    print(json.dumps(json_ready(summary), indent=2), flush=True)


if __name__ == "__main__":
    main()
