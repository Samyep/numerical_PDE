"""Train one SWE HCFL model directly on 512 cells and audit its ringing.

The controlled comparison is the same central/nonnegative-Roe/feasibility
architecture trained on either 64-cell or 512-cell trajectories.  Both are
deployed through the identical hard-safety stack on the same held-out
512-cell periodic and transmissive problems.
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
BASE_EXPERIMENT = HERE.parent / "swe_1d_consistent_hcfl"
if str(BASE_EXPERIMENT) not in sys.path:
    sys.path.insert(0, str(BASE_EXPERIMENT))

import evaluate_deployment as deployment  # noqa: E402
import run_swe_consistent as base  # noqa: E402


REFERENCE_CELLS = 2048
TRAIN_CELLS = 512
RESTRICTION_FACTOR = REFERENCE_CELLS // TRAIN_CELLS
INTERVAL_SUBSTEPS = TRAIN_CELLS // base.NCOARSE
TRAIN_LAMBDA = base.DT_SNAPSHOT * base.NCOARSE
ARM = base.Arm(
    "central_nonnegative_feas_s4_train512",
    "central_nonnegative",
    1.0e-3,
    4,
)
ORIGINAL_ARM = next(
    arm
    for arm in base.ARMS
    if arm.name == "central_nonnegative_feas_s4"
)


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


# ---------------------------------------------------------------------------
# HLL-2048 reference generation, restricted conservatively to 512 cells
# ---------------------------------------------------------------------------


def reference_rhs(state: torch.Tensor, dx: float) -> torch.Tensor:
    flux = base.t_hll(state)
    return -(flux - torch.roll(flux, 1, dims=1)) / dx


def reference_ssprk2(
    state: torch.Tensor,
    dt: float,
    dx: float,
) -> torch.Tensor:
    stage_one = state + dt * reference_rhs(state, dx)
    if not bool(torch.isfinite(stage_one).all()) or not bool(
        (stage_one[..., 0] > 0.0).all()
    ):
        raise RuntimeError("HLL-2048 reference failed at SSP-RK2 stage one")
    stage_two = stage_one + dt * reference_rhs(stage_one, dx)
    output = 0.5 * state + 0.5 * stage_two
    if not bool(torch.isfinite(output).all()) or not bool(
        (output[..., 0] > 0.0).all()
    ):
        raise RuntimeError("HLL-2048 reference failed at SSP-RK2 completion")
    return output


def restrict_reference(state: torch.Tensor) -> torch.Tensor:
    batch = state.shape[0]
    return state.reshape(
        batch, TRAIN_CELLS, RESTRICTION_FACTOR, 2
    ).mean(dim=2)


@torch.no_grad()
def generate_reference_batch(
    initial: np.ndarray,
    device: torch.device,
    dtype: torch.dtype = torch.float32,
) -> torch.Tensor:
    state = torch.as_tensor(initial, dtype=dtype, device=device)
    dx = 1.0 / REFERENCE_CELLS
    saved = [restrict_reference(state).float().cpu()]
    for _ in range(1, base.NSNAP):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            h = state[..., 0]
            u = state[..., 1] / h
            maximum_speed = float(
                (u.abs() + torch.sqrt(base.G * h)).max().item()
            )
            dt = min(remaining, 0.20 * dx / max(maximum_speed, 1.0e-12))
            state = reference_ssprk2(state, dt, dx)
            remaining -= dt
        saved.append(restrict_reference(state).float().cpu())
    return torch.stack(saved, dim=1)


def initial_condition_groups(split: str, seed: int) -> list[tuple[str, np.ndarray]]:
    if split == "train":
        return [
            (
                "ordinary",
                base.legacy.generate_ic(
                    220, REFERENCE_CELLS, 21000 + seed, ood=False
                ),
            ),
            (
                "broad",
                base.legacy.generate_ic(
                    160, REFERENCE_CELLS, 22000 + seed, ood=True
                ),
            ),
            (
                "froude",
                base.generate_froude_coverage_ic(
                    200, REFERENCE_CELLS, 23000 + seed
                ),
            ),
        ]
    if split == "validation":
        return [
            (
                "ordinary",
                base.legacy.generate_ic(
                    44, REFERENCE_CELLS, 24000 + seed, ood=False
                ),
            ),
            (
                "broad",
                base.legacy.generate_ic(
                    32, REFERENCE_CELLS, 25000 + seed, ood=True
                ),
            ),
            (
                "froude",
                base.generate_froude_coverage_ic(
                    60, REFERENCE_CELLS, 26000 + seed
                ),
            ),
        ]
    raise ValueError(split)


def generate_split(
    split: str,
    seed: int,
    device: torch.device,
    reference_batch_size: int,
) -> torch.Tensor:
    pieces: list[torch.Tensor] = []
    for group_name, initial in initial_condition_groups(split, seed):
        for start in range(0, initial.shape[0], reference_batch_size):
            stop = min(start + reference_batch_size, initial.shape[0])
            started = time.perf_counter()
            result = generate_reference_batch(initial[start:stop], device)
            pieces.append(result)
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


def numpy_float64_reference(initial: np.ndarray) -> np.ndarray:
    state = np.asarray(initial, dtype=np.float64).copy()
    trajectories = state.shape[0]
    dx = 1.0 / REFERENCE_CELLS

    def restrict(values: np.ndarray) -> np.ndarray:
        return values.reshape(
            trajectories, TRAIN_CELLS, RESTRICTION_FACTOR, 2
        ).mean(axis=2)

    saved = [restrict(state)]
    for _ in range(1, base.NSNAP):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            h = state[..., 0]
            u = state[..., 1] / h
            speed = float(np.max(np.abs(u) + np.sqrt(base.G * h)))
            dt = min(remaining, 0.20 * dx / max(speed, 1.0e-12))
            state = base.np_periodic_ssprk2(state, dt, dx)
            remaining -= dt
        saved.append(restrict(state))
    return np.stack(saved, axis=1)


def reference_precision_audit(
    seed: int,
    device: torch.device,
) -> dict[str, float]:
    sample = base.legacy.generate_ic(
        2, REFERENCE_CELLS, 24000 + seed, ood=False
    )
    gpu = generate_reference_batch(sample, device).numpy().astype(np.float64)
    cpu = numpy_float64_reference(sample)
    error = gpu - cpu
    scale = np.std(cpu, axis=(0, 1, 2))
    normalized = error / np.maximum(scale, 1.0e-12)
    return {
        "trajectories": 2,
        "maximum_absolute_state_difference": float(np.max(np.abs(error))),
        "rms_state_difference": float(np.sqrt(np.mean(error * error))),
        "normalized_rms_state_difference": float(
            np.sqrt(np.mean(normalized * normalized))
        ),
    }


def load_or_generate_data(
    output: Path,
    seed: int,
    device: torch.device,
    reference_batch_size: int,
) -> tuple[torch.Tensor, torch.Tensor, dict[str, float]]:
    cache = output / f"data_cache_train512_seed{seed}.pt"
    expected_train = (580, base.NSNAP, TRAIN_CELLS, 2)
    expected_validation = (136, base.NSNAP, TRAIN_CELLS, 2)
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
    train = generate_split(
        "train", seed, device, reference_batch_size
    )
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


# ---------------------------------------------------------------------------
# Direct 512-cell training
# ---------------------------------------------------------------------------


def prepare_statistics(
    train_data: torch.Tensor,
) -> tuple[np.ndarray, np.ndarray, torch.Tensor]:
    values = base.primitive(train_data)
    mean = values.mean(dim=(0, 1, 2)).numpy()
    std = values.std(dim=(0, 1, 2)).numpy()
    state_std = train_data.std(dim=(0, 1, 2))
    return mean, std, state_std


def interval_forward(
    model: base.Solver,
    state: torch.Tensor,
    collect_feasibility: bool,
) -> tuple[torch.Tensor, torch.Tensor]:
    feasibility = torch.zeros((), dtype=state.dtype, device=state.device)
    for _ in range(INTERVAL_SUBSTEPS):
        raw_flux = model.raw_flux(state)
        if collect_feasibility:
            raw_residual = base.entropy_residual(raw_flux, state)
            feasibility = feasibility + torch.relu(raw_residual).square().mean()
        projected = base.strict_entropy_projection(raw_flux, state)
        state = base.fv_step(state, projected, TRAIN_LAMBDA)
    return state, feasibility / INTERVAL_SUBSTEPS


@torch.no_grad()
def validation_rollout(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    state = data[:, 0].to(device)
    scale = state_std.to(device)
    squared_errors: list[float] = []
    minimum_depth = float(state[..., 0].min().item())
    for snapshot in range(1, data.shape[1]):
        state, _ = interval_forward(model, state, False)
        target = data[:, snapshot].to(device)
        squared_errors.append(float((((state - target) / scale) ** 2).mean()))
        minimum_depth = min(minimum_depth, float(state[..., 0].min().item()))
    return {
        "rollout_nrmse": float(np.sqrt(np.mean(squared_errors))),
        "minimum_depth": minimum_depth,
    }


@torch.no_grad()
def teacher_forced_validation_loss(
    model: base.Solver,
    validation_data: torch.Tensor,
    state_std: torch.Tensor,
    device: torch.device,
    pairs: int = 256,
    batch_size: int = 32,
) -> float:
    inputs = validation_data[:, :-1].reshape(-1, TRAIN_CELLS, 2)
    targets = validation_data[:, 1:].reshape(-1, TRAIN_CELLS, 2)
    indices = torch.linspace(0, inputs.shape[0] - 1, pairs).long()
    scale = state_std.to(device)
    total = 0.0
    count = 0
    for start in range(0, pairs, batch_size):
        chosen = indices[start : start + batch_size]
        state = inputs[chosen].to(device)
        target = targets[chosen].to(device)
        prediction, _ = interval_forward(model, state, False)
        error = (prediction - target) / scale
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
    torch.manual_seed(41000 + args.seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(41000 + args.seed)
    model = base.Solver(
        ARM.model_name,
        mean,
        std,
        width=args.width,
        stencil_cells=ARM.stencil_cells,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    generator = torch.Generator().manual_seed(42000 + args.seed)
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
        improved = metric < best_metric
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
            "validation_minimum_depth": rollout["minimum_depth"],
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

        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach().item())

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
        raise RuntimeError("No 512-cell validation checkpoint was produced")
    model.load_state_dict(best_state)
    model.eval()
    final_validation = validation_rollout(
        model, validation_data, state_std, device
    )
    summary = {
        "seed": args.seed,
        "arm": ARM.name,
        "model": ARM.model_name,
        "stencil_cells": ARM.stencil_cells,
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


def load_train512_model(
    output: Path,
    seed: int,
    mean: np.ndarray,
    std: np.ndarray,
    width: int,
) -> base.Solver:
    checkpoint = output / f"{ARM.name}_converged_best_seed{seed}.pt"
    model = base.Solver(
        ARM.model_name,
        mean,
        std,
        width=width,
        stencil_cells=ARM.stencil_cells,
    )
    model.load_state_dict(
        torch.load(checkpoint, map_location="cpu", weights_only=True)
    )
    model.eval()
    return model


# ---------------------------------------------------------------------------
# Held-out 512-cell deployment comparison
# ---------------------------------------------------------------------------


def method_row(
    boundary: str,
    case: str,
    method: str,
    training_cells: int | None,
    reference: np.ndarray,
    prediction: np.ndarray,
    state_std: np.ndarray,
    periodic: bool,
    run_stats: dict[str, Any],
) -> dict[str, Any]:
    return {
        "boundary": boundary,
        "case": case,
        "method": method,
        "training_cells": training_cells,
        **deployment.diagnostics(reference, prediction, state_std, periodic),
        **{
            f"run_{key}": value
            for key, value in run_stats.items()
            if key != "seconds"
        },
    }


def evaluate_deployment(
    original_model: base.Solver,
    train512_model: base.Solver,
    classical_roe_model: base.Solver,
    state_std: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    trajectories: dict[str, Any] = {}
    case_groups = (
        ("periodic", True, deployment.PERIODIC_CASES),
        ("transmissive", False, deployment.NONPERIODIC_CASES),
    )
    for boundary, periodic, cases in case_groups:
        for case, display in cases.items():
            print(
                json.dumps(
                    {"stage": "held_out_deployment", "case": case}
                ),
                flush=True,
            )
            raw_reference, reference_stats = deployment.strict_hll_rollout(
                deployment.initial_condition(case, REFERENCE_CELLS), periodic
            )
            reference = deployment.conservative_restrict(
                raw_reference, TRAIN_CELLS
            )
            native, native_stats = deployment.strict_hll_rollout(
                reference[0], periodic
            )
            classical_roe, classical_roe_stats = deployment.learned_rollout(
                classical_roe_model, reference[0], periodic
            )
            trained64, stats64 = deployment.learned_rollout(
                original_model, reference[0], periodic
            )
            trained512, stats512 = deployment.learned_rollout(
                train512_model, reference[0], periodic
            )
            classical_proposal = deployment.proposal_diagnostics(
                classical_roe_model, reference, periodic
            )
            proposal64 = deployment.proposal_diagnostics(
                original_model, reference, periodic
            )
            proposal512 = deployment.proposal_diagnostics(
                train512_model, reference, periodic
            )
            rows.extend(
                [
                    method_row(
                        boundary,
                        case,
                        "native_hll_512",
                        None,
                        reference,
                        native,
                        state_std,
                        periodic,
                        native_stats,
                    ),
                    {
                        **method_row(
                            boundary,
                            case,
                            "classical_roe_512",
                            None,
                            reference,
                            classical_roe,
                            state_std,
                            periodic,
                            classical_roe_stats,
                        ),
                        **classical_proposal,
                    },
                    {
                        **method_row(
                            boundary,
                            case,
                            "hcfl_trained_64",
                            64,
                            reference,
                            trained64,
                            state_std,
                            periodic,
                            stats64,
                        ),
                        **proposal64,
                    },
                    {
                        **method_row(
                            boundary,
                            case,
                            "hcfl_trained_512",
                            512,
                            reference,
                            trained512,
                            state_std,
                            periodic,
                            stats512,
                        ),
                        **proposal512,
                    },
                ]
            )
            trajectories[case] = {
                "boundary": boundary,
                "display": display,
                "raw_reference": raw_reference,
                "native": native,
                "classical_roe": classical_roe,
                "trained64": trained64,
                "trained512": trained512,
                "reference_stats": reference_stats,
            }
    return rows, trajectories


def aggregate(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for boundary in ("periodic", "transmissive"):
        selected_boundary = [row for row in rows if row["boundary"] == boundary]
        output[boundary] = {}
        for method in (
            "native_hll_512",
            "classical_roe_512",
            "hcfl_trained_64",
            "hcfl_trained_512",
        ):
            selected = [
                row for row in selected_boundary if row["method"] == method
            ]
            output[boundary][method] = {
                "mean_rollout_nrmse": float(
                    np.mean([row["rollout_nrmse"] for row in selected])
                ),
                "mean_final_snapshot_nrmse": float(
                    np.mean([row["final_snapshot_nrmse"] for row in selected])
                ),
                "mean_positive_final_tv_excess": float(
                    np.mean(
                        [row["mean_positive_final_tv_excess"] for row in selected]
                    )
                ),
                "total_excess_significant_extrema": int(
                    np.sum(
                        [row["total_excess_significant_extrema"] for row in selected]
                    )
                ),
                "minimum_depth": float(
                    min(row["minimum_depth"] for row in selected)
                ),
            }
            if method != "native_hll_512":
                output[boundary][method][
                    "mean_absolute_roe_multiplier_change_from_one"
                ] = float(
                    np.mean(
                        [
                            row["mean_absolute_roe_multiplier_change_from_one"]
                            for row in selected
                        ]
                    )
                )
    return output


def plot_profiles(
    cases: dict[str, str],
    trajectories: dict[str, Any],
    title: str,
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, len(cases), figsize=(14.8, 6.3), squeeze=False)
    fine_x = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    coarse_x = (np.arange(TRAIN_CELLS) + 0.5) / TRAIN_CELLS
    for column, (case, display) in enumerate(cases.items()):
        item = trajectories[case]
        values = {
            "reference": deployment.primitive_np(item["raw_reference"][-1]),
            "native": deployment.primitive_np(item["native"][-1]),
            "classical_roe": deployment.primitive_np(item["classical_roe"][-1]),
            "trained64": deployment.primitive_np(item["trained64"][-1]),
            "trained512": deployment.primitive_np(item["trained512"][-1]),
        }
        for row, (component, label) in enumerate(((0, "h"), (1, "u"))):
            axis = axes[row, column]
            first = column == 0 and row == 0
            axis.plot(
                fine_x,
                values["reference"][:, component],
                color="#B7B7B7",
                linewidth=1.1,
                label="HLL-2048" if first else None,
            )
            axis.plot(
                coarse_x,
                values["native"][:, component],
                color="#222222",
                linestyle="--",
                linewidth=1.05,
                label="HLL-512" if first else None,
            )
            axis.plot(
                coarse_x,
                values["classical_roe"][:, component],
                color="#7A5195",
                linestyle=":",
                linewidth=1.05,
                label="classical Roe-512" if first else None,
            )
            axis.plot(
                coarse_x,
                values["trained64"][:, component],
                color="#D55E00",
                linewidth=1.0,
                alpha=0.85,
                label="HCFL trained at 64" if first else None,
            )
            axis.plot(
                coarse_x,
                values["trained512"][:, component],
                color="#0072B2",
                linewidth=1.15,
                label="HCFL trained at 512" if first else None,
            )
            axis.grid(alpha=0.18)
            axis.set_xlim(0.0, 1.0)
            if row == 0:
                axis.set_title(display)
            if column == 0:
                axis.set_ylabel(label)
            if row == 1:
                axis.set_xlabel("x")
    figure.suptitle(
        f"{title}, t={base.DT_SNAPSHOT * (deployment.EVAL_SNAPSHOTS - 1):.4f}",
        fontsize=13,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    figure.tight_layout(rect=(0.0, 0.07, 1.0, 0.94))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_combined_profiles(
    trajectories: dict[str, Any],
    output: Path,
) -> None:
    """Put periodic and transmissive profile audits in one figure."""
    figure, axes = plt.subplots(4, 4, figsize=(14.8, 11.0), squeeze=False)
    fine_x = (np.arange(REFERENCE_CELLS) + 0.5) / REFERENCE_CELLS
    coarse_x = (np.arange(TRAIN_CELLS) + 0.5) / TRAIN_CELLS
    blocks = (
        (0, deployment.PERIODIC_CASES),
        (2, deployment.NONPERIODIC_CASES),
    )
    for row_offset, cases in blocks:
        for column, (case, display) in enumerate(cases.items()):
            item = trajectories[case]
            values = {
                "reference": deployment.primitive_np(item["raw_reference"][-1]),
                "classical_roe": deployment.primitive_np(
                    item["classical_roe"][-1]
                ),
                "trained512": deployment.primitive_np(item["trained512"][-1]),
            }
            for local_row, (component, variable) in enumerate(((0, "h"), (1, "u"))):
                row = row_offset + local_row
                axis = axes[row, column]
                first = row == 0 and column == 0
                axis.plot(
                    fine_x,
                    values["reference"][:, component],
                    color="#B7B7B7",
                    linewidth=1.1,
                    label="HLL-2048" if first else None,
                )
                axis.plot(
                    coarse_x,
                    values["classical_roe"][:, component],
                    color="#7A5195",
                    linestyle=":",
                    linewidth=1.05,
                    label="classical Roe-512" if first else None,
                )
                axis.plot(
                    coarse_x,
                    values["trained512"][:, component],
                    color="#0072B2",
                    linewidth=1.15,
                    label="HCFL trained at 512" if first else None,
                )
                axis.grid(alpha=0.18)
                axis.set_xlim(0.0, 1.0)
                if local_row == 0:
                    axis.set_title(display.removeprefix("periodic "))
                if column == 0:
                    axis.set_ylabel(variable)
                if local_row == 1:
                    axis.set_xlabel("x")
    figure.suptitle(
        "512-grid deployment — method comparison, "
        f"t={base.DT_SNAPSHOT * (deployment.EVAL_SNAPSHOTS - 1):.4f}",
        fontsize=13,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    figure.tight_layout(rect=(0.0, 0.045, 1.0, 0.965))
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_summary(rows: list[dict[str, Any]], output: Path) -> None:
    summary = aggregate(rows)
    methods = (
        "native_hll_512",
        "classical_roe_512",
        "hcfl_trained_64",
        "hcfl_trained_512",
    )
    labels = ("HLL-512", "Roe-512", "HCFL train-64", "HCFL train-512")
    colors = ("#222222", "#7A5195", "#D55E00", "#0072B2")
    metrics = (
        ("mean_rollout_nrmse", "mean rollout NRMSE"),
        ("mean_positive_final_tv_excess", "mean positive final TV excess"),
        ("total_excess_significant_extrema", "total excess extrema"),
    )
    figure, axes = plt.subplots(2, 3, figsize=(13.2, 7.4))
    for row, boundary in enumerate(("periodic", "transmissive")):
        for column, (metric, ylabel) in enumerate(metrics):
            axis = axes[row, column]
            values = [summary[boundary][method][metric] for method in methods]
            axis.bar(np.arange(4), values, color=colors, alpha=0.9)
            axis.set_xticks(np.arange(4), labels, rotation=20, ha="right")
            axis.set_ylabel(ylabel)
            axis.grid(axis="y", alpha=0.2)
            if column == 0:
                axis.set_title(boundary)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_training_curve(curve: list[dict[str, Any]], output: Path) -> None:
    updates = np.array([row["update"] for row in curve])
    rollout = np.array([row["validation_rollout_nrmse"] for row in curve])
    teacher = np.sqrt(
        np.array([row["validation_teacher_forced_loss"] for row in curve])
    )
    figure, axis = plt.subplots(figsize=(7.0, 4.2))
    axis.plot(updates, rollout, label="validation rollout NRMSE", linewidth=1.4)
    axis.plot(updates, teacher, label="teacher-forced validation RMSE", linewidth=1.2)
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
    mean = np.array([1.0, 0.0], dtype=np.float32)
    std = np.ones(2, dtype=np.float32)
    model = base.Solver(
        ARM.model_name, mean, std, width=16, stencil_cells=4
    ).to(device)
    x = (torch.arange(TRAIN_CELLS, device=device) + 0.5) / TRAIN_CELLS
    h = 1.0 + 0.1 * torch.sin(2.0 * math.pi * x)
    u = 0.2 * torch.cos(2.0 * math.pi * x)
    state = torch.stack([h, h * u], dim=-1)[None]
    prediction, feasibility = interval_forward(model, state, True)
    loss = prediction.square().mean() + 1.0e-3 * feasibility
    loss.backward()
    if not bool(torch.isfinite(prediction).all()) or float(
        prediction[..., 0].detach().min()
    ) <= 0:
        raise RuntimeError("512-cell differentiable interval smoke test failed")
    if not any(parameter.grad is not None for parameter in model.parameters()):
        raise RuntimeError("No gradient crossed the 512-cell forward map")
    print("train512 self-test passed", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase", choices=("data", "train", "evaluate", "all"), default="all"
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--reference-batch-size", type=int, default=128)
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
                "gpu": torch.cuda.get_device_name(device) if device.type == "cuda" else None,
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
    report_path = output / f"report_{ARM.name}_seed{args.seed}.json"
    curve_path = output / f"training_curve_{ARM.name}_seed{args.seed}.csv"
    checkpoint = output / f"{ARM.name}_converged_best_seed{args.seed}.pt"

    if args.phase == "data":
        return

    model: base.Solver
    training_report: dict[str, Any]
    curve: list[dict[str, Any]]
    can_resume = args.resume and checkpoint.exists() and report_path.exists()
    if args.phase == "evaluate" or can_resume:
        model = load_train512_model(
            output, args.seed, mean, std, args.width
        )
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

    old_results = BASE_EXPERIMENT / "results"
    old_train, _ = base.load_or_make_data(old_results, args.seed)
    old_mean, old_std, old_state_std = base.prepare_statistics(old_train)
    original_model = base.load_model(
        ORIGINAL_ARM,
        base.checkpoint_path(old_results, ORIGINAL_ARM, args.seed),
        old_mean,
        old_std,
        args.width,
    )
    torch.manual_seed(43000 + args.seed)
    classical_roe_model = base.Solver(
        ARM.model_name,
        mean,
        std,
        width=args.width,
        stencil_cells=ARM.stencil_cells,
    )
    classical_roe_model.eval()
    rows, trajectories = evaluate_deployment(
        original_model,
        model,
        classical_roe_model,
        old_state_std.numpy(),
    )
    metrics_path = output / f"deployment_train64_vs_train512_seed{args.seed}.csv"
    write_csv(metrics_path, rows)
    plot_profiles(
        deployment.PERIODIC_CASES,
        trajectories,
        "Periodic: identical HCFL architecture, training resolution ablation",
        output / f"periodic_train64_vs_train512_seed{args.seed}.png",
    )
    plot_profiles(
        deployment.NONPERIODIC_CASES,
        trajectories,
        "Transmissive: identical HCFL architecture, training resolution ablation",
        output / f"nonperiodic_train64_vs_train512_seed{args.seed}.png",
    )
    plot_combined_profiles(
        trajectories,
        output / f"periodic_and_nonperiodic_train64_vs_train512_seed{args.seed}.png",
    )
    plot_summary(
        rows,
        output / f"stability_accuracy_train64_vs_train512_seed{args.seed}.png",
    )
    deployment_summary = {
        "scope": "single-seed direct 512-cell training validation",
        "seed": args.seed,
        "architecture": "central + nonnegative Roe + feasibility, 4-cell stencil",
        "training_report": training_report,
        "reference_precision_audit": precision_audit,
        "held_out_final_time": base.DT_SNAPSHOT * (deployment.EVAL_SNAPSHOTS - 1),
        "aggregate": aggregate(rows),
    }
    (output / f"summary_train64_vs_train512_seed{args.seed}.json").write_text(
        json.dumps(json_ready(deployment_summary), indent=2),
        encoding="utf-8",
    )
    print(json.dumps(json_ready(deployment_summary), indent=2), flush=True)


if __name__ == "__main__":
    main()
