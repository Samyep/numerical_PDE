"""Train and audit the exploratory 2-D SWE radial-dam comparison.

Compared methods:

* HLL-32 finite volume;
* 2-D HCFL with shared orientation-aware four/six-cell flux networks;
* vanilla FNO3d from the official clawNO repository;
* clawFNO3d from the official clawNO repository.

The neural-operator arms reuse the authors' source architecture and published
radial-dam hyperparameters, but train on the local recorded split.  They are
therefore official-architecture adaptations, not published-number
reproductions.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[4]
OFFICIAL_CLAWNO = WORKSPACE / ".hcfl_third_party_audit" / "clawNO"
OFFICIAL_COMMIT = "1c549dbf1d06dc35a8a5df2b897c62aa2f9db186"

sys.path.insert(0, str(HERE))
import swe_radial_common as C  # noqa: E402


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def save_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def official_commit() -> str:
    if not (OFFICIAL_CLAWNO / ".git").exists():
        raise FileNotFoundError(
            f"Official clawNO source is missing at {OFFICIAL_CLAWNO}"
        )
    return subprocess.check_output(
        ["git", "-C", str(OFFICIAL_CLAWNO), "rev-parse", "HEAD"],
        text=True,
    ).strip()


def training_statistics(train: torch.Tensor) -> tuple[np.ndarray, np.ndarray, torch.Tensor, torch.Tensor]:
    values = C.primitive(train)
    mean = values.mean(dim=(0, 1, 2, 3)).numpy()
    std = values.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-4).numpy()
    conserved_std = train.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-4)
    primitive_std = C.primitive_channel_scale(train)
    return mean, std, conserved_std, primitive_std


@torch.no_grad()
def rollout_hcfl(
    model: C.HCFL2D,
    initial: torch.Tensor,
    saved_states: int,
    batch_size: int = 5,
) -> tuple[torch.Tensor, dict[str, float]]:
    outputs: list[torch.Tensor] = []
    aggregate: defaultdict[str, float] = defaultdict(float)
    aggregate["minimum_depth_theta"] = 1.0
    aggregate["minimum_entropy_beta"] = 1.0
    aggregate["maximum_entropy_balance"] = -float("inf")
    aggregate["maximum_interface_residual"] = -float("inf")
    for start in range(0, initial.shape[0], batch_size):
        state = initial[start : start + batch_size].clone()
        sequence = [state.cpu()]
        for _ in range(1, saved_states):
            state, diagnostics = C.advance_hcfl_interval(model, state)
            C.merge_stats(aggregate, diagnostics)
            sequence.append(state.cpu())
        outputs.append(torch.stack(sequence, dim=1))
    return torch.cat(outputs, dim=0), dict(aggregate)


@torch.no_grad()
def rollout_hll(
    initial: torch.Tensor,
    saved_states: int,
    batch_size: int = 10,
) -> tuple[torch.Tensor, dict[str, float]]:
    outputs: list[torch.Tensor] = []
    aggregate: defaultdict[str, float] = defaultdict(float)
    aggregate["maximum_entropy_balance"] = -float("inf")
    for start in range(0, initial.shape[0], batch_size):
        state = initial[start : start + batch_size].clone()
        sequence = [state.cpu()]
        for _ in range(1, saved_states):
            state, diagnostics = C.advance_hll_interval(state)
            C.merge_stats(aggregate, diagnostics)
            sequence.append(state.cpu())
        outputs.append(torch.stack(sequence, dim=1))
    return torch.cat(outputs, dim=0), dict(aggregate)


@torch.no_grad()
def validate_hcfl(
    model: C.HCFL2D,
    validation: torch.Tensor,
    primitive_std: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    prediction, diagnostics = rollout_hcfl(
        model, validation[:, 0].to(device), validation.shape[1]
    )
    finite = torch.isfinite(prediction).all(dim=(-4, -3, -2, -1))
    positive = prediction[..., 0].amin(dim=(-3, -2, -1)) >= C.H_FLOOR
    complete = finite & positive
    if bool(complete.any()):
        error = C.primitive_nrmse(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
    else:
        error = float("inf")
    return {
        "completion": float(complete.float().mean()),
        "rollout_nrmse_completed": error,
        "minimum_depth": float(prediction[..., 0].min()),
        "maximum_entropy_balance": diagnostics["maximum_entropy_balance"],
        "maximum_interface_residual": diagnostics["maximum_interface_residual"],
    }


def train_hcfl_arm(
    stencil: int,
    train: torch.Tensor,
    validation: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    conserved_std: torch.Tensor,
    primitive_std: torch.Tensor,
    output: Path,
    device: torch.device,
    max_updates: int,
    validation_interval: int,
    batch_size: int,
    seed: int,
) -> tuple[Path, dict[str, Any]]:
    seed_everything(71000 + seed + stencil)
    model = C.HCFL2D(mean, std, width=72, stencil_cells=stencil).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=5.0e-4)
    generator = torch.Generator(device=device).manual_seed(72000 + seed + stencil)
    train_device = train.to(device)
    conserved_std = conserved_std.to(device)
    checkpoint = output / f"hcfl_s{stencil}_best_seed{seed}.pt"
    curve: list[dict[str, Any]] = []
    best_completion = -1.0
    best_error = float("inf")
    best_update = -1
    plateau_anchor = float("inf")
    stale = 0
    stopped = max_updates
    reason = "maximum_updates"
    started = time.perf_counter()

    def check(update: int, minibatch_loss: float | None) -> None:
        nonlocal best_completion, best_error, best_update
        model.eval()
        metrics = validate_hcfl(model, validation, primitive_std, device)
        completion = metrics["completion"]
        error = metrics["rollout_nrmse_completed"]
        improved = completion > best_completion or (
            completion == best_completion and error < best_error
        )
        if improved:
            best_completion = completion
            best_error = error
            best_update = update
            torch.save(C.clone_state_dict(model), checkpoint)
        row = {
            "method": f"HCFL-s{stencil}",
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": minibatch_loss,
            "new_best": improved,
            **metrics,
        }
        curve.append(row)
        print(row, flush=True)

    check(0, None)
    last_loss: float | None = None
    for update in range(1, max_updates + 1):
        model.train()
        trajectory_indices = torch.randint(
            0, train_device.shape[0], (batch_size,), generator=generator, device=device
        )
        time_indices = torch.randint(
            0, train_device.shape[1] - 1, (batch_size,), generator=generator, device=device
        )
        state = train_device[trajectory_indices, time_indices]
        target = train_device[trajectory_indices, time_indices + 1]
        feasibility = torch.zeros((), device=device)
        for _ in range(4):
            raw = model.raw_fluxes(state)
            feasibility = feasibility + model.feasibility_loss_from_raw(state, raw)
            state = C.flux_divergence(
                state, *model.project_raw_fluxes_training(state, raw), C.SNAPSHOT_DT / 4.0
            )
        trajectory_loss = (((state - target) / conserved_std) ** 2).mean()
        loss = trajectory_loss + 1.0e-3 * feasibility / 4.0
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach())

        if update % validation_interval and update != max_updates:
            continue
        check(update, last_loss)
        current = curve[-1]["rollout_nrmse_completed"]
        if update < 1500:
            continue
        if plateau_anchor == float("inf"):
            plateau_anchor = best_error
            continue
        if current <= plateau_anchor * 0.998:
            plateau_anchor = current
            stale = 0
        else:
            stale += 1
        if stale < 2:
            continue
        old_lr = float(optimizer.param_groups[0]["lr"])
        if old_lr <= 2.0e-5 * (1.0 + 1.0e-12):
            stopped = update
            reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_lr = max(2.0e-5, old_lr * 0.2)
        for group in optimizer.param_groups:
            group["lr"] = new_lr
        stale = 0
        plateau_anchor = best_error
        print({"method": f"HCFL-s{stencil}", "event": "reduce_lr", "new_lr": new_lr}, flush=True)

    pd.DataFrame(curve).to_csv(output / f"hcfl_s{stencil}_curve_seed{seed}.csv", index=False)
    report = {
        "method": "2-D HLL + Roe correction + hard projection + feasibility",
        "stencil_cells": stencil,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "best_update": best_update,
        "stop_update": stopped,
        "stop_reason": reason,
        "best_validation_completion": best_completion,
        "best_validation_rollout_nrmse": best_error,
        "training_seconds": time.perf_counter() - started,
        "maximum_updates": max_updates,
        "proposal_feasibility_weight": 1.0e-3,
        "uses_low_order_flux_during_training": False,
        "uses_low_order_flux_only_in_deployment_safety_wrapper": True,
    }
    save_json(output / f"hcfl_s{stencil}_report_seed{seed}.json", report)
    return checkpoint, report


def _load_official_models() -> tuple[type[torch.nn.Module], type[torch.nn.Module]]:
    if official_commit() != OFFICIAL_COMMIT:
        raise RuntimeError("The checked-out clawNO commit differs from the design record")
    sys.path.insert(0, str(OFFICIAL_CLAWNO))
    from models.FNO import FNO3d  # type: ignore
    from models.clawFNO import clawFNO3d  # type: ignore

    return FNO3d, clawFNO3d


def build_operator(kind: str, device: torch.device) -> torch.nn.Module:
    FNO3d, clawFNO3d = _load_official_models()
    cls = FNO3d if kind == "fno" else clawFNO3d
    return cls(
        num_channels=3,
        initial_step=1,
        modes1=8,
        modes2=8,
        modes3=6,
        width=20,
        grid_type="symmetric",
        time=True,
        time_pad=False,
    ).to(device)


def operator_xy_time(trajectory: torch.Tensor) -> torch.Tensor:
    # [batch,time,y,x,(h,hu,hv)] -> [batch,x,y,time,(hu,hv,h)]
    return C.conserved_to_claw(trajectory).permute(0, 3, 2, 1, 4).contiguous()


def make_operator_input(initial: torch.Tensor, frames: int = 24) -> torch.Tensor:
    # initial [batch,x,y,channels]
    return initial[:, :, :, None, None, :].repeat(1, 1, 1, frames, 1, 1)


def primitive_claw(value: torch.Tensor) -> torch.Tensor:
    h = value[..., 2]
    safe = torch.where(h.abs() >= 1.0e-6, h, torch.where(h >= 0.0, 1.0e-6, -1.0e-6))
    return torch.stack((value[..., 0] / safe, value[..., 1] / safe, h), dim=-1)


def relative_channel_l2(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return relative_channel_l2_per_sample(prediction, target).mean()


def relative_channel_l2_per_sample(
    prediction: torch.Tensor, target: torch.Tensor
) -> torch.Tensor:
    difference = (prediction - target).reshape(prediction.shape[0], -1, 3)
    target_flat = target.reshape(target.shape[0], -1, 3)
    numerator = torch.linalg.vector_norm(difference, dim=1)
    denominator = torch.linalg.vector_norm(target_flat, dim=1).clamp_min(1.0e-8)
    return (numerator / denominator).mean(dim=-1)


@torch.no_grad()
def operator_validation(
    model: torch.nn.Module,
    validation: torch.Tensor,
    device: torch.device,
    batch_size: int = 5,
) -> dict[str, float | int | None]:
    model.eval()
    converted = operator_xy_time(validation)
    all_losses: list[torch.Tensor] = []
    completed_losses: list[torch.Tensor] = []
    completed_count = 0
    minimum_depth = float("inf")
    for start in range(0, converted.shape[0], batch_size):
        data = converted[start : start + batch_size].to(device)
        prediction = model(make_operator_input(data[..., 0, :])).reshape(
            data.shape[0], C.N_COARSE, C.N_COARSE, 24, 3
        )
        target = data[..., 1:25, :]
        sample_losses = relative_channel_l2_per_sample(
            primitive_claw(prediction), primitive_claw(target)
        ).detach().cpu()
        complete = torch.isfinite(prediction).reshape(prediction.shape[0], -1).all(dim=1)
        complete = complete & (prediction[..., 2] >= C.H_FLOOR).reshape(
            prediction.shape[0], -1
        ).all(dim=1)
        all_losses.append(sample_losses)
        if bool(complete.any()):
            completed_losses.append(sample_losses[complete.cpu()])
        completed_count += int(complete.sum())
        minimum_depth = min(minimum_depth, float(prediction[..., 2].min()))
    diagnostic = float(torch.cat(all_losses).mean())
    completed_error = (
        float(torch.cat(completed_losses).mean()) if completed_losses else None
    )
    return {
        "completion": completed_count / converted.shape[0],
        "completed_trajectories": completed_count,
        "primitive_relative_l2_completed_only": completed_error,
        "diagnostic_primitive_relative_l2_all_predictions": diagnostic,
        "minimum_depth": minimum_depth,
    }


def train_operator(
    kind: str,
    train: torch.Tensor,
    validation: torch.Tensor,
    output: Path,
    device: torch.device,
    epochs: int,
    seed: int,
) -> tuple[Path, dict[str, Any]]:
    seed_everything(73000 + seed + (0 if kind == "fno" else 1000))
    model = build_operator(kind, device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1.0e-2, weight_decay=1.0e-4)
    batch_size = 10
    steps_per_epoch = math.ceil(train.shape[0] / batch_size)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=epochs * steps_per_epoch
    )
    generator = torch.Generator().manual_seed(74000 + seed)
    converted = operator_xy_time(train)
    checkpoint = output / f"{kind}_official_architecture_best_seed{seed}.pt"
    best_completion = -1.0
    best_ranking_error = float("inf")
    best_metrics: dict[str, float | int | None] | None = None
    best_epoch = -1
    stale = 0
    curve: list[dict[str, Any]] = []
    started = time.perf_counter()
    completed_epochs = 0

    for epoch in range(1, epochs + 1):
        model.train()
        permutation = torch.randperm(converted.shape[0], generator=generator)
        training_loss = 0.0
        batches = 0
        for start in range(0, converted.shape[0], batch_size):
            indices = permutation[start : start + batch_size]
            data = converted[indices].to(device)
            initial = data[..., 0, :]
            target = data[..., 1:25, :]
            prediction = model(make_operator_input(initial)).reshape(
                data.shape[0], C.N_COARSE, C.N_COARSE, 24, 3
            )
            loss = relative_channel_l2(prediction, target)
            if not bool(torch.isfinite(loss)):
                raise RuntimeError(f"{kind} training became non-finite at epoch {epoch}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            scheduler.step()
            training_loss += float(loss.detach())
            batches += 1

        validation_metrics = operator_validation(model, validation, device)
        completion = float(validation_metrics["completion"])
        completed_error = validation_metrics["primitive_relative_l2_completed_only"]
        # If no validation trajectory is physical, retain a diagnostic
        # checkpoint only.  It is explicitly marked inadmissible in the
        # report and is never presented as a physically valid winner.
        ranking_error = float(
            completed_error
            if completed_error is not None
            else validation_metrics["diagnostic_primitive_relative_l2_all_predictions"]
        )
        improved = completion > best_completion + 1.0e-12 or (
            abs(completion - best_completion) <= 1.0e-12
            and ranking_error < best_ranking_error
        )
        if improved:
            best_completion = completion
            best_ranking_error = ranking_error
            best_metrics = dict(validation_metrics)
            best_epoch = epoch
            stale = 0
            torch.save(C.clone_state_dict(model), checkpoint)
        else:
            stale += 1
        row = {
            "method": kind,
            "epoch": epoch,
            "updates": epoch * steps_per_epoch,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "training_relative_l2": training_loss / max(1, batches),
            **{f"validation_{key}": value for key, value in validation_metrics.items()},
            "new_best": improved,
        }
        curve.append(row)
        completed_epochs = epoch
        if epoch == 1 or epoch % 10 == 0 or improved:
            print(row, flush=True)
        if stale >= 100:
            break

    pd.DataFrame(curve).to_csv(output / f"{kind}_training_curve_seed{seed}.csv", index=False)
    report = {
        "method": "official FNO3d" if kind == "fno" else "official clawFNO3d",
        "label": "official-architecture adaptation",
        "official_source_commit": OFFICIAL_COMMIT,
        "parameters": sum(parameter.numel() for parameter in model.parameters()),
        "best_epoch": best_epoch,
        "stop_epoch": completed_epochs,
        "best_validation_completion": best_completion,
        "has_fully_admissible_validation_checkpoint": best_completion >= 1.0,
        "best_validation_primitive_relative_l2_completed_only": (
            None if best_metrics is None else best_metrics["primitive_relative_l2_completed_only"]
        ),
        "diagnostic_validation_relative_l2_all_predictions": (
            None
            if best_metrics is None
            else best_metrics["diagnostic_primitive_relative_l2_all_predictions"]
        ),
        "best_validation_minimum_depth": (
            None if best_metrics is None else best_metrics["minimum_depth"]
        ),
        "checkpoint_selection_rule": "lexicographic(completion, -completed_only_error)",
        "training_seconds": time.perf_counter() - started,
        "published_radial_hyperparameters": {
            "modes": 8,
            "time_modes": 6,
            "width": 20,
            "batch_size": 10,
            "learning_rate": 1.0e-2,
            "weight_decay": 1.0e-4,
            "maximum_epochs": epochs,
            "early_stopping_patience": 100,
        },
    }
    save_json(output / f"{kind}_training_report_seed{seed}.json", report)
    return checkpoint, report


@torch.no_grad()
def rollout_operator(
    model: torch.nn.Module,
    initial: torch.Tensor,
    saved_states: int,
    device: torch.device,
    batch_size: int = 5,
) -> torch.Tensor:
    outputs: list[torch.Tensor] = []
    blocks = (saved_states - 1) // 24
    if 1 + 24 * blocks != saved_states:
        raise ValueError("Operator rollout length must be an integer number of 24-frame blocks")
    for start in range(0, initial.shape[0], batch_size):
        current = operator_xy_time(initial[start : start + batch_size, None])[:, :, :, 0].to(device)
        pieces = [current]
        for _ in range(blocks):
            prediction = model(make_operator_input(current)).reshape(
                current.shape[0], C.N_COARSE, C.N_COARSE, 24, 3
            )
            pieces.append(prediction)
            current = prediction[..., -1, :]
        claw_sequence = torch.cat((pieces[0][..., None, :], *pieces[1:]), dim=-2)
        conserved = C.claw_to_conserved(claw_sequence).permute(0, 3, 2, 1, 4).cpu()
        outputs.append(conserved)
    return torch.cat(outputs, dim=0)


def load_hcfl(
    checkpoint: Path,
    stencil: int,
    mean: np.ndarray,
    std: np.ndarray,
    device: torch.device,
) -> C.HCFL2D:
    model = C.HCFL2D(mean, std, width=72, stencil_cells=stencil).to(device)
    model.load_state_dict(torch.load(checkpoint, map_location=device, weights_only=True))
    model.eval()
    return model


def load_operator(kind: str, checkpoint: Path, device: torch.device) -> torch.nn.Module:
    model = build_operator(kind, device)
    model.load_state_dict(torch.load(checkpoint, map_location=device, weights_only=True))
    model.eval()
    return model


def evaluate_all(
    dataset: dict[str, Any],
    output: Path,
    mean: np.ndarray,
    std: np.ndarray,
    primitive_std: torch.Tensor,
    hcfl_checkpoint: Path,
    hcfl_stencil: int,
    fno_checkpoint: Path,
    claw_checkpoint: Path,
    device: torch.device,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    hcfl = load_hcfl(hcfl_checkpoint, hcfl_stencil, mean, std, device)
    fno = load_operator("fno", fno_checkpoint, device)
    claw = load_operator("claw", claw_checkpoint, device)
    methods = ("HLL-32", f"HCFL-s{hcfl_stencil}", "FNO", "clawFNO")
    rows: list[dict[str, Any]] = []
    trajectories: dict[str, dict[str, torch.Tensor]] = {}
    exact_safety: dict[str, dict[str, Any]] = {}

    for split_name in ("test_id", "test_radius_ood", "test_height_ood", "test_long"):
        reference = dataset[split_name]["trajectory"]
        initial = reference[:, 0].to(device)
        print({"stage": "evaluate_swe", "split": split_name}, flush=True)
        hll_prediction, hll_stats = rollout_hll(initial, reference.shape[1])
        hcfl_prediction, hcfl_stats = rollout_hcfl(hcfl, initial, reference.shape[1])
        fno_prediction = rollout_operator(fno, reference[:, 0], reference.shape[1], device)
        claw_prediction = rollout_operator(claw, reference[:, 0], reference.shape[1], device)
        prediction_by_method = {
            methods[0]: hll_prediction,
            methods[1]: hcfl_prediction,
            methods[2]: fno_prediction,
            methods[3]: claw_prediction,
        }
        trajectories[split_name] = {"reference": reference, **prediction_by_method}
        exact_safety[split_name] = {"HLL-32": hll_stats, methods[1]: hcfl_stats}
        for method, prediction in prediction_by_method.items():
            for index in range(reference.shape[0]):
                audit = C.audit_sequence(prediction[index], reference[index], primitive_std)
                rows.append(
                    {
                        "split": split_name,
                        "case_index": index,
                        "radius": float(dataset[split_name]["radius"][index]),
                        "inner_height": float(dataset[split_name]["inner_height"][index]),
                        "method": method,
                        **asdict(audit),
                    }
                )

    frame = pd.DataFrame(rows)
    frame.to_csv(output / "swe_radial_case_metrics.csv", index=False)
    torch.save(trajectories, output / "swe_radial_trajectories.pt")
    summary: dict[str, Any] = {
        "completion_definition": (
            "finite and positive only; this is not a conservation/entropy/oscillation "
            "success criterion"
        ),
        "exact_flux_method_safety": exact_safety,
        "aggregate": {},
    }
    for (split_name, method), group in frame.groupby(["split", "method"], sort=False):
        completed = group["completed"].astype(bool)
        completed_rows = group[completed]
        summary["aggregate"].setdefault(split_name, {})[method] = {
            "completion": int(completed.sum()),
            "total": int(len(group)),
            "mean_rollout_nrmse_completed": None
            if completed_rows.empty
            else float(completed_rows["rollout_nrmse"].mean()),
            "mean_rollout_nmae_completed": None
            if completed_rows.empty
            else float(completed_rows["rollout_nmae"].mean()),
            "mean_final_nrmse_completed": None
            if completed_rows.empty
            else float(completed_rows["final_nrmse"].mean()),
            "mean_final_nmae_completed": None
            if completed_rows.empty
            else float(completed_rows["final_nmae"].mean()),
            "mean_rollout_height_mae_completed": None
            if completed_rows.empty
            else float(completed_rows["rollout_height_mae"].mean()),
            "mean_rollout_x_velocity_mae_completed": None
            if completed_rows.empty
            else float(completed_rows["rollout_x_velocity_mae"].mean()),
            "mean_rollout_y_velocity_mae_completed": None
            if completed_rows.empty
            else float(completed_rows["rollout_y_velocity_mae"].mean()),
            "minimum_depth": float(group["minimum_depth"].min()),
            "maximum_conservation_residual": float(group["maximum_conservation_residual"].max()),
            "maximum_entropy_balance_completed_only": (
                None
                if group["maximum_entropy_balance"].dropna().empty
                else float(group["maximum_entropy_balance"].dropna().max())
            ),
            "mean_entropy_violation_rate_completed_only": (
                None
                if group["entropy_violation_rate"].dropna().empty
                else float(group["entropy_violation_rate"].dropna().mean())
            ),
            "mean_height_tv_excess_completed": None
            if completed_rows.empty
            else float(completed_rows["height_tv_excess"].mean()),
            "mean_height_range_overshoot_completed": None
            if completed_rows.empty
            else float(completed_rows["height_range_overshoot"].mean()),
            "mean_centerline_curvature_ratio_completed": None
            if completed_rows.empty
            else float(completed_rows["centerline_curvature_ratio"].mean()),
        }
    summary["method_qualification"] = {
        "FNO": {
            "status": "not_physics_qualified",
            "lower_nrmse_is_not_solver_success": True,
            "lower_nmae_is_not_solver_success": True,
            "reason": (
                "finite/positive output but no FV conservation or entropy guarantee; "
                "measured conservation and entropy violations plus visible OOD ringing"
            ),
            "by_split": {
                split_name: {
                    "maximum_conservation_residual": values["FNO"][
                        "maximum_conservation_residual"
                    ],
                    "mean_entropy_violation_rate": values["FNO"][
                        "mean_entropy_violation_rate_completed_only"
                    ],
                    "mean_centerline_curvature_ratio": values["FNO"][
                        "mean_centerline_curvature_ratio_completed"
                    ],
                    "mean_height_range_overshoot": values["FNO"][
                        "mean_height_range_overshoot_completed"
                    ],
                }
                for split_name, values in summary["aggregate"].items()
            },
        },
        "clawFNO": {
            "status": "failed",
            "reason": "no finite-positive validation or evaluation trajectory",
        },
    }
    save_json(output / "swe_radial_summary.json", summary)
    make_figures(trajectories, frame, output, hcfl_stencil)
    return frame, summary


def make_figures(
    trajectories: dict[str, dict[str, torch.Tensor]],
    metrics: pd.DataFrame,
    output: Path,
    hcfl_stencil: int,
) -> None:
    method_order = ("reference", "HLL-32", f"HCFL-s{hcfl_stencil}", "FNO", "clawFNO")
    display_names = {
        "reference": "reference",
        "HLL-32": "HLL-32",
        f"HCFL-s{hcfl_stencil}": f"HCFL-s{hcfl_stencil}",
        "FNO": "FNO (non-admissible)",
        "clawFNO": "clawFNO (failed)",
    }
    split_order = ("test_id", "test_radius_ood", "test_height_ood", "test_long")
    labels = ("ID", "radius OOD", "height OOD", "2x horizon")
    fig, axes = plt.subplots(len(split_order), len(method_order), figsize=(16, 11), constrained_layout=True)
    for row, (split_name, label) in enumerate(zip(split_order, labels)):
        values = trajectories[split_name]
        reference_h = values["reference"][0, -1, ..., 0]
        vmin, vmax = float(reference_h.min()), float(reference_h.max())
        for column, method in enumerate(method_order):
            image = values[method][0, -1, ..., 0]
            axis = axes[row, column]
            handle = axis.imshow(
                image,
                origin="lower",
                extent=(C.DOMAIN_LOWER, C.DOMAIN_UPPER, C.DOMAIN_LOWER, C.DOMAIN_UPPER),
                cmap="viridis",
                vmin=vmin,
                vmax=vmax,
            )
            axis.set_xticks([])
            axis.set_yticks([])
            if row == 0:
                axis.set_title(display_names[method])
            if column == 0:
                axis.set_ylabel(label)
            if column == len(method_order) - 1:
                fig.colorbar(handle, ax=axis, fraction=0.046, pad=0.02)
    fig.suptitle("2-D SWE radial dam break: final water depth", fontsize=14)
    fig.savefig(output / "swe_radial_final_depth.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=False)
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.08,
        top=0.84,
        wspace=0.16,
        hspace=0.32,
    )
    colors = {
        "reference": "black",
        "HLL-32": "#8c8c8c",
        f"HCFL-s{hcfl_stencil}": "#0072B2",
        "FNO": "#D55E00",
        "clawFNO": "#009E73",
    }
    styles = {"reference": "-", "HLL-32": "--", f"HCFL-s{hcfl_stencil}": "-", "FNO": "-.", "clawFNO": (0, (5, 2))}
    coordinate = np.linspace(C.DOMAIN_LOWER + C.DOMAIN_LENGTH / 64, C.DOMAIN_UPPER - C.DOMAIN_LENGTH / 64, C.N_COARSE)
    middle = C.N_COARSE // 2
    for axis, split_name, label in zip(axes.flat, split_order, labels):
        for method in method_order:
            profile = trajectories[split_name][method][0, -1, middle, :, 0]
            axis.plot(
                coordinate,
                profile,
                color=colors[method],
                linestyle=styles[method],
                linewidth=1.8,
                label=display_names[method],
            )
        axis.set_title(label)
        axis.set_xlabel("x")
        axis.set_ylabel("h")
        axis.grid(alpha=0.2)
    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=5,
        frameon=False,
    )
    fig.suptitle("Central depth profiles", y=0.985)
    fig.savefig(output / "swe_radial_center_profiles.png", dpi=180)
    plt.close(fig)

    aggregate = (
        metrics.groupby(["split", "method"], sort=False)
        .agg(
            completion=("completed", "mean"),
            nrmse=("rollout_nrmse", "mean"),
            nmae=("rollout_nmae", "mean"),
        )
        .reset_index()
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.20, top=0.76, wspace=0.28)
    plot_methods = ("HLL-32", f"HCFL-s{hcfl_stencil}", "FNO", "clawFNO")
    width = 0.19
    x = np.arange(len(split_order))
    for offset, method in enumerate(plot_methods):
        subset = aggregate[aggregate.method == method].set_index("split").reindex(split_order)
        axes[0].bar(
            x + (offset - 1.5) * width,
            subset.completion,
            width,
            label=display_names[method],
        )
        axes[1].bar(x + (offset - 1.5) * width, subset.nrmse, width, label=method)
        axes[2].bar(x + (offset - 1.5) * width, subset.nmae, width, label=method)
    for axis in axes:
        axis.set_xticks(x, labels, rotation=15)
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("finite + positive rate\n(not entropy-qualified)")
    axes[1].set_ylabel("rollout NRMSE\n(finite + positive only)")
    axes[2].set_ylabel("rollout normalized MAE\n(finite + positive only)")
    axes[0].set_title("Completion screen")
    axes[1].set_title("Normalized RMSE")
    axes[2].set_title("Normalized MAE")
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=4,
        frameon=False,
    )
    fig.suptitle("Low error alone is not physical solver success", y=0.98)
    fig.savefig(output / "swe_radial_aggregate.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def choose_hcfl_report(reports: list[dict[str, Any]]) -> dict[str, Any]:
    return sorted(
        reports,
        key=lambda report: (
            -report["best_validation_completion"],
            report["best_validation_rollout_nrmse"],
            report["stencil_cells"],
        ),
    )[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("data", "hcfl", "operators", "evaluate", "all"), default="all")
    parser.add_argument("--output", type=Path, default=HERE / "results" / "swe_radial")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--hcfl-max-updates", type=int, default=8000)
    parser.add_argument("--hcfl-validation-interval", type=int, default=500)
    parser.add_argument("--hcfl-batch-size", type=int, default=4)
    parser.add_argument("--operator-epochs", type=int, default=500)
    parser.add_argument("--operator-kind", choices=("all", "fno", "claw"), default="all")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print({"device": str(device), "output": str(args.output)}, flush=True)
    dataset = C.load_or_generate_dataset(args.output / "swe_radial_dataset.pt", device)
    train = dataset["train"]["trajectory"]
    validation = dataset["validation"]["trajectory"]
    mean, std, conserved_std, primitive_std = training_statistics(train)
    np.savez(
        args.output / "swe_radial_statistics.npz",
        primitive_mean=mean,
        primitive_std=std,
        conserved_std=conserved_std.numpy(),
        metric_primitive_std=primitive_std.numpy(),
    )
    if args.phase == "data":
        return

    hcfl_reports: list[dict[str, Any]] = []
    if args.phase in ("hcfl", "all"):
        for stencil in (4, 6):
            _, report = train_hcfl_arm(
                stencil,
                train,
                validation,
                mean,
                std,
                conserved_std,
                primitive_std,
                args.output,
                device,
                args.hcfl_max_updates,
                args.hcfl_validation_interval,
                args.hcfl_batch_size,
                args.seed,
            )
            hcfl_reports.append(report)
    else:
        for stencil in (4, 6):
            report_path = args.output / f"hcfl_s{stencil}_report_seed{args.seed}.json"
            if report_path.exists():
                hcfl_reports.append(json.loads(report_path.read_text(encoding="utf-8")))
    if args.phase == "hcfl":
        return

    operator_reports: dict[str, Any] = {}
    if args.phase in ("operators", "all"):
        kinds = ("fno", "claw") if args.operator_kind == "all" else (args.operator_kind,)
        for kind in kinds:
            _, operator_reports[kind] = train_operator(
                kind,
                train,
                validation,
                args.output,
                device,
                args.operator_epochs,
                args.seed,
            )
    else:
        for kind in ("fno", "claw"):
            report_path = args.output / f"{kind}_training_report_seed{args.seed}.json"
            if report_path.exists():
                operator_reports[kind] = json.loads(report_path.read_text(encoding="utf-8"))
    if args.phase == "operators":
        return

    if not hcfl_reports:
        raise FileNotFoundError("No HCFL reports are available")
    selected = choose_hcfl_report(hcfl_reports)
    selected_stencil = int(selected["stencil_cells"])
    checkpoints = {
        "hcfl": args.output / f"hcfl_s{selected_stencil}_best_seed{args.seed}.pt",
        "fno": args.output / f"fno_official_architecture_best_seed{args.seed}.pt",
        "claw": args.output / f"claw_official_architecture_best_seed{args.seed}.pt",
    }
    missing = [str(path) for path in checkpoints.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing checkpoints: {missing}")
    frame, summary = evaluate_all(
        dataset,
        args.output,
        mean,
        std,
        primitive_std,
        checkpoints["hcfl"],
        selected_stencil,
        checkpoints["fno"],
        checkpoints["claw"],
        device,
    )
    ledger = {
        "exploratory_design_record": "PROTOCOL.md",
        "official_clawno_commit": official_commit(),
        "hcfl_attempts": hcfl_reports,
        "selected_hcfl_stencil": selected_stencil,
        "operator_reports": operator_reports,
        "selection_used_validation_only": True,
        "test_sets_used_for_tuning": False,
        "higher_grid_attempted": False,
        "higher_grid_rule": (
            "Only if matched-grid HCFL is materially deficient; changing only HCFL after "
            "viewing results would invalidate the comparison."
        ),
    }
    save_json(args.output / "TUNING_LEDGER.json", ledger)
    print(json.dumps(summary["aggregate"], indent=2), flush=True)


if __name__ == "__main__":
    main()
