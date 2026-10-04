"""Train, select, and audit the prospective 2-D Euler HCFL experiment."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import euler_2d_common as C  # noqa: E402


SNAPSHOT_DT = 0.005
TRAINING_SUBSTEPS = 6
FEASIBILITY_WEIGHT = 1.0e-3
PRIMITIVE_NAMES = ("density", "x_velocity", "y_velocity", "pressure")


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def save_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def save_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot save an empty table")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_npz(name: str) -> dict[str, Any]:
    path = DATA / name
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def training_statistics(
    training: torch.Tensor,
) -> tuple[np.ndarray, np.ndarray, torch.Tensor, torch.Tensor]:
    values = C.primitive(training)
    mean = values.mean(dim=(0, 1, 2, 3)).numpy().astype(np.float32)
    std = (
        values.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-4).numpy().astype(np.float32)
    )
    conserved_std = training.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-4)
    primitive_std = C.primitive_scale(training)
    return mean, std, conserved_std, primitive_std


def _initial_diagnostics() -> defaultdict[str, float]:
    values: defaultdict[str, float] = defaultdict(float)
    values["minimum_positivity_theta"] = 1.0
    values["minimum_entropy_beta"] = 1.0
    values["maximum_entropy_balance"] = -float("inf")
    values["maximum_interface_residual"] = -float("inf")
    values["maximum_conservation_closure"] = 0.0
    return values


@torch.no_grad()
def rollout_hcfl(
    model: C.HCFL2DEuler,
    initial: torch.Tensor,
    times: np.ndarray,
    batch_size: int = 4,
) -> tuple[torch.Tensor, dict[str, float]]:
    model.eval()
    all_sequences: list[torch.Tensor] = []
    aggregate = _initial_diagnostics()
    for start in range(0, initial.shape[0], batch_size):
        state = initial[start : start + batch_size].clone()
        sequence = [state.cpu()]
        for index in range(1, len(times)):
            interval = float(times[index] - times[index - 1])
            state, diagnostics = C.advance_hcfl_interval(model, state, interval)
            C.merge_statistics(aggregate, diagnostics)
            sequence.append(state.cpu())
        all_sequences.append(torch.stack(sequence, dim=1))
    return torch.cat(all_sequences, dim=0), dict(aggregate)


@torch.no_grad()
def rollout_hllc(
    initial: torch.Tensor,
    times: np.ndarray,
    batch_size: int = 6,
) -> tuple[torch.Tensor, dict[str, float]]:
    all_sequences: list[torch.Tensor] = []
    aggregate: defaultdict[str, float] = defaultdict(float)
    aggregate["maximum_conservation_closure"] = 0.0
    for start in range(0, initial.shape[0], batch_size):
        state = initial[start : start + batch_size].clone()
        sequence = [state.cpu()]
        for index in range(1, len(times)):
            interval = float(times[index] - times[index - 1])
            state, diagnostics = C.advance_hllc_interval(state, interval)
            C.merge_statistics(aggregate, diagnostics)
            sequence.append(state.cpu())
        all_sequences.append(torch.stack(sequence, dim=1))
    return torch.cat(all_sequences, dim=0), dict(aggregate)


@torch.no_grad()
def validate_hcfl(
    model: C.HCFL2DEuler,
    validation: torch.Tensor,
    times: np.ndarray,
    primitive_std: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    prediction, diagnostics = rollout_hcfl(
        model, validation[:, 0].to(device), times
    )
    finite = torch.isfinite(prediction).all(dim=(-4, -3, -2, -1))
    positive_density = prediction[..., 0].amin(dim=(-3, -2, -1)) >= C.RHO_FLOOR
    positive_pressure = C.pressure_raw(prediction).amin(dim=(-3, -2, -1)) >= C.PRESSURE_FLOOR
    complete = finite & positive_density & positive_pressure
    if bool(complete.any()):
        nmae = C.primitive_nmae(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
        nrmse = C.primitive_nrmse(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
    else:
        nmae = float("inf")
        nrmse = float("inf")
    return {
        "completion": float(complete.float().mean()),
        "rollout_nmae_completed": nmae,
        "rollout_nrmse_completed": nrmse,
        "minimum_density": float(prediction[..., 0].min()),
        "minimum_pressure": float(C.pressure_raw(prediction).min()),
        "maximum_entropy_balance": diagnostics["maximum_entropy_balance"],
        "maximum_interface_residual": diagnostics["maximum_interface_residual"],
        "minimum_positivity_theta": diagnostics["minimum_positivity_theta"],
        "minimum_entropy_beta": diagnostics["minimum_entropy_beta"],
    }


def train_model(
    stencil: int,
    seed: int,
    max_updates: int,
    validation_interval: int,
    device: torch.device,
) -> dict[str, Any]:
    train_archive = load_npz("train.npz")
    validation_archive = load_npz("validation.npz")
    training = torch.from_numpy(train_archive["q"]).float()
    validation = torch.from_numpy(validation_archive["q"]).float()
    times = validation_archive["times"]
    mean, std, conserved_std, primitive_std = training_statistics(training)

    initialization_seed = 91000 + 100 * seed + stencil
    seed_everything(initialization_seed)
    model = C.HCFL2DEuler(mean, std, width=72, stencil_cells=stencil).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=5.0e-4)
    generator = torch.Generator(device=device).manual_seed(92000 + 100 * seed + stencil)
    training_device = training.to(device)
    conserved_std_device = conserved_std.to(device)

    RESULTS.mkdir(parents=True, exist_ok=True)
    checkpoint_path = RESULTS / f"hcfl_s{stencil}_best_seed{seed}.pt"
    curve_path = RESULTS / f"hcfl_s{stencil}_curve_seed{seed}.csv"
    report_path = RESULTS / f"hcfl_s{stencil}_report_seed{seed}.json"
    curve: list[dict[str, Any]] = []
    best_completion = -1.0
    best_nmae = float("inf")
    best_update = -1
    plateau_anchor = float("inf")
    stale_validations = 0
    stop_update = max_updates
    stop_reason = "maximum_updates"
    last_loss: float | None = None
    started = time.perf_counter()

    def check(update: int) -> None:
        nonlocal best_completion, best_nmae, best_update
        metrics = validate_hcfl(model, validation, times, primitive_std, device)
        completion = metrics["completion"]
        nmae = metrics["rollout_nmae_completed"]
        improved = completion > best_completion or (
            completion == best_completion and nmae < best_nmae
        )
        if improved:
            best_completion = completion
            best_nmae = nmae
            best_update = update
            torch.save(
                {
                    "state_dict": C.clone_state_dict(model),
                    "mean": mean,
                    "std": std,
                    "stencil_cells": stencil,
                    "seed": seed,
                    "update": update,
                    "validation": metrics,
                },
                checkpoint_path,
            )
        row = {
            "stencil_cells": stencil,
            "seed": seed,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": last_loss,
            "new_best": improved,
            **metrics,
        }
        curve.append(row)
        save_rows(curve_path, curve)
        print(row, flush=True)
        model.train()

    check(0)
    for update in range(1, max_updates + 1):
        trajectory_indices = torch.randint(
            training_device.shape[0],
            (4,),
            generator=generator,
            device=device,
        )
        time_indices = torch.randint(
            training_device.shape[1] - 1,
            (4,),
            generator=generator,
            device=device,
        )
        state = training_device[trajectory_indices, time_indices]
        target = training_device[trajectory_indices, time_indices + 1]
        feasibility = torch.zeros((), device=device)
        positivity = torch.zeros((), device=device)
        dt = SNAPSHOT_DT / TRAINING_SUBSTEPS
        for _ in range(TRAINING_SUBSTEPS):
            raw_fluxes = model.raw_fluxes(state)
            feasibility = feasibility + model.feasibility_loss_from_raw(
                state, raw_fluxes
            )
            projected_fluxes = model.project_raw_fluxes_training(state, raw_fluxes)
            state = C.flux_divergence(state, *projected_fluxes, dt)
            density_defect = torch.relu(C.RHO_FLOOR - state[..., 0])
            pressure_defect = torch.relu(C.PRESSURE_FLOOR - C.pressure_raw(state))
            positivity = positivity + density_defect.square().mean() + pressure_defect.square().mean()

        trajectory_loss = (((state - target) / conserved_std_device) ** 2).mean()
        loss = (
            trajectory_loss
            + FEASIBILITY_WEIGHT * feasibility / TRAINING_SUBSTEPS
            + 10.0 * positivity / TRAINING_SUBSTEPS
        )
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach())

        if update % validation_interval and update != max_updates:
            continue
        check(update)
        if update < 2000:
            continue
        current_best = best_nmae
        if not math.isfinite(plateau_anchor):
            plateau_anchor = current_best
            continue
        if current_best <= plateau_anchor * 0.998:
            plateau_anchor = current_best
            stale_validations = 0
            continue
        stale_validations += 1
        if stale_validations < 3:
            continue
        old_learning_rate = float(optimizer.param_groups[0]["lr"])
        if old_learning_rate <= 2.0e-5 * (1.0 + 1.0e-12):
            stop_update = update
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_learning_rate = max(2.0e-5, 0.2 * old_learning_rate)
        for group in optimizer.param_groups:
            group["lr"] = new_learning_rate
        plateau_anchor = current_best
        stale_validations = 0
        print(
            {
                "event": "reduce_learning_rate",
                "stencil": stencil,
                "seed": seed,
                "update": update,
                "new_learning_rate": new_learning_rate,
            },
            flush=True,
        )

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    report = {
        "method": "HLLC + signed Roe correction + hard entropy projection + proposal feasibility",
        "stencil_cells": stencil,
        "seed": seed,
        "parameters": C.parameter_count(model),
        "best_update": best_update,
        "best_validation_completion": best_completion,
        "best_validation_nmae": best_nmae,
        "stop_update": stop_update,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "maximum_updates": max_updates,
        "validation_interval": validation_interval,
        "training_substeps_per_saved_interval": TRAINING_SUBSTEPS,
        "proposal_feasibility_weight": FEASIBILITY_WEIGHT,
        "uses_low_order_flux_during_training": False,
        "uses_low_order_flux_only_in_deployment_safety_wrapper": True,
        "checkpoint": checkpoint_path.name,
    }
    save_json(report_path, report)
    return report


def select_stencil() -> dict[str, Any]:
    reports = [
        json.loads((RESULTS / f"hcfl_s{stencil}_report_seed0.json").read_text(encoding="utf-8"))
        for stencil in (4, 6)
    ]
    reports.sort(
        key=lambda row: (
            -float(row["best_validation_completion"]),
            float(row["best_validation_nmae"]),
        )
    )
    selection = {
        "rule": "completion_then_validation_NMAE",
        "selected_stencil_cells": reports[0]["stencil_cells"],
        "candidates": reports,
        "official_test_was_evaluated": False,
    }
    save_json(RESULTS / "stencil_selection.json", selection)
    print(json.dumps(selection, indent=2), flush=True)
    return selection


def load_model(stencil: int, seed: int, device: torch.device) -> C.HCFL2DEuler:
    checkpoint = torch.load(
        RESULTS / f"hcfl_s{stencil}_best_seed{seed}.pt",
        map_location="cpu",
        weights_only=False,
    )
    model = C.HCFL2DEuler(
        checkpoint["mean"],
        checkpoint["std"],
        width=72,
        stencil_cells=stencil,
    ).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    return model


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _metric_pair(
    first: torch.Tensor, second: torch.Tensor, scale: torch.Tensor
) -> dict[str, float]:
    return {
        "nmae": C.primitive_nmae(first, second, scale),
        "nrmse": C.primitive_nrmse(first, second, scale),
    }


def audit_data() -> dict[str, Any]:
    names = (
        "train.npz",
        "validation.npz",
        "test_id.npz",
        "official_ref256.npz",
        "official_ref512.npz",
        "official_ref1024.npz",
    )
    archives = {name: load_npz(name) for name in names}
    training = torch.from_numpy(archives["train.npz"]["q"]).float()
    scale = C.primitive_scale(training)
    current_generator_hash = _sha256(HERE / "generate_pyclaw_dataset.py")
    files: dict[str, Any] = {}
    for name in names:
        archive = archives[name]
        state = torch.from_numpy(archive["q"]).float()
        metadata = json.loads(str(archive["metadata_json"]))
        files[name] = {
            "sha256": _sha256(DATA / name),
            "shape": list(state.shape),
            "minimum_density": float(state[..., 0].min()),
            "minimum_pressure": float(C.pressure_raw(state).min()),
            "finite": bool(torch.isfinite(state).all()),
            "generator_hash_matches": metadata["generator_sha256"]
            == current_generator_hash,
            "metadata": metadata,
        }

    def parameter_keys(name: str) -> set[bytes]:
        values = archives[name]["primitive_states"]
        splits = archives[name]["split_locations"]
        return {
            np.concatenate((state.reshape(-1), split)).tobytes()
            for state, split in zip(values, splits)
        }

    train_keys = parameter_keys("train.npz")
    validation_keys = parameter_keys("validation.npz")
    test_keys = parameter_keys("test_id.npz")
    official256 = torch.from_numpy(archives["official_ref256.npz"]["q"][0]).float()
    official512 = torch.from_numpy(archives["official_ref512.npz"]["q"][0]).float()
    official1024 = torch.from_numpy(archives["official_ref1024.npz"]["q"][0]).float()
    native = torch.from_numpy(
        archives["official_ref1024.npz"]["native_roe_coarse"][0]
    ).float()
    audit = {
        "files": files,
        "training_primitive_scale": {
            name: float(value) for name, value in zip(PRIMITIVE_NAMES, scale)
        },
        "splits_are_disjoint": {
            "train_validation": not bool(train_keys & validation_keys),
            "train_test_id": not bool(train_keys & test_keys),
            "validation_test_id": not bool(validation_keys & test_keys),
        },
        "reference_convergence": {
            "256_vs_1024_rollout": _metric_pair(official256, official1024, scale),
            "512_vs_1024_rollout": _metric_pair(official512, official1024, scale),
            "256_vs_1024_final": _metric_pair(official256[-1:], official1024[-1:], scale),
            "512_vs_1024_final": _metric_pair(official512[-1:], official1024[-1:], scale),
        },
        "native_roe64_vs_1024": {
            "rollout": _metric_pair(native[1:], official1024[1:], scale),
            "final": _metric_pair(native[-1:], official1024[-1:], scale),
        },
        "common_initial_state_maximum_absolute_difference": float(
            (native[0] - official1024[0]).abs().max()
        ),
    }
    required = [
        all(item["finite"] for item in files.values()),
        all(item["minimum_density"] > 0.0 for item in files.values()),
        all(item["minimum_pressure"] > 0.0 for item in files.values()),
        all(item["generator_hash_matches"] for item in files.values()),
        all(audit["splits_are_disjoint"].values()),
        audit["common_initial_state_maximum_absolute_difference"] == 0.0,
    ]
    audit["passed"] = all(required)
    save_json(RESULTS / "data_and_reference_audit.json", audit)
    print(json.dumps(audit, indent=2), flush=True)
    if not audit["passed"]:
        raise RuntimeError("Data/reference audit failed")
    return audit


def accuracy_rows(
    split: str,
    method: str,
    seed: int | None,
    prediction: torch.Tensor,
    reference: torch.Tensor,
    scale: torch.Tensor,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(reference.shape[0]):
        audit = C.audit_accuracy(prediction[index], reference[index], scale)
        rows.append(
            {
                "split": split,
                "trajectory": index,
                "method": method,
                "seed": seed,
                **asdict(audit),
            }
        )
    return rows


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["split"]), str(row["method"]))].append(row)
    summary: dict[str, Any] = {}
    for (split, method), selected in grouped.items():
        completed = [row for row in selected if bool(row["completed"])]
        key = f"{split}/{method}"
        summary[key] = {
            "rows": len(selected),
            "completion_rate": len(completed) / len(selected),
            "nmae_mean_completed": (
                float(np.mean([row["nmae"] for row in completed]))
                if completed
                else None
            ),
            "nmae_std_completed": (
                float(np.std([row["nmae"] for row in completed]))
                if completed
                else None
            ),
            "nrmse_mean_completed": (
                float(np.mean([row["nrmse"] for row in completed]))
                if completed
                else None
            ),
            "minimum_density": min(float(row["minimum_density"]) for row in selected),
            "minimum_pressure": min(float(row["minimum_pressure"]) for row in selected),
            "density_tv_excess_mean_completed": (
                float(np.mean([row["density_tv_excess"] for row in completed]))
                if completed
                else None
            ),
            "density_range_overshoot_mean_completed": (
                float(np.mean([row["density_range_overshoot"] for row in completed]))
                if completed
                else None
            ),
            "pressure_range_overshoot_mean_completed": (
                float(np.mean([row["pressure_range_overshoot"] for row in completed]))
                if completed
                else None
            ),
        }
        seeded = [row for row in completed if row["seed"] not in (None, "")]
        seed_means: dict[int, list[float]] = defaultdict(list)
        for row in seeded:
            seed_means[int(row["seed"])].append(float(row["nmae"]))
        means = [float(np.mean(values)) for values in seed_means.values()]
        summary[key]["nmae_between_seed_std"] = (
            float(np.std(means)) if len(means) > 1 else None
        )
        for field in (
            "density_mae",
            "x_velocity_mae",
            "y_velocity_mae",
            "pressure_mae",
        ):
            summary[key][f"{field}_mean_completed"] = (
                float(np.mean([row[field] for row in completed])) if completed else None
            )
    return summary


def _safety_summary(diagnostics: dict[str, float]) -> dict[str, float]:
    result = dict(diagnostics)
    denominator = max(result.get("batch_substeps", 0.0), 1.0)
    result["positivity_fallback_rate"] = result.get(
        "positivity_limiter_active", 0.0
    ) / denominator
    result["entropy_fallback_rate"] = result.get(
        "entropy_limiter_active", 0.0
    ) / denominator
    return result


def make_figures(
    reference: torch.Tensor,
    native: torch.Tensor,
    hllc: torch.Tensor,
    hcfl: torch.Tensor,
) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    states = [reference[-1], native[-1], hllc[-1], hcfl[-1]]
    titles = [
        "PyClaw Roe-1024 reference",
        "PyClaw Roe-64",
        "HLLC-64",
        "HCFL-64 (seed 0)",
    ]
    density_min = min(float(state[..., 0].min()) for state in states)
    density_max = max(float(state[..., 0].max()) for state in states)
    figure, axes = plt.subplots(1, 4, figsize=(15.2, 3.7), constrained_layout=True)
    image = None
    for axis, state, title in zip(axes, states, titles):
        image = axis.imshow(
            state[..., 0].numpy(),
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            cmap="viridis",
            vmin=density_min,
            vmax=density_max,
            interpolation="nearest",
        )
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("x")
    axes[0].set_ylabel("y")
    assert image is not None
    figure.colorbar(image, ax=axes, label="density", shrink=0.86)
    figure.savefig(RESULTS / "official_final_density.png", dpi=220)
    plt.close(figure)

    coordinates = (np.arange(C.N_COARSE) + 0.5) / C.N_COARSE
    line_indices = (32, 48)
    figure, axes = plt.subplots(2, 1, figsize=(8.2, 6.3), sharex=True, constrained_layout=True)
    styles = [
        ("black", "-", 2.2),
        ("#d55e00", (0, (6, 3)), 1.8),
        ("#0072b2", (0, (3, 2)), 1.8),
        ("#009e73", (0, (8, 2, 2, 2)), 2.0),
    ]
    for axis, row in zip(axes, line_indices):
        for state, title, style in zip(states, titles, styles):
            axis.plot(
                coordinates,
                state[row, :, 0].numpy(),
                color=style[0],
                linestyle=style[1],
                linewidth=style[2],
                label=title,
            )
        axis.set_ylabel("density")
        axis.set_title(f"horizontal cut y={(row + 0.5) / C.N_COARSE:.3f}")
        axis.grid(alpha=0.2)
    axes[-1].set_xlabel("x")
    axes[0].legend(ncol=2, fontsize=8)
    figure.savefig(RESULTS / "official_density_linecuts.png", dpi=220)
    plt.close(figure)


def make_nmae_figure(rows: list[dict[str, Any]]) -> None:
    methods = ("PyClaw Roe-64", "HLLC-64", "HCFL-64")
    splits = ("test_id", "official_quadrants")
    labels = ("Held-out ID", "Official quadrants OOD")
    values = np.zeros((len(splits), len(methods)), dtype=np.float64)
    errors = np.zeros_like(values)
    for split_index, split in enumerate(splits):
        for method_index, method in enumerate(methods):
            selected = [
                row
                for row in rows
                if row["split"] == split
                and row["method"] == method
                and bool(row["completed"])
            ]
            if method == "HCFL-64":
                by_seed: dict[int, list[float]] = defaultdict(list)
                for row in selected:
                    by_seed[int(row["seed"])].append(float(row["nmae"]))
                seed_means = np.asarray(
                    [np.mean(item) for item in by_seed.values()], dtype=np.float64
                )
                values[split_index, method_index] = seed_means.mean()
                errors[split_index, method_index] = seed_means.std()
            else:
                values[split_index, method_index] = np.mean(
                    [float(row["nmae"]) for row in selected]
                )

    colors = ("#d55e00", "#0072b2", "#009e73")
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
    for split_index, axis in enumerate(axes):
        positions = np.arange(len(methods))
        axis.bar(
            positions,
            values[split_index],
            yerr=errors[split_index],
            color=colors,
            capsize=4,
            edgecolor="black",
            linewidth=0.7,
        )
        axis.set_xticks(positions, methods, rotation=15, ha="right")
        axis.set_ylabel("NMAE")
        axis.set_title(labels[split_index])
        axis.grid(axis="y", alpha=0.25)
        for position, value in zip(positions, values[split_index]):
            axis.text(position, value, f"{value:.4f}", ha="center", va="bottom", fontsize=8)
    figure.savefig(RESULTS / "nmae_comparison.png", dpi=220)
    plt.close(figure)


def evaluate(stencil: int, seeds: list[int], device: torch.device) -> dict[str, Any]:
    training = torch.from_numpy(load_npz("train.npz")["q"]).float()
    scale = C.primitive_scale(training)
    test_archive = load_npz("test_id.npz")
    official_archive = load_npz("official_ref1024.npz")
    test_reference = torch.from_numpy(test_archive["q"]).float()
    official_reference = torch.from_numpy(official_archive["q"]).float()
    test_native = torch.from_numpy(test_archive["native_roe_coarse"]).float()
    official_native = torch.from_numpy(official_archive["native_roe_coarse"]).float()
    rows: list[dict[str, Any]] = []
    rows += accuracy_rows("test_id", "PyClaw Roe-64", None, test_native, test_reference, scale)
    rows += accuracy_rows(
        "official_quadrants", "PyClaw Roe-64", None, official_native, official_reference, scale
    )

    started = time.perf_counter()
    hllc_test, hllc_test_diagnostics = rollout_hllc(
        test_reference[:, 0].to(device), test_archive["times"]
    )
    hllc_official, hllc_official_diagnostics = rollout_hllc(
        official_reference[:, 0].to(device), official_archive["times"]
    )
    hllc_seconds = time.perf_counter() - started
    rows += accuracy_rows("test_id", "HLLC-64", None, hllc_test, test_reference, scale)
    rows += accuracy_rows(
        "official_quadrants", "HLLC-64", None, hllc_official, official_reference, scale
    )

    safety: dict[str, Any] = {
        "HLLC-64/test_id": hllc_test_diagnostics,
        "HLLC-64/official_quadrants": hllc_official_diagnostics,
    }
    timing: dict[str, float] = {"HLLC-64": hllc_seconds}
    plot_hcfl: torch.Tensor | None = None
    for seed in seeds:
        model = load_model(stencil, seed, device)
        started = time.perf_counter()
        prediction_test, diagnostics_test = rollout_hcfl(
            model, test_reference[:, 0].to(device), test_archive["times"]
        )
        prediction_official, diagnostics_official = rollout_hcfl(
            model,
            official_reference[:, 0].to(device),
            official_archive["times"],
            batch_size=1,
        )
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        timing[f"HCFL-64 seed {seed}"] = time.perf_counter() - started
        rows += accuracy_rows(
            "test_id", "HCFL-64", seed, prediction_test, test_reference, scale
        )
        rows += accuracy_rows(
            "official_quadrants",
            "HCFL-64",
            seed,
            prediction_official,
            official_reference,
            scale,
        )
        safety[f"HCFL-64 seed {seed}/test_id"] = _safety_summary(diagnostics_test)
        safety[f"HCFL-64 seed {seed}/official_quadrants"] = _safety_summary(
            diagnostics_official
        )
        np.savez_compressed(
            RESULTS / f"predictions_s{stencil}_seed{seed}.npz",
            test_id=prediction_test.numpy(),
            official=prediction_official.numpy(),
        )
        if plot_hcfl is None:
            plot_hcfl = prediction_official[0]

    save_rows(RESULTS / "case_metrics.csv", rows)
    summary = {
        "selected_stencil_cells": stencil,
        "seeds": seeds,
        "accuracy": summarize_rows(rows),
        "safety": safety,
        "timing_seconds": timing,
        "training_primitive_scale": {
            name: float(value) for name, value in zip(PRIMITIVE_NAMES, scale)
        },
    }
    save_json(RESULTS / "summary.json", summary)
    if plot_hcfl is None:
        raise RuntimeError("No HCFL seed was evaluated")
    make_figures(
        official_reference[0], official_native[0], hllc_official[0], plot_hcfl
    )
    make_nmae_figure(rows)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("audit-data")

    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--stencil", type=int, choices=(4, 6), required=True)
    train_parser.add_argument("--seed", type=int, required=True)
    train_parser.add_argument("--max-updates", type=int, default=50000)
    train_parser.add_argument("--validation-interval", type=int, default=500)

    subparsers.add_parser("select")
    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--stencil", type=int, choices=(4, 6), required=True)
    evaluate_parser.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print({"command": args.command, "device": str(device)}, flush=True)
    if args.command == "audit-data":
        audit_data()
    elif args.command == "train":
        train_model(
            args.stencil,
            args.seed,
            args.max_updates,
            args.validation_interval,
            device,
        )
    elif args.command == "select":
        select_stencil()
    elif args.command == "evaluate":
        evaluate(args.stencil, args.seeds, device)


if __name__ == "__main__":
    main()
