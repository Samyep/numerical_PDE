"""Train and evaluate the controlled 18-cell transverse HCFL ablation."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "euler_2d_clawpack_quadrants"
RESULTS = HERE / "results"
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(HERE))
import euler_2d_common as C  # noqa: E402
import run_experiment as B  # noqa: E402
import transverse_models as M  # noqa: E402


SNAPSHOT_DT = B.SNAPSHOT_DT
TRAINING_SUBSTEPS = B.TRAINING_SUBSTEPS
FEASIBILITY_WEIGHT = B.FEASIBILITY_WEIGHT
BATCH_SIZE = 4
INITIAL_LEARNING_RATE = 5.0e-4
MINIMUM_LEARNING_RATE = 2.0e-5


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


def result_path(variant: str, kind: str, seed: int, extension: str) -> Path:
    return RESULTS / f"{variant}_{kind}_seed{seed}.{extension}"


@torch.no_grad()
def transverse_response(
    model: torch.nn.Module,
    states: torch.Tensor,
    batch_size: int = 4,
) -> dict[str, float]:
    """Measure how strongly learned coefficients use transverse rows."""

    absolute_sum = 0.0
    coefficient_sum = 0.0
    gate_sum = 0.0
    count = 0
    model.eval()
    if not hasattr(model.flux_net, "coefficient_details"):
        return {
            "mean_absolute_transverse_coefficient_response": 0.0,
            "relative_transverse_coefficient_response": 0.0,
            "mean_transverse_gate": 0.0,
        }
    for start in range(0, states.shape[0], batch_size):
        batch = states[start : start + batch_size]
        for direction in ("x", "y"):
            oriented = C.orient_state(batch, direction)
            coefficients, centre_only, gate = model.flux_net.coefficient_details(
                oriented
            )
            absolute_sum += float((coefficients - centre_only).abs().sum())
            coefficient_sum += float(coefficients.abs().sum())
            gate_sum += float(gate.sum()) * coefficients.shape[-1]
            count += coefficients.numel()
    mean_response = absolute_sum / max(count, 1)
    return {
        "mean_absolute_transverse_coefficient_response": mean_response,
        "relative_transverse_coefficient_response": absolute_sum
        / max(coefficient_sum, 1.0e-12),
        "mean_transverse_gate": gate_sum / max(count, 1),
    }


@torch.no_grad()
def validation_tv_excess(
    prediction: torch.Tensor, reference: torch.Tensor
) -> float:
    pred_tv = C.total_variation_density(prediction[:, -1])
    ref_tv = C.total_variation_density(reference[:, -1])
    return float(
        (torch.relu(pred_tv - ref_tv) / ref_tv.clamp_min(1.0e-12)).double().mean()
    )


@torch.no_grad()
def validate(
    model: torch.nn.Module,
    validation: torch.Tensor,
    times: np.ndarray,
    primitive_std: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    prediction, diagnostics = B.rollout_hcfl(
        model, validation[:, 0].to(device), times
    )
    finite = torch.isfinite(prediction).all(dim=(-4, -3, -2, -1))
    positive_density = prediction[..., 0].amin(dim=(-3, -2, -1)) >= C.RHO_FLOOR
    positive_pressure = (
        C.pressure_raw(prediction).amin(dim=(-3, -2, -1)) >= C.PRESSURE_FLOOR
    )
    complete = finite & positive_density & positive_pressure
    if bool(complete.any()):
        nmae = C.primitive_nmae(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
        nrmse = C.primitive_nrmse(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
        tv_excess = validation_tv_excess(
            prediction[complete], validation[complete]
        )
    else:
        nmae = float("inf")
        nrmse = float("inf")
        tv_excess = float("inf")
    # Four saved times per trajectory are enough to track architectural use
    # without materially increasing validation cost.
    response_states = validation[:, ::5].reshape(-1, *validation.shape[2:]).to(device)
    response = transverse_response(model, response_states)
    return {
        "completion": float(complete.float().mean()),
        "rollout_nmae_completed": nmae,
        "rollout_nrmse_completed": nrmse,
        "final_density_tv_excess_completed": tv_excess,
        "minimum_density": float(prediction[..., 0].min()),
        "minimum_pressure": float(C.pressure_raw(prediction).min()),
        "maximum_entropy_balance": diagnostics["maximum_entropy_balance"],
        "maximum_interface_residual": diagnostics["maximum_interface_residual"],
        "minimum_positivity_theta": diagnostics["minimum_positivity_theta"],
        "minimum_entropy_beta": diagnostics["minimum_entropy_beta"],
        **response,
    }


def train_model(
    variant: str,
    seed: int,
    max_updates: int,
    validation_interval: int,
    device: torch.device,
) -> dict[str, Any]:
    if variant not in M.MODEL_TYPES:
        raise ValueError(variant)
    cells_per_face_input = 6 if variant == "normal6wide" else 18
    oriented_patch = [1, 6] if variant == "normal6wide" else [3, 6]
    train_archive = B.load_npz("train.npz")
    validation_archive = B.load_npz("validation.npz")
    training = torch.from_numpy(train_archive["q"]).float()
    validation = torch.from_numpy(validation_archive["q"]).float()
    times = validation_archive["times"]
    mean, std, conserved_std, primitive_std = B.training_statistics(training)

    # Data sampling is deliberately identical across variants for each seed.
    B.seed_everything(93100 + 100 * seed)
    model = M.make_model(mean, std, variant).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=INITIAL_LEARNING_RATE)
    generator = torch.Generator(device=device).manual_seed(93200 + 100 * seed)
    training_device = training.to(device)
    conserved_std_device = conserved_std.to(device)

    RESULTS.mkdir(parents=True, exist_ok=True)
    checkpoint_path = result_path(variant, "best", seed, "pt")
    curve_path = result_path(variant, "curve", seed, "csv")
    report_path = result_path(variant, "report", seed, "json")
    curve: list[dict[str, Any]] = []
    best_completion = -1.0
    best_nmae = float("inf")
    best_update = -1
    best_metrics: dict[str, float] = {}
    plateau_anchor = float("inf")
    stale_validations = 0
    stop_update = max_updates
    stop_reason = "maximum_updates"
    last_loss: float | None = None
    started = time.perf_counter()

    def check(update: int) -> None:
        nonlocal best_completion, best_nmae, best_update, best_metrics
        metrics = validate(model, validation, times, primitive_std, device)
        completion = metrics["completion"]
        nmae = metrics["rollout_nmae_completed"]
        improved = completion > best_completion or (
            completion == best_completion and nmae < best_nmae
        )
        if improved:
            best_completion = completion
            best_nmae = nmae
            best_update = update
            best_metrics = metrics
            torch.save(
                {
                    "state_dict": C.clone_state_dict(model),
                    "mean": mean,
                    "std": std,
                    "variant": variant,
                    "cells_per_face_input": cells_per_face_input,
                    "seed": seed,
                    "update": update,
                    "validation": metrics,
                },
                checkpoint_path,
            )
        row = {
            "variant": variant,
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
            (BATCH_SIZE,),
            generator=generator,
            device=device,
        )
        time_indices = torch.randint(
            training_device.shape[1] - 1,
            (BATCH_SIZE,),
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
            positivity = positivity + density_defect.square().mean()
            positivity = positivity + pressure_defect.square().mean()

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
        if old_learning_rate <= MINIMUM_LEARNING_RATE * (1.0 + 1.0e-12):
            stop_update = update
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_learning_rate = max(MINIMUM_LEARNING_RATE, 0.2 * old_learning_rate)
        for group in optimizer.param_groups:
            group["lr"] = new_learning_rate
        plateau_anchor = current_best
        stale_validations = 0
        print(
            {
                "event": "reduce_learning_rate",
                "variant": variant,
                "seed": seed,
                "update": update,
                "new_learning_rate": new_learning_rate,
            },
            flush=True,
        )

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    report = {
        "variant": variant,
        "method": (
            "normal-only capacity control"
            if variant == "normal6wide"
            else "18-cell transverse HLLC + signed Roe correction + hard entropy projection + proposal feasibility"
        ),
        "cells_per_face_input": cells_per_face_input,
        "oriented_patch": oriented_patch,
        "seed": seed,
        "parameters": C.parameter_count(model),
        "best_update": best_update,
        "best_validation_completion": best_completion,
        "best_validation_nmae": best_nmae,
        "best_validation_metrics": best_metrics,
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


def load_model(
    variant: str, seed: int, device: torch.device
) -> torch.nn.Module:
    checkpoint = torch.load(
        result_path(variant, "best", seed, "pt"),
        map_location=device,
        weights_only=False,
    )
    model = M.make_model(checkpoint["mean"], checkpoint["std"], variant).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    return model


def accuracy_rows(
    split: str,
    method: str,
    seed: int | None,
    prediction: torch.Tensor,
    reference: torch.Tensor,
    scale: torch.Tensor,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(prediction.shape[0]):
        audit = C.audit_accuracy(prediction[index], reference[index], scale)
        pred_density = prediction[index, -1, ..., 0].double()
        ref_density = reference[index, -1, ..., 0].double()
        pred_tv = C.total_variation_density(prediction[index, -1:]).double()[0]
        ref_tv = C.total_variation_density(reference[index, -1:]).double()[0]
        pred_curvature = (
            (
                pred_density[..., 2:]
                - 2.0 * pred_density[..., 1:-1]
                + pred_density[..., :-2]
            )
            .abs()
            .sum()
            + (
                pred_density[..., 2:, :]
                - 2.0 * pred_density[..., 1:-1, :]
                + pred_density[..., :-2, :]
            )
            .abs()
            .sum()
        )
        ref_curvature = (
            (
                ref_density[..., 2:]
                - 2.0 * ref_density[..., 1:-1]
                + ref_density[..., :-2]
            )
            .abs()
            .sum()
            + (
                ref_density[..., 2:, :]
                - 2.0 * ref_density[..., 1:-1, :]
                + ref_density[..., :-2, :]
            )
            .abs()
            .sum()
        )
        signed_tv_error = float((pred_tv - ref_tv) / ref_tv.clamp_min(1.0e-12))
        signed_curvature_error = float(
            (pred_curvature - ref_curvature) / ref_curvature.clamp_min(1.0e-12)
        )
        rows.append(
            {
                "split": split,
                "method": method,
                "seed": seed,
                "case": index,
                **asdict(audit),
                "density_tv_relative_error_signed": signed_tv_error,
                "density_curvature_relative_error_signed": signed_curvature_error,
                "density_curvature_excess": max(0.0, signed_curvature_error),
            }
        )
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    groups: defaultdict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(row["split"], row["method"])].append(row)
    metric_names = (
        "nmae",
        "nrmse",
        "density_mae",
        "x_velocity_mae",
        "y_velocity_mae",
        "pressure_mae",
        "density_tv_excess",
        "density_tv_relative_error_signed",
        "density_curvature_relative_error_signed",
        "density_curvature_excess",
        "density_range_overshoot",
        "pressure_range_overshoot",
    )
    for (split, method), group in groups.items():
        completed = [row for row in group if row["completed"]]
        summary: dict[str, Any] = {
            "rows": len(group),
            "completion_rate": len(completed) / len(group),
            "minimum_density": min(float(row["minimum_density"]) for row in group),
            "minimum_pressure": min(float(row["minimum_pressure"]) for row in group),
        }
        for name in metric_names:
            values = [float(row[name]) for row in completed if row[name] is not None]
            summary[f"{name}_mean_completed"] = (
                float(np.mean(values)) if values else None
            )
            summary[f"{name}_std_completed"] = (
                float(np.std(values)) if values else None
            )
        result[f"{split}/{method}"] = summary
    return result


def summarize_between_seeds(
    rows: list[dict[str, Any]], method: str
) -> dict[str, Any]:
    """Separate case variation from variation of per-seed mean accuracy."""

    metrics = (
        "nmae",
        "nrmse",
        "density_tv_excess",
        "density_tv_relative_error_signed",
        "density_curvature_relative_error_signed",
    )
    result: dict[str, Any] = {}
    for split in sorted({row["split"] for row in rows if row["method"] == method}):
        selected = [
            row
            for row in rows
            if row["split"] == split and row["method"] == method and row["completed"]
        ]
        seeds = sorted({int(row["seed"]) for row in selected})
        seed_means: dict[str, dict[str, float]] = {}
        for seed in seeds:
            seed_rows = [row for row in selected if int(row["seed"]) == seed]
            seed_means[str(seed)] = {
                metric: float(np.mean([float(row[metric]) for row in seed_rows]))
                for metric in metrics
            }
        aggregate: dict[str, Any] = {
            "seeds": seeds,
            "per_seed_means": seed_means,
        }
        for metric in metrics:
            values = [seed_means[str(seed)][metric] for seed in seeds]
            aggregate[f"{metric}_mean_between_seeds"] = float(np.mean(values))
            aggregate[f"{metric}_std_between_seeds"] = float(np.std(values))
        result[split] = aggregate
    return result


def safety_summary(diagnostics: dict[str, float]) -> dict[str, float]:
    batch_substeps = diagnostics.get("batch_substeps", 0.0)
    result = dict(diagnostics)
    result["positivity_fallback_rate"] = diagnostics.get(
        "positivity_limiter_active", 0.0
    ) / max(batch_substeps, 1.0)
    result["entropy_fallback_rate"] = diagnostics.get(
        "entropy_limiter_active", 0.0
    ) / max(batch_substeps, 1.0)
    return result


def make_official_figures(
    reference: torch.Tensor,
    native: torch.Tensor,
    normal: torch.Tensor,
    predictions: dict[str, torch.Tensor],
) -> None:
    methods = {
        "Reference-1024": reference,
        "PyClaw Roe-64": native,
        "Normal 6-cell HCFL (seed 0)": normal,
        "Normal 6-cell wide (seed 0)": predictions["normal6wide"],
        "Flat 18-cell HCFL (seed 0)": predictions["flat18"],
        "Gated 18-cell HCFL (seed 0)": predictions["gated18"],
    }
    density_arrays = {
        name: C.primitive(sequence[-1])[..., 0].numpy()
        for name, sequence in methods.items()
    }
    lower = min(float(value.min()) for value in density_arrays.values())
    upper = max(float(value.max()) for value in density_arrays.values())
    figure, axes_grid = plt.subplots(2, 3, figsize=(10.2, 6.7), constrained_layout=True)
    axes = axes_grid.reshape(-1)
    image = None
    for axis, (name, density) in zip(axes, density_arrays.items()):
        image = axis.imshow(
            density,
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            vmin=lower,
            vmax=upper,
            cmap="viridis",
            interpolation="nearest",
        )
        axis.set_title(name, fontsize=9)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    if image is not None:
        figure.colorbar(image, ax=axes, shrink=0.82, label="density")
    figure.savefig(RESULTS / "official_final_density.png", dpi=220)
    plt.close(figure)

    coordinates = (np.arange(reference.shape[-2]) + 0.5) / reference.shape[-2]
    middle = reference.shape[-2] // 2
    figure, axes = plt.subplots(2, 1, figsize=(8.5, 6.4), constrained_layout=True)
    styles = ("-", "--", "-.", (0, (6, 2)), (0, (3, 1)), (0, (1, 1)))
    for style, (name, density) in zip(styles, density_arrays.items()):
        axes[0].plot(coordinates, density[middle], linestyle=style, label=name)
        axes[1].plot(coordinates, density[:, middle], linestyle=style, label=name)
    axes[0].set_title("Final density: horizontal centre line")
    axes[1].set_title("Final density: vertical centre line")
    for axis in axes:
        axis.set_xlabel("coordinate")
        axis.set_ylabel("density")
        axis.grid(alpha=0.22)
    axes[0].legend(ncol=2, fontsize=8)
    figure.savefig(RESULTS / "official_density_linecuts.png", dpi=220)
    plt.close(figure)


def make_validation_curve_figure(variants: list[str], seed: int = 0) -> None:
    labels = {
        "normal6wide": "Normal 6-cell wide",
        "flat18": "Flat 18-cell",
        "gated18": "Gated 18-cell",
    }
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 3.8), constrained_layout=True)
    for variant in variants:
        path = result_path(variant, "curve", seed, "csv")
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        updates = [int(row["update"]) for row in rows]
        nmae = [float(row["rollout_nmae_completed"]) for row in rows]
        tv = [float(row["final_density_tv_excess_completed"]) for row in rows]
        axes[0].plot(updates, nmae, marker="o", markersize=2.5, label=labels[variant])
        axes[1].plot(updates, tv, marker="o", markersize=2.5, label=labels[variant])
    axes[0].axhline(0.029313357129620066, color="0.4", linestyle=":", label="Old normal-6 best")
    axes[0].set_ylabel("Validation rollout NMAE")
    axes[1].set_ylabel("Final density TV excess")
    for axis in axes:
        axis.set_xlabel("Optimizer update")
        axis.grid(alpha=0.22)
        axis.legend(fontsize=8)
    figure.savefig(RESULTS / "validation_curves.png", dpi=220)
    plt.close(figure)


def evaluate(
    variants: list[str], seeds: list[int], device: torch.device
) -> dict[str, Any]:
    training = torch.from_numpy(B.load_npz("train.npz")["q"]).float()
    scale = C.primitive_scale(training)
    test_archive = B.load_npz("test_id.npz")
    official_archive = B.load_npz("official_ref1024.npz")
    test_reference = torch.from_numpy(test_archive["q"]).float()
    official_reference = torch.from_numpy(official_archive["q"]).float()
    test_native = torch.from_numpy(test_archive["native_roe_coarse"]).float()
    official_native = torch.from_numpy(official_archive["native_roe_coarse"]).float()

    rows = accuracy_rows(
        "test_id", "PyClaw Roe-64", None, test_native, test_reference, scale
    )
    rows += accuracy_rows(
        "official_quadrants",
        "PyClaw Roe-64",
        None,
        official_native,
        official_reference,
        scale,
    )
    safety: dict[str, Any] = {}
    response: dict[str, Any] = {}
    timings: dict[str, float] = {}
    plot_predictions: dict[str, torch.Tensor] = {}

    for variant in variants:
        if variant not in M.MODEL_TYPES:
            raise ValueError(variant)
        variant_seeds = seeds if variant == "flat18" else [seeds[0]]
        for seed in variant_seeds:
            model = load_model(variant, seed, device)
            started = time.perf_counter()
            prediction_test, diagnostics_test = B.rollout_hcfl(
                model, test_reference[:, 0].to(device), test_archive["times"]
            )
            prediction_official, diagnostics_official = B.rollout_hcfl(
                model,
                official_reference[:, 0].to(device),
                official_archive["times"],
                batch_size=1,
            )
            if device.type == "cuda":
                torch.cuda.synchronize(device)
            label = f"{variant} seed {seed}"
            timings[label] = time.perf_counter() - started
            rows += accuracy_rows(
                "test_id", variant, seed, prediction_test, test_reference, scale
            )
            rows += accuracy_rows(
                "official_quadrants",
                variant,
                seed,
                prediction_official,
                official_reference,
                scale,
            )
            safety[f"{label}/test_id"] = safety_summary(diagnostics_test)
            safety[f"{label}/official_quadrants"] = safety_summary(
                diagnostics_official
            )
            response[label] = transverse_response(
                model,
                test_reference[:, ::5]
                .reshape(-1, *test_reference.shape[2:])
                .to(device),
            )
            np.savez_compressed(
                RESULTS / f"{variant}_predictions_seed{seed}.npz",
                test_id=prediction_test.numpy(),
                official=prediction_official.numpy(),
            )
            if seed == variant_seeds[0]:
                plot_predictions[variant] = prediction_official[0]

    normal_archive = np.load(
        BASE / "results" / "predictions_s6_seed0.npz", allow_pickle=False
    )
    normal_official = torch.from_numpy(normal_archive["official"])[0]
    summary = {
        "status": "exploratory; official test is post-hoc only",
        "variants": variants,
        "flat18_seeds": seeds,
        "control_seed": seeds[0],
        "accuracy": summarize(rows),
        "flat18_between_seed_accuracy": summarize_between_seeds(rows, "flat18"),
        "safety": safety,
        "transverse_response_on_test_id": response,
        "timing_seconds": timings,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    save_rows(RESULTS / "case_metrics.csv", rows)
    save_json(RESULTS / "summary.json", summary)
    make_official_figures(
        official_reference[0],
        official_native[0],
        normal_official,
        plot_predictions,
    )
    make_validation_curve_figure(variants, seeds[0])
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--variant", choices=M.MODEL_TYPES, required=True)
    train_parser.add_argument("--seed", type=int, required=True)
    train_parser.add_argument("--max-updates", type=int, default=50000)
    train_parser.add_argument("--validation-interval", type=int, default=500)
    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument(
        "--variants", choices=M.MODEL_TYPES, nargs="+", default=list(M.MODEL_TYPES)
    )
    evaluate_parser.add_argument("--seeds", type=int, nargs="+", default=[0])
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print({"command": args.command, "device": str(device)}, flush=True)
    if args.command == "train":
        report = train_model(
            args.variant,
            args.seed,
            args.max_updates,
            args.validation_interval,
            device,
        )
        print(json.dumps(report, indent=2), flush=True)
    elif args.command == "evaluate":
        evaluate(args.variants, args.seeds, device)


if __name__ == "__main__":
    main()
