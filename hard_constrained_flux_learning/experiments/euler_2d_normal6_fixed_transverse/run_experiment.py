"""Train and audit six-cell normal HCFL with fixed transverse Roe transport.

The only learned object is the shared six-cell normal flux used in both
coordinate directions.  A parameter-free Clawpack-style transverse
increment-wave split supplies the genuinely two-dimensional corner transport.
Training, validation stopping, data splits, and deployment safety match the
fixed-64 diverse Euler benchmark.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
DIVERSE = HERE.parent / "euler_2d_diverse64"
TRANSVERSE = HERE.parent / "euler_2d_transverse_ablation"
sys.path.insert(0, str(DIVERSE))
sys.path.insert(0, str(TRANSVERSE))
import run_diverse64 as R  # noqa: E402
import run_consistent64 as S  # noqa: E402
import fixed_transverse as T  # noqa: E402


C = R.C
B = R.B
METHOD_ID = T.METHOD_ID
METHOD_LABEL = "HCFL-6 + fixed transverse Roe"
COMPARATOR_LABEL = "HCFL-18 learned transverse"
INFERENCE_ABLATION_LABEL = "HCFL-6 transverse disabled after training"
RESULTS = HERE / "results"
DEFAULT_DATA = Path(
    os.environ.get("HCFL_DIVERSE64_DATA", str(DIVERSE / "data"))
)


def result_path(kind: str, seed: int, extension: str) -> Path:
    return RESULTS / f"{METHOD_ID}_{kind}_seed{seed}.{extension}"


def train_model(
    seed: int,
    max_updates: int,
    validation_interval: int,
    device: torch.device,
    data: Path,
) -> dict[str, Any]:
    train_archive = R.load_npz(data / "train.npz")
    validation_archive = R.load_npz(data / "validation.npz")
    training = torch.from_numpy(train_archive["q"]).float()
    validation = torch.from_numpy(validation_archive["q"]).float()
    family_names = validation_archive["family_names"]
    times = validation_archive["times"]
    mean, std, conserved_std, primitive_std = B.training_statistics(training)

    B.seed_everything(95100 + 100 * seed)
    model = T.make_model(mean, std).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=R.INITIAL_LEARNING_RATE)
    generator = torch.Generator(device=device).manual_seed(95200 + 100 * seed)
    training_device = training.to(device)
    conserved_std_device = conserved_std.to(device)

    RESULTS.mkdir(parents=True, exist_ok=True)
    checkpoint_path = result_path("best", seed, "pt")
    curve_path = result_path("curve", seed, "csv")
    report_path = result_path("report", seed, "json")
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
        metrics = R.validate(
            model,
            validation,
            family_names,
            times,
            primitive_std,
            device,
        )
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
                    "variant": METHOD_ID,
                    "seed": seed,
                    "update": update,
                    "validation": metrics,
                    "transverse_level": T.TRANSVERSE_LEVEL,
                    "train_data_sha256": R._sha256(data / "train.npz"),
                    "validation_data_sha256": R._sha256(
                        data / "validation.npz"
                    ),
                },
                checkpoint_path,
            )
        row = {
            "seed": seed,
            "update": update,
            "learning_rate": float(optimizer.param_groups[0]["lr"]),
            "last_minibatch_loss": last_loss,
            "new_best": improved,
            **metrics,
        }
        curve.append(row)
        R.save_rows(curve_path, curve)
        print(row, flush=True)

    check(0)
    for update in range(1, max_updates + 1):
        trajectory_indices = torch.randint(
            training_device.shape[0],
            (R.BATCH_SIZE,),
            generator=generator,
            device=device,
        )
        time_indices = torch.randint(
            training_device.shape[1] - 1,
            (R.BATCH_SIZE,),
            generator=generator,
            device=device,
        )
        state = training_device[trajectory_indices, time_indices]
        target = training_device[trajectory_indices, time_indices + 1]
        feasibility = torch.zeros((), device=device)
        positivity = torch.zeros((), device=device)
        dt = B.SNAPSHOT_DT / B.TRAINING_SUBSTEPS
        for _ in range(B.TRAINING_SUBSTEPS):
            raw_fluxes = model.raw_fluxes(state, dt=dt)
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
            + B.FEASIBILITY_WEIGHT * feasibility / B.TRAINING_SUBSTEPS
            + 10.0 * positivity / B.TRAINING_SUBSTEPS
        )
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        last_loss = float(loss.detach())

        if update % validation_interval and update != max_updates:
            continue
        check(update)
        if update < 3000:
            continue
        current_best = best_nmae
        if not math.isfinite(plateau_anchor):
            plateau_anchor = current_best
            continue
        if current_best <= plateau_anchor * (1.0 - R.RELATIVE_IMPROVEMENT):
            plateau_anchor = current_best
            stale_validations = 0
            continue
        stale_validations += 1
        if stale_validations < R.PLATEAU_VALIDATIONS:
            continue
        old_learning_rate = float(optimizer.param_groups[0]["lr"])
        if old_learning_rate <= R.MINIMUM_LEARNING_RATE * (1.0 + 1.0e-12):
            stop_update = update
            stop_reason = "validation_plateau_at_minimum_learning_rate"
            break
        new_learning_rate = max(
            R.MINIMUM_LEARNING_RATE, 0.2 * old_learning_rate
        )
        for group in optimizer.param_groups:
            group["lr"] = new_learning_rate
        plateau_anchor = current_best
        stale_validations = 0
        print(
            {
                "event": "reduce_learning_rate",
                "seed": seed,
                "update": update,
                "new_learning_rate": new_learning_rate,
            },
            flush=True,
        )

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    report = {
        "method": (
            "shared normal six-cell central physical flux + nonnegative "
            "entropy-fixed Roe dissipation + fixed Roe transverse increment "
            "transport + proposal feasibility + hard interface entropy projection"
        ),
        "transverse_transport": (
            "Clawpack classic transverse increment-wave semantics "
            "(transverse_waves=1); no second-order correction waves"
        ),
        "why_not_transverse_waves_2": (
            "level 2 also propagates an explicit second-order normal correction "
            "that is absent from the retained 1-D HCFL method"
        ),
        "deployment_grid": [64, 64],
        "cross_grid_claim": False,
        "reference_fine_cells": R._metadata(train_archive)["fine_cells"],
        "train_trajectories": int(training.shape[0]),
        "validation_trajectories": int(validation.shape[0]),
        "families": list(R.EXPECTED_FAMILIES),
        "seed": seed,
        "parameters": C.parameter_count(model),
        "transverse_trainable_parameters": 0,
        "best_update": best_update,
        "best_validation_completion": best_completion,
        "best_validation_nmae": best_nmae,
        "best_validation_metrics": best_metrics,
        "stop_update": stop_update,
        "stop_reason": stop_reason,
        "training_seconds": time.perf_counter() - started,
        "maximum_updates": max_updates,
        "validation_interval": validation_interval,
        "plateau_validations": R.PLATEAU_VALIDATIONS,
        "training_substeps_per_saved_interval": B.TRAINING_SUBSTEPS,
        "proposal_feasibility_weight": B.FEASIBILITY_WEIGHT,
        "roe_multiplier_parameterization": "1 + tanh(logit), range [0, 2]",
        "uses_low_order_flux_during_training": False,
        "uses_low_order_flux_only_in_deployment_safety_wrapper": True,
        "train_data_sha256": R._sha256(data / "train.npz"),
        "validation_data_sha256": R._sha256(data / "validation.npz"),
        "checkpoint": checkpoint_path.name,
    }
    R.save_json(report_path, report)
    print(json.dumps(report, indent=2), flush=True)
    return report


def load_model(seed: int, device: torch.device) -> torch.nn.Module:
    checkpoint = torch.load(
        result_path("best", seed, "pt"), map_location=device, weights_only=False
    )
    if checkpoint["variant"] != METHOD_ID:
        raise RuntimeError(f"Unexpected checkpoint variant: {checkpoint['variant']}")
    model = T.make_model(checkpoint["mean"], checkpoint["std"]).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    return model


@torch.no_grad()
def learned_flux_audit(
    model: torch.nn.Module,
    states: torch.Tensor,
    dt: float,
    batch_size: int = 4,
) -> dict[str, float]:
    multiplier_sum = 0.0
    multiplier_change = 0.0
    multiplier_changed = 0.0
    proposal_change = 0.0
    standard_flux_sum = 0.0
    transverse_change = 0.0
    normal_flux_sum = 0.0
    violations = 0.0
    interfaces = 0
    count = 0
    minimum = float("inf")
    maximum = -float("inf")
    model.eval()
    for start in range(0, states.shape[0], batch_size):
        batch = states[start : start + batch_size]
        for direction in ("x", "y"):
            oriented = C.orient_state(batch, direction)
            multipliers = model.flux_net.multipliers(oriented)
            standard = model.flux_net.standard_roe_faces_oriented(oriented)
            proposal = model.flux_net.forward_oriented(oriented)
            change = (multipliers - 1.0).abs()
            multiplier_sum += float(multipliers.sum())
            multiplier_change += float(change.sum())
            multiplier_changed += float((change > 1.0e-3).sum())
            proposal_change += float((proposal - standard).abs().sum())
            standard_flux_sum += float(standard.abs().sum())
            count += multipliers.numel()
            minimum = min(minimum, float(multipliers.min()))
            maximum = max(maximum, float(multipliers.max()))

        normal = model.normal_raw_fluxes(batch)
        full = model.raw_fluxes(batch, dt=dt)
        transverse_change += sum(
            float((with_transverse - base).abs().sum())
            for with_transverse, base in zip(full, normal)
        )
        normal_flux_sum += sum(float(value.abs().sum()) for value in normal)
        residual_x = C.interface_entropy_residual_oriented(
            full[0], C.orient_state(batch, "x")
        )
        residual_y = C.interface_entropy_residual_oriented(
            C.orient_state(full[1], "y"), C.orient_state(batch, "y")
        )
        violations += float((residual_x > 0.0).sum() + (residual_y > 0.0).sum())
        interfaces += residual_x.numel() + residual_y.numel()
    return {
        "mean_roe_multiplier": multiplier_sum / max(count, 1),
        "minimum_roe_multiplier": minimum,
        "maximum_roe_multiplier": maximum,
        "mean_absolute_multiplier_change_from_one": multiplier_change
        / max(count, 1),
        "multiplier_change_fraction_above_1e-3": multiplier_changed
        / max(count, 1),
        "normal_proposal_change_from_standard_roe_l1_ratio": proposal_change
        / max(standard_flux_sum, 1.0e-12),
        "fixed_transverse_to_normal_flux_l1_ratio": transverse_change
        / max(normal_flux_sum, 1.0e-12),
        "raw_full_proposal_entropy_violation_fraction": violations
        / max(interfaces, 1),
        "transverse_trainable_parameters": 0.0,
    }


@torch.no_grad()
def evaluate(
    seeds: list[int], device: torch.device, data: Path
) -> dict[str, Any]:
    training = torch.from_numpy(R.load_npz(data / "train.npz")["q"]).float()
    scale = C.primitive_scale(training)
    archive = R.load_npz(data / "test.npz")
    reference = torch.from_numpy(archive["q"]).float()
    native = torch.from_numpy(archive["native_roe_coarse"]).float()
    families = archive["family_names"]
    times = archive["times"]
    rows = R.accuracy_rows(
        "PyClaw Roe-64", None, native, reference, scale, families
    )

    hllc, hllc_diagnostics = B.rollout_hllc(
        reference[:, 0].to(device), times
    )
    rows += R.accuracy_rows("HLLC-64", None, hllc, reference, scale, families)
    safety: dict[str, Any] = {"HLLC-64": hllc_diagnostics}
    learned_flux: dict[str, Any] = {}
    timing: dict[str, float] = {}
    audit_states = reference[:, ::5].reshape(-1, *reference.shape[2:]).to(device)
    audit_dt = B.SNAPSHOT_DT / B.TRAINING_SUBSTEPS

    for seed in seeds:
        comparator_checkpoint = (
            DIVERSE / "results" / f"{S.METHOD_ID}_best_seed{seed}.pt"
        )
        if comparator_checkpoint.exists():
            comparator = S.load_model(seed, device)
            comparator_prediction, comparator_diagnostics = B.rollout_hcfl(
                comparator, reference[:, 0].to(device), times
            )
            rows += R.accuracy_rows(
                COMPARATOR_LABEL,
                seed,
                comparator_prediction,
                reference,
                scale,
                families,
            )
            safety[f"{COMPARATOR_LABEL} seed {seed}"] = comparator_diagnostics

        model = load_model(seed, device)
        model.transverse_enabled = False
        ablation_prediction, ablation_diagnostics = B.rollout_hcfl(
            model, reference[:, 0].to(device), times
        )
        rows += R.accuracy_rows(
            INFERENCE_ABLATION_LABEL,
            seed,
            ablation_prediction,
            reference,
            scale,
            families,
        )
        safety[f"{INFERENCE_ABLATION_LABEL} seed {seed}"] = ablation_diagnostics
        model.transverse_enabled = True
        started = time.perf_counter()
        prediction, diagnostics = B.rollout_hcfl(
            model, reference[:, 0].to(device), times
        )
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        key = f"{METHOD_LABEL} seed {seed}"
        timing[key] = time.perf_counter() - started
        rows += R.accuracy_rows(
            METHOD_LABEL, seed, prediction, reference, scale, families
        )
        safety[key] = diagnostics
        learned_flux[key] = learned_flux_audit(model, audit_states, audit_dt)

    summary = {
        "status": "held-out fixed-64 dimension-consistent transverse audit",
        "seeds": seeds,
        "learned_stencil_cells": 6,
        "shared_xy_network": True,
        "fixed_transverse_level": T.TRANSVERSE_LEVEL,
        "fixed_transverse_trainable_parameters": 0,
        "inference_ablation": (
            "the same trained checkpoint is also evaluated with only the fixed "
            "transverse term disabled; this is not a separately trained model"
        ),
        "reference_fine_cells": R._metadata(archive)["fine_cells"],
        "test_data_sha256": R._sha256(data / "test.npz"),
        "accuracy": R.summarize(rows),
        "safety": safety,
        "learned_flux_audit": learned_flux,
        "timing_seconds": timing,
    }
    R.save_rows(RESULTS / "test_case_metrics.csv", rows)
    R.save_json(RESULTS / "test_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    commands = parser.add_subparsers(dest="command", required=True)
    train_parser = commands.add_parser("train")
    train_parser.add_argument("--seed", type=int, required=True)
    train_parser.add_argument("--max-updates", type=int, default=50000)
    train_parser.add_argument("--validation-interval", type=int, default=500)
    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    args = parser.parse_args()

    data = args.data_dir.resolve()
    for name in ("train.npz", "validation.npz", "test.npz"):
        if not (data / name).is_file():
            raise FileNotFoundError(data / name)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(
        {"command": args.command, "device": str(device), "data": str(data)},
        flush=True,
    )
    if args.command == "train":
        train_model(
            args.seed, args.max_updates, args.validation_interval, device, data
        )
    elif args.command == "evaluate":
        evaluate(args.seeds, device, data)


if __name__ == "__main__":
    main()
