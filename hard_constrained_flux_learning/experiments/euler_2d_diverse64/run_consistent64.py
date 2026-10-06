"""Train the method-consistent fixed-64 2-D Euler HCFL model.

This is the direct 2-D analogue of the retained 1-D Euler method:

    central physical flux
    + nonnegative entropy-fixed Roe dissipation multipliers
    + raw-proposal feasibility loss
    + hard interface entropy projection

The x and y evaluations share exactly one network.  All data, optimization,
validation stopping, deployment safety, and test cases match ``run_diverse64``
so that the flux parameterization is the controlled experimental factor.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

import run_diverse64 as R


C = R.C
B = R.B
M = R.M
METHOD_ID = M.CONSISTENT_VARIANT
METHOD_LABEL = "Central + nonnegative Roe-18"
SIGNED_LABEL = "HLLC + signed Roe-18"


def result_path(kind: str, seed: int, extension: str) -> Path:
    return R.RESULTS / f"{METHOD_ID}_{kind}_seed{seed}.{extension}"


def train_model(
    seed: int,
    max_updates: int,
    validation_interval: int,
    device: torch.device,
) -> dict[str, Any]:
    train_archive = R.load_npz(R.DATA / "train.npz")
    validation_archive = R.load_npz(R.DATA / "validation.npz")
    training = torch.from_numpy(train_archive["q"]).float()
    validation = torch.from_numpy(validation_archive["q"]).float()
    family_names = validation_archive["family_names"]
    times = validation_archive["times"]
    mean, std, conserved_std, primitive_std = B.training_statistics(training)

    # Match the signed-model initialization and minibatch schedule exactly.
    B.seed_everything(95100 + 100 * seed)
    model = M.make_model(mean, std, METHOD_ID).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=R.INITIAL_LEARNING_RATE)
    generator = torch.Generator(device=device).manual_seed(95200 + 100 * seed)
    training_device = training.to(device)
    conserved_std_device = conserved_std.to(device)

    R.RESULTS.mkdir(parents=True, exist_ok=True)
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
                    "train_data_sha256": R._sha256(R.DATA / "train.npz"),
                    "validation_data_sha256": R._sha256(
                        R.DATA / "validation.npz"
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
            "flat18 central physical flux + nonnegative entropy-fixed Roe "
            "dissipation + hard entropy projection + proposal feasibility"
        ),
        "controlled_comparator": (
            "flat18 HLLC + signed Roe correction on identical data and schedule"
        ),
        "deployment_grid": [64, 64],
        "cross_grid_claim": False,
        "reference_fine_cells": R._metadata(train_archive)["fine_cells"],
        "train_trajectories": int(training.shape[0]),
        "validation_trajectories": int(validation.shape[0]),
        "families": list(R.EXPECTED_FAMILIES),
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
        "plateau_validations": R.PLATEAU_VALIDATIONS,
        "training_substeps_per_saved_interval": B.TRAINING_SUBSTEPS,
        "proposal_feasibility_weight": B.FEASIBILITY_WEIGHT,
        "roe_multiplier_parameterization": "1 + tanh(logit), range [0, 2]",
        "uses_hllc_inside_learned_proposal": False,
        "uses_low_order_flux_during_training": False,
        "uses_low_order_flux_only_in_deployment_safety_wrapper": True,
        "train_data_sha256": R._sha256(R.DATA / "train.npz"),
        "validation_data_sha256": R._sha256(R.DATA / "validation.npz"),
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
    model = M.make_model(checkpoint["mean"], checkpoint["std"], METHOD_ID).to(
        device
    )
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    return model


@torch.no_grad()
def multiplier_audit(
    model: torch.nn.Module,
    states: torch.Tensor,
    batch_size: int = 4,
) -> dict[str, float]:
    multiplier_sum = 0.0
    change_sum = 0.0
    changed = 0.0
    transverse_sum = 0.0
    proposal_change_sum = 0.0
    baseline_flux_sum = 0.0
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
            multipliers, centre_only, _ = model.flux_net.multiplier_details(
                oriented
            )
            baseline = model.flux_net.standard_roe_faces_oriented(oriented)
            proposal = model.flux_net.forward_oriented(oriented)
            residual = C.interface_entropy_residual_oriented(proposal, oriented)
            change = (multipliers - 1.0).abs()
            multiplier_sum += float(multipliers.sum())
            change_sum += float(change.sum())
            changed += float((change > 1.0e-3).sum())
            transverse_sum += float((multipliers - centre_only).abs().sum())
            proposal_change_sum += float((proposal - baseline).abs().sum())
            baseline_flux_sum += float(baseline.abs().sum())
            violations += float((residual > 0.0).sum())
            interfaces += residual.numel()
            count += multipliers.numel()
            minimum = min(minimum, float(multipliers.min()))
            maximum = max(maximum, float(multipliers.max()))
    return {
        "mean_roe_multiplier": multiplier_sum / max(count, 1),
        "minimum_roe_multiplier": minimum,
        "maximum_roe_multiplier": maximum,
        "mean_absolute_multiplier_change_from_one": change_sum / max(count, 1),
        "multiplier_change_fraction_above_1e-3": changed / max(count, 1),
        "relative_transverse_multiplier_response": transverse_sum
        / max(change_sum, 1.0e-12),
        "raw_proposal_change_from_standard_roe_l1_ratio": proposal_change_sum
        / max(baseline_flux_sum, 1.0e-12),
        "raw_proposal_entropy_violation_fraction": violations
        / max(interfaces, 1),
    }


@torch.no_grad()
def evaluate(seeds: list[int], device: torch.device) -> dict[str, Any]:
    training = torch.from_numpy(R.load_npz(R.DATA / "train.npz")["q"]).float()
    scale = C.primitive_scale(training)
    archive = R.load_npz(R.DATA / "test.npz")
    reference = torch.from_numpy(archive["q"]).float()
    native = torch.from_numpy(archive["native_roe_coarse"]).float()
    families = archive["family_names"]
    times = archive["times"]
    rows = R.accuracy_rows(
        "PyClaw Roe-64", None, native, reference, scale, families
    )

    hllc, hllc_diagnostics = B.rollout_hllc(reference[:, 0].to(device), times)
    rows += R.accuracy_rows(
        "HLLC-64", None, hllc, reference, scale, families
    )
    safety: dict[str, Any] = {"HLLC-64": hllc_diagnostics}
    learned_flux: dict[str, Any] = {}
    timing: dict[str, float] = {}
    audit_states = reference[:, ::5].reshape(-1, *reference.shape[2:]).to(device)

    for seed in seeds:
        signed_model = R.load_model(seed, device)
        started = time.perf_counter()
        signed_prediction, signed_diagnostics = B.rollout_hcfl(
            signed_model, reference[:, 0].to(device), times
        )
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        signed_key = f"{SIGNED_LABEL} seed {seed}"
        timing[signed_key] = time.perf_counter() - started
        rows += R.accuracy_rows(
            SIGNED_LABEL,
            seed,
            signed_prediction,
            reference,
            scale,
            families,
        )
        safety[signed_key] = signed_diagnostics

        model = load_model(seed, device)
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
        learned_flux[key] = multiplier_audit(model, audit_states)

    summary = {
        "status": "held-out method-consistent fixed-64 comparison",
        "seeds": seeds,
        "controlled_factor": (
            "central+nonnegative Roe versus HLLC+signed Roe; identical flat18 "
            "network budget, data, optimizer schedule, projection, and safety"
        ),
        "reference_fine_cells": R._metadata(archive)["fine_cells"],
        "test_data_sha256": R._sha256(R.DATA / "test.npz"),
        "accuracy": R.summarize(rows),
        "safety": safety,
        "learned_flux_audit": learned_flux,
        "timing_seconds": timing,
    }
    R.save_rows(R.RESULTS / "consistent_test_case_metrics.csv", rows)
    R.save_json(R.RESULTS / "consistent_test_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("audit-data")
    commands.add_parser("reference-audit")
    train_parser = commands.add_parser("train")
    train_parser.add_argument("--seed", type=int, required=True)
    train_parser.add_argument("--max-updates", type=int, default=50000)
    train_parser.add_argument("--validation-interval", type=int, default=500)
    evaluate_parser = commands.add_parser("evaluate")
    evaluate_parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print({"command": args.command, "device": str(device)}, flush=True)
    if args.command == "audit-data":
        R.audit_data()
    elif args.command == "reference-audit":
        R.reference_audit()
    elif args.command == "train":
        train_model(args.seed, args.max_updates, args.validation_interval, device)
    elif args.command == "evaluate":
        evaluate(args.seeds, device)


if __name__ == "__main__":
    main()
