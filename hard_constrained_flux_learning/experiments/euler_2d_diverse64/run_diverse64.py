"""Train and evaluate the fixed-grid diverse-data 2-D Euler HCFL experiment."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
BASE = HERE.parent / "euler_2d_clawpack_quadrants"
TRANSVERSE = HERE.parent / "euler_2d_transverse_ablation"
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(TRANSVERSE))
import euler_2d_common as C  # noqa: E402
import run_experiment as B  # noqa: E402
import transverse_models as M  # noqa: E402


BATCH_SIZE = 4
INITIAL_LEARNING_RATE = 5.0e-4
MINIMUM_LEARNING_RATE = 2.0e-5
PLATEAU_VALIDATIONS = 6
RELATIVE_IMPROVEMENT = 0.001
EXPECTED_FAMILIES = (
    "oblique_riemann",
    "contact_shear",
    "oblique_quadrant",
    "radial_interface",
    "colliding_waves",
    "smooth_packet",
)


def load_npz(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def _metadata(archive: dict[str, Any]) -> dict[str, Any]:
    return json.loads(str(archive["metadata_json"].item()))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def save_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Cannot save an empty table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def result_path(kind: str, seed: int, extension: str) -> Path:
    return RESULTS / f"flat18_{kind}_seed{seed}.{extension}"


def audit_data() -> dict[str, Any]:
    names = ("train.npz", "validation.npz", "test.npz")
    archives = {name: load_npz(DATA / name) for name in names}
    files: dict[str, Any] = {}
    all_specs: dict[str, set[str]] = {}
    passed = True
    expected_times: np.ndarray | None = None

    for name, archive in archives.items():
        q = torch.from_numpy(archive["q"]).float()
        times = archive["times"]
        metadata = _metadata(archive)
        families = [str(value) for value in archive["family_names"].tolist()]
        family_counts = {family: families.count(family) for family in EXPECTED_FAMILIES}
        specs = json.loads(str(archive["specs_json"].item()))
        spec_keys = {json.dumps(spec, sort_keys=True) for spec in specs}
        all_specs[name] = spec_keys
        if expected_times is None:
            expected_times = times
        times_match = np.array_equal(times, expected_times)
        pressure = C.pressure_raw(q)
        item = {
            "shape": list(q.shape),
            "sha256": _sha256(DATA / name),
            "metadata": metadata,
            "family_counts": family_counts,
            "finite": bool(torch.isfinite(q).all()),
            "minimum_density": float(q[..., 0].min()),
            "minimum_pressure": float(pressure.min()),
            "times_match": times_match,
            "has_native_roe_coarse": "native_roe_coarse" in archive,
        }
        item_passed = (
            item["finite"]
            and item["minimum_density"] > 0.0
            and item["minimum_pressure"] > 0.0
            and times_match
            and q.shape[-3:-1] == (64, 64)
            and set(families) == set(EXPECTED_FAMILIES)
            and max(family_counts.values()) - min(family_counts.values()) <= 1
            and metadata["coarse_cells"] == 64
            and metadata["restriction"]
            == "conservative_block_average_of_conserved_variables"
            and (name != "test.npz" or item["has_native_roe_coarse"])
        )
        item["passed"] = item_passed
        passed = passed and item_passed
        files[name] = item

    disjoint = {
        "train_validation": not bool(all_specs["train.npz"] & all_specs["validation.npz"]),
        "train_test": not bool(all_specs["train.npz"] & all_specs["test.npz"]),
        "validation_test": not bool(
            all_specs["validation.npz"] & all_specs["test.npz"]
        ),
    }
    passed = passed and all(disjoint.values())
    report = {"passed": passed, "files": files, "splits_are_disjoint": disjoint}
    save_json(RESULTS / "data_audit.json", report)
    print(json.dumps(report, indent=2), flush=True)
    if not passed:
        raise RuntimeError("Diverse64 data audit failed")
    return report


def reference_audit() -> dict[str, Any]:
    archives = {
        cells: load_npz(DATA / f"reference_audit{cells}.npz")
        for cells in (256, 512, 1024)
    }
    reference = torch.from_numpy(archives[1024]["q"]).float()
    primitive_reference = C.primitive(reference)
    scale = primitive_reference.std(dim=(0, 1, 2, 3)).clamp_min(1.0e-8)
    families = [str(value) for value in archives[1024]["family_names"].tolist()]
    report: dict[str, Any] = {
        "reference": 1024,
        "scale": scale.tolist(),
        "common_initial_state_maximum_absolute_difference": {},
        "comparisons": {},
    }
    for cells in (256, 512):
        candidate = torch.from_numpy(archives[cells]["q"]).float()
        report["common_initial_state_maximum_absolute_difference"][str(cells)] = float(
            (candidate[:, 0] - reference[:, 0]).abs().max()
        )
        rows: dict[str, Any] = {}
        for family in EXPECTED_FAMILIES:
            indices = [index for index, value in enumerate(families) if value == family]
            rows[family] = {
                "rollout_nmae": C.primitive_nmae(
                    candidate[indices, 1:], reference[indices, 1:], scale
                ),
                "final_nmae": C.primitive_nmae(
                    candidate[indices, -1:], reference[indices, -1:], scale
                ),
            }
        rows["all"] = {
            "rollout_nmae": C.primitive_nmae(candidate[:, 1:], reference[:, 1:], scale),
            "final_nmae": C.primitive_nmae(candidate[:, -1:], reference[:, -1:], scale),
            "rollout_nrmse": C.primitive_nrmse(
                candidate[:, 1:], reference[:, 1:], scale
            ),
        }
        report["comparisons"][f"{cells}_vs_1024"] = rows
    save_json(RESULTS / "reference_audit.json", report)
    print(json.dumps(report, indent=2), flush=True)
    return report


@torch.no_grad()
def validate(
    model: torch.nn.Module,
    validation: torch.Tensor,
    family_names: np.ndarray,
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
    result: dict[str, float] = {
        "completion": float(complete.float().mean()),
        "minimum_density": float(prediction[..., 0].min()),
        "minimum_pressure": float(C.pressure_raw(prediction).min()),
        "maximum_entropy_balance": diagnostics["maximum_entropy_balance"],
        "maximum_interface_residual": diagnostics["maximum_interface_residual"],
        "minimum_positivity_theta": diagnostics["minimum_positivity_theta"],
        "minimum_entropy_beta": diagnostics["minimum_entropy_beta"],
    }
    if bool(complete.any()):
        result["rollout_nmae_completed"] = C.primitive_nmae(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
        result["rollout_nrmse_completed"] = C.primitive_nrmse(
            prediction[complete, 1:], validation[complete, 1:], primitive_std
        )
    else:
        result["rollout_nmae_completed"] = float("inf")
        result["rollout_nrmse_completed"] = float("inf")

    family_values: list[float] = []
    for family in EXPECTED_FAMILIES:
        mask = torch.tensor(family_names == family) & complete
        value = (
            C.primitive_nmae(
                prediction[mask, 1:], validation[mask, 1:], primitive_std
            )
            if bool(mask.any())
            else float("inf")
        )
        result[f"nmae_{family}"] = value
        family_values.append(value)
    result["worst_family_nmae"] = max(family_values)
    model.train()
    return result


def train_model(
    seed: int,
    max_updates: int,
    validation_interval: int,
    device: torch.device,
) -> dict[str, Any]:
    train_archive = load_npz(DATA / "train.npz")
    validation_archive = load_npz(DATA / "validation.npz")
    training = torch.from_numpy(train_archive["q"]).float()
    validation = torch.from_numpy(validation_archive["q"]).float()
    family_names = validation_archive["family_names"]
    times = validation_archive["times"]
    mean, std, conserved_std, primitive_std = B.training_statistics(training)

    B.seed_everything(95100 + 100 * seed)
    model = M.make_model(mean, std, "flat18").to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=INITIAL_LEARNING_RATE)
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
        metrics = validate(
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
                    "variant": "flat18",
                    "seed": seed,
                    "update": update,
                    "validation": metrics,
                    "train_data_sha256": _sha256(DATA / "train.npz"),
                    "validation_data_sha256": _sha256(DATA / "validation.npz"),
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
        save_rows(curve_path, curve)
        print(row, flush=True)

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
        if current_best <= plateau_anchor * (1.0 - RELATIVE_IMPROVEMENT):
            plateau_anchor = current_best
            stale_validations = 0
            continue
        stale_validations += 1
        if stale_validations < PLATEAU_VALIDATIONS:
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
                "seed": seed,
                "update": update,
                "new_learning_rate": new_learning_rate,
            },
            flush=True,
        )

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    report = {
        "method": "flat18 HLLC + signed Roe correction + hard entropy projection + proposal feasibility",
        "deployment_grid": [64, 64],
        "cross_grid_claim": False,
        "reference_fine_cells": _metadata(train_archive)["fine_cells"],
        "train_trajectories": int(training.shape[0]),
        "validation_trajectories": int(validation.shape[0]),
        "families": list(EXPECTED_FAMILIES),
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
        "plateau_validations": PLATEAU_VALIDATIONS,
        "training_substeps_per_saved_interval": B.TRAINING_SUBSTEPS,
        "proposal_feasibility_weight": B.FEASIBILITY_WEIGHT,
        "uses_low_order_flux_during_training": False,
        "uses_low_order_flux_only_in_deployment_safety_wrapper": True,
        "train_data_sha256": _sha256(DATA / "train.npz"),
        "validation_data_sha256": _sha256(DATA / "validation.npz"),
        "checkpoint": checkpoint_path.name,
    }
    save_json(report_path, report)
    print(json.dumps(report, indent=2), flush=True)
    return report


def load_model(seed: int, device: torch.device) -> torch.nn.Module:
    checkpoint = torch.load(
        result_path("best", seed, "pt"), map_location=device, weights_only=False
    )
    model = M.make_model(checkpoint["mean"], checkpoint["std"], "flat18").to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    return model


def accuracy_rows(
    method: str,
    seed: int | None,
    prediction: torch.Tensor,
    reference: torch.Tensor,
    scale: torch.Tensor,
    families: np.ndarray,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(prediction.shape[0]):
        audit = C.audit_accuracy(prediction[index], reference[index], scale)
        pred_tv = C.total_variation_density(prediction[index, -1:]).double()[0]
        ref_tv = C.total_variation_density(reference[index, -1:]).double()[0]
        rows.append(
            {
                "method": method,
                "seed": seed,
                "case": index,
                "family": str(families[index]),
                **asdict(audit),
                "density_tv_relative_error_signed": float(
                    (pred_tv - ref_tv) / ref_tv.clamp_min(1.0e-12)
                ),
            }
        )
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: defaultdict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["method"]), str(row["family"]))].append(row)
        groups[(str(row["method"]), "all")].append(row)
    result: dict[str, Any] = {}
    metrics = (
        "nmae",
        "nrmse",
        "density_mae",
        "x_velocity_mae",
        "y_velocity_mae",
        "pressure_mae",
        "density_tv_excess",
        "density_tv_relative_error_signed",
        "density_range_overshoot",
        "pressure_range_overshoot",
    )
    for (method, family), group in groups.items():
        completed = [row for row in group if row["completed"]]
        item: dict[str, Any] = {
            "rows": len(group),
            "completion_rate": len(completed) / len(group),
            "minimum_density": min(float(row["minimum_density"]) for row in group),
            "minimum_pressure": min(float(row["minimum_pressure"]) for row in group),
        }
        for metric in metrics:
            values = [float(row[metric]) for row in completed if row[metric] is not None]
            item[f"{metric}_mean_completed"] = float(np.mean(values)) if values else None
            item[f"{metric}_std_completed"] = float(np.std(values)) if values else None
        result[f"{method}/{family}"] = item
    return result


@torch.no_grad()
def learned_flux_audit(
    model: torch.nn.Module,
    states: torch.Tensor,
    batch_size: int = 4,
) -> dict[str, float]:
    """Verify that the selected network is not the zero-correction solution."""

    coefficient_absolute_sum = 0.0
    coefficient_nonzero = 0.0
    transverse_response_sum = 0.0
    raw_correction_sum = 0.0
    base_flux_sum = 0.0
    count = 0
    model.eval()
    for start in range(0, states.shape[0], batch_size):
        batch = states[start : start + batch_size]
        for direction in ("x", "y"):
            oriented = C.orient_state(batch, direction)
            coefficients, centre_only, _ = model.flux_net.coefficient_details(oriented)
            base = C.hllc_faces_oriented(oriented)
            proposal = model.flux_net.forward_oriented(oriented)
            coefficient_absolute_sum += float(coefficients.abs().sum())
            coefficient_nonzero += float((coefficients.abs() > 1.0e-3).sum())
            transverse_response_sum += float((coefficients - centre_only).abs().sum())
            raw_correction_sum += float((proposal - base).abs().sum())
            base_flux_sum += float(base.abs().sum())
            count += coefficients.numel()
    return {
        "mean_absolute_roe_coefficient": coefficient_absolute_sum / max(count, 1),
        "coefficient_fraction_above_1e-3": coefficient_nonzero / max(count, 1),
        "relative_transverse_coefficient_response": transverse_response_sum
        / max(coefficient_absolute_sum, 1.0e-12),
        "raw_correction_to_hllc_flux_l1_ratio": raw_correction_sum
        / max(base_flux_sum, 1.0e-12),
    }


@torch.no_grad()
def evaluate(seeds: list[int], device: torch.device) -> dict[str, Any]:
    training = torch.from_numpy(load_npz(DATA / "train.npz")["q"]).float()
    scale = C.primitive_scale(training)
    archive = load_npz(DATA / "test.npz")
    reference = torch.from_numpy(archive["q"]).float()
    native = torch.from_numpy(archive["native_roe_coarse"]).float()
    families = archive["family_names"]
    times = archive["times"]
    rows = accuracy_rows("PyClaw Roe-64", None, native, reference, scale, families)

    hllc, hllc_diagnostics = B.rollout_hllc(
        reference[:, 0].to(device), times
    )
    rows += accuracy_rows("HLLC-64", None, hllc, reference, scale, families)
    safety: dict[str, Any] = {"HLLC-64": hllc_diagnostics}
    learned_flux: dict[str, Any] = {}
    timing: dict[str, float] = {}
    for seed in seeds:
        model = load_model(seed, device)
        started = time.perf_counter()
        prediction, diagnostics = B.rollout_hcfl(
            model, reference[:, 0].to(device), times
        )
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        timing[f"flat18 seed {seed}"] = time.perf_counter() - started
        rows += accuracy_rows("flat18", seed, prediction, reference, scale, families)
        safety[f"flat18 seed {seed}"] = diagnostics
        audit_states = reference[:, ::5].reshape(-1, *reference.shape[2:]).to(device)
        learned_flux[f"flat18 seed {seed}"] = learned_flux_audit(
            model, audit_states
        )

    summary = {
        "status": "held-out diverse fixed-64 evaluation",
        "seeds": seeds,
        "reference_fine_cells": _metadata(archive)["fine_cells"],
        "test_data_sha256": _sha256(DATA / "test.npz"),
        "accuracy": summarize(rows),
        "safety": safety,
        "learned_flux_nonzero_audit": learned_flux,
        "timing_seconds": timing,
    }
    save_rows(RESULTS / "test_case_metrics.csv", rows)
    save_json(RESULTS / "test_summary.json", summary)
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
        audit_data()
    elif args.command == "reference-audit":
        reference_audit()
    elif args.command == "train":
        train_model(args.seed, args.max_updates, args.validation_interval, device)
    elif args.command == "evaluate":
        evaluate(args.seeds, device)


if __name__ == "__main__":
    main()
