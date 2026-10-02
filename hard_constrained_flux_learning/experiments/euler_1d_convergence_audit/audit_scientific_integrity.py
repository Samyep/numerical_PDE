"""Reproduce the seed-0 data-leakage and matched-control integrity checks.

This audit deliberately does not train or select a checkpoint.  It verifies
the already frozen convergence artifacts, regenerates the train/validation
initial states to check exact split separation, and inspects the learned flux
against its zero-correction HLLC control.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

import plot_best_vs_fvm as comparison
import run_convergence_audit as convergence


HERE = Path(__file__).resolve().parent
base = convergence.base
shared = convergence.shared
SELECTED_ARM = "dissipation_broad"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_artifact(path: Path) -> str:
    """Hash text with canonical newlines and binary checkpoints bytewise."""
    if path.suffix.lower() in {".py", ".csv", ".json", ".md"}:
        content = path.read_text(encoding="utf-8")
        canonical = content.replace("\r\n", "\n").replace("\r", "\n")
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return sha256_file(path)


def state_fingerprints(states: torch.Tensor) -> set[str]:
    initial = states[:, 0].detach().cpu().numpy()
    return {
        hashlib.sha256(np.ascontiguousarray(row).tobytes()).hexdigest()
        for row in initial
    }


def canonical_fingerprints() -> dict[str, str]:
    output: dict[str, str] = {}
    for name, state in shared.canonical_initial_conditions().items():
        restricted = state.reshape(
            1,
            base.NCOARSE,
            base.FACTOR,
            3,
        ).mean(axis=2).astype(np.float32)
        output[name] = hashlib.sha256(
            np.ascontiguousarray(restricted[0]).tobytes()
        ).hexdigest()
    return output


def audit_dataset_separation(seed: int) -> dict[str, Any]:
    broad = shared.make_baseline_training_data(seed)
    wave = shared.make_wave_coverage_training_data(broad, seed)
    validation = convergence.make_validation_data(seed)

    broad_hashes = state_fingerprints(broad)
    wave_hashes = state_fingerprints(wave)
    training_hashes = broad_hashes | wave_hashes
    validation_hashes = state_fingerprints(validation)
    canonical = canonical_fingerprints()
    canonical_hashes = set(canonical.values())

    train_validation_overlap = training_hashes & validation_hashes
    train_canonical_overlap = training_hashes & canonical_hashes
    validation_canonical_overlap = validation_hashes & canonical_hashes
    if train_validation_overlap:
        raise AssertionError("Exact train/validation initial-state overlap")
    if train_canonical_overlap:
        raise AssertionError("A canonical initial state occurs in training")
    if validation_canonical_overlap:
        raise AssertionError("A canonical initial state occurs in validation")

    return {
        "training_tensor_shape": list(broad.shape),
        "wave_training_tensor_shape": list(wave.shape),
        "validation_tensor_shape": list(validation.shape),
        "training_cells_seen_by_network": int(broad.shape[-2]),
        "fine_cells_used_to_generate_training_targets": base.NREF,
        "restriction_factor": base.FACTOR,
        "training_generator_seeds": {
            "ordinary": 6000 + seed,
            "broad_random": 7000 + seed,
            "random_extreme": 8000 + seed,
            "structured_wave": 9000 + seed,
        },
        "validation_generator_seeds": {
            "ordinary": 14000 + seed,
            "broad_random": 15000 + seed,
            "random_extreme": 16000 + seed,
            "structured_wave": 17000 + seed,
        },
        "exact_initial_state_overlap_counts": {
            "training_validation": len(train_validation_overlap),
            "training_canonical": len(train_canonical_overlap),
            "validation_canonical": len(validation_canonical_overlap),
        },
        "canonical_initial_state_sha256": canonical,
    }


def audit_checkpoint_selection(results_dir: Path, seed: int) -> dict[str, Any]:
    curve = read_csv(results_dir / f"training_curve_seed{seed}.csv")
    convergence_rows = read_csv(
        results_dir / f"convergence_seed{seed}.csv"
    )
    summaries = {row["arm"]: row for row in convergence_rows}
    expected_arms = set(convergence.ARMS)
    if set(summaries) != expected_arms:
        raise AssertionError("Convergence summary does not contain every arm")

    arms: dict[str, Any] = {}
    for arm in convergence.ARMS:
        candidates = [row for row in curve if row["arm"] == arm]
        if not candidates:
            raise AssertionError(f"Missing validation curve for {arm}")
        selected = min(
            candidates,
            key=lambda row: float(row["validation_rollout_nrmse"]),
        )
        summary = summaries[arm]
        selected_update = int(selected["update"])
        selected_metric = float(selected["validation_rollout_nrmse"])
        recorded_update = int(summary["best_update"])
        recorded_metric = float(summary["best_validation_rollout_nrmse"])
        if selected_update != recorded_update or not np.isclose(
            selected_metric,
            recorded_metric,
            rtol=0.0,
            atol=1.0e-12,
        ):
            raise AssertionError(
                f"{arm} checkpoint is not the validation-curve minimum"
            )
        if summary["converged"].lower() != "true":
            raise AssertionError(f"{arm} did not meet the convergence rule")
        if summary["stop_reason"] != (
            "validation_plateau_at_minimum_learning_rate"
        ):
            raise AssertionError(f"Unexpected stop reason for {arm}")
        checkpoint = results_dir / (
            f"{arm}_converged_best_seed{seed}.pt"
        )
        if not checkpoint.is_file():
            raise AssertionError(f"Missing converged checkpoint: {checkpoint}")
        arms[arm] = {
            "selected_update": selected_update,
            "validation_rollout_nrmse": selected_metric,
            "stop_update": int(summary["stop_update"]),
            "stop_reason": summary["stop_reason"],
            "checkpoint_sha256": sha256_file(checkpoint),
        }

    winner = min(
        arms,
        key=lambda arm: arms[arm]["validation_rollout_nrmse"],
    )
    if winner != SELECTED_ARM:
        raise AssertionError(
            f"Expected validation winner {SELECTED_ARM}, found {winner}"
        )

    retained_checkpoints = sorted(HERE.parents[1].glob("**/*.pt"))
    nonconverged_checkpoints = [
        path
        for path in retained_checkpoints
        if "_converged_best_" not in path.name
    ]
    if nonconverged_checkpoints:
        raise AssertionError(
            "Weights without a validation-converged label are retained: "
            + ", ".join(path.name for path in nonconverged_checkpoints)
        )
    return {
        "selection_signal": "independent validation rollout NRMSE only",
        "canonical_cases_used_for_selection": False,
        "selected_arm": winner,
        "all_arms_met_declared_convergence_rule": True,
        "retained_nonconverged_checkpoints": [],
        "retained_converged_checkpoint_count": len(retained_checkpoints),
        "arms": arms,
    }


@torch.no_grad()
def audit_nonzero_correction(
    results_dir: Path,
    seed: int,
    width: int,
) -> dict[str, float]:
    model = comparison.load_model(
        results_dir,
        SELECTED_ARM,
        seed,
        width,
    )
    correction_squared = 0.0
    hllc_squared = 0.0
    correction_divergence_squared = 0.0
    hllc_divergence_squared = 0.0
    maximum_correction = 0.0

    for name in comparison.CASES:
        state = torch.from_numpy(
            comparison.initial_condition(name, 512)
        ).float()
        hllc = base.t_hllc(state)
        learned = model.flux_net(state)
        correction = learned - hllc
        hllc_divergence = hllc - torch.roll(hllc, 1, dims=-2)
        correction_divergence = correction - torch.roll(
            correction,
            1,
            dims=-2,
        )
        correction_squared += float((correction.double() ** 2).sum())
        hllc_squared += float((hllc.double() ** 2).sum())
        correction_divergence_squared += float(
            (correction_divergence.double() ** 2).sum()
        )
        hllc_divergence_squared += float(
            (hllc_divergence.double() ** 2).sum()
        )
        maximum_correction = max(
            maximum_correction,
            float(correction.abs().max()),
        )

    constant = torch.tensor(
        [[[1.0, 0.2, 2.52]]],
        dtype=torch.float32,
    ).repeat(1, 512, 1)
    constant_correction = model.flux_net(constant) - base.t_hllc(constant)
    constant_error = float(constant_correction.abs().max())
    final_layer = model.flux_net.net[-1]
    parameter_l2 = float(
        torch.sqrt(
            (final_layer.weight.double() ** 2).sum()
            + (final_layer.bias.double() ** 2).sum()
        )
    )
    flux_ratio = float(
        np.sqrt(correction_squared / max(hllc_squared, 1.0e-30))
    )
    update_ratio = float(
        np.sqrt(
            correction_divergence_squared
            / max(hllc_divergence_squared, 1.0e-30)
        )
    )
    if parameter_l2 <= 0.0 or maximum_correction <= 0.0 or update_ratio <= 0.0:
        raise AssertionError("The selected neural correction is identically zero")
    if constant_error > 1.0e-7:
        raise AssertionError("Dissipation correction is inconsistent on constants")

    return {
        "final_layer_parameter_l2": parameter_l2,
        "maximum_initial_flux_correction": maximum_correction,
        "initial_flux_correction_relative_rms": flux_ratio,
        "initial_update_correction_relative_rms": update_ratio,
        "maximum_constant_state_correction": constant_error,
    }


def audit_matched_hllc(results_dir: Path, seed: int) -> dict[str, Any]:
    path = results_dir / (
        f"hllc512_vs_hcfl512_with_hllc2048_seed{seed}.json"
    )
    result = json.loads(path.read_text(encoding="utf-8"))
    if result["checkpoint_training_cells"] != base.NCOARSE:
        raise AssertionError("Comparison metadata misstates the training grid")
    if result["deployment_cells"] != 512:
        raise AssertionError("Comparison was not run on 512 cells")
    if "2048" not in result["reference"] or "512" not in result["reference"]:
        raise AssertionError("Comparison reference is not HLLC-2048 to 512")

    per_case: dict[str, Any] = {}
    maximum_conservation_drift = 0.0
    for name, case in result["cases"].items():
        matched = case["matched_safe_hllc_512_zero_nn"]
        learned = case["learned_hcfl_512"]
        matched_error = float(matched["rollout_nrmse"])
        learned_error = float(learned["rollout_nrmse"])
        if learned_error >= matched_error:
            raise AssertionError(
                f"Learned HCFL did not beat the matched control on {name}"
            )
        if learned["local_limiter_intervention_rate"] != 0.0:
            raise AssertionError(f"Positivity fallback was active on {name}")
        if learned["fd_entropy_intervention_rate"] != 0.0:
            raise AssertionError(f"Entropy fallback was active on {name}")
        case_drift = max(
            float(value)
            for value in learned["max_conservation_drift"].values()
        )
        maximum_conservation_drift = max(
            maximum_conservation_drift,
            case_drift,
        )
        per_case[name] = {
            "matched_zero_nn_rollout_nrmse": matched_error,
            "learned_rollout_nrmse": learned_error,
            "learned_reduction_percent": 100.0
            * (matched_error - learned_error)
            / matched_error,
        }
    if maximum_conservation_drift >= 1.0e-6:
        raise AssertionError("Learned rollout conservation drift exceeds 1e-6")

    mean = result["mean_metrics"]
    return {
        "reference": result["reference"],
        "training_cells": result["checkpoint_training_cells"],
        "deployment_cells": result["deployment_cells"],
        "matched_zero_nn_mean_rollout_nrmse": mean[
            "matched_safe_hllc_512_zero_nn"
        ]["mean_rollout_nrmse"],
        "learned_mean_rollout_nrmse": mean["learned_hcfl_512"][
            "mean_rollout_nrmse"
        ],
        "learned_reduction_vs_matched_zero_nn_percent": mean[
            "learned_hcfl_512"
        ]["reduction_vs_matched_zero_nn_percent"],
        "learned_better_on_every_canonical_case": True,
        "positivity_or_fully_discrete_entropy_fallback_used": False,
        "maximum_conservation_drift": maximum_conservation_drift,
        "per_case": per_case,
    }


def artifact_hashes(results_dir: Path, seed: int) -> dict[str, str]:
    paths = (
        HERE / "run_convergence_audit.py",
        HERE / "plot_hllc512_vs_hcfl512.py",
        results_dir / f"training_curve_seed{seed}.csv",
        results_dir / f"convergence_seed{seed}.csv",
        results_dir / f"dissipation_broad_converged_best_seed{seed}.pt",
        results_dir
        / f"hllc512_vs_hcfl512_with_hllc2048_seed{seed}.json",
        results_dir
        / f"hllc512_vs_hcfl512_with_hllc2048_seed{seed}.csv",
    )
    return {
        path.relative_to(HERE).as_posix(): sha256_artifact(path)
        for path in paths
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=HERE / "results",
    )
    args = parser.parse_args()
    results_dir = args.results_dir.resolve()

    dataset = audit_dataset_separation(args.seed)
    selection = audit_checkpoint_selection(results_dir, args.seed)
    nonzero = audit_nonzero_correction(
        results_dir,
        args.seed,
        args.width,
    )
    matched = audit_matched_hllc(results_dir, args.seed)
    report = {
        "seed": args.seed,
        "verdict": (
            "PASS for the recorded seed-0 controlled experiment; this is "
            "evidence of utility without detected leakage, not a universal "
            "or multi-seed superiority claim"
        ),
        "checks": {
            "exact_split_overlap": "pass",
            "validation_only_checkpoint_selection": "pass",
            "all_retained_training_checkpoints_converged": "pass",
            "learned_correction_nonzero": "pass",
            "constant_state_consistency": "pass",
            "matched_zero_nn_control": "pass",
            "admissibility_and_conservation": "pass",
        },
        "dataset_separation": dataset,
        "checkpoint_selection": selection,
        "learned_correction": nonzero,
        "matched_hllc_512_control": matched,
        "artifact_sha256": artifact_hashes(results_dir, args.seed),
        "limitations": [
            "one training seed",
            "five named canonical problems",
            "HLLC-2048 is a finite-resolution reference rather than an exact solution",
            "the 512-cell result is zero-shot resolution transfer from 64-cell states",
            "MUSCL-HLLC and WENO baselines have not yet been tested",
        ],
    }
    output = results_dir / f"scientific_integrity_audit_seed{args.seed}.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(output)
    print(json.dumps(report["checks"], indent=2))


if __name__ == "__main__":
    main()
