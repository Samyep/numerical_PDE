"""Controlled extensions of the frozen a3 Terminal-Aware Tangent-Actor."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    make_targets,
    tangent_labels,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.models import (
    ConstantTangentActor,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    train_actor,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    evaluate_tangent_actor,
    load_tangent_actor,
)

from .common import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    RAW_DIR,
    STRONG_RESULTS_DIR,
    atomic_json,
    atomic_torch_save,
    driver_for,
    ensure_dirs,
    evaluation_data,
    parameter_count,
    stable_seed,
    state_fingerprint,
    synchronize,
)


BRANCHES: dict[str, dict[str, int]] = {
    "iterations": {"training_points": 2048, "paths_per_label": 128, "phases": 3},
    "points4096": {"training_points": 4096, "paths_per_label": 128, "phases": 1},
    "points8192": {"training_points": 8192, "paths_per_label": 128, "phases": 1},
    "paths256": {"training_points": 2048, "paths_per_label": 256, "phases": 1},
    "paths512": {"training_points": 2048, "paths_per_label": 512, "phases": 1},
}
CONTINUATION_BRANCHES = ("points8192_continued",)


def _strong_a3_path(benchmark: str, seed: int) -> Path:
    return (
        STRONG_RESULTS_DIR / "checkpoints"
        / f"tangent_actor_{benchmark}_s{seed}_a3.pt"
    )


def _base_training_seconds(benchmark: str, seed: int) -> float:
    path = STRONG_RESULTS_DIR / "raw" / f"tangent_actor_{benchmark}_s{seed}_main.json"
    with path.open() as handle:
        result = json.load(handle)
    return float(result["phases"][2]["cumulative_training_seconds"])


def _state_bank(
    benchmark: str,
    seed: int,
    extension_index: int,
    *,
    maximum_points: int = 8192,
    validation_points: int = 256,
) -> dict[str, np.ndarray]:
    """Create exact nested point prefixes shared by all controlled branches."""
    ensure_dirs()
    path = CACHE_DIR / f"nested_points_{benchmark}_s{seed}_e{extension_index}.npz"
    if path.exists():
        with np.load(path) as archive:
            return {name: archive[name] for name in archive.files}
    driver = driver_for(benchmark)
    train_rng = np.random.default_rng(
        stable_seed("high_budget_nested_train_points", benchmark, seed, extension_index)
    )
    valid_rng = np.random.default_rng(
        stable_seed("high_budget_fixed_validation_points", benchmark, seed)
    )
    payload = {
        "train_t": train_rng.uniform(
            0.0, np.nextafter(driver.T, 0.0), maximum_points
        ).astype(np.float64),
        "train_x": train_rng.standard_normal((maximum_points, driver.d)).astype(np.float64),
        "valid_t": valid_rng.uniform(
            0.0, np.nextafter(driver.T, 0.0), validation_points
        ).astype(np.float64),
        "valid_x": valid_rng.standard_normal((validation_points, driver.d)).astype(np.float64),
    }
    np.savez_compressed(path, **payload)
    return payload


def nested_bank_audit(benchmark: str, seed: int) -> dict[str, Any]:
    bank = _state_bank(benchmark, seed, 1)
    fingerprints = {
        str(count): state_fingerprint(bank["train_t"][:count], bank["train_x"][:count])
        for count in (2048, 4096, 8192)
    }
    return {
        "benchmark": benchmark, "seed": seed,
        "counts": [2048, 4096, 8192], "prefix_fingerprints": fingerprints,
        "prefix_equality_verified": all(
            np.array_equal(bank["train_x"][:smaller], bank["train_x"][:larger][:smaller])
            and np.array_equal(bank["train_t"][:smaller], bank["train_t"][:larger][:smaller])
            for smaller, larger in ((2048, 4096), (4096, 8192))
        ),
        "validation_fingerprint": state_fingerprint(bank["valid_t"], bank["valid_x"]),
    }


def _generate_labels(
    benchmark: str,
    seed: int,
    extension_index: int,
    source: nn.Module,
    config: FitConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    driver = driver_for(benchmark)
    bank = _state_bank(benchmark, seed, extension_index)
    train_t = bank["train_t"][: config.training_points]
    train_x = bank["train_x"][: config.training_points]
    valid_t, valid_x = bank["valid_t"], bank["valid_x"]
    common = {
        "paths": config.paths_per_label,
        "time_steps": config.label_time_steps,
        "device": device,
        "parent_batch": config.parent_batch,
        "path_chunk": config.path_chunk,
    }
    # Paths use common seeds across M, so the antithetic stream is shared as far
    # as the batched traversal permits. Point banks are exact nested prefixes.
    train_value, train_q, train_meta = tangent_labels(
        driver, source, train_t, train_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "train"),
        **common,
    )
    valid_value, valid_q, valid_meta = tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "valid"),
        **common,
    )
    replicate_value, replicate_q, replicate_meta = tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "replicate"),
        **common,
    )
    train_r, train_q, train_projection = make_targets(driver, train_value, train_q)
    valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_q)
    value_replication = float(np.sqrt(np.mean((valid_value - replicate_value) ** 2)))
    if valid_q is None or replicate_q is None:
        gradient_component = gradient_vector = None
    else:
        difference = valid_q - replicate_q
        gradient_component = float(np.sqrt(np.mean(difference**2)))
        gradient_vector = float(np.sqrt(np.mean(np.sum(difference**2, axis=1))))
    arrays: dict[str, np.ndarray] = {
        "train_t": train_t, "train_x": train_x, "train_r": train_r,
        "valid_t": valid_t, "valid_x": valid_x, "valid_r": valid_r,
    }
    if train_q is not None and valid_q is not None:
        arrays.update(train_q=train_q, valid_q=valid_q)
    metadata = {
        "train": train_meta, "validation": valid_meta, "replication": replicate_meta,
        "train_projection": train_projection, "validation_projection": valid_projection,
        "value_label_replication_rmse": value_replication,
        "gradient_label_replication_component_rmse": gradient_component,
        "gradient_label_replication_vector_rmse": gradient_vector,
        "wall_seconds": float(
            train_meta["wall_seconds"] + valid_meta["wall_seconds"]
            + replicate_meta["wall_seconds"]
        ),
        "nested_state_bank_fingerprint": state_fingerprint(train_t, train_x),
        "fixed_validation_bank_fingerprint": state_fingerprint(valid_t, valid_x),
        "ground_truth_used": False,
    }
    return arrays, metadata


def _evaluate_replicates(
    benchmark: str,
    seed: int,
    branch: str,
    extension_index: int,
    actor: nn.Module,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, list[float], list[dict[str, Any]]]:
    predictions: list[np.ndarray] = []
    metadata: list[dict[str, Any]] = []
    for replicate in range(4):
        prediction, _, meta = evaluate_tangent_actor(
            benchmark, actor, x, device, paths=128, time_steps=32,
            seed=stable_seed(
                "high_budget_policy_validation", benchmark, seed, branch,
                extension_index, replicate,
            ),
        )
        predictions.append(prediction)
        metadata.append(meta)
    matrix = np.asarray(predictions)
    return matrix.mean(axis=0), list(matrix.mean(axis=1)), metadata


def reference_free_policy_validation(
    benchmark: str,
    seed: int,
    branch: str,
    extension_index: int,
    previous: nn.Module,
    actor: nn.Module,
    device: torch.device,
) -> dict[str, Any]:
    """Paired policy-value comparisons; no numerical reference is accessed."""
    driver = driver_for(benchmark)
    # Read x only. The independent reference fields in this common bank are not
    # accessed by training or selection code.
    data = evaluation_data(benchmark)
    x = data["x"][data["validation_indices"][:128]]
    initial = ConstantTangentActor(driver).to(device)
    new_mean, new_rep, new_meta = _evaluate_replicates(
        benchmark, seed, branch, extension_index, actor, x, device
    )
    old_mean, old_rep, old_meta = _evaluate_replicates(
        benchmark, seed, branch, extension_index, previous, x, device
    )
    a0_mean, a0_rep, a0_meta = _evaluate_replicates(
        benchmark, seed, branch, extension_index, initial, x, device
    )
    paired_previous = np.asarray(old_rep) - np.asarray(new_rep)
    paired_initial = np.asarray(a0_rep) - np.asarray(new_rep)

    def summary(replicates: np.ndarray, state_difference: np.ndarray) -> dict[str, float]:
        mean = float(np.mean(replicates))
        mc_se = float(np.std(replicates, ddof=1) / np.sqrt(len(replicates)))
        state_se = float(np.std(state_difference, ddof=1) / np.sqrt(len(state_difference)))
        # t_0.975,3 for four independent MC replications.
        return {
            "mean": mean, "mc_standard_error": mc_se,
            "mc_ci95_low": mean - 3.182446 * mc_se,
            "mc_ci95_high": mean + 3.182446 * mc_se,
            "state_standard_error": state_se,
            "state_ci95_low": mean - 1.959964 * state_se,
            "state_ci95_high": mean + 1.959964 * state_se,
        }

    previous_summary = summary(paired_previous, old_mean - new_mean)
    initial_summary = summary(paired_initial, a0_mean - new_mean)
    all_meta = new_meta + old_meta + a0_meta
    return {
        "interpretation": "positive means the new Actor has lower tangent-policy value",
        "previous_minus_new": previous_summary,
        "initial_a0_minus_new": initial_summary,
        "states": len(x), "mc_replicates": 4, "paths_per_replicate": 128,
        "time_steps": 32,
        "mean_weight_ess_fraction": float(
            np.mean([item["weight_ess_fraction"] for item in all_meta])
        ),
        "minimum_weight_ess_fraction": float(
            np.min([item["weight_ess_fraction"] for item in all_meta])
        ),
        "nonfinite_count": int(sum(item["nonfinite_count"] for item in all_meta)),
        "ground_truth_used": False,
    }


def _config(training_points: int, paths_per_label: int) -> FitConfig:
    return FitConfig(
        training_points=training_points, validation_points=256,
        paths_per_label=paths_per_label, label_time_steps=16,
        optimizer_steps=2500, width=192, depth=3, learning_rate=2e-4,
        batch_size=512, parent_batch=8, path_chunk=32,
        checkpoint_interval=125, gradient_clip=10.0,
    )


def train_branch(
    benchmark: str,
    seed: int,
    branch: str,
    device: torch.device,
) -> dict[str, Any]:
    if branch not in BRANCHES:
        raise KeyError(branch)
    ensure_dirs()
    driver = driver_for(benchmark)
    spec = BRANCHES[branch]
    config = _config(spec["training_points"], spec["paths_per_label"])
    source, _ = load_tangent_actor(_strong_a3_path(benchmark, seed), device)
    base_seconds = _base_training_seconds(benchmark, seed)
    cumulative = base_seconds
    phases: list[dict[str, Any]] = []
    no_improvement_streak = 0
    synchronize(device)
    resource_started = time.perf_counter()
    for extension_index in range(1, spec["phases"] + 1):
        synchronize(device)
        started = time.perf_counter()
        arrays, label_meta = _generate_labels(
            benchmark, seed, extension_index, source, config, device
        )
        actor, training, trace = train_actor(
            driver, config, "terminal_aware",
            arrays["train_t"], arrays["train_x"], arrays["train_r"], arrays.get("train_q"),
            arrays["valid_t"], arrays["valid_x"], arrays["valid_r"], arrays.get("valid_q"),
            seed=stable_seed("high_budget_fit", benchmark, seed, extension_index),
            device=device, initial=source,
        )
        policy = reference_free_policy_validation(
            benchmark, seed, branch, extension_index, source, actor, device
        )
        synchronize(device)
        phase_seconds_before_save = time.perf_counter() - started
        checkpoint_label = (
            f"a{3 + extension_index}" if branch == "iterations"
            else f"a4_{branch}"
        )
        checkpoint_path = (
            CHECKPOINT_DIR / f"tangent_actor_{benchmark}_s{seed}_{checkpoint_label}.pt"
        )
        atomic_torch_save(
            checkpoint_path,
            {
                "state_dict": {
                    name: value.detach().cpu().clone()
                    for name, value in actor.state_dict().items()
                },
                "actor_kind": "terminal_aware", "width": config.width,
                "depth": config.depth, "benchmark": benchmark, "seed": seed,
                "branch": branch, "extension_index": extension_index,
                "checkpoint_label": checkpoint_label,
            },
        )
        synchronize(device)
        phase_seconds = time.perf_counter() - started
        cumulative += phase_seconds
        improvement = policy["previous_minus_new"]
        meaningful = bool(
            improvement["mc_ci95_low"] > 0.0 and improvement["mean"] >= 1e-4
        )
        no_improvement_streak = 0 if meaningful else no_improvement_streak + 1
        phase = {
            "method": "tangent_actor_high_budget", "benchmark": benchmark,
            "seed": seed, "branch": branch, "extension_index": extension_index,
            "checkpoint_label": checkpoint_label, "checkpoint": str(checkpoint_path),
            "source_checkpoint": "a3" if extension_index == 1 else phases[-1]["checkpoint_label"],
            "phase_seconds": phase_seconds,
            "phase_seconds_before_checkpoint_io": phase_seconds_before_save,
            "base_a3_training_seconds": base_seconds,
            "cumulative_training_seconds": cumulative,
            "labels": label_meta, "training": training, "trace": trace,
            "reference_free_policy_validation": policy,
            "meaningful_reference_free_improvement": meaningful,
            "parameter_count": parameter_count(actor),
            "selection_criterion": "held-out tangent labels and paired reference-free policy value",
            "ground_truth_used_for_selection": False,
        }
        phases.append(phase)
        source = actor
        atomic_json(
            RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json",
            {
                "method": "tangent_actor_high_budget", "benchmark": benchmark,
                "seed": seed, "branch": branch, "configuration": asdict(config),
                "base_checkpoint": str(_strong_a3_path(benchmark, seed)),
                "base_a3_training_seconds": base_seconds, "phases": phases,
                "resource_seconds_this_run": time.perf_counter() - resource_started,
                "nested_bank_audit": nested_bank_audit(benchmark, seed),
                "ground_truth_used_for_training_or_selection": False,
            },
        )
        print(
            f"{benchmark} s{seed} {branch} {checkpoint_label}: "
            f"cum={cumulative:.2f}s paired={improvement['mean']:.6g} "
            f"CI=[{improvement['mc_ci95_low']:.6g},{improvement['mc_ci95_high']:.6g}]",
            flush=True,
        )
        if branch == "iterations" and extension_index >= 2 and no_improvement_streak >= 2:
            break
    result = {
        "method": "tangent_actor_high_budget", "benchmark": benchmark,
        "seed": seed, "branch": branch, "configuration": asdict(config),
        "base_checkpoint": str(_strong_a3_path(benchmark, seed)),
        "base_a3_training_seconds": base_seconds, "phases": phases,
        "resource_seconds_this_run": time.perf_counter() - resource_started,
        "nested_bank_audit": nested_bank_audit(benchmark, seed),
        "ground_truth_used_for_training_or_selection": False,
    }
    atomic_json(
        RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json", result
    )
    return result


def train_pilot_suite(
    benchmark: str, seed: int, device: torch.device
) -> list[dict[str, Any]]:
    results = []
    for branch in BRANCHES:
        results.append(train_branch(benchmark, seed, branch, device))
    return results


def continue_points8192(
    benchmark: str,
    seed: int,
    device: torch.device,
    *,
    maximum_additional_phases: int = 4,
    maximum_cumulative_seconds: float = 1200.0,
) -> dict[str, Any]:
    """Continue only the data-favored 8192-point branch after staged review."""
    ensure_dirs()
    pilot_path = RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_points8192.json"
    pilot = json.loads(pilot_path.read_text())
    pilot_phase = pilot["phases"][-1]
    config = _config(training_points=8192, paths_per_label=128)
    driver = driver_for(benchmark)
    source, _ = load_tangent_actor(Path(pilot_phase["checkpoint"]), device)
    cumulative = float(pilot_phase["cumulative_training_seconds"])
    previous_phase_seconds = float(pilot_phase["phase_seconds"])
    phases: list[dict[str, Any]] = []
    no_improvement_streak = 0
    synchronize(device)
    resource_started = time.perf_counter()
    for extension_index in range(2, 2 + maximum_additional_phases):
        # Do not begin a phase that is likely to exceed the +25% B5 ceiling.
        if cumulative + previous_phase_seconds > maximum_cumulative_seconds:
            break
        synchronize(device)
        started = time.perf_counter()
        arrays, label_meta = _generate_labels(
            benchmark, seed, extension_index, source, config, device
        )
        actor, training, trace = train_actor(
            driver, config, "terminal_aware",
            arrays["train_t"], arrays["train_x"], arrays["train_r"], arrays.get("train_q"),
            arrays["valid_t"], arrays["valid_x"], arrays["valid_r"], arrays.get("valid_q"),
            seed=stable_seed("high_budget_fit", benchmark, seed, extension_index),
            device=device, initial=source,
        )
        policy = reference_free_policy_validation(
            benchmark, seed, "points8192_continued", extension_index,
            source, actor, device,
        )
        synchronize(device)
        before_save = time.perf_counter() - started
        label = f"a{3 + extension_index}_points8192"
        checkpoint_path = CHECKPOINT_DIR / f"tangent_actor_{benchmark}_s{seed}_{label}.pt"
        atomic_torch_save(
            checkpoint_path,
            {
                "state_dict": {
                    name: value.detach().cpu().clone()
                    for name, value in actor.state_dict().items()
                },
                "actor_kind": "terminal_aware", "width": config.width,
                "depth": config.depth, "benchmark": benchmark, "seed": seed,
                "branch": "points8192_continued", "extension_index": extension_index,
                "checkpoint_label": label,
            },
        )
        synchronize(device)
        phase_seconds = time.perf_counter() - started
        cumulative += phase_seconds
        improvement = policy["previous_minus_new"]
        meaningful = bool(
            improvement["mc_ci95_low"] > 0.0 and improvement["mean"] >= 1e-4
        )
        no_improvement_streak = 0 if meaningful else no_improvement_streak + 1
        phase = {
            "method": "tangent_actor_high_budget", "benchmark": benchmark,
            "seed": seed, "branch": "points8192_continued",
            "extension_index": extension_index, "checkpoint_label": label,
            "checkpoint": str(checkpoint_path),
            "source_checkpoint": pilot_phase["checkpoint_label"]
            if extension_index == 2 else phases[-1]["checkpoint_label"],
            "phase_seconds": phase_seconds,
            "phase_seconds_before_checkpoint_io": before_save,
            "base_a3_training_seconds": pilot["base_a3_training_seconds"],
            "cumulative_training_seconds": cumulative,
            "labels": label_meta, "training": training, "trace": trace,
            "reference_free_policy_validation": policy,
            "meaningful_reference_free_improvement": meaningful,
            "parameter_count": parameter_count(actor),
            "selection_criterion": "held-out tangent labels and paired reference-free policy value",
            "ground_truth_used_for_selection": False,
        }
        phases.append(phase)
        source = actor
        previous_phase_seconds = phase_seconds
        partial = {
            "method": "tangent_actor_high_budget", "benchmark": benchmark,
            "seed": seed, "branch": "points8192_continued",
            "configuration": asdict(config), "source_pilot": str(pilot_path),
            "source_checkpoint": pilot_phase["checkpoint"],
            "pilot_a4_cumulative_training_seconds": pilot_phase["cumulative_training_seconds"],
            "phases": phases,
            "resource_seconds_this_run": time.perf_counter() - resource_started,
            "maximum_cumulative_seconds": maximum_cumulative_seconds,
            "staged_decision_basis": "8192 points had the largest paired improvement in all four pilot tasks",
            "ground_truth_used_for_training_or_selection": False,
        }
        atomic_json(
            RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_points8192_continued.json",
            partial,
        )
        print(
            f"{benchmark} s{seed} points8192_continued {label}: "
            f"cum={cumulative:.2f}s paired={improvement['mean']:.6g} "
            f"CI=[{improvement['mc_ci95_low']:.6g},{improvement['mc_ci95_high']:.6g}]",
            flush=True,
        )
        if no_improvement_streak >= 2:
            break
    result = {
        "method": "tangent_actor_high_budget", "benchmark": benchmark,
        "seed": seed, "branch": "points8192_continued",
        "configuration": asdict(config), "source_pilot": str(pilot_path),
        "source_checkpoint": pilot_phase["checkpoint"],
        "pilot_a4_cumulative_training_seconds": pilot_phase["cumulative_training_seconds"],
        "phases": phases,
        "resource_seconds_this_run": time.perf_counter() - resource_started,
        "maximum_cumulative_seconds": maximum_cumulative_seconds,
        "staged_decision_basis": "8192 points had the largest paired improvement in all four pilot tasks",
        "ground_truth_used_for_training_or_selection": False,
    }
    atomic_json(
        RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_points8192_continued.json",
        result,
    )
    return result


__all__ = [
    "BRANCHES", "CONTINUATION_BRANCHES", "nested_bank_audit", "reference_free_policy_validation",
    "continue_points8192", "train_branch", "train_pilot_suite",
]
