"""Staged MCV-trained Tangent-Actor continuations from the frozen a3 Actor."""

from __future__ import annotations

import argparse
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
)
from invariant_region_mlp.experiments.tangent_actor_high_budget.actor_scaling import (
    _base_training_seconds,
    _config,
    _state_bank,
    reference_free_policy_validation,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
    train_actor,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    load_tangent_actor,
)

from .common import (
    CHECKPOINT_DIR,
    LABEL_DIR,
    RAW_DIR,
    a3_checkpoint,
    atomic_json,
    atomic_torch_save,
    driver_for,
    ensure_dirs,
    stable_seed,
    state_fingerprint,
    synchronize,
)
from .labels import mcv_tangent_labels


def _label_cache_paths(
    benchmark: str, seed: int, points: int, extension_index: int
) -> tuple[Path, Path]:
    stem = f"mcv_{benchmark}_s{seed}_points{points}_e{extension_index}"
    return LABEL_DIR / f"{stem}.npz", LABEL_DIR / f"{stem}.json"


def generate_mcv_labels(
    benchmark: str,
    seed: int,
    extension_index: int,
    source: nn.Module,
    config: FitConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Generate paired MCV labels with the exact existing state/RNG protocol."""
    ensure_dirs()
    cache_path, metadata_path = _label_cache_paths(
        benchmark, seed, config.training_points, extension_index
    )
    if cache_path.exists() and metadata_path.exists():
        with np.load(cache_path) as archive:
            arrays = {key: archive[key] for key in archive.files}
        metadata = json.loads(metadata_path.read_text())
        metadata = {**metadata, "loaded_from_cache": True}
        return arrays, metadata

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
    train_value, train_q, train_meta = mcv_tangent_labels(
        driver, source, train_t, train_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "train"),
        **common,
    )
    valid_value, valid_q, valid_meta = mcv_tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "valid"),
        **common,
    )
    replicate_value, replicate_q, replicate_meta = mcv_tangent_labels(
        driver, source, valid_t, valid_x,
        seed=stable_seed("high_budget_labels", benchmark, seed, extension_index, "replicate"),
        **common,
    )
    train_r, train_q, train_projection = make_targets(driver, train_value, train_q)
    valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_q)
    assert train_q is not None and valid_q is not None and replicate_q is not None
    value_replication = float(np.sqrt(np.mean((valid_value - replicate_value) ** 2)))
    difference = valid_q - replicate_q
    gradient_component = float(np.sqrt(np.mean(difference**2)))
    gradient_vector = float(np.sqrt(np.mean(np.sum(difference**2, axis=1))))
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_r": train_r,
        "train_q": train_q,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_r": valid_r,
        "valid_q": valid_q,
    }
    metadata: dict[str, Any] = {
        "estimator": "direct_alpha1_mcv_full_corrected_payoff_gradient",
        "alpha": 1.0,
        "train": train_meta,
        "validation": valid_meta,
        "replication": replicate_meta,
        "train_projection": train_projection,
        "validation_projection": valid_projection,
        "value_label_replication_rmse": value_replication,
        "gradient_label_replication_component_rmse": gradient_component,
        "gradient_label_replication_vector_rmse": gradient_vector,
        "wall_seconds": float(
            train_meta["wall_seconds"] + valid_meta["wall_seconds"]
            + replicate_meta["wall_seconds"]
        ),
        "nested_state_bank_fingerprint": state_fingerprint(train_t, train_x),
        "fixed_validation_bank_fingerprint": state_fingerprint(valid_t, valid_x),
        "plain_protocol_rng_names_reused": True,
        "complete_corrected_payoff_differentiated": True,
        "source_actor_parameters_frozen": True,
        "source_actor_state_graph_preserved": True,
        "ground_truth_used": False,
        "loaded_from_cache": False,
    }
    temporary = cache_path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(cache_path)
    atomic_json(metadata_path, metadata)
    return arrays, metadata


def train_mcv_branch(
    benchmark: str,
    seed: int,
    *,
    training_points: int,
    maximum_phases: int,
    maximum_cumulative_seconds: float,
    device: torch.device,
) -> dict[str, Any]:
    """Continue MCV policy improvement from a3 using complete updates only."""
    ensure_dirs()
    branch = f"mcv_points{training_points}"
    output_path = RAW_DIR / f"{branch}_{benchmark}_s{seed}.json"
    config = _config(training_points, 128)
    driver = driver_for(benchmark)
    source, _ = load_tangent_actor(a3_checkpoint(benchmark, seed), device)
    base_seconds = _base_training_seconds(benchmark, seed)
    cumulative = base_seconds
    phases: list[dict[str, Any]] = []
    previous_phase_seconds: float | None = None

    # Resume only complete phases. Each cached label file is tied to the source
    # checkpoint sequence for this unique branch.
    if output_path.exists():
        existing = json.loads(output_path.read_text())
        phases = existing.get("phases", [])
        if phases:
            last = phases[-1]
            source, _ = load_tangent_actor(Path(last["checkpoint"]), device)
            cumulative = float(last["cumulative_training_seconds"])
            previous_phase_seconds = float(last["phase_seconds"])

    synchronize(device)
    resource_started = time.perf_counter()
    for extension_index in range(len(phases) + 1, maximum_phases + 1):
        if (
            previous_phase_seconds is not None
            and cumulative + previous_phase_seconds > maximum_cumulative_seconds
        ):
            break
        synchronize(device)
        started = time.perf_counter()
        arrays, label_metadata = generate_mcv_labels(
            benchmark, seed, extension_index, source, config, device
        )
        actor, training, trace = train_actor(
            driver,
            config,
            "terminal_aware",
            arrays["train_t"],
            arrays["train_x"],
            arrays["train_r"],
            arrays["train_q"],
            arrays["valid_t"],
            arrays["valid_x"],
            arrays["valid_r"],
            arrays["valid_q"],
            seed=stable_seed("high_budget_fit", benchmark, seed, extension_index),
            device=device,
            initial=source,
        )
        policy = reference_free_policy_validation(
            benchmark, seed, branch, extension_index, source, actor, device
        )
        synchronize(device)
        before_save = time.perf_counter() - started
        label = f"a{3 + extension_index}_mcv_points{training_points}"
        checkpoint_path = CHECKPOINT_DIR / f"tangent_actor_{benchmark}_s{seed}_{label}.pt"
        atomic_torch_save(
            checkpoint_path,
            {
                "state_dict": {
                    name: value.detach().cpu().clone()
                    for name, value in actor.state_dict().items()
                },
                "actor_kind": "terminal_aware",
                "width": config.width,
                "depth": config.depth,
                "benchmark": benchmark,
                "seed": seed,
                "branch": branch,
                "extension_index": extension_index,
                "checkpoint_label": label,
                "training_label_estimator": "direct_alpha1_mcv",
            },
        )
        synchronize(device)
        phase_seconds = time.perf_counter() - started
        cumulative += phase_seconds
        phase = {
            "method": "tangent_actor_mcv_training",
            "benchmark": benchmark,
            "seed": seed,
            "branch": branch,
            "extension_index": extension_index,
            "iteration": 3 + extension_index,
            "checkpoint_label": label,
            "checkpoint": str(checkpoint_path),
            "source_checkpoint": "a3" if extension_index == 1 else phases[-1]["checkpoint_label"],
            "phase_seconds": phase_seconds,
            "phase_seconds_before_checkpoint_io": before_save,
            "base_a3_training_seconds": base_seconds,
            "additional_continuation_seconds": cumulative - base_seconds,
            "cumulative_training_seconds": cumulative,
            "labels": label_metadata,
            "training": training,
            "trace": trace,
            "reference_free_policy_validation": policy,
            "selection_criterion": "minimum normalized held-out MCV tangent-label loss",
            "ground_truth_used_for_selection": False,
        }
        phases.append(phase)
        source = actor
        previous_phase_seconds = phase_seconds
        partial = {
            "method": "tangent_actor_mcv_training",
            "benchmark": benchmark,
            "seed": seed,
            "branch": branch,
            "configuration": asdict(config),
            "base_checkpoint": str(a3_checkpoint(benchmark, seed)),
            "base_a3_training_seconds": base_seconds,
            "maximum_cumulative_seconds": maximum_cumulative_seconds,
            "phases": phases,
            "resource_seconds_this_run": time.perf_counter() - resource_started,
            "ground_truth_used_for_training_or_selection": False,
        }
        atomic_json(output_path, partial)
        improvement = policy["previous_minus_new"]
        print(
            f"{benchmark} s{seed} {branch} {label}: cumulative={cumulative:.2f}s "
            f"paired_previous_minus_new={improvement['mean']:.6g} "
            f"CI=[{improvement['mc_ci95_low']:.6g},{improvement['mc_ci95_high']:.6g}]",
            flush=True,
        )
    result = {
        "method": "tangent_actor_mcv_training",
        "benchmark": benchmark,
        "seed": seed,
        "branch": branch,
        "configuration": asdict(config),
        "base_checkpoint": str(a3_checkpoint(benchmark, seed)),
        "base_a3_training_seconds": base_seconds,
        "maximum_cumulative_seconds": maximum_cumulative_seconds,
        "phases": phases,
        "resource_seconds_this_run": time.perf_counter() - resource_started,
        "ground_truth_used_for_training_or_selection": False,
    }
    atomic_json(output_path, result)
    return result


def run(args: argparse.Namespace) -> None:
    decision_path = Path(__file__).resolve().parents[2] / "results" / "tangent_actor_mcv_training" / "stage0_decision.json"
    if not decision_path.exists() or not json.loads(decision_path.read_text())["passed"]:
        raise RuntimeError("Stage 0 did not pass; refusing MCV retraining")
    device = torch.device(args.device)
    if args.stage == 1:
        for benchmark in ("nonseparable_d100_k2", "joint_concave_d400_k2"):
            train_mcv_branch(
                benchmark,
                0,
                training_points=2048,
                maximum_phases=1,
                maximum_cumulative_seconds=500.0,
                device=device,
            )
    else:
        stage1_path = decision_path.with_name("stage1_decision.json")
        if not stage1_path.exists() or not json.loads(stage1_path.read_text())["passed"]:
            raise RuntimeError("Stage 1 did not pass; refusing high-budget Stage 2")
        for benchmark in ("nonseparable_d100_k2", "joint_concave_d400_k2"):
            for seed in (0, 1):
                train_mcv_branch(
                    benchmark,
                    seed,
                    training_points=8192,
                    maximum_phases=5,
                    maximum_cumulative_seconds=1200.0,
                    device=device,
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=int, choices=(1, 2), required=True)
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
