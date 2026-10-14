"""Run the three-problem teacher-free tangent-Actor feasibility pilot."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import platform
import subprocess
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from .common import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    atomic_torch_save,
    ensure_dirs,
    fixed_points,
    model_metrics,
    stable_seed,
)
from .core import (
    ActorConfig,
    actor_tangents,
    evaluate_tangent_mc,
    make_targets,
    tangent_labels,
    train_actor,
)
from .diagnostics import bregman_identity_audit
from .drivers import ConcaveSemilinearDriver, make_driver
from .models import ConstantTangentActor, TangentActor
from .reference import PeriodicReference, ReferenceTangentActor, solve_periodic_reference


PROBLEMS = ("allen_cahn", "joint_concave", "joint_convex")


def load_protocol() -> dict[str, Any]:
    with (RESULTS_DIR / "experiment_protocol.json").open() as handle:
        return json.load(handle)


def _load_or_build_labels(
    driver: ConcaveSemilinearDriver,
    source: nn.Module,
    iteration: int,
    config: ActorConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache_path = CACHE_DIR / f"{driver.name}_a{iteration}_labels.npz"
    metadata_path = CACHE_DIR / f"{driver.name}_a{iteration}_labels.json"
    if cache_path.exists() and metadata_path.exists():
        with np.load(cache_path) as data:
            arrays = {key: data[key] for key in data.files}
        with metadata_path.open() as handle:
            return arrays, json.load(handle)

    train_t, train_x = fixed_points(
        driver.name, driver.d, f"iteration_{iteration}_train", config.training_points,
        T=driver.T, include_time=True,
    )
    valid_t, valid_x = fixed_points(
        driver.name, driver.d, f"iteration_{iteration}_label_valid", config.validation_points,
        T=driver.T, include_time=True,
    )
    train_value, train_z, train_meta = tangent_labels(
        driver, source, train_t, train_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(driver.name, iteration, "train_labels"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    valid_value, valid_z, valid_meta = tangent_labels(
        driver, source, valid_t, valid_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(driver.name, iteration, "valid_labels"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    replicate_value, replicate_z, replicate_meta = tangent_labels(
        driver, source, valid_t, valid_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(driver.name, iteration, "valid_labels_replication"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    train_r, train_q, train_projection = make_targets(driver, train_value, train_z)
    valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_z)
    arrays: dict[str, np.ndarray] = {
        "train_t": train_t,
        "train_x": train_x,
        "train_raw_value": train_value,
        "train_r": train_r,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_raw_value": valid_value,
        "valid_r": valid_r,
    }
    if train_q is not None and valid_q is not None:
        arrays["train_q"] = train_q
        arrays["valid_q"] = valid_q
    value_replication_rmse = float(np.sqrt(np.mean((valid_value - replicate_value) ** 2)))
    if valid_z is None:
        gradient_replication_rmse = None
        gradient_replication_relative_l2 = None
    else:
        assert replicate_z is not None
        gradient_replication_rmse = float(np.sqrt(np.mean((valid_z - replicate_z) ** 2)))
        gradient_replication_relative_l2 = float(
            np.linalg.norm(valid_z - replicate_z)
            / max(np.linalg.norm(0.5 * (valid_z + replicate_z)), np.finfo(float).tiny)
        )
    metadata = {
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "train": train_meta,
        "valid": valid_meta,
        "replication": replicate_meta,
        "train_projection": train_projection,
        "valid_projection": valid_projection,
        "value_replication_rmse": value_replication_rmse,
        "gradient_replication_rmse": gradient_replication_rmse,
        "gradient_replication_relative_l2": gradient_replication_relative_l2,
        "train_bank_sha256": array_sha256(train_t, train_x),
        "valid_bank_sha256": array_sha256(valid_t, valid_x),
        "reference_information_used": False,
    }
    atomic_npz(cache_path, **arrays)
    atomic_json(metadata_path, metadata)
    return arrays, metadata


def _load_or_train_actor(
    driver: ConcaveSemilinearDriver,
    iteration: int,
    source: nn.Module,
    arrays: dict[str, np.ndarray],
    label_metadata: dict[str, Any],
    config: ActorConfig,
    device: torch.device,
) -> tuple[TangentActor, dict[str, Any]]:
    checkpoint_path = CHECKPOINT_DIR / f"{driver.name}_a{iteration}.pt"
    metadata_path = CHECKPOINT_DIR / f"{driver.name}_a{iteration}.json"
    actor = TangentActor(driver, width=config.width, depth=config.depth).to(device)
    if checkpoint_path.exists() and metadata_path.exists():
        actor.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
        actor.eval()
        with metadata_path.open() as handle:
            return actor, json.load(handle)
    initial = source if isinstance(source, TangentActor) else None
    actor, training, trace = train_actor(
        driver, config,
        arrays["train_t"], arrays["train_x"], arrays["train_r"], arrays.get("train_q"),
        arrays["valid_t"], arrays["valid_x"], arrays["valid_r"], arrays.get("valid_q"),
        seed=stable_seed(driver.name, iteration, "actor_training"),
        device=device,
        initial=initial,
    )
    metadata = {
        "problem": driver.name,
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "label": label_metadata,
        "training": training,
        "trace": trace,
        "checkpoint_selected_without_exact_reference": True,
    }
    atomic_torch_save(checkpoint_path, actor.state_dict())
    atomic_json(metadata_path, metadata)
    return actor, metadata


def _exact_working_values(
    driver: ConcaveSemilinearDriver,
    reference: PeriodicReference,
    t: np.ndarray,
    x: np.ndarray,
) -> np.ndarray:
    return reference.value(t, x @ driver.beta)


def _posthoc_tangent_metrics(
    driver: ConcaveSemilinearDriver,
    reference: PeriodicReference,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray | None, np.ndarray, np.ndarray | None]:
    predicted_r, predicted_q = actor_tangents(driver, actor, t, x, device)
    s = x @ driver.beta
    exact_r = reference.value(t, s)
    r_error = predicted_r - exact_r
    metrics: dict[str, Any] = {
        "tangent_r_rmse": float(np.sqrt(np.mean(r_error**2))),
        "tangent_r_relative_l2": float(
            np.linalg.norm(r_error) / max(np.linalg.norm(exact_r), np.finfo(float).tiny)
        ),
        "tangent_r_min": float(np.min(predicted_r)),
        "tangent_r_max": float(np.max(predicted_r)),
        "tangent_r_domain_violation_rate": float(
            np.mean((predicted_r < driver.r_low - 1e-7) | (predicted_r > driver.r_high + 1e-7))
        ),
    }
    if driver.has_gradient:
        assert predicted_q is not None
        exact_q = (
            np.sqrt(2.0) * reference.derivative(t, s)[:, None] * driver.beta[None, :]
        )
        q_error = predicted_q - exact_q
        metrics.update({
            "tangent_q_component_rmse": float(np.sqrt(np.mean(q_error**2))),
            "tangent_q_vector_rmse": float(np.sqrt(np.mean(np.sum(q_error**2, axis=1)))),
            "tangent_q_relative_l2": float(
                np.linalg.norm(q_error) / max(np.linalg.norm(exact_q), np.finfo(float).tiny)
            ),
            "tangent_q_norm_mean": float(np.mean(np.linalg.norm(predicted_q, axis=1))),
            "exact_q_norm_mean": float(np.mean(np.linalg.norm(exact_q, axis=1))),
        })
    else:
        exact_q = None
        metrics.update({
            "tangent_q_component_rmse": None,
            "tangent_q_vector_rmse": None,
            "tangent_q_relative_l2": None,
            "tangent_q_norm_mean": None,
            "exact_q_norm_mean": None,
        })
    return metrics, predicted_r, predicted_q, exact_r, exact_q


def run_problem(name: str, device: torch.device, protocol: dict[str, Any]) -> None:
    driver = make_driver(name, d=int(protocol["dimension"]))
    output_path = RESULTS_DIR / f"raw_{name}.json"
    if output_path.exists():
        print(f"SKIP completed {name}: {output_path}", flush=True)
        return
    config = ActorConfig(**protocol["actor"])
    actors: list[nn.Module] = [ConstantTangentActor(driver).to(device)]
    phase_metadata: list[dict[str, Any]] = []
    for iteration in (1, 2):
        arrays, label_metadata = _load_or_build_labels(
            driver, actors[-1], iteration, config, device
        )
        actor, metadata = _load_or_train_actor(
            driver, iteration, actors[-1], arrays, label_metadata, config, device
        )
        actors.append(actor)
        phase_metadata.append(metadata)
        print(
            f"{name} a{iteration}: labels="
            f"{label_metadata['train']['wall_seconds'] + label_metadata['valid']['wall_seconds']:.1f}s "
            f"train={metadata['training']['wall_seconds']:.1f}s "
            f"fit={metadata['training']['best_validation_normalized_loss']:.5g}",
            flush=True,
        )

    evaluation = protocol["evaluation"]
    validation_t, validation_x = fixed_points(
        name, driver.d, "policy_valid", evaluation["policy_validation_points"],
        T=driver.T, include_time=False,
    )
    validation_outputs: list[dict[str, np.ndarray]] = []
    validation_metadata: list[dict[str, Any]] = []
    for actor in actors:
        output, metadata = evaluate_tangent_mc(
            driver, actor, validation_t, validation_x,
            paths=evaluation["policy_validation_paths"],
            time_steps=evaluation["policy_validation_steps"],
            seed=stable_seed(name, "policy_validation_crn"),
            device=device,
        )
        validation_outputs.append(output)
        validation_metadata.append(metadata)

    # Reference construction begins only after all actor checkpoints and the
    # reference-free iteration comparison have been frozen.
    reference = solve_periodic_reference(driver)
    exact_actor = ReferenceTangentActor(driver, reference).to(device)
    actor_names = ["a0", "a1", "a2", "exact_tangent"]
    all_actors = actors + [exact_actor]
    test_t, test_x = fixed_points(
        name, driver.d, "test", evaluation["test_points"],
        T=driver.T, include_time=False,
    )
    exact_working = _exact_working_values(driver, reference, test_t, test_x)
    exact_original = driver.working_to_original(exact_working)
    audit_t, audit_x = fixed_points(
        name, driver.d, "tangent_audit", evaluation["audit_points"],
        T=driver.T, include_time=True,
    )
    rows: list[dict[str, Any]] = []
    saved: dict[str, np.ndarray] = {
        "test_t": test_t,
        "test_x": test_x,
        "reference_working_value": exact_working,
        "reference_original_value": exact_original,
        "audit_t": audit_t,
        "audit_x": audit_x,
    }
    main_variant = evaluation["mc_floor_variants"][0]
    main_outputs: dict[str, dict[str, np.ndarray]] = {}
    main_metadata: dict[str, dict[str, Any]] = {}
    for actor_name, actor in zip(actor_names, all_actors, strict=True):
        output, metadata = evaluate_tangent_mc(
            driver, actor, test_t, test_x,
            paths=main_variant["paths"],
            time_steps=main_variant["steps"],
            seed=stable_seed(name, "test_crn", main_variant["name"]),
            device=device,
        )
        prediction = driver.working_to_original(output["working_value"])
        metrics = model_metrics(prediction, exact_original)
        tangent_metrics, predicted_r, predicted_q, exact_r, exact_q = _posthoc_tangent_metrics(
            driver, reference, actor, audit_t, audit_x, device
        )
        phase_seconds = metadata["wall_seconds"]
        if actor_name in ("a1", "a2"):
            phase = phase_metadata[int(actor_name[1:]) - 1]
            phase_seconds += phase["label"]["train"]["wall_seconds"]
            phase_seconds += phase["label"]["valid"]["wall_seconds"]
            phase_seconds += phase["label"]["replication"]["wall_seconds"]
            phase_seconds += phase["training"]["wall_seconds"]
        rows.append({
            "pde": name,
            "d": driver.d,
            "iteration": actor_name,
            **metrics,
            **tangent_metrics,
            "final_mc_mean_se": metadata["mean_value_mc_se"],
            "final_mc_paths": main_variant["paths"],
            "final_mc_steps": main_variant["steps"],
            "inference_seconds": metadata["wall_seconds"],
            "phase_gpu_seconds": phase_seconds,
            "c_min": metadata["c_min"],
            "c_max": metadata["c_max"],
            "c_mean": metadata["c_mean"],
            "c_std": metadata["c_std"],
            "terminal_weight_min": metadata["weight_min"],
            "terminal_weight_max": metadata["weight_max"],
            "terminal_weight_mean": metadata["weight_mean"],
            "terminal_weight_std": metadata["weight_std"],
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "nonfinite_count": metadata["nonfinite_count"],
        })
        main_outputs[actor_name] = output
        main_metadata[actor_name] = metadata
        saved[f"{actor_name}_working_value"] = output["working_value"]
        saved[f"{actor_name}_original_value"] = prediction
        saved[f"{actor_name}_mc_variance_of_mean"] = output["mc_variance_of_mean"]
        saved[f"{actor_name}_audit_r"] = predicted_r
        if actor_name == "a0":
            saved["exact_audit_r"] = exact_r
            if exact_q is not None:
                saved["exact_audit_q"] = exact_q
        if predicted_q is not None:
            saved[f"{actor_name}_audit_q"] = predicted_q
        print(
            f"{name} {actor_name}: RMSE={metrics['rmse']:.6g} "
            f"r={tangent_metrics['tangent_r_rmse']:.6g} "
            f"q={tangent_metrics['tangent_q_vector_rmse']}",
            flush=True,
        )

    iteration_rows: list[dict[str, Any]] = []
    for index, actor_name in enumerate(("a0", "a1", "a2")):
        values = validation_outputs[index]["working_value"]
        if index:
            previous = validation_outputs[index - 1]["working_value"]
            delta = values - previous
            paired_delta = float(np.mean(delta))
            paired_se = float(np.std(delta, ddof=1) / np.sqrt(len(delta)))
            phase = phase_metadata[index - 1]
        else:
            paired_delta = paired_se = None
            phase = None
        iteration_rows.append({
            "pde": name,
            "iteration": actor_name,
            "reference_free_weighted_value": float(np.mean(values)),
            "reference_free_value_mc_se": validation_metadata[index]["mean_value_mc_se"],
            "paired_delta_vs_previous": paired_delta,
            "paired_delta_se_across_states": paired_se,
            "value_label_replication_rmse": (
                phase["label"]["value_replication_rmse"] if phase else None
            ),
            "gradient_label_replication_rmse": (
                phase["label"]["gradient_replication_rmse"] if phase else None
            ),
            "actor_validation_normalized_loss": (
                phase["training"]["best_validation_normalized_loss"] if phase else None
            ),
            "actor_validation_value_mse": (
                phase["training"]["best_validation_value_mse"] if phase else None
            ),
            "actor_validation_gradient_mse": (
                phase["training"]["best_validation_gradient_mse"] if phase else None
            ),
            "selected_optimizer_step": phase["training"]["selected_step"] if phase else None,
            "label_seconds": (
                phase["label"]["train"]["wall_seconds"]
                + phase["label"]["valid"]["wall_seconds"]
                + phase["label"]["replication"]["wall_seconds"] if phase else 0.0
            ),
            "training_seconds": phase["training"]["wall_seconds"] if phase else 0.0,
            "max_label_weight": phase["label"]["train"]["weight_max"] if phase else None,
            "label_weight_ess_fraction": (
                phase["label"]["train"]["weight_ess_fraction"] if phase else None
            ),
            "label_gradient_norm_mean": (
                phase["label"]["train"]["gradient_norm_mean"] if phase else None
            ),
            "label_nonfinite_count": phase["label"]["train"]["nonfinite_count"] if phase else 0,
        })
        saved[f"{actor_name}_validation_working_value"] = values

    floor_rows: list[dict[str, Any]] = []
    # The main configuration is already evaluated above. Add higher-path and
    # higher-step comparisons for the final learned and exact tangent Actors.
    for variant in evaluation["mc_floor_variants"]:
        for actor_name, actor in (("a2", actors[2]), ("exact_tangent", exact_actor)):
            if variant["name"] == main_variant["name"]:
                output = main_outputs[actor_name]
                metadata = main_metadata[actor_name]
            else:
                output, metadata = evaluate_tangent_mc(
                    driver, actor, test_t, test_x,
                    paths=variant["paths"], time_steps=variant["steps"],
                    seed=stable_seed(name, "test_crn", variant["name"]),
                    device=device,
                )
            prediction = driver.working_to_original(output["working_value"])
            metrics = model_metrics(prediction, exact_original)
            floor_rows.append({
                "pde": name,
                "actor": actor_name,
                "variant": variant["name"],
                "paths": variant["paths"],
                "steps": variant["steps"],
                **metrics,
                "mean_mc_se": metadata["mean_value_mc_se"],
                "wall_seconds": metadata["wall_seconds"],
                "terminal_weight_max": metadata["weight_max"],
                "weight_ess_fraction": metadata["weight_ess_fraction"],
                "nonfinite_count": metadata["nonfinite_count"],
            })
            saved[f"floor_{actor_name}_{variant['name']}_original_value"] = prediction
            saved[f"floor_{actor_name}_{variant['name']}_mc_variance_of_mean"] = output[
                "mc_variance_of_mean"
            ]

    bregman_t, bregman_x = fixed_points(
        name, driver.d, "bregman_audit", evaluation["bregman_points"],
        T=driver.T, include_time=False,
    )
    bregman_summary, bregman_arrays = bregman_identity_audit(
        driver, actors[2], reference, bregman_t, bregman_x,
        paths=evaluation["bregman_paths"],
        time_steps=evaluation["bregman_steps"],
        seed=stable_seed(name, "bregman_audit"),
        device=device,
    )
    for key, value in bregman_arrays.items():
        saved[f"bregman_{key}"] = value
    prediction_path = PREDICTION_DIR / f"{name}.npz"
    atomic_npz(prediction_path, **saved)
    raw = {
        "problem": driver.metadata(),
        "protocol": protocol,
        "rows": rows,
        "tangent_iteration_rows": iteration_rows,
        "mc_floor_rows": floor_rows,
        "bregman": bregman_summary,
        "phase_metadata": phase_metadata,
        "prediction_path": str(prediction_path.relative_to(RESULTS_DIR.parent.parent)),
        "test_bank_sha256": array_sha256(test_t, test_x),
        "reference_created_after_checkpoint_and_validation_freeze": True,
        "exact_reference_used_for_selection": False,
        "a3_run": False,
        "a3_reason": "rapid primary feasibility protocol stopped after the required two updates",
    }
    atomic_json(output_path, raw)


def environment_payload() -> dict[str, Any]:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    return {
        "host": platform.node(),
        "python": platform.python_version(),
        "torch": torch.__version__,
        "cuda": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "git_commit_at_run_start": commit,
        "started_unix": time.time(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--problems", nargs="+", default=list(PROBLEMS), choices=PROBLEMS)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required for the GPU pilot")
    protocol = load_protocol()
    atomic_json(RESULTS_DIR / "runtime_environment.json", environment_payload())
    for problem in args.problems:
        run_problem(problem, device, protocol)


if __name__ == "__main__":
    main()
