"""Run the nonseparable and effective-dimension Tangent-Actor study."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
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
from .drivers import ConcaveSemilinearDriver, make_driver, terminal_signal_diagnostics
from .models import ConstantTangentActor, TangentActor
from .reference2d import (
    PeriodicReference2D,
    ReferenceTangentActor2D,
    validate_periodic_reference_2d,
)


FAMILIES = ("allen_cahn", "joint_concave", "joint_convex")


@dataclass(frozen=True)
class ExperimentSpec:
    family: str
    d: int
    k: int
    training_seed: int

    @property
    def base_id(self) -> str:
        return f"{self.family}_d{self.d}_k{self.k}"

    @property
    def run_id(self) -> str:
        return f"{self.base_id}_s{self.training_seed}"


def initial_specs() -> list[ExperimentSpec]:
    result = [ExperimentSpec("nonseparable", 100, 2, seed) for seed in (0, 1, 2)]
    result.extend(
        ExperimentSpec(family, 100, k, 0)
        for family in FAMILIES
        for k in (2, 5, 10)
    )
    return result


def expansion_specs() -> list[ExperimentSpec]:
    return [
        ExperimentSpec(family, d, k, 0)
        for family in FAMILIES
        for d, k in ((100, 20), (400, 2), (400, 5), (400, 10), (400, 20))
    ]


def load_protocol() -> dict[str, Any]:
    with (RESULTS_DIR / "experiment_protocol.json").open() as handle:
        return json.load(handle)


def _driver(spec: ExperimentSpec) -> ConcaveSemilinearDriver:
    return make_driver(spec.family, d=spec.d, k=spec.k)


def _load_or_build_labels(
    spec: ExperimentSpec,
    driver: ConcaveSemilinearDriver,
    source: nn.Module,
    iteration: int,
    config: ActorConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    stem = f"{spec.run_id}_a{iteration}_labels"
    cache_path, metadata_path = CACHE_DIR / f"{stem}.npz", CACHE_DIR / f"{stem}.json"
    if cache_path.exists() and metadata_path.exists():
        with np.load(cache_path) as data:
            arrays = {key: data[key] for key in data.files}
        with metadata_path.open() as handle:
            return arrays, json.load(handle)

    train_t, train_x = fixed_points(
        spec.base_id,
        driver.d,
        f"iteration_{iteration}_train",
        config.training_points,
        T=driver.T,
        include_time=True,
        training_seed=spec.training_seed,
    )
    valid_t, valid_x = fixed_points(
        spec.base_id,
        driver.d,
        f"iteration_{iteration}_label_valid",
        config.validation_points,
        T=driver.T,
        include_time=True,
        training_seed=spec.training_seed,
    )
    label_kwargs = dict(
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    train_value, train_z, train_meta = tangent_labels(
        driver,
        source,
        train_t,
        train_x,
        seed=stable_seed(spec.run_id, iteration, "train_labels"),
        **label_kwargs,
    )
    valid_value, valid_z, valid_meta = tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(spec.run_id, iteration, "valid_labels"),
        **label_kwargs,
    )
    replicate_value, replicate_z, replicate_meta = tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        seed=stable_seed(spec.run_id, iteration, "valid_labels_replication"),
        **label_kwargs,
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
        arrays["train_q"], arrays["valid_q"] = train_q, valid_q
    value_replication_rmse = float(np.sqrt(np.mean((valid_value - replicate_value) ** 2)))
    if valid_z is None:
        gradient_component_replication_rmse = None
        gradient_vector_replication_rmse = None
        gradient_replication_relative_l2 = None
    else:
        assert replicate_z is not None
        difference = valid_z - replicate_z
        gradient_component_replication_rmse = float(np.sqrt(np.mean(difference**2)))
        gradient_vector_replication_rmse = float(
            np.sqrt(np.mean(np.square(difference).sum(axis=1)))
        )
        gradient_replication_relative_l2 = float(
            np.linalg.norm(difference)
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
        "gradient_component_replication_rmse": gradient_component_replication_rmse,
        "gradient_vector_replication_rmse": gradient_vector_replication_rmse,
        "gradient_replication_relative_l2": gradient_replication_relative_l2,
        "train_bank_sha256": array_sha256(train_t, train_x),
        "valid_bank_sha256": array_sha256(valid_t, valid_x),
        "analytic_nonterminal_reference_used": False,
        "oracle_subspace_projection_used": False,
    }
    atomic_npz(cache_path, **arrays)
    atomic_json(metadata_path, metadata)
    return arrays, metadata


def _load_or_train_actor(
    spec: ExperimentSpec,
    driver: ConcaveSemilinearDriver,
    iteration: int,
    source: nn.Module,
    arrays: dict[str, np.ndarray],
    label_metadata: dict[str, Any],
    config: ActorConfig,
    device: torch.device,
) -> tuple[TangentActor, dict[str, Any]]:
    stem = f"{spec.run_id}_a{iteration}"
    checkpoint_path, metadata_path = CHECKPOINT_DIR / f"{stem}.pt", CHECKPOINT_DIR / f"{stem}.json"
    actor = TangentActor(driver, width=config.width, depth=config.depth).to(device)
    if checkpoint_path.exists() and metadata_path.exists():
        actor.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
        actor.eval()
        with metadata_path.open() as handle:
            return actor, json.load(handle)
    initial = source if isinstance(source, TangentActor) else None
    actor, training, trace = train_actor(
        driver,
        config,
        arrays["train_t"],
        arrays["train_x"],
        arrays["train_r"],
        arrays.get("train_q"),
        arrays["valid_t"],
        arrays["valid_x"],
        arrays["valid_r"],
        arrays.get("valid_q"),
        seed=stable_seed(spec.run_id, iteration, "actor_training"),
        device=device,
        initial=initial,
    )
    metadata = {
        "run_id": spec.run_id,
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "label": label_metadata,
        "training": training,
        "trace": trace,
        "checkpoint_selected_by": "held-out noisy tangent-label loss only",
        "ground_truth_used": False,
    }
    atomic_torch_save(checkpoint_path, actor.state_dict())
    atomic_json(metadata_path, metadata)
    return actor, metadata


def _load_actors(
    spec: ExperimentSpec, driver: ConcaveSemilinearDriver, config: ActorConfig, device: torch.device
) -> list[nn.Module]:
    actors: list[nn.Module] = [ConstantTangentActor(driver).to(device)]
    for iteration in (1, 2, 3):
        actor = TangentActor(driver, width=config.width, depth=config.depth).to(device)
        actor.load_state_dict(
            torch.load(
                CHECKPOINT_DIR / f"{spec.run_id}_a{iteration}.pt",
                map_location=device,
                weights_only=True,
            )
        )
        actor.eval()
        actors.append(actor)
    return actors


def run_training_and_reference_free(
    spec: ExperimentSpec, device: torch.device, protocol: dict[str, Any]
) -> None:
    output_path = RESULTS_DIR / f"training_{spec.run_id}.json"
    if output_path.exists():
        print(f"SKIP frozen training {spec.run_id}", flush=True)
        return
    driver = _driver(spec)
    config = ActorConfig(**protocol["actor"])
    actors: list[nn.Module] = [ConstantTangentActor(driver).to(device)]
    phases: list[dict[str, Any]] = []
    for iteration in (1, 2, 3):
        arrays, label_metadata = _load_or_build_labels(
            spec, driver, actors[-1], iteration, config, device
        )
        actor, metadata = _load_or_train_actor(
            spec, driver, iteration, actors[-1], arrays, label_metadata, config, device
        )
        actors.append(actor)
        phases.append(metadata)
        print(
            f"{spec.run_id} a{iteration}: label="
            f"{label_metadata['train']['wall_seconds'] + label_metadata['valid']['wall_seconds']:.1f}s "
            f"train={metadata['training']['wall_seconds']:.1f}s "
            f"fit={metadata['training']['best_validation_normalized_loss']:.5g}",
            flush=True,
        )

    evaluation = protocol["evaluation"]
    validation_t, validation_x = fixed_points(
        spec.base_id,
        driver.d,
        "policy_validation",
        evaluation["policy_validation_points"],
        T=driver.T,
        include_time=False,
        training_seed=protocol["fixed_evaluation_seed"],
    )
    outputs, metadata = [], []
    for actor in actors:
        value, details = evaluate_tangent_mc(
            driver,
            actor,
            validation_t,
            validation_x,
            paths=evaluation["policy_validation_paths"],
            time_steps=evaluation["policy_validation_steps"],
            seed=stable_seed(spec.base_id, "policy_validation_crn"),
            device=device,
            parent_batch=evaluation["evaluation_parent_batch"],
            path_chunk=evaluation["evaluation_path_chunk"],
        )
        outputs.append(value)
        metadata.append(details)
    rows: list[dict[str, Any]] = []
    for index, actor_name in enumerate(("a0", "a1", "a2", "a3")):
        values = outputs[index]["working_value"]
        phase = None if index == 0 else phases[index - 1]
        delta = None if index == 0 else values - outputs[index - 1]["working_value"]
        details = metadata[index]
        gradient_mse = None if phase is None else phase["training"]["best_validation_gradient_mse"]
        rows.append(
            {
                "run_id": spec.run_id,
                "family": spec.family,
                "d": spec.d,
                "k": spec.k,
                "training_seed": spec.training_seed,
                "iteration": actor_name,
                "reference_free_working_value_mean": float(np.mean(values)),
                "reference_free_value_mc_se": details["mean_value_mc_se"],
                "paired_delta_vs_previous": None if delta is None else float(np.mean(delta)),
                "paired_delta_se_across_states": None
                if delta is None
                else float(np.std(delta, ddof=1) / np.sqrt(len(delta))),
                "actor_validation_value_rmse": None
                if phase is None
                else float(np.sqrt(phase["training"]["best_validation_value_mse"])),
                "actor_validation_gradient_component_rmse": None
                if gradient_mse is None
                else float(np.sqrt(gradient_mse)),
                "actor_validation_gradient_vector_rmse": None
                if gradient_mse is None
                else float(np.sqrt(driver.d * gradient_mse)),
                "actor_validation_normalized_loss": None
                if phase is None
                else phase["training"]["best_validation_normalized_loss"],
                "selected_optimizer_step": None
                if phase is None
                else phase["training"]["selected_step"],
                "value_label_replication_rmse": None
                if phase is None
                else phase["label"]["value_replication_rmse"],
                "gradient_label_replication_component_rmse": None
                if phase is None
                else phase["label"]["gradient_component_replication_rmse"],
                "gradient_label_replication_vector_rmse": None
                if phase is None
                else phase["label"]["gradient_vector_replication_rmse"],
                "label_weight_ess_fraction": None
                if phase is None
                else phase["label"]["train"]["weight_ess_fraction"],
                "label_weight_max": None
                if phase is None
                else phase["label"]["train"]["weight_max"],
                "evaluation_weight_ess_fraction": details["weight_ess_fraction"],
                "evaluation_weight_min": details["weight_min"],
                "evaluation_weight_max": details["weight_max"],
                "evaluation_nonfinite_count": details["nonfinite_count"],
                "label_seconds": 0.0
                if phase is None
                else sum(
                    phase["label"][split]["wall_seconds"]
                    for split in ("train", "valid", "replication")
                ),
                "training_seconds": 0.0
                if phase is None
                else phase["training"]["wall_seconds"],
                "inference_seconds": details["wall_seconds"],
                "actor_parameters": 0
                if phase is None
                else phase["training"]["parameters"],
            }
        )
    saved = {
        "validation_t": validation_t,
        "validation_x": validation_x,
        **{
            f"a{index}_working_value": outputs[index]["working_value"]
            for index in range(4)
        },
        **{
            f"a{index}_mc_variance_of_mean": outputs[index]["mc_variance_of_mean"]
            for index in range(4)
        },
    }
    prediction_path = PREDICTION_DIR / f"reference_free_{spec.run_id}.npz"
    atomic_npz(prediction_path, **saved)
    atomic_json(
        output_path,
        {
            "spec": spec.__dict__,
            "problem": driver.metadata(),
            "terminal_signal": terminal_signal_diagnostics(driver),
            "protocol": protocol,
            "iteration_rows": rows,
            "phase_metadata": phases,
            "prediction_path": str(prediction_path.relative_to(driver_path_root())),
            "policy_validation_bank_sha256": array_sha256(validation_t, validation_x),
            "ground_truth_available_during_checkpoint_selection": False,
            "ground_truth_used_for_checkpoint_or_iteration_selection": False,
        },
    )


def driver_path_root():
    return RESULTS_DIR.parent.parent


def _posthoc_actor_metrics(
    driver: ConcaveSemilinearDriver,
    reference: PeriodicReference2D,
    actor: nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    predicted_r, predicted_q = actor_tangents(driver, actor, t, x, device)
    coordinates = x @ driver.B[:, :2]
    exact_r = reference.value(t, coordinates)
    r_error = predicted_r - exact_r
    result: dict[str, Any] = {
        "actor_exact_value_rmse": float(np.sqrt(np.mean(r_error**2))),
        "actor_exact_value_relative_l2": float(
            np.linalg.norm(r_error) / max(np.linalg.norm(exact_r), np.finfo(float).tiny)
        ),
        "actor_value_min": float(np.min(predicted_r)),
        "actor_value_max": float(np.max(predicted_r)),
    }
    arrays: dict[str, np.ndarray] = {"predicted_r": predicted_r, "exact_r": exact_r}
    if driver.has_gradient:
        assert predicted_q is not None
        reduced = reference.gradient(t, coordinates)
        exact_q = np.sqrt(2.0) * (reduced @ driver.B[:, :2].T)
        error = predicted_q - exact_q
        result.update(
            {
                "actor_exact_gradient_component_rmse": float(np.sqrt(np.mean(error**2))),
                "actor_exact_gradient_vector_rmse": float(
                    np.sqrt(np.mean(np.square(error).sum(axis=1)))
                ),
                "actor_exact_gradient_relative_l2": float(
                    np.linalg.norm(error) / max(np.linalg.norm(exact_q), np.finfo(float).tiny)
                ),
                "actor_predicted_gradient_norm_mean": float(
                    np.mean(np.linalg.norm(predicted_q, axis=1))
                ),
                "actor_exact_gradient_norm_mean": float(np.mean(np.linalg.norm(exact_q, axis=1))),
            }
        )
        arrays.update({"predicted_q": predicted_q, "exact_q": exact_q})
    else:
        result.update(
            {
                "actor_exact_gradient_component_rmse": None,
                "actor_exact_gradient_vector_rmse": None,
                "actor_exact_gradient_relative_l2": None,
                "actor_predicted_gradient_norm_mean": None,
                "actor_exact_gradient_norm_mean": None,
            }
        )
    return result, arrays


def run_ground_truth_evaluation(
    spec: ExperimentSpec,
    device: torch.device,
    protocol: dict[str, Any],
    reference: PeriodicReference2D,
) -> None:
    if spec.k != 2:
        return
    output_path = RESULTS_DIR / f"ground_truth_{spec.run_id}.json"
    if output_path.exists():
        print(f"SKIP ground truth {spec.run_id}", flush=True)
        return
    driver = _driver(spec)
    config = ActorConfig(**protocol["actor"])
    actors = _load_actors(spec, driver, config, device)
    exact_actor = ReferenceTangentActor2D(driver, reference).to(device)
    actor_names = ("a0", "a1", "a2", "a3", "exact_tangent")
    all_actors = actors + [exact_actor]
    evaluation = protocol["evaluation"]
    test_t, test_x = fixed_points(
        spec.base_id,
        driver.d,
        "ground_truth_test",
        evaluation["ground_truth_test_points"],
        T=driver.T,
        include_time=False,
        training_seed=protocol["fixed_evaluation_seed"],
    )
    coordinates = test_x @ driver.B[:, :2]
    exact_working = reference.value(test_t, coordinates)
    exact_original = driver.working_to_original(exact_working)
    audit_t, audit_x = fixed_points(
        spec.base_id,
        driver.d,
        "exact_actor_audit",
        evaluation["actor_audit_points"],
        T=driver.T,
        include_time=True,
        training_seed=protocol["fixed_evaluation_seed"],
    )
    rows, saved = [], {
        "test_t": test_t,
        "test_x": test_x,
        "reference_working_value": exact_working,
        "reference_original_value": exact_original,
        "audit_t": audit_t,
        "audit_x": audit_x,
    }
    for actor_name, actor in zip(actor_names, all_actors, strict=True):
        output, details = evaluate_tangent_mc(
            driver,
            actor,
            test_t,
            test_x,
            paths=evaluation["ground_truth_paths"],
            time_steps=evaluation["ground_truth_steps"],
            seed=stable_seed(spec.base_id, "ground_truth_crn"),
            device=device,
            parent_batch=evaluation["evaluation_parent_batch"],
            path_chunk=evaluation["evaluation_path_chunk"],
        )
        prediction = driver.working_to_original(output["working_value"])
        accuracy = model_metrics(prediction, exact_original)
        actor_metrics, actor_arrays = _posthoc_actor_metrics(
            driver, reference, actor, audit_t, audit_x, device
        )
        rows.append(
            {
                "run_id": spec.run_id,
                "family": spec.family,
                "d": spec.d,
                "k": spec.k,
                "training_seed": spec.training_seed,
                "iteration": actor_name,
                **accuracy,
                **actor_metrics,
                "mc_mean_se": details["mean_value_mc_se"],
                "paths": evaluation["ground_truth_paths"],
                "steps": evaluation["ground_truth_steps"],
                "weight_ess_fraction": details["weight_ess_fraction"],
                "weight_min": details["weight_min"],
                "weight_max": details["weight_max"],
                "nonfinite_count": details["nonfinite_count"],
                "inference_seconds": details["wall_seconds"],
            }
        )
        saved[f"{actor_name}_working_value"] = output["working_value"]
        saved[f"{actor_name}_original_value"] = prediction
        saved[f"{actor_name}_mc_variance_of_mean"] = output["mc_variance_of_mean"]
        for key, value in actor_arrays.items():
            saved[f"{actor_name}_{key}"] = value
        print(
            f"{spec.run_id} {actor_name}: RMSE={accuracy['rmse']:.6g} "
            f"r={actor_metrics['actor_exact_value_rmse']:.6g} "
            f"q={actor_metrics['actor_exact_gradient_vector_rmse']}",
            flush=True,
        )

    bregman_t, bregman_x = fixed_points(
        spec.base_id,
        driver.d,
        f"bregman_audit_seed_{spec.training_seed}",
        evaluation["bregman_points"],
        T=driver.T,
        include_time=False,
        training_seed=protocol["fixed_evaluation_seed"],
    )
    bregman, arrays = bregman_identity_audit(
        driver,
        actors[-1],
        reference,
        bregman_t,
        bregman_x,
        paths=evaluation["bregman_paths"],
        time_steps=evaluation["bregman_steps"],
        seed=stable_seed(spec.run_id, "bregman"),
        device=device,
    )
    for key, value in arrays.items():
        saved[f"bregman_{key}"] = value
    prediction_path = PREDICTION_DIR / f"ground_truth_{spec.run_id}.npz"
    atomic_npz(prediction_path, **saved)
    atomic_json(
        output_path,
        {
            "spec": spec.__dict__,
            "rows": rows,
            "bregman": bregman,
            "prediction_path": str(prediction_path.relative_to(driver_path_root())),
            "test_bank_sha256": array_sha256(test_t, test_x),
            "checkpoint_and_reference_free_evaluation_frozen_before_this_step": True,
            "ground_truth_used_for_selection": False,
        },
    )


def environment_payload(stage: str) -> dict[str, Any]:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "unknown"
    return {
        "stage": stage,
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
    parser.add_argument("--stage", choices=("initial", "expansion", "all"), default="initial")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--run-ids", nargs="*")
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required")
    protocol = load_protocol()
    specs = initial_specs() if args.stage == "initial" else expansion_specs()
    if args.stage == "all":
        specs = initial_specs() + expansion_specs()
    if args.run_ids:
        requested = set(args.run_ids)
        specs = [spec for spec in specs if spec.run_id in requested]
        missing = requested - {spec.run_id for spec in specs}
        if missing:
            raise ValueError(f"unknown run IDs: {sorted(missing)}")
    atomic_json(RESULTS_DIR / f"runtime_environment_{args.stage}.json", environment_payload(args.stage))
    for spec in specs:
        run_training_and_reference_free(spec, device, protocol)

    # Ground truth is created only after every selected training run and its
    # reference-free comparisons above have been frozen to disk.
    references: dict[str, PeriodicReference2D] = {}
    validations: dict[str, dict[str, Any]] = {}
    for spec in specs:
        if spec.k != 2:
            continue
        if spec.family not in references:
            reference, validation = validate_periodic_reference_2d(_driver(spec))
            references[spec.family] = reference
            validations[spec.family] = validation
        run_ground_truth_evaluation(spec, device, protocol, references[spec.family])
    if validations:
        existing_path = RESULTS_DIR / "reference_validation.json"
        if existing_path.exists():
            with existing_path.open() as handle:
                existing = json.load(handle)
        else:
            existing = {}
        existing.update(validations)
        atomic_json(existing_path, existing)


if __name__ == "__main__":
    main()
