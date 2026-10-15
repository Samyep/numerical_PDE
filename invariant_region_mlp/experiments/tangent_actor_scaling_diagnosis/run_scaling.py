"""Run the staged high-dimensional Tangent-Actor scaling diagnosis."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import platform
import subprocess
import time
from typing import Any

import numpy as np
import torch
from torch import nn

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    evaluate_tangent_mc,
    make_targets,
    tangent_labels,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
    make_driver,
    terminal_signal_diagnostics,
)
from invariant_region_mlp.experiments.convex_concave_semilinear_extended.models import (
    ConstantTangentActor,
)

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
    nested_points,
    stable_seed,
)
from .models import ActorKind, make_actor
from .training import FitConfig, train_actor


@dataclass(frozen=True)
class Variant:
    name: str
    training_points: int
    width: int
    actor_kind: ActorKind


@dataclass(frozen=True)
class ExperimentSpec:
    family: str
    d: int
    k: int
    training_seed: int
    variant: str

    @property
    def problem_id(self) -> str:
        return f"{self.family}_d{self.d}_k{self.k}"

    @property
    def run_id(self) -> str:
        return f"{self.problem_id}_{self.variant}_s{self.training_seed}"


def load_protocol() -> dict[str, Any]:
    with (RESULTS_DIR / "experiment_protocol.json").open() as handle:
        return json.load(handle)


def variants(protocol: dict[str, Any]) -> dict[str, Variant]:
    return {
        name: Variant(name=name, **values)
        for name, values in protocol["variants"].items()
    }


def fit_config(protocol: dict[str, Any], variant: Variant) -> FitConfig:
    values = dict(protocol["actor_base"])
    values.update(training_points=variant.training_points, width=variant.width)
    return FitConfig(**values)


def stage_specs(
    stage: str, protocol: dict[str, Any], selected_variant: str | None
) -> list[ExperimentSpec]:
    if stage == "stage-a":
        return [
            ExperimentSpec("joint_concave", 400, 5, seed, variant)
            for seed in protocol["training_seeds"]
            for variant in protocol["stage_a_variants"]
        ]
    if selected_variant is None:
        selection_path = RESULTS_DIR / "stage_b_selection.json"
        if not selection_path.exists():
            raise ValueError("Stage B/C needs --variant or stage_b_selection.json")
        with selection_path.open() as handle:
            selected_variant = json.load(handle)["selected_variant"]
    if stage == "stage-b":
        problems = [
            ("joint_concave", 400, 10),
            ("allen_cahn", 400, 5),
        ]
        if protocol["stage_b_optional_allen_k10"]:
            problems.append(("allen_cahn", 400, 10))
        return [
            ExperimentSpec(family, d, k, seed, selected_variant)
            for family, d, k in problems
            for seed in protocol["training_seeds"]
        ]
    if stage == "stage-c":
        return [
            ExperimentSpec("joint_concave", 100, 5, seed, selected_variant)
            for seed in protocol["training_seeds"]
        ]
    raise KeyError(stage)


def _driver(spec: ExperimentSpec) -> ConcaveSemilinearDriver:
    return make_driver(spec.family, d=spec.d, k=spec.k)


def _label_stem(
    spec: ExperimentSpec, iteration: int, generated_training_points: int
) -> str:
    if iteration == 1:
        return (
            f"{spec.problem_id}_s{spec.training_seed}_a1_shared_"
            f"n{generated_training_points}"
        )
    return f"{spec.run_id}_a{iteration}_labels"


def _build_or_load_labels(
    spec: ExperimentSpec,
    driver: ConcaveSemilinearDriver,
    source: nn.Module,
    iteration: int,
    config: FitConfig,
    protocol: dict[str, Any],
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    is_stage_a_problem = (
        spec.family == "joint_concave" and spec.d == 400 and spec.k == 5
    )
    generated_points = (
        protocol["nested_bank_max_points"]
        if iteration == 1 and is_stage_a_problem
        else config.training_points
    )
    stem = _label_stem(spec, iteration, generated_points)
    arrays_path, metadata_path = CACHE_DIR / f"{stem}.npz", CACHE_DIR / f"{stem}.json"
    if arrays_path.exists() and metadata_path.exists():
        with np.load(arrays_path) as archive:
            arrays = {name: archive[name] for name in archive.files}
        with metadata_path.open() as handle:
            metadata = json.load(handle)
        metadata = dict(metadata)
        metadata["loaded_from_cache"] = True
    else:
        train_t, train_x = nested_points(
            spec.problem_id,
            driver.d,
            f"iteration_{iteration}_train",
            generated_points,
            T=driver.T,
            include_time=True,
            training_seed=spec.training_seed,
        )
        valid_t, valid_x = nested_points(
            spec.problem_id,
            driver.d,
            f"iteration_{iteration}_validation",
            config.validation_points,
            T=driver.T,
            include_time=True,
            training_seed=spec.training_seed,
        )
        kwargs = {
            "paths": config.paths_per_label,
            "time_steps": config.label_time_steps,
            "device": device,
            "parent_batch": config.parent_batch,
            "path_chunk": config.path_chunk,
        }
        common_seed = stable_seed(
            "scaling_labels", spec.problem_id, spec.training_seed, iteration
        )
        train_value, train_q, train_meta = tangent_labels(
            driver,
            source,
            train_t,
            train_x,
            seed=common_seed,
            **kwargs,
        )
        valid_value, valid_q, valid_meta = tangent_labels(
            driver,
            source,
            valid_t,
            valid_x,
            seed=stable_seed(common_seed, "validation"),
            **kwargs,
        )
        replicate_value, replicate_q, replicate_meta = tangent_labels(
            driver,
            source,
            valid_t,
            valid_x,
            seed=stable_seed(common_seed, "replication"),
            **kwargs,
        )
        train_r, train_q, train_projection = make_targets(driver, train_value, train_q)
        valid_r, valid_q, valid_projection = make_targets(driver, valid_value, valid_q)
        arrays = {
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
        difference = valid_value - replicate_value
        value_replication_rmse = float(np.sqrt(np.mean(difference**2)))
        if valid_q is None:
            gradient_component_replication_rmse = None
            gradient_vector_replication_rmse = None
        else:
            assert replicate_q is not None
            q_difference = valid_q - replicate_q
            gradient_component_replication_rmse = float(
                np.sqrt(np.mean(q_difference**2))
            )
            gradient_vector_replication_rmse = float(
                np.sqrt(np.mean(np.square(q_difference).sum(axis=1)))
            )
        metadata = {
            "source": f"a{iteration - 1}",
            "iteration": iteration,
            "generated_training_points": generated_points,
            "train": train_meta,
            "validation": valid_meta,
            "replication": replicate_meta,
            "train_projection": train_projection,
            "validation_projection": valid_projection,
            "value_replication_rmse": value_replication_rmse,
            "gradient_component_replication_rmse": gradient_component_replication_rmse,
            "gradient_vector_replication_rmse": gradient_vector_replication_rmse,
            "train_bank_sha256": array_sha256(train_t, train_x),
            "validation_bank_sha256": array_sha256(valid_t, valid_x),
            "label_random_seed": common_seed,
            "analytic_nonterminal_reference_used": False,
            "oracle_subspace_projection_used": False,
            "loaded_from_cache": False,
        }
        atomic_npz(arrays_path, **arrays)
        atomic_json(metadata_path, metadata)
    used = config.training_points
    sliced = {
        name: (value[:used] if name.startswith("train_") else value)
        for name, value in arrays.items()
    }
    metadata["used_training_points"] = used
    metadata["used_train_bank_sha256"] = array_sha256(
        sliced["train_t"], sliced["train_x"]
    )
    metadata["nested_prefix_fraction"] = used / metadata["generated_training_points"]
    return sliced, metadata


def _train_or_load_actor(
    spec: ExperimentSpec,
    driver: ConcaveSemilinearDriver,
    variant: Variant,
    iteration: int,
    source: nn.Module,
    arrays: dict[str, np.ndarray],
    labels: dict[str, Any],
    config: FitConfig,
    device: torch.device,
) -> tuple[nn.Module, dict[str, Any]]:
    stem = f"{spec.run_id}_a{iteration}"
    checkpoint_path = CHECKPOINT_DIR / f"{stem}.pt"
    metadata_path = CHECKPOINT_DIR / f"{stem}.json"
    actor = make_actor(
        variant.actor_kind, driver, width=config.width, depth=config.depth
    ).to(device)
    if checkpoint_path.exists() and metadata_path.exists():
        actor.load_state_dict(
            torch.load(checkpoint_path, map_location=device, weights_only=True)
        )
        actor.eval()
        with metadata_path.open() as handle:
            metadata = json.load(handle)
        metadata["loaded_from_cache"] = True
        return actor, metadata

    initial = (
        source
        if getattr(source, "actor_kind", None) == variant.actor_kind
        else None
    )
    actor, training, trace = train_actor(
        driver,
        config,
        variant.actor_kind,
        arrays["train_t"],
        arrays["train_x"],
        arrays["train_r"],
        arrays.get("train_q"),
        arrays["valid_t"],
        arrays["valid_x"],
        arrays["valid_r"],
        arrays.get("valid_q"),
        seed=stable_seed("scaling_actor", spec.run_id, iteration),
        device=device,
        initial=initial,
    )
    metadata = {
        "run_id": spec.run_id,
        "iteration": iteration,
        "variant": asdict(variant),
        "labels": labels,
        "training": training,
        "trace": trace,
        "checkpoint_selected_by": "held-out noisy tangent-label loss only",
        "ground_truth_used": False,
        "loaded_from_cache": False,
    }
    atomic_torch_save(checkpoint_path, actor.state_dict())
    atomic_json(metadata_path, metadata)
    return actor, metadata


def _phase_row(phase: dict[str, Any] | None) -> dict[str, Any]:
    if phase is None:
        return {
            "selected_step": None,
            "actor_parameters": 0,
            "train_value_rmse": None,
            "validation_value_rmse": None,
            "value_generalization_gap_rmse": None,
            "train_gradient_component_rmse": None,
            "validation_gradient_component_rmse": None,
            "train_gradient_vector_rmse": None,
            "validation_gradient_vector_rmse": None,
            "gradient_generalization_gap_vector_rmse": None,
            "value_label_replication_rmse": None,
            "gradient_label_replication_component_rmse": None,
            "gradient_label_replication_vector_rmse": None,
            "label_weight_ess_fraction": None,
            "label_seconds": 0.0,
            "training_seconds": 0.0,
            "warm_started": None,
        }
    training, labels = phase["training"], phase["labels"]
    label_seconds = sum(
        labels[name]["wall_seconds"] for name in ("train", "validation", "replication")
    )
    if labels["iteration"] == 1:
        label_seconds *= labels["nested_prefix_fraction"]
    return {
        "selected_step": training["selected_step"],
        "actor_parameters": training["parameters"],
        "train_value_rmse": training["selected_train_value_rmse"],
        "validation_value_rmse": training["selected_validation_value_rmse"],
        "value_generalization_gap_rmse": training[
            "selected_value_generalization_gap_rmse"
        ],
        "train_gradient_component_rmse": training[
            "selected_train_gradient_component_rmse"
        ],
        "validation_gradient_component_rmse": training[
            "selected_validation_gradient_component_rmse"
        ],
        "train_gradient_vector_rmse": training[
            "selected_train_gradient_vector_rmse"
        ],
        "validation_gradient_vector_rmse": training[
            "selected_validation_gradient_vector_rmse"
        ],
        "gradient_generalization_gap_vector_rmse": training[
            "selected_gradient_generalization_gap_vector_rmse"
        ],
        "value_label_replication_rmse": labels["value_replication_rmse"],
        "gradient_label_replication_component_rmse": labels[
            "gradient_component_replication_rmse"
        ],
        "gradient_label_replication_vector_rmse": labels[
            "gradient_vector_replication_rmse"
        ],
        "label_weight_ess_fraction": labels["train"]["weight_ess_fraction"],
        "label_seconds": label_seconds,
        "training_seconds": training["wall_seconds"],
        "warm_started": training["warm_started"],
    }


def run_spec(
    spec: ExperimentSpec,
    variant: Variant,
    protocol: dict[str, Any],
    device: torch.device,
) -> None:
    output_path = RESULTS_DIR / f"run_{spec.run_id}.json"
    if output_path.exists():
        print(f"SKIP complete {spec.run_id}", flush=True)
        return
    started = time.perf_counter()
    driver = _driver(spec)
    config = fit_config(protocol, variant)
    actors: list[nn.Module] = [ConstantTangentActor(driver).to(device)]
    phases: list[dict[str, Any]] = []
    for iteration in protocol["iterations"]:
        arrays, labels = _build_or_load_labels(
            spec, driver, actors[-1], iteration, config, protocol, device
        )
        actor, phase = _train_or_load_actor(
            spec,
            driver,
            variant,
            iteration,
            actors[-1],
            arrays,
            labels,
            config,
            device,
        )
        actors.append(actor)
        phases.append(phase)
        print(
            f"{spec.run_id} a{iteration}: checkpoint="
            f"{phase['training']['selected_step']} "
            f"val={phase['training']['best_validation_normalized_loss']:.4g} "
            f"fit={phase['training']['selected_validation_value_rmse']:.4g}",
            flush=True,
        )

    evaluation = protocol["evaluation"]
    test_t, test_x = nested_points(
        spec.problem_id,
        driver.d,
        "policy_test",
        evaluation["points"],
        T=driver.T,
        include_time=False,
        training_seed=evaluation["fixed_test_seed"],
    )
    outputs: list[dict[str, np.ndarray]] = []
    details: list[dict[str, Any]] = []
    for actor in actors:
        output, detail = evaluate_tangent_mc(
            driver,
            actor,
            test_t,
            test_x,
            paths=evaluation["paths"],
            time_steps=evaluation["time_steps"],
            seed=stable_seed("scaling_policy_crn", spec.problem_id),
            device=device,
            parent_batch=evaluation["parent_batch"],
            path_chunk=evaluation["path_chunk"],
        )
        outputs.append(output)
        details.append(detail)

    rows = []
    for index, (output, detail) in enumerate(zip(outputs, details, strict=True)):
        phase = None if index == 0 else phases[index - 1]
        values = output["working_value"]
        previous_delta = (
            np.zeros_like(values)
            if index == 0
            else values - outputs[index - 1]["working_value"]
        )
        initial_delta = values - outputs[0]["working_value"]
        rows.append(
            {
                "run_id": spec.run_id,
                **asdict(spec),
                "actor_kind": variant.actor_kind,
                "training_points": variant.training_points,
                "width": variant.width,
                "iteration": f"a{index}",
                "working_value_mean": float(np.mean(values)),
                "working_value_std_across_states": float(np.std(values, ddof=1)),
                "delta_from_previous_mean": float(np.mean(previous_delta)),
                "delta_from_previous_se_across_states": float(
                    np.std(previous_delta, ddof=1) / np.sqrt(len(previous_delta))
                ),
                "improvement_from_previous_mean": float(-np.mean(previous_delta)),
                "delta_from_a0_mean": float(np.mean(initial_delta)),
                "delta_from_a0_se_across_states": float(
                    np.std(initial_delta, ddof=1) / np.sqrt(len(initial_delta))
                ),
                "improvement_from_a0_mean": float(-np.mean(initial_delta)),
                "mc_mean_se": detail["mean_value_mc_se"],
                "evaluation_weight_ess_fraction": detail["weight_ess_fraction"],
                "evaluation_weight_min": detail["weight_min"],
                "evaluation_weight_max": detail["weight_max"],
                "evaluation_nonfinite_count": detail["nonfinite_count"],
                "inference_seconds": detail["wall_seconds"],
                **_phase_row(phase),
            }
        )

    prediction_path = PREDICTION_DIR / f"{spec.run_id}.npz"
    atomic_npz(
        prediction_path,
        test_t=test_t,
        test_x=test_x,
        **{
            f"a{index}_working_value": output["working_value"]
            for index, output in enumerate(outputs)
        },
        **{
            f"a{index}_mc_variance_of_mean": output["mc_variance_of_mean"]
            for index, output in enumerate(outputs)
        },
    )
    payload = {
        "spec": asdict(spec),
        "variant": asdict(variant),
        "problem": driver.metadata(),
        "terminal_signal": terminal_signal_diagnostics(driver),
        "actor_config": asdict(config),
        "iteration_rows": rows,
        "phases": phases,
        "test_bank_sha256": array_sha256(test_t, test_x),
        "prediction_path": str(prediction_path.relative_to(RESULTS_DIR.parent.parent)),
        "run_wall_seconds": time.perf_counter() - started,
        "selection_used_exact_test": False,
        "analytic_nonterminal_training_data_used": False,
        "projected_actor_input_used": False,
    }
    atomic_json(output_path, payload)
    print(
        f"DONE {spec.run_id}: a3-a0={rows[-1]['delta_from_a0_mean']:+.6g} "
        f"wall={payload['run_wall_seconds']:.1f}s",
        flush=True,
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
    parser.add_argument("--stage", choices=("stage-a", "stage-b", "stage-c"), required=True)
    parser.add_argument("--variant")
    parser.add_argument("--run-ids", nargs="*")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    ensure_dirs()
    protocol = load_protocol()
    all_variants = variants(protocol)
    specs = stage_specs(args.stage, protocol, args.variant)
    if args.run_ids:
        requested = set(args.run_ids)
        specs = [spec for spec in specs if spec.run_id in requested]
        missing = requested - {spec.run_id for spec in specs}
        if missing:
            raise ValueError(f"unknown run IDs: {sorted(missing)}")
    for spec in specs:
        if spec.variant not in all_variants:
            raise KeyError(f"unknown variant {spec.variant}")
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required")
    atomic_json(
        RESULTS_DIR / f"runtime_environment_{args.stage}.json",
        environment_payload(args.stage),
    )
    for spec in specs:
        run_spec(spec, all_variants[spec.variant], protocol, device)


if __name__ == "__main__":
    main()

