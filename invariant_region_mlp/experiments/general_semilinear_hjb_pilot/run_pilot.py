"""Execute the three-problem, two-update general semilinear HJB pilot."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import subprocess
import time
from typing import Any

import numpy as np
import torch

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
    actor_actions,
    evaluate_policy_plain_mc,
    greedy_targets,
    policy_gradient_labels,
    train_actor,
)
from .models import ActorNet, ZeroPolicy
from .problems import BoundedControlProblem, ControlledHJB, DirectionalProblem, LQGProblem, make_problem
from .reference import RiccatiReference, ScalarReference, solve_scalar_reference
from .reference_policies import LQGReferencePolicy, ScalarReferencePolicy


def load_protocol() -> dict[str, Any]:
    with (RESULTS_DIR / "experiment_protocol.json").open() as handle:
        return json.load(handle)


def actor_config(protocol: dict[str, Any]) -> ActorConfig:
    return ActorConfig(**protocol["actor"])


def make_reference(
    problem: ControlledHJB,
) -> tuple[RiccatiReference | ScalarReference, torch.nn.Module]:
    if isinstance(problem, LQGProblem):
        reference = RiccatiReference(problem)
        return reference, LQGReferencePolicy(problem, reference)
    if not isinstance(problem, DirectionalProblem):
        raise TypeError(problem)
    reference = solve_scalar_reference(problem, dx=0.025, dt=0.0005)
    return reference, ScalarReferencePolicy(problem, reference)


def exact_values(
    problem: ControlledHJB,
    reference: RiccatiReference | ScalarReference,
    t: np.ndarray,
    x: np.ndarray,
) -> np.ndarray:
    if isinstance(reference, RiccatiReference):
        return reference.value(t, x)
    assert isinstance(problem, DirectionalProblem)
    return reference.value(t, x @ problem.beta)


def exact_actions(
    problem: ControlledHJB,
    reference: RiccatiReference | ScalarReference,
    t: np.ndarray,
    x: np.ndarray,
) -> np.ndarray:
    if isinstance(reference, RiccatiReference):
        return reference.action(t, x)
    assert isinstance(problem, DirectionalProblem)
    return reference.action(t, x @ problem.beta)[:, None]


def _load_or_build_labels(
    problem: ControlledHJB,
    source: torch.nn.Module,
    iteration: int,
    config: ActorConfig,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    cache_path = CACHE_DIR / f"{problem.name}_a{iteration}_labels.npz"
    metadata_path = CACHE_DIR / f"{problem.name}_a{iteration}_labels.json"
    if cache_path.exists() and metadata_path.exists():
        with np.load(cache_path) as data:
            arrays = {key: data[key] for key in data.files}
        with metadata_path.open() as handle:
            return arrays, json.load(handle)

    train_t, train_x = fixed_points(
        problem.name, problem.d, f"iteration_{iteration}_train", config.training_points,
        T=problem.T, include_time=True,
    )
    valid_t, valid_x = fixed_points(
        problem.name, problem.d, f"iteration_{iteration}_label_valid", config.validation_points,
        T=problem.T, include_time=True,
    )
    train_value, train_gradient, train_meta = policy_gradient_labels(
        problem, source, train_t, train_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(problem.name, iteration, "train_labels"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    valid_value, valid_gradient, valid_meta = policy_gradient_labels(
        problem, source, valid_t, valid_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(problem.name, iteration, "valid_labels"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    _, replicate_gradient, replicate_meta = policy_gradient_labels(
        problem, source, valid_t, valid_x,
        paths=config.paths_per_label,
        time_steps=config.label_time_steps,
        seed=stable_seed(problem.name, iteration, "valid_labels_replication"),
        device=device,
        parent_batch=config.parent_batch,
        path_chunk=config.path_chunk,
    )
    train_target = greedy_targets(problem, train_t, train_x, train_gradient, device)
    valid_target = greedy_targets(problem, valid_t, valid_x, valid_gradient, device)
    arrays = {
        "train_t": train_t,
        "train_x": train_x,
        "train_value": train_value,
        "train_gradient": train_gradient,
        "train_target": train_target,
        "valid_t": valid_t,
        "valid_x": valid_x,
        "valid_value": valid_value,
        "valid_gradient": valid_gradient,
        "valid_target": valid_target,
    }
    metadata = {
        "train": train_meta,
        "valid": valid_meta,
        "replication": replicate_meta,
        "gradient_replication_rmse": float(np.sqrt(np.mean((valid_gradient - replicate_gradient) ** 2))),
        "gradient_replication_relative_l2": float(
            np.linalg.norm(valid_gradient - replicate_gradient)
            / max(np.linalg.norm(0.5 * (valid_gradient + replicate_gradient)), np.finfo(float).tiny)
        ),
        "train_bank_sha256": array_sha256(train_t, train_x),
        "valid_bank_sha256": array_sha256(valid_t, valid_x),
        "reference_information_used": False,
    }
    atomic_npz(cache_path, **arrays)
    atomic_json(metadata_path, metadata)
    return arrays, metadata


def _load_or_train_actor(
    problem: ControlledHJB,
    iteration: int,
    initial: ActorNet | None,
    arrays: dict[str, np.ndarray],
    label_metadata: dict[str, Any],
    config: ActorConfig,
    device: torch.device,
) -> tuple[ActorNet, dict[str, Any]]:
    checkpoint_path = CHECKPOINT_DIR / f"{problem.name}_a{iteration}.pt"
    metadata_path = CHECKPOINT_DIR / f"{problem.name}_a{iteration}.json"
    actor = ActorNet(problem, width=config.width, depth=config.depth).to(device)
    if checkpoint_path.exists() and metadata_path.exists():
        actor.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
        actor.eval()
        with metadata_path.open() as handle:
            return actor, json.load(handle)
    actor, train_metadata, trace = train_actor(
        problem,
        config,
        arrays["train_t"], arrays["train_x"], arrays["train_target"],
        arrays["valid_t"], arrays["valid_x"], arrays["valid_target"],
        seed=stable_seed(problem.name, iteration, "actor_training"),
        device=device,
        initial=initial,
    )
    metadata = {
        "problem": problem.name,
        "iteration": iteration,
        "source": f"a{iteration - 1}",
        "label": label_metadata,
        "training": train_metadata,
        "trace": trace,
        "checkpoint_selected_without_exact_reference": True,
    }
    atomic_torch_save(checkpoint_path, actor.state_dict())
    atomic_json(metadata_path, metadata)
    return actor, metadata


def _posthoc_action_metrics(
    problem: ControlledHJB,
    reference: RiccatiReference | ScalarReference,
    policy: torch.nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    prediction = actor_actions(problem, policy, t, x, device)
    truth = exact_actions(problem, reference, t, x)
    difference = prediction - truth
    metrics = {
        "action_rmse": float(np.sqrt(np.mean(difference**2))),
        "action_relative_l2": float(np.linalg.norm(difference) / max(np.linalg.norm(truth), np.finfo(float).tiny)),
        "max_abs_action": float(np.max(np.abs(prediction))),
        "feasibility_violation_rate": 0.0,
    }
    if hasattr(problem, "A"):
        metrics["feasibility_violation_rate"] = float(np.mean(np.abs(prediction) > problem.A + 1e-7))
        metrics["actor_saturation_rate"] = float(np.mean(np.abs(prediction) >= 0.99 * problem.A))
        metrics["optimal_saturation_rate"] = float(np.mean(np.abs(truth) >= 0.99 * problem.A))
    if isinstance(reference, ScalarReference):
        assert isinstance(problem, DirectionalProblem)
        p = reference.gradient(t, x @ problem.beta)
        optimal_h = problem.hamiltonian_numpy(truth[:, 0], p)
        actor_h = problem.hamiltonian_numpy(prediction[:, 0], p)
        residual = actor_h - optimal_h
        metrics.update(
            hamiltonian_residual_mean=float(np.mean(residual)),
            hamiltonian_residual_rmse=float(np.sqrt(np.mean(residual**2))),
            hamiltonian_residual_max=float(np.max(residual)),
            hamiltonian_negative_fraction=float(np.mean(residual < -1e-6)),
        )
    return metrics, prediction, truth


def run_problem(problem_name: str, device: torch.device, protocol: dict[str, Any]) -> Path:
    ensure_dirs()
    output_path = RESULTS_DIR / f"raw_{problem_name}.json"
    if output_path.exists():
        print(f"SKIP completed {problem_name}: {output_path}", flush=True)
        return output_path
    problem = make_problem(problem_name)
    config = actor_config(protocol)
    policies: list[torch.nn.Module] = [ZeroPolicy(problem).to(device)]
    phase_metadata: list[dict[str, Any]] = []
    for iteration in (1, 2):
        arrays, label_metadata = _load_or_build_labels(
            problem, policies[-1], iteration, config, device
        )
        initial = policies[-1] if isinstance(policies[-1], ActorNet) else None
        actor, metadata = _load_or_train_actor(
            problem, iteration, initial, arrays, label_metadata, config, device
        )
        policies.append(actor)
        phase_metadata.append(metadata)
        print(
            f"{problem.name} a{iteration}: label={label_metadata['train']['wall_seconds']:.1f}s "
            f"train={metadata['training']['wall_seconds']:.1f}s "
            f"valid_action_mse={metadata['training']['best_validation_action_mse']:.6g}",
            flush=True,
        )

    evaluation = protocol["evaluation"]
    validation_t, validation_x = fixed_points(
        problem.name, problem.d, "policy_valid", evaluation["policy_validation_points"],
        T=problem.T, include_time=False,
    )
    validation_outputs: list[dict[str, np.ndarray]] = []
    validation_metadata: list[dict[str, Any]] = []
    for policy in policies:
        output, metadata = evaluate_policy_plain_mc(
            problem, policy, validation_t, validation_x,
            paths=evaluation["policy_validation_paths"],
            time_steps=evaluation["policy_validation_steps"],
            seed=stable_seed(problem.name, "policy_validation_crn"),
            device=device,
        )
        validation_outputs.append(output)
        validation_metadata.append(metadata)

    # Exact/numerical references are intentionally constructed only after all
    # actor checkpoints and reference-free policy validation values are frozen.
    reference, reference_policy = make_reference(problem)
    reference_policy = reference_policy.to(device)
    all_policies = policies + [reference_policy]
    all_names = ["a0", "a1", "a2", "optimal_reference_policy"]
    test_t, test_x = fixed_points(
        problem.name, problem.d, "test", evaluation["test_points"],
        T=problem.T, include_time=False,
    )
    truth = exact_values(problem, reference, test_t, test_x)
    audit_t, audit_x = fixed_points(
        problem.name, problem.d, "action_audit", 1024, T=problem.T, include_time=True
    )
    rows: list[dict[str, Any]] = []
    prediction_arrays: dict[str, np.ndarray] = {
        "test_t": test_t,
        "test_x": test_x,
        "reference_value": truth,
        "audit_t": audit_t,
        "audit_x": audit_x,
    }
    total_test_seconds = 0.0
    for index, (name, policy) in enumerate(zip(all_names, all_policies, strict=True)):
        output, mc_metadata = evaluate_policy_plain_mc(
            problem, policy, test_t, test_x,
            paths=evaluation["test_paths"],
            time_steps=evaluation["test_steps"],
            seed=stable_seed(problem.name, "test_crn"),
            device=device,
        )
        total_test_seconds += mc_metadata["wall_seconds"]
        metrics = model_metrics(output["value"], truth)
        action_metrics, action, optimal_action = _posthoc_action_metrics(
            problem, reference, policy, audit_t, audit_x, device
        )
        phase_gpu = mc_metadata["wall_seconds"]
        if 0 < index <= 2:
            phase = phase_metadata[index - 1]
            phase_gpu += phase["label"]["train"]["wall_seconds"]
            phase_gpu += phase["label"]["valid"]["wall_seconds"]
            phase_gpu += phase["label"]["replication"]["wall_seconds"]
            phase_gpu += phase["training"]["wall_seconds"]
        row = {
            "pde": problem.name,
            "d": problem.d,
            "m": problem.m,
            "iteration": name,
            **metrics,
            **action_metrics,
            "value_mean_mc_se": mc_metadata["mean_value_mc_se"],
            "validation_policy_cost": (
                float(np.mean(validation_outputs[index]["value"])) if index < 3 else None
            ),
            "validation_policy_cost_mc_se": (
                validation_metadata[index]["mean_value_mc_se"] if index < 3 else None
            ),
            "mean_running_cost": float(np.mean(output["running_cost"])),
            "mean_abs_action_along_paths": float(np.mean(output["mean_abs_action"])),
            "path_saturation_rate": float(np.mean(output["saturation_rate"])),
            "gpu_seconds": phase_gpu,
            "final_mc_seconds": mc_metadata["wall_seconds"],
            "final_mc_paths": evaluation["test_paths"],
            "final_mc_steps": evaluation["test_steps"],
        }
        rows.append(row)
        prediction_arrays[f"{name}_value"] = output["value"]
        prediction_arrays[f"{name}_mc_variance_of_mean"] = output["mc_variance_of_mean"]
        prediction_arrays[f"{name}_audit_action"] = action
        if index == 0:
            prediction_arrays["optimal_audit_action"] = optimal_action
        print(
            f"{problem.name} {name}: RMSE={metrics['rmse']:.6g}, "
            f"regret={metrics['policy_regret']:.6g}, action={action_metrics['action_rmse']:.6g}",
            flush=True,
        )

    validation_rows = []
    for index, name in enumerate(("a0", "a1", "a2")):
        value = validation_outputs[index]["value"]
        previous = validation_outputs[index - 1]["value"] if index else None
        difference = value - previous if previous is not None else np.full_like(value, np.nan)
        phase = phase_metadata[index - 1] if index else None
        validation_rows.append({
            "pde": problem.name,
            "iteration": name,
            "policy_cost": float(np.mean(value)),
            "policy_cost_mc_se": validation_metadata[index]["mean_value_mc_se"],
            "paired_delta_vs_previous": float(np.mean(difference)) if index else None,
            "paired_delta_se_across_states": float(np.std(difference, ddof=1) / np.sqrt(len(difference))) if index else None,
            "gradient_replication_rmse": phase["label"]["gradient_replication_rmse"] if phase else None,
            "gradient_replication_relative_l2": phase["label"]["gradient_replication_relative_l2"] if phase else None,
            "actor_validation_action_mse": phase["training"]["best_validation_action_mse"] if phase else None,
            "label_seconds": (
                phase["label"]["train"]["wall_seconds"] + phase["label"]["valid"]["wall_seconds"]
                + phase["label"]["replication"]["wall_seconds"] if phase else 0.0
            ),
            "training_seconds": phase["training"]["wall_seconds"] if phase else 0.0,
        })
        prediction_arrays[f"{name}_validation_value"] = value
        prediction_arrays[f"{name}_validation_mc_variance_of_mean"] = validation_outputs[index]["mc_variance_of_mean"]

    if isinstance(problem, DirectionalProblem):
        control_s = np.linspace(-4.0, 4.0, 401)
        control_t = np.concatenate([np.full(401, value) for value in (0.0, 0.25, 0.5)])
        control_s_tiled = np.tile(control_s, 3)
        control_x = control_s_tiled[:, None] * problem.beta[None, :]
        prediction_arrays["control_plot_s"] = control_s_tiled
        prediction_arrays["control_plot_t"] = control_t
        prediction_arrays["control_plot_optimal"] = exact_actions(
            problem, reference, control_t, control_x
        )[:, 0]
        for name, policy in zip(all_names[:3], policies, strict=True):
            prediction_arrays[f"control_plot_{name}"] = actor_actions(
                problem, policy, control_t, control_x, device
            )[:, 0]

    prediction_path = PREDICTION_DIR / f"{problem.name}.npz"
    atomic_npz(prediction_path, **prediction_arrays)
    raw = {
        "problem": problem.metadata(),
        "protocol": protocol,
        "rows": rows,
        "policy_iteration_rows": validation_rows,
        "phase_metadata": phase_metadata,
        "prediction_path": str(prediction_path.relative_to(RESULTS_DIR.parent.parent)),
        "test_bank_sha256": array_sha256(test_t, test_x),
        "test_evaluation_total_seconds": total_test_seconds,
        "reference_created_after_checkpoint_freeze": True,
        "exact_reference_used_for_selection": False,
    }
    atomic_json(output_path, raw)
    return output_path


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
    parser.add_argument("--problems", nargs="+", default=["lqg", "quartic", "bounded"])
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required for the GPU pilot")
    protocol = load_protocol()
    atomic_json(RESULTS_DIR / "runtime_environment.json", environment_payload())
    for name in args.problems:
        run_problem(name, device, protocol)


if __name__ == "__main__":
    main()

