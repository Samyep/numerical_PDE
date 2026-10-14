"""Reference-free conditional third update for cases with useful a1->a2 gain."""

from __future__ import annotations

import argparse
import json

import numpy as np
import torch

from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RESULTS_DIR,
    atomic_json,
    atomic_npz,
    fixed_points,
    model_metrics,
    stable_seed,
)
from .core import ActorConfig, evaluate_policy_plain_mc
from .models import ActorNet
from .problems import DirectionalProblem, make_problem
from .run_pilot import (
    _load_or_build_labels,
    _load_or_train_actor,
    _posthoc_action_metrics,
    exact_values,
    make_reference,
)


def run_a3(problem_name: str, device: torch.device) -> None:
    raw_path = RESULTS_DIR / f"raw_{problem_name}.json"
    with raw_path.open() as handle:
        raw = json.load(handle)
    if any(row["iteration"] == "a3" for row in raw["rows"]):
        print(f"SKIP completed {problem_name} a3", flush=True)
        return
    protocol = raw["protocol"]
    config = ActorConfig(**protocol["actor"])
    problem = make_problem(problem_name)
    by_iteration = {row["iteration"]: row for row in raw["policy_iteration_rows"]}
    a2_validation = by_iteration["a2"]
    relative_reduction = -a2_validation["paired_delta_vs_previous"] / abs(
        by_iteration["a1"]["policy_cost"]
    )
    z_score = -a2_validation["paired_delta_vs_previous"] / max(
        a2_validation["paired_delta_se_across_states"], np.finfo(float).tiny
    )
    if not (relative_reduction > 0.001 and z_score > 1.96):
        raise RuntimeError(
            f"{problem_name} does not meet reference-free a3 criterion: "
            f"relative={relative_reduction}, z={z_score}"
        )

    a2 = ActorNet(problem, width=config.width, depth=config.depth).to(device)
    a2.load_state_dict(
        torch.load(CHECKPOINT_DIR / f"{problem.name}_a2.pt", map_location=device, weights_only=True)
    )
    a2.eval()
    arrays, label_metadata = _load_or_build_labels(problem, a2, 3, config, device)
    a3, phase = _load_or_train_actor(
        problem, 3, a2, arrays, label_metadata, config, device
    )
    evaluation = protocol["evaluation"]
    validation_t, validation_x = fixed_points(
        problem.name, problem.d, "policy_valid", evaluation["policy_validation_points"],
        T=problem.T, include_time=False,
    )
    validation = []
    validation_meta = []
    for policy in (a2, a3):
        output, metadata = evaluate_policy_plain_mc(
            problem, policy, validation_t, validation_x,
            paths=evaluation["policy_validation_paths"],
            time_steps=evaluation["policy_validation_steps"],
            seed=stable_seed(problem.name, "policy_validation_crn"),
            device=device,
        )
        validation.append(output)
        validation_meta.append(metadata)
    difference = validation[1]["value"] - validation[0]["value"]
    validation_row = {
        "pde": problem.name,
        "iteration": "a3",
        "policy_cost": float(np.mean(validation[1]["value"])),
        "policy_cost_mc_se": validation_meta[1]["mean_value_mc_se"],
        "paired_delta_vs_previous": float(np.mean(difference)),
        "paired_delta_se_across_states": float(np.std(difference, ddof=1) / np.sqrt(len(difference))),
        "gradient_replication_rmse": phase["label"]["gradient_replication_rmse"],
        "gradient_replication_relative_l2": phase["label"]["gradient_replication_relative_l2"],
        "actor_validation_action_mse": phase["training"]["best_validation_action_mse"],
        "label_seconds": (
            phase["label"]["train"]["wall_seconds"]
            + phase["label"]["valid"]["wall_seconds"]
            + phase["label"]["replication"]["wall_seconds"]
        ),
        "training_seconds": phase["training"]["wall_seconds"],
    }

    # Only after checkpointing and reference-free validation is the reference queried.
    reference, _ = make_reference(problem)
    test_t, test_x = fixed_points(
        problem.name, problem.d, "test", evaluation["test_points"],
        T=problem.T, include_time=False,
    )
    truth = exact_values(problem, reference, test_t, test_x)
    output, mc_metadata = evaluate_policy_plain_mc(
        problem, a3, test_t, test_x,
        paths=evaluation["test_paths"],
        time_steps=evaluation["test_steps"],
        seed=stable_seed(problem.name, "test_crn"),
        device=device,
    )
    audit_t, audit_x = fixed_points(
        problem.name, problem.d, "action_audit", 1024, T=problem.T, include_time=True
    )
    action_metrics, action, _ = _posthoc_action_metrics(
        problem, reference, a3, audit_t, audit_x, device
    )
    metrics = model_metrics(output["value"], truth)
    phase_gpu = (
        phase["label"]["train"]["wall_seconds"]
        + phase["label"]["valid"]["wall_seconds"]
        + phase["label"]["replication"]["wall_seconds"]
        + phase["training"]["wall_seconds"]
        + mc_metadata["wall_seconds"]
    )
    row = {
        "pde": problem.name,
        "d": problem.d,
        "m": problem.m,
        "iteration": "a3",
        **metrics,
        **action_metrics,
        "value_mean_mc_se": mc_metadata["mean_value_mc_se"],
        "validation_policy_cost": validation_row["policy_cost"],
        "validation_policy_cost_mc_se": validation_row["policy_cost_mc_se"],
        "mean_running_cost": float(np.mean(output["running_cost"])),
        "mean_abs_action_along_paths": float(np.mean(output["mean_abs_action"])),
        "path_saturation_rate": float(np.mean(output["saturation_rate"])),
        "gpu_seconds": phase_gpu,
        "final_mc_seconds": mc_metadata["wall_seconds"],
        "final_mc_paths": evaluation["test_paths"],
        "final_mc_steps": evaluation["test_steps"],
    }
    with np.load(PREDICTION_DIR / f"{problem.name}.npz") as saved:
        prediction_arrays = {key: saved[key] for key in saved.files}
    prediction_arrays.update({
        "a3_value": output["value"],
        "a3_mc_variance_of_mean": output["mc_variance_of_mean"],
        "a3_audit_action": action,
        "a3_validation_value": validation[1]["value"],
        "a3_validation_mc_variance_of_mean": validation[1]["mc_variance_of_mean"],
    })
    if isinstance(problem, DirectionalProblem):
        control_s = prediction_arrays["control_plot_s"]
        control_t = prediction_arrays["control_plot_t"]
        control_x = control_s[:, None] * problem.beta[None, :]
        from .core import actor_actions
        prediction_arrays["control_plot_a3"] = actor_actions(
            problem, a3, control_t, control_x, device
        )[:, 0]
    atomic_npz(PREDICTION_DIR / f"{problem.name}.npz", **prediction_arrays)
    raw["rows"].append(row)
    raw["policy_iteration_rows"].append(validation_row)
    raw["phase_metadata"].append(phase)
    raw["conditional_a3_decision"] = {
        "criterion": "a2 paired validation z > 1.96 and relative cost reduction > 0.1%",
        "a2_relative_validation_reduction": relative_reduction,
        "a2_paired_z": z_score,
        "exact_reference_used": False,
    }
    atomic_json(raw_path, raw)
    print(
        f"{problem.name} a3: delta={validation_row['paired_delta_vs_previous']:.6g}, "
        f"RMSE={row['rmse']:.6g}, action={row['action_rmse']:.6g}",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--problems", nargs="+", default=["lqg", "quartic"])
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required")
    for problem in args.problems:
        run_a3(problem, torch.device("cuda"))


if __name__ == "__main__":
    main()

