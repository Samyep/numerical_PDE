"""Run the pre-authorized a3 extension only after reference-free a2 improvement."""

from __future__ import annotations

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
from .core import ActorConfig, evaluate_tangent_mc
from .diagnostics import bregman_identity_audit
from .drivers import make_driver
from .models import TangentActor
from .reference import ReferenceTangentActor, solve_periodic_reference
from .run_pilot import (
    PROBLEMS,
    _load_or_build_labels,
    _load_or_train_actor,
    _posthoc_tangent_metrics,
    load_protocol,
)


def _load_actor(name: str, iteration: int, config: ActorConfig, device: torch.device) -> TangentActor:
    driver = make_driver(name)
    actor = TangentActor(driver, width=config.width, depth=config.depth).to(device)
    actor.load_state_dict(
        torch.load(
            CHECKPOINT_DIR / f"{name}_a{iteration}.pt", map_location=device, weights_only=True
        )
    )
    actor.eval()
    return actor


def run_extension(name: str, device: torch.device, protocol: dict) -> None:
    raw_path = RESULTS_DIR / f"raw_{name}.json"
    with raw_path.open() as handle:
        raw = json.load(handle)
    if raw.get("a3_run"):
        print(f"SKIP completed {name} a3", flush=True)
        return
    a2_iteration = next(
        row for row in raw["tangent_iteration_rows"] if row["iteration"] == "a2"
    )
    upper_95 = (
        a2_iteration["paired_delta_vs_previous"]
        + 1.96 * a2_iteration["paired_delta_se_across_states"]
    )
    if upper_95 >= 0.0:
        raw["a3_run"] = False
        raw["a3_reason"] = f"a2 reference-free paired 95% upper bound {upper_95} was not negative"
        atomic_json(raw_path, raw)
        print(f"SKIP ineligible {name}: upper95={upper_95:.6g}", flush=True)
        return

    driver = make_driver(name, d=int(protocol["dimension"]))
    config = ActorConfig(**protocol["actor"])
    a2 = _load_actor(name, 2, config, device)
    arrays, label_metadata = _load_or_build_labels(driver, a2, 3, config, device)
    a3, phase = _load_or_train_actor(
        driver, 3, a2, arrays, label_metadata, config, device
    )

    evaluation = protocol["evaluation"]
    validation_t, validation_x = fixed_points(
        name, driver.d, "policy_valid", evaluation["policy_validation_points"],
        T=driver.T, include_time=False,
    )
    validation_output, validation_metadata = evaluate_tangent_mc(
        driver, a3, validation_t, validation_x,
        paths=evaluation["policy_validation_paths"],
        time_steps=evaluation["policy_validation_steps"],
        seed=stable_seed(name, "policy_validation_crn"),
        device=device,
    )
    prediction_path = PREDICTION_DIR / f"{name}.npz"
    with np.load(prediction_path) as stored:
        saved = {key: stored[key] for key in stored.files}
    previous_validation = saved["a2_validation_working_value"]
    validation_value = validation_output["working_value"]
    delta = validation_value - previous_validation
    iteration_row = {
        "pde": name,
        "iteration": "a3",
        "reference_free_weighted_value": float(np.mean(validation_value)),
        "reference_free_value_mc_se": validation_metadata["mean_value_mc_se"],
        "paired_delta_vs_previous": float(np.mean(delta)),
        "paired_delta_se_across_states": float(np.std(delta, ddof=1) / np.sqrt(len(delta))),
        "value_label_replication_rmse": phase["label"]["value_replication_rmse"],
        "gradient_label_replication_rmse": phase["label"]["gradient_replication_rmse"],
        "actor_validation_normalized_loss": phase["training"]["best_validation_normalized_loss"],
        "actor_validation_value_mse": phase["training"]["best_validation_value_mse"],
        "actor_validation_gradient_mse": phase["training"]["best_validation_gradient_mse"],
        "selected_optimizer_step": phase["training"]["selected_step"],
        "label_seconds": (
            phase["label"]["train"]["wall_seconds"]
            + phase["label"]["valid"]["wall_seconds"]
            + phase["label"]["replication"]["wall_seconds"]
        ),
        "training_seconds": phase["training"]["wall_seconds"],
        "max_label_weight": phase["label"]["train"]["weight_max"],
        "label_weight_ess_fraction": phase["label"]["train"]["weight_ess_fraction"],
        "label_gradient_norm_mean": phase["label"]["train"]["gradient_norm_mean"],
        "label_nonfinite_count": phase["label"]["train"]["nonfinite_count"],
    }
    saved["a3_validation_working_value"] = validation_value

    # Construct the post-hoc reference only after the a3 checkpoint and its
    # reference-free comparison are fixed.
    reference = solve_periodic_reference(driver)
    exact_actor = ReferenceTangentActor(driver, reference).to(device)
    test_t, test_x = saved["test_t"], saved["test_x"]
    exact_original = saved["reference_original_value"]
    audit_t, audit_x = saved["audit_t"], saved["audit_x"]
    main_variant = evaluation["mc_floor_variants"][0]
    main_output, main_metadata = evaluate_tangent_mc(
        driver, a3, test_t, test_x,
        paths=main_variant["paths"], time_steps=main_variant["steps"],
        seed=stable_seed(name, "test_crn", main_variant["name"]), device=device,
    )
    main_prediction = driver.working_to_original(main_output["working_value"])
    metrics = model_metrics(main_prediction, exact_original)
    tangent_metrics, predicted_r, predicted_q, _, _ = _posthoc_tangent_metrics(
        driver, reference, a3, audit_t, audit_x, device
    )
    phase_seconds = (
        phase["label"]["train"]["wall_seconds"]
        + phase["label"]["valid"]["wall_seconds"]
        + phase["label"]["replication"]["wall_seconds"]
        + phase["training"]["wall_seconds"]
        + main_metadata["wall_seconds"]
    )
    result_row = {
        "pde": name,
        "d": driver.d,
        "iteration": "a3",
        **metrics,
        **tangent_metrics,
        "final_mc_mean_se": main_metadata["mean_value_mc_se"],
        "final_mc_paths": main_variant["paths"],
        "final_mc_steps": main_variant["steps"],
        "inference_seconds": main_metadata["wall_seconds"],
        "phase_gpu_seconds": phase_seconds,
        "c_min": main_metadata["c_min"],
        "c_max": main_metadata["c_max"],
        "c_mean": main_metadata["c_mean"],
        "c_std": main_metadata["c_std"],
        "terminal_weight_min": main_metadata["weight_min"],
        "terminal_weight_max": main_metadata["weight_max"],
        "terminal_weight_mean": main_metadata["weight_mean"],
        "terminal_weight_std": main_metadata["weight_std"],
        "weight_ess_fraction": main_metadata["weight_ess_fraction"],
        "nonfinite_count": main_metadata["nonfinite_count"],
    }
    saved["a3_working_value"] = main_output["working_value"]
    saved["a3_original_value"] = main_prediction
    saved["a3_mc_variance_of_mean"] = main_output["mc_variance_of_mean"]
    saved["a3_audit_r"] = predicted_r
    if predicted_q is not None:
        saved["a3_audit_q"] = predicted_q

    floor_rows = []
    for variant in evaluation["mc_floor_variants"]:
        if variant["name"] == main_variant["name"]:
            output, metadata = main_output, main_metadata
        else:
            output, metadata = evaluate_tangent_mc(
                driver, a3, test_t, test_x,
                paths=variant["paths"], time_steps=variant["steps"],
                seed=stable_seed(name, "test_crn", variant["name"]), device=device,
            )
        prediction = driver.working_to_original(output["working_value"])
        floor_rows.append({
            "pde": name,
            "actor": "a3",
            "variant": variant["name"],
            "paths": variant["paths"],
            "steps": variant["steps"],
            **model_metrics(prediction, exact_original),
            "mean_mc_se": metadata["mean_value_mc_se"],
            "wall_seconds": metadata["wall_seconds"],
            "terminal_weight_max": metadata["weight_max"],
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "nonfinite_count": metadata["nonfinite_count"],
        })
        saved[f"floor_a3_{variant['name']}_original_value"] = prediction
        saved[f"floor_a3_{variant['name']}_mc_variance_of_mean"] = output[
            "mc_variance_of_mean"
        ]

    bregman_t, bregman_x = fixed_points(
        name, driver.d, "bregman_audit", evaluation["bregman_points"],
        T=driver.T, include_time=False,
    )
    bregman_summary, bregman_arrays = bregman_identity_audit(
        driver, a3, reference, bregman_t, bregman_x,
        paths=evaluation["bregman_paths"], time_steps=evaluation["bregman_steps"],
        seed=stable_seed(name, "bregman_audit"), device=device,
    )
    for key, value in bregman_arrays.items():
        saved[f"bregman_a3_{key}"] = value
    atomic_npz(prediction_path, **saved)

    raw["rows"] = [row for row in raw["rows"] if row["iteration"] != "a3"] + [result_row]
    raw["tangent_iteration_rows"] = [
        row for row in raw["tangent_iteration_rows"] if row["iteration"] != "a3"
    ] + [iteration_row]
    raw["mc_floor_rows"] = [
        row for row in raw["mc_floor_rows"] if row["actor"] != "a3"
    ] + floor_rows
    raw["phase_metadata"] = [
        item for item in raw["phase_metadata"] if item["iteration"] != 3
    ] + [phase]
    raw["bregman_a2"] = raw["bregman"]
    raw["bregman"] = bregman_summary
    raw["bregman_actor"] = "a3"
    raw["a3_run"] = True
    raw["a3_reason"] = (
        "a2 reference-free paired improvement had a negative 95% upper bound: "
        f"{upper_95:.8g}"
    )
    raw["a3_exact_reference_used_for_eligibility_or_selection"] = False
    atomic_json(raw_path, raw)
    print(
        f"{name} a3: delta={iteration_row['paired_delta_vs_previous']:.6g} "
        f"se={iteration_row['paired_delta_se_across_states']:.6g} "
        f"RMSE={metrics['rmse']:.6g}", flush=True,
    )


def main() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required")
    protocol = load_protocol()
    device = torch.device("cuda")
    extension = {
        "stage": "conditional_a3",
        "eligibility": "a2 paired_delta + 1.96 * paired_SE < 0 using reference-free validation",
        "eligible_problems": [],
        "test_ground_truth_used_for_eligibility": False,
    }
    for name in PROBLEMS:
        with (RESULTS_DIR / f"raw_{name}.json").open() as handle:
            raw = json.load(handle)
        row = next(item for item in raw["tangent_iteration_rows"] if item["iteration"] == "a2")
        if row["paired_delta_vs_previous"] + 1.96 * row["paired_delta_se_across_states"] < 0:
            extension["eligible_problems"].append(name)
    atomic_json(RESULTS_DIR / "conditional_a3_protocol.json", extension)
    for name in PROBLEMS:
        run_extension(name, device, protocol)


if __name__ == "__main__":
    main()
