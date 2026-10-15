"""Posthoc common-bank evaluation and MC/systematic-error decomposition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    evaluate_tangent_actor,
    load_tangent_actor,
)

from .actor_scaling import BRANCHES, CONTINUATION_BRANCHES
from .common import (
    CHECKPOINT_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    STRONG_RESULTS_DIR,
    atomic_json,
    driver_for,
    ensure_dirs,
    evaluation_data,
    metrics,
    stable_seed,
)
from .original_baselines import (
    load_deep_bsde_original,
    load_sdgd_pinn_original,
    predict_deep_bsde_original,
    predict_sdgd_pinn_original,
)
from .scasml_original import predict_scasml_original_level1


def _json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def _actor_candidates(benchmark: str, seed: int) -> list[dict[str, Any]]:
    base_raw = _json(
        STRONG_RESULTS_DIR / "raw" / f"tangent_actor_{benchmark}_s{seed}_main.json"
    )
    base = base_raw["phases"][2]
    candidates: list[dict[str, Any]] = [
        {
            "evaluation_id": "frozen_original__a3_frozen",
            "checkpoint": "a3_frozen", "branch": "frozen_original",
            "path": base["checkpoint"],
            "training_seconds": base["cumulative_training_seconds"],
            "validation_objective": base["training"]["best_validation_normalized_loss"],
            "phase": base,
        }
    ]
    for branch in BRANCHES:
        path = RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json"
        result = _json(path)
        for phase in result["phases"]:
            candidates.append(
                {
                    "evaluation_id": f"{branch}__{phase['checkpoint_label']}",
                    "checkpoint": phase["checkpoint_label"], "branch": branch,
                    "path": phase["checkpoint"],
                    "training_seconds": phase["cumulative_training_seconds"],
                    "validation_objective": phase["training"]["best_validation_normalized_loss"],
                    "phase": phase,
                }
            )
    for branch in CONTINUATION_BRANCHES:
        path = RAW_DIR / f"tangent_actor_high_budget_{benchmark}_s{seed}_{branch}.json"
        result = _json(path)
        for phase in result["phases"]:
            candidates.append(
                {
                    "evaluation_id": f"{branch}__{phase['checkpoint_label']}",
                    "checkpoint": phase["checkpoint_label"], "branch": branch,
                    "path": phase["checkpoint"],
                    "training_seconds": phase["cumulative_training_seconds"],
                    "validation_objective": phase["training"]["best_validation_normalized_loss"],
                    "phase": phase,
                }
            )
    return candidates


def _decomposition(
    prediction: np.ndarray,
    reference: np.ndarray,
    variance_of_mean: np.ndarray,
    *,
    seed: int,
    bootstrap_repetitions: int = 2000,
) -> dict[str, float]:
    squared = np.square(np.asarray(prediction) - np.asarray(reference))
    variance = np.asarray(variance_of_mean)
    observed = float(np.mean(squared))
    mc = float(np.mean(variance))
    systematic = max(observed - mc, 0.0)
    rng = np.random.default_rng(seed)
    count = len(squared)
    observed_boot = np.empty(bootstrap_repetitions)
    mc_boot = np.empty(bootstrap_repetitions)
    systematic_boot = np.empty(bootstrap_repetitions)
    for start in range(0, bootstrap_repetitions, 100):
        stop = min(start + 100, bootstrap_repetitions)
        indices = rng.integers(0, count, size=(stop - start, count))
        local_observed = np.mean(squared[indices], axis=1)
        local_mc = np.mean(variance[indices], axis=1)
        observed_boot[start:stop] = local_observed
        mc_boot[start:stop] = local_mc
        systematic_boot[start:stop] = np.maximum(local_observed - local_mc, 0.0)

    def interval(values: np.ndarray) -> tuple[float, float]:
        return tuple(float(value) for value in np.quantile(values, (0.025, 0.975)))

    obs_low, obs_high = interval(observed_boot)
    mc_low, mc_high = interval(mc_boot)
    sys_low, sys_high = interval(systematic_boot)
    return {
        "observed_mse": observed, "observed_mse_ci95_low": obs_low,
        "observed_mse_ci95_high": obs_high,
        "mc_variance_contribution": mc, "mc_variance_ci95_low": mc_low,
        "mc_variance_ci95_high": mc_high,
        "systematic_mse_proxy": systematic, "systematic_mse_ci95_low": sys_low,
        "systematic_mse_ci95_high": sys_high,
        "systematic_rmse_proxy": float(np.sqrt(systematic)),
        "systematic_rmse_ci95_low": float(np.sqrt(sys_low)),
        "systematic_rmse_ci95_high": float(np.sqrt(sys_high)),
        "mc_fraction_observed_mse": mc / observed if observed > 0 else 0.0,
        "bootstrap_repetitions": bootstrap_repetitions,
    }


def evaluate_actors(benchmark: str, seed: int, device: torch.device) -> dict[str, Any]:
    ensure_dirs()
    data = evaluation_data(benchmark)
    validation_indices = data["validation_indices"]
    test_indices = data["test_indices"]
    validation_x = data["x"][validation_indices]
    test_x = data["x"][test_indices]
    candidates = _actor_candidates(benchmark, seed)
    selection: list[dict[str, Any]] = []
    actors: dict[str, torch.nn.Module] = {}
    # This pass does not read either reference array.
    common_validation_seed = stable_seed(
        "high_budget_common_actor_selection", benchmark, seed
    )
    for candidate in candidates:
        actor, _ = load_tangent_actor(Path(candidate["path"]), device)
        evaluation_id = candidate["evaluation_id"]
        actors[evaluation_id] = actor
        prediction, variance, metadata = evaluate_tangent_actor(
            benchmark, actor, validation_x, device, paths=512, time_steps=32,
            seed=common_validation_seed,
        )
        selection.append(
            {
                "evaluation_id": evaluation_id,
                "checkpoint": candidate["checkpoint"], "branch": candidate["branch"],
                "mean_policy_value": float(np.mean(prediction)),
                "mean_policy_value_mc_se": float(np.sqrt(np.sum(variance)) / len(variance)),
                "selection_inference_seconds": metadata["wall_seconds"],
                "weight_ess_fraction": metadata["weight_ess_fraction"],
                "nonfinite_count": metadata["nonfinite_count"],
                "ground_truth_used": False,
            }
        )
    selected = min(selection, key=lambda item: item["mean_policy_value"])
    selected_id = selected["evaluation_id"]

    # Only after checkpoint selection is frozen do we access ground truth.
    test_reference = data["reference_original"][test_indices]
    rows: list[dict[str, Any]] = []
    decompositions: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {"test_reference": test_reference}
    for candidate in candidates:
        label = candidate["checkpoint"]
        evaluation_id = candidate["evaluation_id"]
        actor = actors[evaluation_id]
        prediction, variance, metadata = evaluate_tangent_actor(
            benchmark, actor, test_x, device, paths=512, time_steps=32,
            seed=stable_seed("high_budget_test_common", benchmark, seed, 512, 32),
        )
        key = f"actor_{evaluation_id}_p512_s32"
        row = {
            "method": "tangent_actor_high_budget",
            "benchmark": benchmark, "seed": seed, "branch": candidate["branch"],
            "checkpoint": label, "evaluation_id": evaluation_id,
            "checkpoint_path": candidate["path"],
            "validation_selected": evaluation_id == selected_id,
            "training_seconds": candidate["training_seconds"],
            "selection_inference_seconds": next(
                item["selection_inference_seconds"] for item in selection
                if item["evaluation_id"] == evaluation_id
            ),
            "inference_seconds_N960": metadata["wall_seconds"],
            "total_seconds_N960": candidate["training_seconds"] + metadata["wall_seconds"],
            "inference_paths": 512, "inference_time_steps": 32,
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "mc_mean_se": metadata["mean_value_mc_se"],
            "prediction_key": key, "validation_objective": candidate["validation_objective"],
            "final_test_used_for_selection": False,
            **metrics(prediction, test_reference),
        }
        rows.append(row)
        decomposition = _decomposition(
            prediction, test_reference, variance,
            seed=stable_seed("high_budget_bootstrap", benchmark, seed, evaluation_id, 512, 32),
        )
        decompositions.append({**{key: row[key] for key in (
            "method", "benchmark", "seed", "branch", "checkpoint", "evaluation_id", "training_seconds",
            "inference_paths", "inference_time_steps", "validation_selected",
        )}, **decomposition})
        predictions[key] = prediction

    selected_candidate = next(
        item for item in candidates if item["evaluation_id"] == selected_id
    )
    selected_label = selected_candidate["checkpoint"]
    selected_actor = actors[selected_id]
    # High-quality inference for the validation-selected checkpoint only.
    for paths in (2048, 8192):
        prediction, variance, metadata = evaluate_tangent_actor(
            benchmark, selected_actor, test_x, device, paths=paths, time_steps=32,
            seed=stable_seed("high_budget_selected_test", benchmark, seed, paths, 32),
        )
        key = f"actor_{selected_id}_p{paths}_s32"
        row = {
            "method": "tangent_actor_high_budget", "benchmark": benchmark, "seed": seed,
            "branch": selected_candidate["branch"], "checkpoint": selected_label,
            "evaluation_id": selected_id,
            "checkpoint_path": selected_candidate["path"], "validation_selected": True,
            "training_seconds": selected_candidate["training_seconds"],
            "selection_inference_seconds": selected["selection_inference_seconds"],
            "inference_seconds_N960": metadata["wall_seconds"],
            "total_seconds_N960": selected_candidate["training_seconds"] + metadata["wall_seconds"],
            "inference_paths": paths, "inference_time_steps": 32,
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "mc_mean_se": metadata["mean_value_mc_se"], "prediction_key": key,
            "validation_objective": selected_candidate["validation_objective"],
            "final_test_used_for_selection": False,
            **metrics(prediction, test_reference),
        }
        rows.append(row)
        decomposition = _decomposition(
            prediction, test_reference, variance,
            seed=stable_seed("high_budget_bootstrap", benchmark, seed, selected_id, paths, 32),
        )
        decompositions.append({**{key: row[key] for key in (
            "method", "benchmark", "seed", "branch", "checkpoint", "evaluation_id", "training_seconds",
            "inference_paths", "inference_time_steps", "validation_selected",
        )}, **decomposition})
        predictions[key] = prediction

    # Euler sensitivity at fixed M=2048. Step 32 is already included above.
    for steps in (16, 64):
        prediction, variance, metadata = evaluate_tangent_actor(
            benchmark, selected_actor, test_x, device, paths=2048, time_steps=steps,
            seed=stable_seed("high_budget_euler_test", benchmark, seed, 2048, steps),
        )
        key = f"actor_{selected_id}_p2048_s{steps}"
        row = {
            "method": "tangent_actor_high_budget", "benchmark": benchmark, "seed": seed,
            "branch": selected_candidate["branch"], "checkpoint": selected_label,
            "evaluation_id": selected_id,
            "checkpoint_path": selected_candidate["path"], "validation_selected": True,
            "training_seconds": selected_candidate["training_seconds"],
            "selection_inference_seconds": selected["selection_inference_seconds"],
            "inference_seconds_N960": metadata["wall_seconds"],
            "total_seconds_N960": selected_candidate["training_seconds"] + metadata["wall_seconds"],
            "inference_paths": 2048, "inference_time_steps": steps,
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "mc_mean_se": metadata["mean_value_mc_se"], "prediction_key": key,
            "validation_objective": selected_candidate["validation_objective"],
            "final_test_used_for_selection": False,
            **metrics(prediction, test_reference),
        }
        rows.append(row)
        decomposition = _decomposition(
            prediction, test_reference, variance,
            seed=stable_seed("high_budget_bootstrap", benchmark, seed, selected_id, 2048, steps),
        )
        decompositions.append({**{key: row[key] for key in (
            "method", "benchmark", "seed", "branch", "checkpoint", "evaluation_id", "training_seconds",
            "inference_paths", "inference_time_steps", "validation_selected",
        )}, **decomposition})
        predictions[key] = prediction

    prediction_path = PREDICTION_DIR / f"actor_evaluation_{benchmark}_s{seed}.npz"
    np.savez_compressed(prediction_path, **predictions)
    result = {
        "benchmark": benchmark, "seed": seed, "selection": selection,
        "selected_checkpoint": selected_label, "selected_evaluation_id": selected_id,
        "selection_rule": "minimum mean fixed-bank tangent-policy value",
        "ground_truth_used_for_selection": False, "rows": rows,
        "decompositions": decompositions, "prediction_path": str(prediction_path),
    }
    atomic_json(RAW_DIR / f"actor_evaluation_{benchmark}_s{seed}.json", result)
    return result


def evaluate_original_baselines(
    benchmark: str, seed: int, device: torch.device
) -> dict[str, Any]:
    ensure_dirs()
    data = evaluation_data(benchmark)
    validation_indices, test_indices = data["validation_indices"], data["test_indices"]
    validation_x = data["x"][validation_indices]
    test_x = data["x"][test_indices]
    test_reference = data["reference_original"][test_indices]
    validation_reference = data["reference_original"][validation_indices]
    driver = driver_for(benchmark)
    rows: list[dict[str, Any]] = []
    predictions: dict[str, np.ndarray] = {"test_reference": test_reference}
    for method in ("deep_bsde_original", "sdgd_pinn_original"):
        result = _json(RAW_DIR / f"{method}_{benchmark}_s{seed}.json")
        for index, checkpoint in enumerate(result["checkpoints"], 1):
            if method == "deep_bsde_original":
                model, _ = load_deep_bsde_original(Path(checkpoint["checkpoint"]), device)
                prediction, inference = predict_deep_bsde_original(model, driver, test_x, device)
            else:
                model, _ = load_sdgd_pinn_original(Path(checkpoint["checkpoint"]), device)
                prediction, inference = predict_sdgd_pinn_original(model, driver, test_x, device)
            key = f"{method}_B{index}"
            rows.append(
                {
                    "method": method, "benchmark": benchmark, "seed": seed,
                    "checkpoint": f"B{index}", "checkpoint_path": checkpoint["checkpoint"],
                    "training_seconds": checkpoint["training_seconds"],
                    "inference_seconds_N960": inference,
                    "total_seconds_N960": checkpoint["training_seconds"] + inference,
                    "target_budget_seconds": checkpoint["target_budget_seconds"],
                    "optimizer_step": checkpoint["step"],
                    "validation_objective": checkpoint["validation_objective"],
                    "prediction_key": key, "terminal_aware_features": False,
                    "final_test_used_for_selection": False,
                    **metrics(prediction, test_reference),
                }
            )
            predictions[key] = prediction

    # Corrected SCaSML uses the final ordinary PINN, never the hard-terminal one.
    pinn_result = _json(RAW_DIR / f"sdgd_pinn_original_{benchmark}_s{seed}.json")
    final_pinn = pinn_result["checkpoints"][-1]
    pinn, _ = load_sdgd_pinn_original(Path(final_pinn["checkpoint"]), device)
    for paths in (8, 32):
        val_prediction, _, val_meta = predict_scasml_original_level1(
            pinn, benchmark, validation_x, device, paths=paths, probes=4,
            seed=stable_seed("high_budget_scasml", benchmark, seed, paths, "validation"),
        )
        prediction, variance, metadata = predict_scasml_original_level1(
            pinn, benchmark, test_x, device, paths=paths, probes=4,
            seed=stable_seed("high_budget_scasml", benchmark, seed, paths, "test"),
        )
        key = f"scasml_original_p{paths}"
        rows.append(
            {
                "method": "scasml_original", "benchmark": benchmark, "seed": seed,
                "checkpoint": "PINN_B5", "checkpoint_path": final_pinn["checkpoint"],
                "training_seconds": final_pinn["training_seconds"],
                "inference_seconds_N960": metadata["inference_seconds"],
                "total_seconds_N960": final_pinn["training_seconds"] + metadata["inference_seconds"],
                "inference_paths": paths, "hutchinson_probes": 4,
                "validation_relative_l2": metrics(val_prediction, validation_reference)["relative_l2"],
                "validation_inference_seconds": val_meta["inference_seconds"],
                "mc_mean_se": metadata["mean_mc_se_across_workload"],
                "mc_rmse_standard_error": float(np.sqrt(np.mean(variance))),
                "mean_terminal_defect_rms": metadata["mean_terminal_defect_rms"],
                "terminal_defect_included": True,
                "validation_objective": final_pinn["validation_objective"],
                "prediction_key": key, "terminal_aware_features": False,
                "final_test_used_for_selection": False,
                **metrics(prediction, test_reference),
            }
        )
        predictions[key] = prediction
    prediction_path = PREDICTION_DIR / f"original_baselines_{benchmark}_s{seed}.npz"
    np.savez_compressed(prediction_path, **predictions)
    result = {
        "benchmark": benchmark, "seed": seed, "rows": rows,
        "prediction_path": str(prediction_path), "ground_truth_used_for_selection": False,
    }
    atomic_json(RAW_DIR / f"original_baseline_evaluation_{benchmark}_s{seed}.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("actors", "baselines", "all"), default="all")
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--seed", type=int, choices=(0, 1), required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    started = time.perf_counter()
    completed = []
    if args.kind in ("actors", "all"):
        evaluate_actors(args.benchmark, args.seed, device)
        completed.append("actors")
    if args.kind in ("baselines", "all"):
        evaluate_original_baselines(args.benchmark, args.seed, device)
        completed.append("baselines")
    print(json.dumps({"status": "complete", "kind": completed, "wall_seconds": time.perf_counter() - started}), flush=True)


if __name__ == "__main__":
    main()
