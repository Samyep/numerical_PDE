"""Common-inference posthoc evaluation for Plain- and MCV-trained Actors."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_high_budget.evaluate import (
    _decomposition,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_quick.mcv import (
    evaluate_pair_samples,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    evaluate_tangent_actor,
    load_tangent_actor,
)

from .common import (
    HIGH_BUDGET_RESULTS_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    STRONG_RESULTS_DIR,
    atomic_json,
    driver_for,
    ensure_dirs,
    metrics,
    stable_seed,
    write_csv,
)


BUDGETS = (240.0, 480.0, 960.0)


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _mcv_candidates(benchmark: str, seed: int) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for points in (2048, 8192):
        path = RAW_DIR / f"mcv_points{points}_{benchmark}_s{seed}.json"
        if not path.exists():
            continue
        result = _json(path)
        for phase in result["phases"]:
            candidates.append({
                "evaluation_id": phase["checkpoint_label"],
                "checkpoint": phase["checkpoint_label"],
                "checkpoint_path": phase["checkpoint"],
                "branch": phase["branch"],
                "training_points": points,
                "paths_per_label": result["configuration"]["paths_per_label"],
                "iteration": phase["iteration"],
                "training_seconds": phase["cumulative_training_seconds"],
                "additional_continuation_seconds": phase["additional_continuation_seconds"],
                "validation_objective": phase["training"][
                    "best_validation_normalized_loss"
                ],
            })
    return candidates


def _plain_result(benchmark: str, seed: int) -> dict[str, Any]:
    return _json(
        HIGH_BUDGET_RESULTS_DIR / "raw" / f"actor_evaluation_{benchmark}_s{seed}.json"
    )


def _plain_catalog(benchmark: str, seed: int) -> list[dict[str, Any]]:
    result = _plain_result(benchmark, seed)
    decompositions = {
        (row["evaluation_id"], row["inference_paths"], row["inference_time_steps"]): row
        for row in result["decompositions"]
    }
    validation = {row["evaluation_id"]: row for row in result["selection"]}
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in result["rows"]:
        if row["inference_paths"] != 512 or row["inference_time_steps"] != 32:
            continue
        evaluation_id = row["evaluation_id"]
        if evaluation_id in seen:
            continue
        seen.add(evaluation_id)
        decomposition = decompositions[(evaluation_id, 512, 32)]
        branch = row["branch"]
        points = 2048
        if "points4096" in branch:
            points = 4096
        elif "points8192" in branch:
            points = 8192
        paths_per_label = 128
        if "paths256" in branch:
            paths_per_label = 256
        elif "paths512" in branch:
            paths_per_label = 512
        rows.append({
            "training_variant": "plain",
            "benchmark": benchmark,
            "seed": seed,
            "evaluation_id": evaluation_id,
            "checkpoint": row["checkpoint"],
            "checkpoint_path": row["checkpoint_path"],
            "branch": branch,
            "training_points": points,
            "paths_per_label": paths_per_label,
            "training_seconds": row["training_seconds"],
            "validation_policy_mean": validation[evaluation_id]["mean_policy_value"],
            "validation_policy_mc_se": validation[evaluation_id][
                "mean_policy_value_mc_se"
            ],
            "inference_variant": "plain",
            "inference_paths": 512,
            "inference_time_steps": 32,
            "inference_seconds": row["inference_seconds_N960"],
            "rmse": row["rmse"],
            "relative_l2": row["relative_l2"],
            "bias": row["bias"],
            "weight_ess_fraction": row["weight_ess_fraction"],
            "systematic_rmse_proxy": decomposition["systematic_rmse_proxy"],
            "systematic_rmse_ci95_low": decomposition["systematic_rmse_ci95_low"],
            "systematic_rmse_ci95_high": decomposition["systematic_rmse_ci95_high"],
            "mc_variance_contribution": decomposition["mc_variance_contribution"],
            "mc_fraction_observed_mse": decomposition["mc_fraction_observed_mse"],
            "final_test_used_for_selection": False,
        })
    return rows


def _evaluate_mcv_training_candidates(
    benchmark: str, seed: int, device: torch.device
) -> list[dict[str, Any]]:
    candidates = _mcv_candidates(benchmark, seed)
    output_path = RAW_DIR / f"common_evaluation_{benchmark}_s{seed}.json"
    result: dict[str, Any] = (
        _json(output_path) if output_path.exists()
        else {"benchmark": benchmark, "seed": seed, "rows": []}
    )
    completed = {row["evaluation_id"] for row in result["rows"]}
    data_path = STRONG_RESULTS_DIR / f"evaluation_points_{benchmark}.npz"
    with np.load(data_path) as archive:
        x = archive["x"].copy()
        validation_indices = archive["validation_indices"].copy()
        test_indices = archive["test_indices"].copy()
        test_reference = archive["reference_original"][test_indices].copy()
    validation_x = x[validation_indices]
    test_x = x[test_indices]
    common_validation_seed = stable_seed(
        "mcv_training_common_validation", benchmark, seed, 512, 32
    )
    common_test_seed = stable_seed("high_budget_test_common", benchmark, seed, 512, 32)
    prediction_payload: dict[str, np.ndarray] = {"test_reference": test_reference}
    for candidate in candidates:
        if candidate["evaluation_id"] in completed:
            continue
        actor, _ = load_tangent_actor(Path(candidate["checkpoint_path"]), device)
        validation_prediction, validation_variance, validation_meta = evaluate_tangent_actor(
            benchmark,
            actor,
            validation_x,
            device,
            paths=512,
            time_steps=32,
            seed=common_validation_seed,
        )
        prediction, variance, metadata = evaluate_tangent_actor(
            benchmark,
            actor,
            test_x,
            device,
            paths=512,
            time_steps=32,
            seed=common_test_seed,
        )
        decomposition = _decomposition(
            prediction,
            test_reference,
            variance,
            seed=stable_seed(
                "mcv_training_systematic_bootstrap",
                benchmark,
                seed,
                candidate["evaluation_id"],
            ),
        )
        key = f"{candidate['evaluation_id']}_plain_p512_s32"
        prediction_payload[key] = prediction
        row = {
            "training_variant": "mcv",
            "benchmark": benchmark,
            "seed": seed,
            **candidate,
            "validation_policy_mean": float(np.mean(validation_prediction)),
            "validation_policy_mc_se": float(
                np.sqrt(np.sum(validation_variance)) / len(validation_variance)
            ),
            "validation_inference_seconds": validation_meta["wall_seconds"],
            "inference_variant": "plain",
            "inference_paths": 512,
            "inference_time_steps": 32,
            "inference_seconds": metadata["wall_seconds"],
            "weight_ess_fraction": metadata["weight_ess_fraction"],
            "nonfinite_count": metadata["nonfinite_count"],
            **metrics(prediction, test_reference),
            **decomposition,
            "prediction_key": key,
            "final_test_used_for_selection": False,
        }
        result["rows"].append(row)
        atomic_json(output_path, result)
        print(
            f"evaluate {benchmark} s{seed} {candidate['evaluation_id']}: "
            f"rel={row['relative_l2']:.6g} sys={row['systematic_rmse_proxy']:.6g}",
            flush=True,
        )
    prediction_path = PREDICTION_DIR / f"common_evaluation_{benchmark}_s{seed}.npz"
    if prediction_payload:
        np.savez_compressed(prediction_path, **prediction_payload)
    return result["rows"]


def _selected_at_budget(
    rows: list[dict[str, Any]], budget: float
) -> dict[str, Any]:
    available = [row for row in rows if row["training_seconds"] <= budget]
    if not available:
        raise RuntimeError(f"no complete checkpoint available by {budget}s")
    return min(available, key=lambda row: row["validation_policy_mean"])


def _evaluate_with_mcv_inference(
    actor_path: str,
    benchmark: str,
    seed: int,
    *,
    checkpoint: str,
    training_variant: str,
    training_seconds: float,
    steps: int,
    device: torch.device,
) -> dict[str, Any]:
    driver = driver_for(benchmark)
    actor, _ = load_tangent_actor(Path(actor_path), device)
    data_path = STRONG_RESULTS_DIR / f"evaluation_points_{benchmark}.npz"
    with np.load(data_path) as archive:
        test_indices = archive["test_indices"].copy()
        test_x = archive["x"][test_indices].copy()
        reference = archive["reference_original"][test_indices].copy()
    paths = 512
    samples, metadata = evaluate_pair_samples(
        driver,
        actor,
        np.zeros(len(test_x), dtype=np.float64),
        test_x,
        paths=paths,
        time_steps=steps,
        seed=stable_seed(
            "mcv_training_final_mcv_inference", benchmark, seed, checkpoint, paths, steps
        ),
        device=device,
        parent_batch=8,
        path_chunk=64,
        collect_martingale=True,
    )
    corrected_pairs = samples["payoff_pairs"] - samples["martingale_pairs"]
    working_prediction = corrected_pairs.mean(axis=1)
    prediction = np.asarray(driver.working_to_original(working_prediction))
    variance = corrected_pairs.var(axis=1, ddof=1) / corrected_pairs.shape[1]
    decomposition = _decomposition(
        prediction,
        reference,
        variance,
        seed=stable_seed(
            "mcv_training_final_mcv_decomposition",
            benchmark,
            seed,
            checkpoint,
            steps,
        ),
    )
    return {
        "training_variant": training_variant,
        "inference_variant": "mcv",
        "benchmark": benchmark,
        "seed": seed,
        "checkpoint": checkpoint,
        "checkpoint_path": actor_path,
        "training_seconds": training_seconds,
        "inference_paths": paths,
        "inference_time_steps": steps,
        "inference_seconds": metadata["wall_seconds"],
        "weight_ess_fraction": metadata["weight_ess_fraction"],
        "nonfinite_count": metadata["nonfinite_count"],
        **metrics(prediction, reference),
        **decomposition,
        "final_test_used_for_selection": False,
    }


def run(args: argparse.Namespace) -> None:
    ensure_dirs()
    device = torch.device(args.device)
    all_rows: list[dict[str, Any]] = []
    paired_final: list[dict[str, Any]] = []
    euler_rows: list[dict[str, Any]] = []
    budget_rows: list[dict[str, Any]] = []
    for benchmark in ("nonseparable_d100_k2", "joint_concave_d400_k2"):
        for seed in (0, 1):
            plain_rows = _plain_catalog(benchmark, seed)
            mcv_rows = _evaluate_mcv_training_candidates(benchmark, seed, device)
            all_rows.extend(plain_rows)
            all_rows.extend(mcv_rows)
            for budget in BUDGETS:
                for variant, rows in (("plain", plain_rows), ("mcv", mcv_rows)):
                    # MCV and Plain share a3 before the first continuation. Add
                    # a3 to the MCV available set for the early B3 boundary.
                    available = rows
                    if variant == "mcv":
                        base = next(
                            row for row in plain_rows
                            if row["branch"] == "frozen_original"
                        )
                        available = [{**base, "training_variant": "mcv"}, *rows]
                    selected = _selected_at_budget(available, budget)
                    budget_rows.append({
                        "benchmark": benchmark,
                        "seed": seed,
                        "budget_seconds": budget,
                        "training_variant": variant,
                        "selected_checkpoint": selected["checkpoint"],
                        "selected_branch": selected["branch"],
                        "actual_training_seconds": selected["training_seconds"],
                        "relative_l2": selected["relative_l2"],
                        "rmse": selected["rmse"],
                        "systematic_rmse_proxy": selected["systematic_rmse_proxy"],
                        "validation_policy_mean": selected["validation_policy_mean"],
                        "selection_uses_ground_truth": False,
                    })

            selected_plain = _selected_at_budget(plain_rows, 960.0)
            selected_mcv = _selected_at_budget(mcv_rows, 960.0)
            # Plain inference rows already use exactly the same N and time grid.
            paired_final.extend([
                {**selected_plain, "b5_validation_selected": True},
                {**selected_mcv, "b5_validation_selected": True},
            ])
            for selected in (selected_plain, selected_mcv):
                mcv_inference = _evaluate_with_mcv_inference(
                    selected["checkpoint_path"],
                    benchmark,
                    seed,
                    checkpoint=selected["checkpoint"],
                    training_variant=selected["training_variant"],
                    training_seconds=selected["training_seconds"],
                    steps=32,
                    device=device,
                )
                mcv_inference["b5_validation_selected"] = True
                paired_final.append(mcv_inference)
                if selected["training_variant"] == "mcv":
                    euler_rows.append(mcv_inference)
                    for steps in (16, 64):
                        euler_rows.append(_evaluate_with_mcv_inference(
                            selected["checkpoint_path"],
                            benchmark,
                            seed,
                            checkpoint=selected["checkpoint"],
                            training_variant="mcv",
                            training_seconds=selected["training_seconds"],
                            steps=steps,
                            device=device,
                        ))
    write_csv(RESULTS_DIR / "actor_training_frontier.csv", all_rows)
    write_csv(RESULTS_DIR / "systematic_error_comparison.csv", paired_final)
    write_csv(RESULTS_DIR / "budget_operating_points.csv", budget_rows)
    write_csv(RESULTS_DIR / "euler_sensitivity.csv", euler_rows)
    atomic_json(RAW_DIR / "evaluation_summary.json", {
        "frontier_rows": len(all_rows),
        "paired_final_rows": paired_final,
        "budget_operating_points": budget_rows,
        "euler_sensitivity": euler_rows,
        "primary_actor_quality_inference": {
            "variant": "plain weighted tangent FK",
            "paths": 512,
            "time_steps": 32,
            "common_across_training_variants": True,
        },
        "selection": "minimum common validation tangent-policy mean within available budget",
        "ground_truth_used_for_selection": False,
    })


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
