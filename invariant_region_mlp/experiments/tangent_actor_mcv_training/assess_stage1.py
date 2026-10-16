"""Reference-free Stage 1 gate with common high-quality Actor diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    make_targets,
)
from invariant_region_mlp.experiments.tangent_actor_high_budget.actor_scaling import (
    _state_bank,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.common import (
    evaluation_data,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    evaluate_tangent_actor,
    load_tangent_actor,
)

from .common import (
    HIGH_BUDGET_RESULTS_DIR,
    RAW_DIR,
    RESULTS_DIR,
    a3_checkpoint,
    atomic_json,
    driver_for,
    ensure_dirs,
    metrics,
    stable_seed,
    write_csv,
)
from .labels import mcv_tangent_labels


def _actor_outputs(
    actor: torch.nn.Module,
    driver: Any,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    values: list[np.ndarray] = []
    gradients: list[np.ndarray] = []
    actor.eval()
    with torch.no_grad():
        for begin in range(0, len(t), 512):
            end = min(begin + 512, len(t))
            bt = torch.as_tensor(t[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(x[begin:end], dtype=torch.float32, device=device)
            r, q = actor.predict(bt, bx, driver.T)
            assert q is not None
            values.append(r.double().cpu().numpy())
            gradients.append(q.double().cpu().numpy())
    return np.concatenate(values), np.concatenate(gradients)


def _plain_checkpoint(benchmark: str) -> Path:
    path = (
        HIGH_BUDGET_RESULTS_DIR / "checkpoints"
        / f"tangent_actor_{benchmark}_s0_a4.pt"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def _mcv_checkpoint(benchmark: str) -> Path:
    path = (
        RESULTS_DIR / "checkpoints"
        / f"tangent_actor_{benchmark}_s0_a4_mcv_points2048.pt"
    )
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def assess_benchmark(
    benchmark: str, device: torch.device
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    seed = 0
    driver = driver_for(benchmark)
    source, _ = load_tangent_actor(a3_checkpoint(benchmark, seed), device)
    actors = {
        "plain": load_tangent_actor(_plain_checkpoint(benchmark), device)[0],
        "mcv": load_tangent_actor(_mcv_checkpoint(benchmark), device)[0],
    }

    bank = _state_bank(benchmark, seed, 1)
    valid_t, valid_x = bank["valid_t"], bank["valid_x"]
    target_value, target_q, target_meta = mcv_tangent_labels(
        driver,
        source,
        valid_t,
        valid_x,
        paths=512,
        time_steps=16,
        seed=stable_seed("mcv_training_stage1_common_actor_target", benchmark, seed),
        device=device,
        parent_batch=8,
        path_chunk=32,
    )
    target_value, target_q, projection = make_targets(driver, target_value, target_q)
    assert target_q is not None

    data = evaluation_data(benchmark)
    validation_indices = data["validation_indices"]
    policy_x = data["x"][validation_indices]
    # References are read only after the checkpoint pair and the reference-free
    # gate protocol have been frozen; they do not enter the gate below.
    posthoc_reference = data["reference_original"][validation_indices]
    rows: list[dict[str, Any]] = []
    policy_predictions: dict[str, list[np.ndarray]] = {name: [] for name in actors}
    policy_metadata: dict[str, list[dict[str, Any]]] = {name: [] for name in actors}
    for replicate in range(4):
        inference_seed = stable_seed(
            "mcv_training_stage1_common_policy", benchmark, seed, replicate
        )
        for variant, actor in actors.items():
            prediction, _, metadata = evaluate_tangent_actor(
                benchmark,
                actor,
                policy_x,
                device,
                paths=512,
                time_steps=32,
                seed=inference_seed,
            )
            policy_predictions[variant].append(prediction)
            policy_metadata[variant].append(metadata)

    for variant, actor in actors.items():
        predicted_value, predicted_q = _actor_outputs(
            actor, driver, valid_t, valid_x, device
        )
        value_error = predicted_value - target_value
        gradient_error = predicted_q - target_q
        policy_matrix = np.asarray(policy_predictions[variant])
        policy_mean = policy_matrix.mean(axis=0)
        rows.append({
            "stage": 1,
            "benchmark": benchmark,
            "seed": seed,
            "variant": variant,
            "training_points": 2048,
            "common_actor_target_paths": 512,
            "actor_validation_value_rmse": float(np.sqrt(np.mean(value_error**2))),
            "actor_validation_gradient_component_rmse": float(
                np.sqrt(np.mean(gradient_error**2))
            ),
            "actor_validation_gradient_vector_rmse": float(
                np.sqrt(np.mean(np.sum(gradient_error**2, axis=1)))
            ),
            "common_policy_value_mean": float(np.mean(policy_mean)),
            "common_policy_value_mc_se": float(
                np.std(policy_matrix.mean(axis=1), ddof=1) / 2.0
            ),
            "posthoc_validation_value_rmse": metrics(
                policy_mean, posthoc_reference
            )["rmse"],
            "posthoc_validation_relative_l2": metrics(
                policy_mean, posthoc_reference
            )["relative_l2"],
            "policy_inference_seconds_mean": float(np.mean([
                item["wall_seconds"] for item in policy_metadata[variant]
            ])),
            "weight_ess_fraction_mean": float(np.mean([
                item["weight_ess_fraction"] for item in policy_metadata[variant]
            ])),
            "nonfinite_count": int(sum(
                item["nonfinite_count"] for item in policy_metadata[variant]
            )),
            "ground_truth_used_for_gate": False,
        })

    plain_policy = np.asarray(policy_predictions["plain"])
    mcv_policy = np.asarray(policy_predictions["mcv"])
    replicate_difference = plain_policy.mean(axis=1) - mcv_policy.mean(axis=1)
    policy_delta = float(np.mean(replicate_difference))
    policy_se = float(np.std(replicate_difference, ddof=1) / 2.0)
    plain = next(row for row in rows if row["variant"] == "plain")
    mcv = next(row for row in rows if row["variant"] == "mcv")
    comparison = {
        "benchmark": benchmark,
        "plain_minus_mcv_policy_value": policy_delta,
        "plain_minus_mcv_policy_value_mc_se": policy_se,
        "plain_minus_mcv_policy_value_ci95_low": policy_delta - 3.182446 * policy_se,
        "plain_minus_mcv_policy_value_ci95_high": policy_delta + 3.182446 * policy_se,
        "actor_value_error_ratio_mcv_over_plain": (
            mcv["actor_validation_value_rmse"] / plain["actor_validation_value_rmse"]
        ),
        "actor_gradient_vector_error_ratio_mcv_over_plain": (
            mcv["actor_validation_gradient_vector_rmse"]
            / plain["actor_validation_gradient_vector_rmse"]
        ),
        "posthoc_rmse_ratio_mcv_over_plain": (
            mcv["posthoc_validation_value_rmse"]
            / plain["posthoc_validation_value_rmse"]
        ),
        "common_actor_target_metadata": target_meta,
        "target_projection": projection,
    }
    return rows, comparison


def run(args: argparse.Namespace) -> None:
    ensure_dirs()
    device = torch.device(args.device)
    all_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for benchmark in ("nonseparable_d100_k2", "joint_concave_d400_k2"):
        rows, comparison = assess_benchmark(benchmark, device)
        all_rows.extend(rows)
        comparisons.append(comparison)
        print(
            f"stage1 {benchmark}: policy plain-mcv="
            f"{comparison['plain_minus_mcv_policy_value']:.6g}, "
            f"value ratio={comparison['actor_value_error_ratio_mcv_over_plain']:.3f}, "
            f"gradient ratio={comparison['actor_gradient_vector_error_ratio_mcv_over_plain']:.3f}",
            flush=True,
        )
    write_csv(RESULTS_DIR / "actor_validation_diagnostics.csv", all_rows)

    value_ratios = np.asarray([
        item["actor_value_error_ratio_mcv_over_plain"] for item in comparisons
    ])
    gradient_ratios = np.asarray([
        item["actor_gradient_vector_error_ratio_mcv_over_plain"] for item in comparisons
    ])
    policy_deltas = np.asarray([
        item["plain_minus_mcv_policy_value"] for item in comparisons
    ])
    # Predeclared practical gate: lower geometric-mean Actor error in both heads,
    # no >10% regression on either benchmark/head, and no average degradation in
    # the reference-free policy value. Ground Truth metrics remain descriptive.
    checks = {
        "geometric_mean_value_error_ratio_below_one": float(
            np.exp(np.mean(np.log(value_ratios)))
        ) < 1.0,
        "geometric_mean_gradient_error_ratio_below_one": float(
            np.exp(np.mean(np.log(gradient_ratios)))
        ) < 1.0,
        "no_actor_head_ratio_above_1p10": bool(
            max(np.max(value_ratios), np.max(gradient_ratios)) <= 1.10
        ),
        "mean_reference_free_policy_not_worse": float(np.mean(policy_deltas)) >= 0.0,
        "nonfinite_count_zero": all(row["nonfinite_count"] == 0 for row in all_rows),
    }
    passed = all(checks.values())
    decision = {
        "stage": 1,
        "passed": passed,
        "decision": "PROCEED_STAGE2" if passed else "STOP_AFTER_STAGE1",
        "gate_uses_ground_truth": False,
        "criteria_predeclared_in_code_before_results": {
            "geometric_mean_value_error_ratio_mcv_over_plain": "< 1.0",
            "geometric_mean_gradient_vector_error_ratio_mcv_over_plain": "< 1.0",
            "maximum_single_actor_error_ratio": "<= 1.10",
            "mean_plain_minus_mcv_reference_free_policy_value": ">= 0",
            "nonfinite_count": 0,
        },
        "checks": checks,
        "comparisons": comparisons,
    }
    atomic_json(RESULTS_DIR / "stage1_decision.json", decision)
    atomic_json(RAW_DIR / "stage1_assessment.json", {
        "rows": all_rows, "decision": decision,
    })
    print(f"STAGE1_DECISION={decision['decision']}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
