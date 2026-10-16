"""Run frozen-Actor MCV calibration and independent replicated evaluation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    evaluate_tangent_mc,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    load_tangent_actor,
)

from .common import (
    BENCHMARKS,
    CALIBRATION_PATHS,
    INFERENCE_PATHS,
    MC_REPETITIONS,
    PREDICTION_DIR,
    RAW_DIR,
    TIME_STEPS,
    atomic_json,
    checkpoint_specs,
    driver_for,
    ensure_dirs,
    metrics,
    posthoc_reference,
    stable_seed,
    state_bank,
    synchronize,
)
from .mcv import calibrated_alpha, evaluate_pair_samples


T_CRITICAL_975_DF7 = 2.364624251


def _paired_state_bootstrap_ratio(
    numerator_state_variance: np.ndarray,
    denominator_state_variance: np.ndarray,
    *,
    seed: int,
    repetitions: int = 4000,
) -> tuple[float, float]:
    numerator = np.asarray(numerator_state_variance, dtype=np.float64)
    denominator = np.asarray(denominator_state_variance, dtype=np.float64)
    if numerator.shape != denominator.shape or numerator.ndim != 1:
        raise ValueError((numerator.shape, denominator.shape))
    rng = np.random.default_rng(seed)
    ratios = np.empty(repetitions, dtype=np.float64)
    for begin in range(0, repetitions, 200):
        end = min(begin + 200, repetitions)
        indices = rng.integers(0, len(numerator), size=(end - begin, len(numerator)))
        ratios[begin:end] = np.mean(numerator[indices], axis=1) / np.maximum(
            np.mean(denominator[indices], axis=1), np.finfo(float).tiny
        )
    low, high = np.quantile(ratios, (0.025, 0.975))
    return float(low), float(high)


def _paired_mean_summary(difference: np.ndarray) -> dict[str, float | bool]:
    """Use eight independent MC replicates as the uncertainty units."""
    array = np.asarray(difference, dtype=np.float64)
    replicate_means = array.mean(axis=1)
    mean = float(replicate_means.mean())
    standard_error = float(replicate_means.std(ddof=1) / np.sqrt(len(replicate_means)))
    radius = T_CRITICAL_975_DF7 * standard_error
    return {
        "mean": mean,
        "standard_error": standard_error,
        "ci95_low": mean - radius,
        "ci95_high": mean + radius,
        "ci95_contains_zero": mean - radius <= 0.0 <= mean + radius,
        "absolute_t_statistic": abs(mean) / standard_error if standard_error > 0 else 0.0,
    }


def _variant_row(
    *,
    benchmark: str,
    seed: int,
    checkpoint: dict[str, Any],
    paths: int,
    variant: str,
    alpha: float,
    predictions: np.ndarray,
    plain_predictions: np.ndarray,
    martingales: np.ndarray,
    pair_variances: np.ndarray,
    plain_pair_variances: np.ndarray,
    reference: np.ndarray,
    plain_seconds: np.ndarray,
    augmented_seconds: np.ndarray,
    calibration_seconds: float,
    weight_ess: np.ndarray,
    payoff_antithetic_ratio: np.ndarray,
    martingale_antithetic_ratio: np.ndarray,
    pair_correlation: np.ndarray,
) -> dict[str, Any]:
    replicate_state_variance = np.var(predictions, axis=0, ddof=1)
    plain_state_variance = np.var(plain_predictions, axis=0, ddof=1)
    denominator = max(float(np.mean(plain_state_variance)), np.finfo(float).tiny)
    variance_ratio = float(np.mean(replicate_state_variance) / denominator)
    ci_low, ci_high = _paired_state_bootstrap_ratio(
        replicate_state_variance,
        plain_state_variance,
        seed=stable_seed(
            "mcv_ratio_bootstrap", benchmark, seed, checkpoint["checkpoint_key"],
            paths, variant,
        ),
    )
    pair_ratio = float(
        np.mean(pair_variances) /
        max(float(np.mean(plain_pair_variances)), np.finfo(float).tiny)
    )
    mean_difference = _paired_mean_summary(predictions - plain_predictions)
    martingale_summary = _paired_mean_summary(martingales)
    replicate_metrics = [metrics(item, reference) for item in predictions]
    mean_prediction = predictions.mean(axis=0)
    mean_metrics = metrics(mean_prediction, reference)
    plain_time = float(np.mean(plain_seconds))
    if variant == "no_mcv":
        evaluation_time = plain_time
        amortized_time = plain_time
        peak_kind = "plain"
    else:
        evaluation_time = float(np.mean(augmented_seconds))
        amortized_time = evaluation_time + (
            calibration_seconds / len(augmented_seconds) if variant == "calibrated" else 0.0
        )
        peak_kind = "augmented"
    runtime_ratio = amortized_time / max(plain_time, np.finfo(float).tiny)
    return {
        "benchmark": benchmark,
        "seed": seed,
        "checkpoint_key": checkpoint["checkpoint_key"],
        "iteration": checkpoint["iteration"],
        "checkpoint_path": str(checkpoint["path"]),
        "source_study": checkpoint["source_study"],
        "variant": variant,
        "alpha": alpha,
        "states": predictions.shape[1],
        "paths": paths,
        "antithetic_pairs": paths // 2,
        "time_steps": TIME_STEPS,
        "mc_repetitions": predictions.shape[0],
        "variance_ratio_replicates": variance_ratio,
        "variance_ratio_ci95_low": ci_low,
        "variance_ratio_ci95_high": ci_high,
        "variance_ratio_pair_estimator": pair_ratio,
        "effective_variance_reduction": 1.0 / max(variance_ratio, np.finfo(float).tiny),
        "paths_for_same_standard_error": paths * variance_ratio,
        "paths_saved_for_same_standard_error": paths * (1.0 - variance_ratio),
        "path_saving_fraction": 1.0 - variance_ratio,
        "mean_state_mc_standard_error": float(np.sqrt(np.mean(replicate_state_variance))),
        "bank_mean_mc_standard_error": float(
            np.sqrt(np.sum(replicate_state_variance)) / predictions.shape[1]
        ),
        "pair_estimated_mean_state_mc_standard_error": float(
            np.sqrt(np.mean(pair_variances))
        ),
        "prediction_difference_mean": mean_difference["mean"],
        "prediction_difference_standard_error": mean_difference["standard_error"],
        "prediction_difference_ci95_low": mean_difference["ci95_low"],
        "prediction_difference_ci95_high": mean_difference["ci95_high"],
        "prediction_difference_ci95_contains_zero": mean_difference["ci95_contains_zero"],
        "prediction_difference_absolute_t_statistic": mean_difference["absolute_t_statistic"],
        "martingale_mean": martingale_summary["mean"],
        "martingale_mean_standard_error": martingale_summary["standard_error"],
        "martingale_mean_ci95_low": martingale_summary["ci95_low"],
        "martingale_mean_ci95_high": martingale_summary["ci95_high"],
        "martingale_ci95_contains_zero": martingale_summary["ci95_contains_zero"],
        "mean_replicate_rmse": float(np.mean([item["rmse"] for item in replicate_metrics])),
        "mean_replicate_relative_l2": float(
            np.mean([item["relative_l2"] for item in replicate_metrics])
        ),
        "rmse_of_replicate_mean": mean_metrics["rmse"],
        "relative_l2_of_replicate_mean": mean_metrics["relative_l2"],
        "bias_of_replicate_mean": mean_metrics["bias"],
        "plain_mean_seconds": plain_time,
        "evaluation_mean_seconds": evaluation_time,
        "calibration_seconds": calibration_seconds if variant == "calibrated" else 0.0,
        "amortized_seconds_per_replicate": amortized_time,
        "runtime_ratio_including_amortized_calibration": runtime_ratio,
        "effective_variance_reduction_per_gpu_second": (
            1.0 / max(variance_ratio * runtime_ratio, np.finfo(float).tiny)
        ),
        "runtime_peak_kind": peak_kind,
        "weight_ess_fraction": float(np.mean(weight_ess)),
        "payoff_antithetic_variance_ratio": float(np.mean(payoff_antithetic_ratio)),
        "martingale_antithetic_variance_ratio": float(
            np.mean(martingale_antithetic_ratio)
        ),
        "payoff_martingale_pair_correlation": float(np.mean(pair_correlation)),
        "ground_truth_used_for_alpha_or_selection": False,
        "ground_truth_role": "posthoc accuracy only",
    }


def run(
    benchmark: str,
    seed: int,
    checkpoint_keys: list[str],
    run_tag: str,
    device: torch.device,
) -> dict[str, Any]:
    ensure_dirs()
    driver = driver_for(benchmark)
    x, indices, bank_fingerprint = state_bank(benchmark)
    t = np.zeros(len(x), dtype=np.float64)
    available = {item["checkpoint_key"]: item for item in checkpoint_specs(benchmark, seed)}
    selected_specs = [available[key] for key in checkpoint_keys]
    rows: list[dict[str, Any]] = []
    runtime_rows: list[dict[str, Any]] = []
    calibration_rows: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    checkpoint_payloads: dict[str, dict[str, Any]] = {}
    prediction_payload: dict[str, np.ndarray] = {
        "state_indices": indices,
        "states": x,
    }
    synchronize(device)
    overall_started = time.perf_counter()
    for checkpoint in selected_specs:
        actor, payload = load_tangent_actor(Path(checkpoint["path"]), device)
        checkpoint_payloads[checkpoint["checkpoint_key"]] = {
            "actor_kind": str(payload.get("actor_kind", "unknown")),
            "width": int(payload.get("width", 0)),
            "depth": int(payload.get("depth", 0)),
        }
        calibration_seed = stable_seed(
            "mcv_calibration", benchmark, seed, checkpoint["checkpoint_key"]
        )
        calibration, calibration_meta = evaluate_pair_samples(
            driver,
            actor,
            t,
            x,
            paths=CALIBRATION_PATHS,
            time_steps=TIME_STEPS,
            seed=calibration_seed,
            device=device,
        )
        alpha, alpha_meta = calibrated_alpha(
            calibration["payoff_pairs"], calibration["martingale_pairs"]
        )
        calibration_rows.append(
            {
                "benchmark": benchmark,
                "seed": seed,
                "checkpoint_key": checkpoint["checkpoint_key"],
                "iteration": checkpoint["iteration"],
                "alpha": alpha,
                "calibration_paths": CALIBRATION_PATHS,
                "calibration_antithetic_pairs": CALIBRATION_PATHS // 2,
                "calibration_states": len(x),
                "calibration_seed": calibration_seed,
                "calibration_seconds": calibration_meta["wall_seconds"],
                "weight_ess_fraction": calibration_meta["weight_ess_fraction"],
                "payoff_martingale_pair_correlation": calibration_meta[
                    "payoff_martingale_pair_correlation"
                ],
                "ground_truth_used": False,
                **alpha_meta,
            }
        )
        for paths in INFERENCE_PATHS:
            plain_predictions = []
            direct_predictions = []
            calibrated_predictions = []
            martingale_predictions = []
            plain_pair_variances = []
            direct_pair_variances = []
            calibrated_pair_variances = []
            plain_seconds = []
            augmented_seconds = []
            plain_peak = []
            augmented_peak = []
            weight_ess = []
            payoff_antithetic_ratio = []
            martingale_antithetic_ratio = []
            pair_correlation = []
            ordinary_match_errors = []
            ordinary_estimator_mcv_flags = []
            for replicate in range(MC_REPETITIONS):
                evaluation_seed = stable_seed(
                    "mcv_evaluation", benchmark, seed, checkpoint["checkpoint_key"],
                    paths, replicate,
                )

                def ordinary_call() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
                    return evaluate_tangent_mc(
                        driver,
                        actor,
                        t,
                        x,
                        paths=paths,
                        time_steps=TIME_STEPS,
                        seed=evaluation_seed,
                        device=device,
                        parent_batch=8,
                        path_chunk=64,
                    )

                def augmented_call() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
                    return evaluate_pair_samples(
                        driver,
                        actor,
                        t,
                        x,
                        paths=paths,
                        time_steps=TIME_STEPS,
                        seed=evaluation_seed,
                        device=device,
                        parent_batch=8,
                        path_chunk=64,
                    )

                if replicate % 2 == 0:
                    ordinary, ordinary_meta = ordinary_call()
                    augmented, augmented_meta = augmented_call()
                else:
                    augmented, augmented_meta = augmented_call()
                    ordinary, ordinary_meta = ordinary_call()
                payoff_pairs = augmented["payoff_pairs"]
                martingale_pairs = augmented["martingale_pairs"]
                plain_working = payoff_pairs.mean(axis=1)
                ordinary_working = ordinary["working_value"]
                ordinary_match_errors.append(
                    float(np.max(np.abs(plain_working - ordinary_working)))
                )
                if ordinary_match_errors[-1] > 2e-6:
                    raise AssertionError(
                        f"ordinary FK mismatch {ordinary_match_errors[-1]:.3g}"
                    )
                direct_working = (payoff_pairs - martingale_pairs).mean(axis=1)
                calibrated_working = (
                    payoff_pairs - alpha * martingale_pairs
                ).mean(axis=1)
                plain_predictions.append(driver.working_to_original(plain_working))
                direct_predictions.append(driver.working_to_original(direct_working))
                calibrated_predictions.append(driver.working_to_original(calibrated_working))
                martingale_predictions.append(
                    driver.reconstruction_sign * martingale_pairs.mean(axis=1)
                )
                pair_count = paths // 2
                plain_pair_variances.append(
                    np.var(payoff_pairs, axis=1, ddof=1) / pair_count
                )
                direct_pair_variances.append(
                    np.var(payoff_pairs - martingale_pairs, axis=1, ddof=1) / pair_count
                )
                calibrated_pair_variances.append(
                    np.var(payoff_pairs - alpha * martingale_pairs, axis=1, ddof=1)
                    / pair_count
                )
                plain_seconds.append(ordinary_meta["wall_seconds"])
                augmented_seconds.append(augmented_meta["wall_seconds"])
                plain_peak.append(ordinary_meta["peak_gpu_bytes"] or 0)
                augmented_peak.append(augmented_meta["peak_gpu_bytes"] or 0)
                weight_ess.append(augmented_meta["weight_ess_fraction"])
                payoff_antithetic_ratio.append(
                    augmented_meta["payoff_antithetic_variance_ratio"]
                )
                martingale_antithetic_ratio.append(
                    augmented_meta["martingale_antithetic_variance_ratio"]
                )
                pair_correlation.append(
                    augmented_meta["payoff_martingale_pair_correlation"]
                )
                ordinary_estimator_mcv_flags.append(
                    bool(ordinary_meta["martingale_control_variate"])
                )
            arrays = {
                "no_mcv": np.asarray(plain_predictions),
                "direct": np.asarray(direct_predictions),
                "calibrated": np.asarray(calibrated_predictions),
            }
            pair_variances = {
                "no_mcv": np.asarray(plain_pair_variances),
                "direct": np.asarray(direct_pair_variances),
                "calibrated": np.asarray(calibrated_pair_variances),
            }
            martingales = np.asarray(martingale_predictions)
            runtime_rows.append(
                {
                    "benchmark": benchmark,
                    "seed": seed,
                    "checkpoint_key": checkpoint["checkpoint_key"],
                    "iteration": checkpoint["iteration"],
                    "paths": paths,
                    "states": len(x),
                    "mc_repetitions": MC_REPETITIONS,
                    "plain_mean_seconds": float(np.mean(plain_seconds)),
                    "plain_seconds_std": float(np.std(plain_seconds, ddof=1)),
                    "mcv_mean_seconds": float(np.mean(augmented_seconds)),
                    "mcv_seconds_std": float(np.std(augmented_seconds, ddof=1)),
                    "mcv_runtime_ratio": float(
                        np.mean(augmented_seconds) / np.mean(plain_seconds)
                    ),
                    "mcv_runtime_overhead_fraction": float(
                        np.mean(augmented_seconds) / np.mean(plain_seconds) - 1.0
                    ),
                    "plain_peak_gpu_bytes": int(max(plain_peak)),
                    "mcv_peak_gpu_bytes": int(max(augmented_peak)),
                    "additional_peak_gpu_bytes": int(max(augmented_peak) - max(plain_peak)),
                    "calibration_seconds": calibration_meta["wall_seconds"],
                    "ordinary_fk_already_uses_mcv": any(ordinary_estimator_mcv_flags),
                    "ordinary_fk_match_max_abs": max(ordinary_match_errors),
                }
            )
            pending.append(
                {
                    "checkpoint": checkpoint,
                    "paths": paths,
                    "alpha": alpha,
                    "arrays": arrays,
                    "pair_variances": pair_variances,
                    "martingales": martingales,
                    "plain_seconds": np.asarray(plain_seconds),
                    "augmented_seconds": np.asarray(augmented_seconds),
                    "calibration_seconds": float(calibration_meta["wall_seconds"]),
                    "weight_ess": np.asarray(weight_ess),
                    "payoff_antithetic_ratio": np.asarray(payoff_antithetic_ratio),
                    "martingale_antithetic_ratio": np.asarray(martingale_antithetic_ratio),
                    "pair_correlation": np.asarray(pair_correlation),
                }
            )
            print(
                f"{benchmark} s{seed} {checkpoint['checkpoint_key']} M={paths} "
                f"alpha={alpha:.4f} plain={np.mean(plain_seconds):.3f}s "
                f"mcv={np.mean(augmented_seconds):.3f}s",
                flush=True,
            )
        del actor
        if device.type == "cuda":
            torch.cuda.empty_cache()

    # The numerical reference is opened only after alpha and all predictions are frozen.
    reference = posthoc_reference(benchmark, indices)
    prediction_payload["reference"] = reference
    for item in pending:
        checkpoint = item["checkpoint"]
        paths = item["paths"]
        plain = item["arrays"]["no_mcv"]
        for variant, alpha in (
            ("no_mcv", 0.0), ("direct", 1.0), ("calibrated", item["alpha"])
        ):
            row = _variant_row(
                benchmark=benchmark,
                seed=seed,
                checkpoint=checkpoint,
                paths=paths,
                variant=variant,
                alpha=alpha,
                predictions=item["arrays"][variant],
                plain_predictions=plain,
                martingales=item["martingales"],
                pair_variances=item["pair_variances"][variant],
                plain_pair_variances=item["pair_variances"]["no_mcv"],
                reference=reference,
                plain_seconds=item["plain_seconds"],
                augmented_seconds=item["augmented_seconds"],
                calibration_seconds=item["calibration_seconds"],
                weight_ess=item["weight_ess"],
                payoff_antithetic_ratio=item["payoff_antithetic_ratio"],
                martingale_antithetic_ratio=item["martingale_antithetic_ratio"],
                pair_correlation=item["pair_correlation"],
            )
            key = f"{checkpoint['checkpoint_key']}_M{paths}_{variant}"
            prediction_payload[f"{key}_replicates"] = item["arrays"][variant]
            rows.append({**row, "prediction_key": f"{key}_replicates"})

    prediction_path = PREDICTION_DIR / f"mcv_{run_tag}_{benchmark}_s{seed}.npz"
    np.savez_compressed(prediction_path, **prediction_payload)
    result = {
        "experiment": "tangent_actor_mcv_quick",
        "run_tag": run_tag,
        "benchmark": benchmark,
        "seed": seed,
        "checkpoint_keys": checkpoint_keys,
        "checkpoint_payloads": checkpoint_payloads,
        "state_bank_fingerprint": bank_fingerprint,
        "state_count": len(x),
        "reference_opened_after_all_stochastic_evaluation": True,
        "ground_truth_used_for_alpha_or_selection": False,
        "ordinary_solver_modified": False,
        "training_performed": False,
        "calibration": calibration_rows,
        "variance_rows": rows,
        "runtime_rows": runtime_rows,
        "prediction_path": str(prediction_path),
        "wall_seconds": time.perf_counter() - overall_started,
    }
    output = RAW_DIR / f"mcv_{run_tag}_{benchmark}_s{seed}.json"
    atomic_json(output, result)
    print(json.dumps({"status": "complete", "output": str(output), "rows": len(rows)}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", choices=BENCHMARKS, required=True)
    parser.add_argument("--seed", type=int, choices=(0, 1), required=True)
    parser.add_argument(
        "--checkpoints", default="all",
        help="Comma-separated subset of a1,a3,late or all",
    )
    parser.add_argument("--run-tag", required=True)
    args = parser.parse_args()
    keys = ["a1", "a3", "late"] if args.checkpoints == "all" else args.checkpoints.split(",")
    if not keys or any(key not in {"a1", "a3", "late"} for key in keys):
        raise ValueError(keys)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    run(args.benchmark, args.seed, keys, args.run_tag, device)


if __name__ == "__main__":
    main()
