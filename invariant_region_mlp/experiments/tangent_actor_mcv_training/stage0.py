"""Mandatory label-noise and differentiability gate for MCV training."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    tangent_labels,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.tangent_baseline import (
    load_tangent_actor,
)

from .common import (
    LABEL_DIR,
    RAW_DIR,
    RESULTS_DIR,
    a3_checkpoint,
    array_sha256,
    atomic_json,
    driver_for,
    ensure_dirs,
    late_checkpoint,
    stable_seed,
    validation_bank,
    write_csv,
)
from .labels import mcv_tangent_labels


T975_DF7 = 2.364624251


def _bootstrap_ratio(
    plain: np.ndarray,
    mcv: np.ndarray,
    *,
    seed: int,
    draws: int = 2000,
) -> tuple[float, float]:
    """Bootstrap states; repetitions stay intact as independent MC draws."""
    if plain.shape != mcv.shape or plain.ndim < 2:
        raise ValueError((plain.shape, mcv.shape))
    rng = np.random.default_rng(seed)
    state_count = plain.shape[1]
    ratios = np.empty(draws, dtype=np.float64)
    for draw in range(draws):
        index = rng.integers(0, state_count, state_count)
        p = np.var(plain[:, index, ...], axis=0, ddof=1).mean()
        c = np.var(mcv[:, index, ...], axis=0, ddof=1).mean()
        ratios[draw] = c / max(p, np.finfo(float).tiny)
    return tuple(float(value) for value in np.quantile(ratios, [0.025, 0.975]))


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    a = left.reshape(-1)
    b = right.reshape(-1)
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _paired_interval(differences: np.ndarray) -> tuple[float, float, float, float]:
    """t interval over independent Brownian repetitions after bank averaging."""
    replicate_means = differences.reshape(differences.shape[0], -1).mean(axis=1)
    mean = float(np.mean(replicate_means))
    se = float(np.std(replicate_means, ddof=1) / np.sqrt(len(replicate_means)))
    return mean, se, mean - T975_DF7 * se, mean + T975_DF7 * se


def _summary_row(
    checkpoint: str,
    paths: int,
    plain_value: np.ndarray,
    mcv_value: np.ndarray,
    plain_z: np.ndarray,
    mcv_z: np.ndarray,
    plain_meta: list[dict[str, Any]],
    mcv_meta: list[dict[str, Any]],
    *,
    seed: int,
) -> dict[str, Any]:
    p_value_var = float(np.var(plain_value, axis=0, ddof=1).mean())
    c_value_var = float(np.var(mcv_value, axis=0, ddof=1).mean())
    p_grad_var = float(np.var(plain_z, axis=0, ddof=1).mean())
    c_grad_var = float(np.var(mcv_z, axis=0, ddof=1).mean())
    p_vector_var = float(np.var(plain_z, axis=0, ddof=1).sum(axis=1).mean())
    c_vector_var = float(np.var(mcv_z, axis=0, ddof=1).sum(axis=1).mean())
    value_ci = _bootstrap_ratio(plain_value, mcv_value, seed=seed)
    grad_ci = _bootstrap_ratio(plain_z, mcv_z, seed=seed + 1)

    pair = np.arange(0, len(plain_value), 2)
    p_value_rep = float(np.sqrt(np.mean(
        (plain_value[pair] - plain_value[pair + 1]) ** 2 / 2.0
    )))
    c_value_rep = float(np.sqrt(np.mean(
        (mcv_value[pair] - mcv_value[pair + 1]) ** 2 / 2.0
    )))
    p_grad_rep = float(np.sqrt(np.mean(
        (plain_z[pair] - plain_z[pair + 1]) ** 2 / 2.0
    )))
    c_grad_rep = float(np.sqrt(np.mean(
        (mcv_z[pair] - mcv_z[pair + 1]) ** 2 / 2.0
    )))
    p_vector_rep = float(np.sqrt(np.mean(np.sum(
        (plain_z[pair] - plain_z[pair + 1]) ** 2 / 2.0, axis=2
    ))))
    c_vector_rep = float(np.sqrt(np.mean(np.sum(
        (mcv_z[pair] - mcv_z[pair + 1]) ** 2 / 2.0, axis=2
    ))))

    # Test differentiated unbiasedness in stable random directions rather than
    # averaging signed coordinates, which can conceal coordinate-level bias.
    rng = np.random.default_rng(seed + 2)
    directions = rng.standard_normal(plain_z.shape[1:])
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    projected_difference = np.sum((mcv_z - plain_z) * directions[None], axis=2)
    value_difference = _paired_interval(mcv_value - plain_value)
    gradient_difference = _paired_interval(projected_difference)

    plain_noise = plain_z - plain_z.mean(axis=0, keepdims=True)
    mcv_noise = mcv_z - mcv_z.mean(axis=0, keepdims=True)
    plain_seconds = np.asarray([item["wall_seconds"] for item in plain_meta])
    mcv_seconds = np.asarray([item["wall_seconds"] for item in mcv_meta])
    plain_peak = [item.get("peak_gpu_bytes") or 0 for item in plain_meta]
    mcv_peak = [item.get("peak_gpu_bytes") or 0 for item in mcv_meta]
    return {
        "checkpoint": checkpoint,
        "paths": paths,
        "replicates": len(plain_value),
        "states": plain_value.shape[1],
        "value_variance_plain": p_value_var,
        "value_variance_mcv": c_value_var,
        "value_variance_ratio_mcv_over_plain": c_value_var / p_value_var,
        "value_variance_ratio_ci95_low": value_ci[0],
        "value_variance_ratio_ci95_high": value_ci[1],
        "gradient_component_variance_plain": p_grad_var,
        "gradient_component_variance_mcv": c_grad_var,
        "gradient_component_variance_ratio_mcv_over_plain": c_grad_var / p_grad_var,
        "gradient_component_variance_ratio_ci95_low": grad_ci[0],
        "gradient_component_variance_ratio_ci95_high": grad_ci[1],
        "gradient_vector_variance_plain": p_vector_var,
        "gradient_vector_variance_mcv": c_vector_var,
        "gradient_vector_variance_ratio_mcv_over_plain": c_vector_var / p_vector_var,
        "value_replication_rmse_plain": p_value_rep,
        "value_replication_rmse_mcv": c_value_rep,
        "gradient_component_replication_rmse_plain": p_grad_rep,
        "gradient_component_replication_rmse_mcv": c_grad_rep,
        "gradient_vector_replication_rmse_plain": p_vector_rep,
        "gradient_vector_replication_rmse_mcv": c_vector_rep,
        "gradient_raw_correlation_plain_mcv": _correlation(plain_z, mcv_z),
        "gradient_noise_correlation_plain_mcv": _correlation(plain_noise, mcv_noise),
        "paired_value_difference_mcv_minus_plain": value_difference[0],
        "paired_value_difference_se": value_difference[1],
        "paired_value_difference_ci95_low": value_difference[2],
        "paired_value_difference_ci95_high": value_difference[3],
        "paired_gradient_directional_difference_mcv_minus_plain": gradient_difference[0],
        "paired_gradient_directional_difference_se": gradient_difference[1],
        "paired_gradient_directional_difference_ci95_low": gradient_difference[2],
        "paired_gradient_directional_difference_ci95_high": gradient_difference[3],
        "plain_seconds_mean": float(np.mean(plain_seconds)),
        "mcv_seconds_mean": float(np.mean(mcv_seconds)),
        "runtime_ratio_mcv_over_plain": float(np.mean(mcv_seconds / plain_seconds)),
        "plain_peak_gpu_bytes": int(max(plain_peak)),
        "mcv_peak_gpu_bytes": int(max(mcv_peak)),
        "memory_ratio_mcv_over_plain": float(max(mcv_peak) / max(max(plain_peak), 1)),
        "weight_ess_fraction_plain_mean": float(np.mean([
            item["weight_ess_fraction"] for item in plain_meta
        ])),
        "weight_ess_fraction_mcv_mean": float(np.mean([
            item["weight_ess_fraction"] for item in mcv_meta
        ])),
        "martingale_path_mean_mean": float(np.mean([
            item["martingale_path_mean"] for item in mcv_meta
        ])),
        "nonfinite_plain": int(sum(item["nonfinite_count"] for item in plain_meta)),
        "nonfinite_mcv": int(sum(item["nonfinite_count"] for item in mcv_meta)),
    }


def _finite_difference_check(
    checkpoint: str,
    actor: torch.nn.Module,
    t: np.ndarray,
    x: np.ndarray,
    *,
    driver: Any,
    device: torch.device,
    seed: int,
) -> dict[str, Any]:
    count = min(8, len(t))
    local_t = t[:count]
    local_x = x[:count]
    rng = np.random.default_rng(seed)
    direction = rng.standard_normal(local_x.shape)
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    epsilon = 2e-3
    common = {
        "paths": 512,
        "time_steps": 16,
        "seed": seed + 1,
        "device": device,
        "parent_batch": 8,
        "path_chunk": 32,
    }
    _, z, metadata = mcv_tangent_labels(
        driver, actor, local_t, local_x, need_z=True, **common
    )
    plus, _, _ = mcv_tangent_labels(
        driver, actor, local_t, local_x + epsilon * direction,
        need_z=False, **common,
    )
    minus, _, _ = mcv_tangent_labels(
        driver, actor, local_t, local_x - epsilon * direction,
        need_z=False, **common,
    )
    assert z is not None
    automatic = np.sum((z / driver.sigma) * direction, axis=1)
    finite_difference = (plus - minus) / (2.0 * epsilon)
    absolute = np.abs(automatic - finite_difference)
    relative = absolute / np.maximum(
        np.maximum(np.abs(automatic), np.abs(finite_difference)), 1e-3
    )
    return {
        "checkpoint": checkpoint,
        "states": count,
        "paths": common["paths"],
        "time_steps": common["time_steps"],
        "epsilon": epsilon,
        "crn": True,
        "autodiff": automatic.tolist(),
        "finite_difference": finite_difference.tolist(),
        "absolute_error": absolute.tolist(),
        "relative_error": relative.tolist(),
        "median_absolute_error": float(np.median(absolute)),
        "maximum_absolute_error": float(np.max(absolute)),
        "median_relative_error": float(np.median(relative)),
        "maximum_relative_error": float(np.max(relative)),
        "complete_corrected_payoff_differentiated": metadata[
            "complete_corrected_payoff_differentiated"
        ],
    }


def run(args: argparse.Namespace) -> None:
    ensure_dirs()
    device = torch.device(args.device)
    benchmark = "nonseparable_d100_k2"
    seed = 0
    driver = driver_for(benchmark)
    t, x, fingerprint = validation_bank(benchmark, seed, count=args.states)
    checkpoint_paths = {
        "a3": a3_checkpoint(benchmark, seed),
        "late": late_checkpoint(benchmark, seed)[0],
    }
    rows: list[dict[str, Any]] = []
    finite_differences: list[dict[str, Any]] = []
    raw_manifest: list[dict[str, Any]] = []
    for checkpoint, path in checkpoint_paths.items():
        actor, payload = load_tangent_actor(path, device)
        for paths in args.paths:
            plain_values: list[np.ndarray] = []
            mcv_values: list[np.ndarray] = []
            plain_gradients: list[np.ndarray] = []
            mcv_gradients: list[np.ndarray] = []
            plain_metadata: list[dict[str, Any]] = []
            mcv_metadata: list[dict[str, Any]] = []
            for replicate in range(args.replicates):
                label_seed = stable_seed(
                    "mcv_training_stage0", benchmark, seed, checkpoint, paths, replicate
                )
                common = {
                    "paths": paths,
                    "time_steps": args.time_steps,
                    "seed": label_seed,
                    "device": device,
                    "parent_batch": args.parent_batch,
                    "path_chunk": min(args.path_chunk, paths),
                }
                # Alternate call order to avoid assigning all warm-cache effects
                # to one estimator. Brownian streams still match by seed.
                if replicate % 2 == 0:
                    pv, pz, pm = tangent_labels(driver, actor, t, x, **common)
                    cv, cz, cm = mcv_tangent_labels(driver, actor, t, x, **common)
                else:
                    cv, cz, cm = mcv_tangent_labels(driver, actor, t, x, **common)
                    pv, pz, pm = tangent_labels(driver, actor, t, x, **common)
                assert pz is not None and cz is not None
                plain_values.append(pv)
                mcv_values.append(cv)
                plain_gradients.append(pz)
                mcv_gradients.append(cz)
                plain_metadata.append(pm)
                mcv_metadata.append(cm)
                print(
                    f"stage0 {checkpoint} M={paths} rep={replicate + 1}/{args.replicates} "
                    f"plain={pm['wall_seconds']:.2f}s mcv={cm['wall_seconds']:.2f}s",
                    flush=True,
                )
            arrays = {
                "plain_value": np.asarray(plain_values),
                "mcv_value": np.asarray(mcv_values),
                "plain_z": np.asarray(plain_gradients),
                "mcv_z": np.asarray(mcv_gradients),
            }
            archive_path = LABEL_DIR / f"stage0_{checkpoint}_paths{paths}.npz"
            np.savez_compressed(archive_path, **arrays)
            summary = _summary_row(
                checkpoint, paths, arrays["plain_value"], arrays["mcv_value"],
                arrays["plain_z"], arrays["mcv_z"], plain_metadata, mcv_metadata,
                seed=stable_seed("mcv_training_bootstrap", checkpoint, paths),
            )
            rows.append({"benchmark": benchmark, "seed": seed, **summary})
            raw_manifest.append({
                "checkpoint": checkpoint,
                "checkpoint_path": str(path),
                "checkpoint_iteration": payload.get("iteration"),
                "paths": paths,
                "archive": str(archive_path),
                "arrays_sha256": array_sha256(*arrays.values()),
                "plain_metadata": plain_metadata,
                "mcv_metadata": mcv_metadata,
            })
        finite_differences.append(_finite_difference_check(
            checkpoint, actor, t, x, driver=driver, device=device,
            seed=stable_seed("mcv_training_finite_difference", checkpoint),
        ))

    write_csv(RESULTS_DIR / "label_variance_comparison.csv", rows)
    atomic_json(RAW_DIR / "stage0_raw_manifest.json", {
        "benchmark": benchmark,
        "seed": seed,
        "state_bank_fingerprint": fingerprint,
        "state_count": len(t),
        "time_steps": args.time_steps,
        "replicates": args.replicates,
        "entries": raw_manifest,
        "finite_difference": finite_differences,
    })
    decisive = [row for row in rows if row["paths"] == max(args.paths)]
    gradient_pass = all(
        row["gradient_component_variance_ratio_mcv_over_plain"] < 0.95
        and row["gradient_component_variance_ratio_ci95_high"] < 1.10
        for row in decisive
    )
    value_pass = all(
        row["value_variance_ratio_mcv_over_plain"] < 1.0 for row in decisive
    )
    stability_pass = all(
        row["nonfinite_plain"] == 0 and row["nonfinite_mcv"] == 0
        for row in decisive
    )
    derivative_pass = all(
        item["median_relative_error"] < 0.10
        and item["maximum_relative_error"] < 0.30
        for item in finite_differences
    )
    overhead_pass = all(row["runtime_ratio_mcv_over_plain"] < 1.50 for row in decisive)
    passed = gradient_pass and value_pass and stability_pass and derivative_pass and overhead_pass
    decision = {
        "stage": 0,
        "passed": passed,
        "decision": "PROCEED_STAGE1" if passed else "STOP_AFTER_STAGE0",
        "criteria_predeclared": {
            "gradient_variance_ratio_point_below": 0.95,
            "gradient_variance_ratio_ci95_high_below": 1.10,
            "value_variance_ratio_below": 1.0,
            "median_fd_relative_error_below": 0.10,
            "maximum_fd_relative_error_below": 0.30,
            "runtime_ratio_below": 1.50,
            "nonfinite_count": 0,
        },
        "checks": {
            "gradient_pass": gradient_pass,
            "value_pass": value_pass,
            "stability_pass": stability_pass,
            "derivative_pass": derivative_pass,
            "overhead_pass": overhead_pass,
        },
        "decisive_rows": decisive,
        "finite_difference": finite_differences,
    }
    atomic_json(RESULTS_DIR / "stage0_decision.json", decision)
    print(f"STAGE0_DECISION={decision['decision']}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--states", type=int, default=128)
    parser.add_argument("--paths", type=int, nargs="+", default=[32, 128])
    parser.add_argument("--time-steps", type=int, default=16)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--parent-batch", type=int, default=8)
    parser.add_argument("--path-chunk", type=int, default=32)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
