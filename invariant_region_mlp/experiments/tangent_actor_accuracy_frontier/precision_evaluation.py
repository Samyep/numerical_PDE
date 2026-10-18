"""Common-random-number Euler audit of frozen Tangent-Actor checkpoints."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.core import (
    antithetic_normals,
    frozen,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    _weighted_rollout_chunk_mcv_labels,
)
from invariant_region_mlp.experiments.tangent_actor_scaling_diagnosis.training import (
    FitConfig,
)

from .common import (
    OLD_CHECKPOINT_DIR,
    OLD_PREDICTION_DIR,
    OLD_REFERENCE_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    REFERENCE_DIR,
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    atomic_npz,
    ensure_dirs,
    environment_payload,
    stable_seed,
    synchronize,
)
from .drivers import SineQuadraticDriver
from .reference2d import load_reference
from .solver import load_actor


STEPS = (64, 128, 256, 512)


def _main_config() -> FitConfig:
    source = Path(__file__).parents[1] / "nonconvex_tangent_actor" / "protocol.json"
    payload = json.loads(source.read_text(encoding="utf-8"))["main"]
    allowed = set(FitConfig.__dataclass_fields__)
    return FitConfig(**{key: value for key, value in payload.items() if key in allowed})


def _old_run_id(seed: int) -> str:
    return f"primary_d100_k2_lam0p2_eta0_main_mcv_s{seed}"


def _load_frozen_primary(seed: int, device: torch.device):
    driver = SineQuadraticDriver(d=100, lam=0.2, eta=0.0)
    path = OLD_CHECKPOINT_DIR / f"{_old_run_id(seed)}_mcv_a4.pt"
    if not path.exists():
        raise FileNotFoundError(path)
    return driver, load_actor(path, driver, _main_config(), device), path


def _test_bank(driver: SineQuadraticDriver) -> tuple[np.ndarray, np.ndarray]:
    path = OLD_PREDICTION_DIR / f"post_freeze_{_old_run_id(0)}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path) as archive:
        t = archive["test_t"].astype(np.float64)
        x = archive["test_x"].astype(np.float64)
    if x.shape != (512, driver.d):
        raise RuntimeError(f"unexpected frozen test bank {x.shape}")
    return t, x


def _reference(driver: SineQuadraticDriver):
    refined = REFERENCE_DIR / f"{driver.name}_refined.npz"
    path = refined if refined.exists() else OLD_REFERENCE_DIR / f"{driver.name}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    metadata_path = REFERENCE_DIR / "primary_refinement_audit.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        uncertainty = float(metadata["reference_uncertainty_value_rmse"])
    else:
        preflight = json.loads(
            (OLD_REFERENCE_DIR.parent / "preflight_audit.json").read_text(encoding="utf-8")
        )
        uncertainty = float(
            preflight["drivers"]["primary"]["reference"][
                "reference_error_proxy_value_rmse"
            ]
        )
    return load_reference(path), path, uncertainty


def _aggregate_normals(normals: torch.Tensor, steps: int) -> torch.Tensor:
    fine_steps = normals.shape[2]
    if fine_steps % steps:
        raise ValueError((fine_steps, steps))
    ratio = fine_steps // steps
    if ratio == 1:
        return normals
    return normals.reshape(
        normals.shape[0], normals.shape[1], steps, ratio, normals.shape[3]
    ).sum(dim=3) / math.sqrt(ratio)


def _evaluate_level(
    driver,
    actor,
    t: np.ndarray,
    x: np.ndarray,
    *,
    paths: int,
    steps: int,
    fine_steps: int,
    seed: int,
    parent_batch: int,
    path_chunk: int,
    device: torch.device,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    if fine_steps % steps or paths % path_chunk or path_chunk % 2:
        raise ValueError("invalid coupled-refinement dimensions")
    pairs = paths // 2
    mcv_pairs = np.empty((len(t), pairs), dtype=np.float64)
    plain_pairs = np.empty_like(mcv_pairs)
    generator = torch.Generator(device=device).manual_seed(seed)
    weight_sum = weight_square_sum = 0.0
    weight_min, weight_max = float("inf"), -float("inf")
    nonfinite = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with frozen(actor):
        for begin in range(0, len(t), parent_batch):
            end = min(begin + parent_batch, len(t))
            bt = torch.as_tensor(t[begin:end], dtype=torch.float32, device=device)
            bx = torch.as_tensor(x[begin:end], dtype=torch.float32, device=device)
            for path_begin in range(0, paths, path_chunk):
                fine_normals = antithetic_normals(
                    (end - begin, path_chunk, fine_steps, driver.d),
                    generator=generator,
                    device=device,
                )
                normals = _aggregate_normals(fine_normals, steps)
                corrected, _, weight, martingale, diagnostics = (
                    _weighted_rollout_chunk_mcv_labels(
                        driver, actor, bt, bx, normals, need_z=False
                    )
                )
                plain = corrected + martingale
                half = path_chunk // 2
                pair_slice = slice(path_begin // 2, (path_begin + path_chunk) // 2)
                mcv_pairs[begin:end, pair_slice] = (
                    0.5 * (corrected[:, :half] + corrected[:, half:])
                ).double().cpu().numpy()
                plain_pairs[begin:end, pair_slice] = (
                    0.5 * (plain[:, :half] + plain[:, half:])
                ).double().cpu().numpy()
                weight_sum += float(weight.double().sum().item())
                weight_square_sum += float(weight.double().square().sum().item())
                weight_min = min(weight_min, float(weight.min().item()))
                weight_max = max(weight_max, float(weight.max().item()))
                nonfinite += int(diagnostics["nonfinite_count"])
                del fine_normals, normals, corrected, martingale, plain, weight
    synchronize(device)
    wall = time.perf_counter() - started
    if nonfinite or not np.all(np.isfinite(mcv_pairs)):
        raise FloatingPointError("nonfinite precision evaluation")
    count = len(t) * paths
    return {
        "mcv_pairs": mcv_pairs,
        "plain_pairs": plain_pairs,
        "mcv": mcv_pairs.mean(axis=1),
        "plain": plain_pairs.mean(axis=1),
        "mcv_pair_variance": mcv_pairs.var(axis=1, ddof=1),
        "plain_pair_variance": plain_pairs.var(axis=1, ddof=1),
    }, {
        "paths": paths,
        "antithetic_pairs": pairs,
        "steps": steps,
        "fine_steps_used_for_crn": fine_steps,
        "seed": seed,
        "wall_seconds": wall,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda" else None,
        "weight_ess_fraction": weight_sum**2 / max(
            count * weight_square_sum, np.finfo(float).tiny
        ),
        "weight_min": weight_min,
        "weight_max": weight_max,
        "nonfinite_count": nonfinite,
    }


def _bootstrap_ci(
    error: np.ndarray,
    mc_variance: np.ndarray,
    reference_uncertainty: float,
    *,
    seed: int,
    replicates: int,
) -> dict[str, list[float]]:
    rng = np.random.default_rng(seed)
    n = len(error)
    rmse = np.empty(replicates, dtype=np.float64)
    residual = np.empty(replicates, dtype=np.float64)
    bias = np.empty(replicates, dtype=np.float64)
    corrected_square = np.square(error) - mc_variance - reference_uncertainty**2
    for start in range(0, replicates, 250):
        count = min(250, replicates - start)
        index = rng.integers(0, n, (count, n))
        sampled = error[index]
        rmse[start:start + count] = np.sqrt(np.mean(np.square(sampled), axis=1))
        residual[start:start + count] = np.sqrt(
            np.maximum(np.mean(corrected_square[index], axis=1), 0.0)
        )
        bias[start:start + count] = np.mean(sampled, axis=1)
    point_rmse = float(np.sqrt(np.mean(np.square(error))))
    mse_mc_se = float(
        np.sqrt(np.sum(4.0 * np.square(error) * mc_variance + 2.0 * np.square(mc_variance)))
        / n
    )
    rmse_mc_delta_se = mse_mc_se / max(2.0 * point_rmse, np.finfo(float).tiny)
    raw_interval = np.quantile(rmse, [0.025, 0.975])
    raw_interval[0] = max(
        0.0, raw_interval[0] - 1.96 * rmse_mc_delta_se - reference_uncertainty
    )
    raw_interval[1] += 1.96 * rmse_mc_delta_se + reference_uncertainty
    bias_mc_se = float(np.sqrt(np.sum(mc_variance)) / n)
    bias_interval = np.quantile(bias, [0.025, 0.975])
    bias_interval[0] -= 1.96 * bias_mc_se + reference_uncertainty
    bias_interval[1] += 1.96 * bias_mc_se + reference_uncertainty
    return {
        "rmse_95": raw_interval.tolist(),
        "residual_rmse_95": np.quantile(residual, [0.025, 0.975]).tolist(),
        "bias_95": bias_interval.tolist(),
        "rmse_mc_delta_se": rmse_mc_delta_se,
    }


def _row(
    prediction: dict[str, np.ndarray],
    truth: np.ndarray,
    details: dict[str, Any],
    *,
    seed: int,
    reference_uncertainty: float,
    bootstrap_replicates: int,
) -> dict[str, Any]:
    error = prediction["mcv"] - truth
    mse = float(np.mean(np.square(error)))
    pairs = int(details["antithetic_pairs"])
    mc_variance = prediction["mcv_pair_variance"] / pairs
    mc_mse = float(np.mean(mc_variance))
    ref_mse = reference_uncertainty**2
    residual = max(mse - mc_mse - ref_mse, 0.0)
    plain_error = prediction["plain"] - truth
    ci = _bootstrap_ci(
        error,
        mc_variance,
        reference_uncertainty,
        seed=seed,
        replicates=bootstrap_replicates,
    )
    result = {
        **details,
        "value_rmse": math.sqrt(mse),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(truth)),
        "mean_bias": float(np.mean(error)),
        "mc_rmse_contribution": math.sqrt(mc_mse),
        "reference_uncertainty_rmse": reference_uncertainty,
        "residual_non_mc_non_reference_rmse_proxy": math.sqrt(residual),
        "plain_value_rmse": float(np.sqrt(np.mean(np.square(plain_error)))),
        "mcv_variance_reduction": float(
            np.mean(prediction["plain_pair_variance"])
            / max(np.mean(prediction["mcv_pair_variance"]), 1e-300)
        ),
        "rmse_ci95_low": ci["rmse_95"][0],
        "rmse_ci95_high": ci["rmse_95"][1],
        "bias_ci95_low": ci["bias_95"][0],
        "bias_ci95_high": ci["bias_95"][1],
        "residual_rmse_ci95_low": ci["residual_rmse_95"][0],
        "residual_rmse_ci95_high": ci["residual_rmse_95"][1],
        "rmse_mc_delta_se": ci["rmse_mc_delta_se"],
        "confidence_method": "state bootstrap widened by delta-method MC SE and conservative reference RMSE",
    }
    return result


def run(
    training_seed: int,
    *,
    paths: int,
    levels: tuple[int, ...],
    fine_steps: int,
    parent_batch: int,
    path_chunk: int,
    bootstrap_replicates: int,
    device: torch.device,
) -> dict[str, Any]:
    ensure_dirs()
    driver, actor, checkpoint = _load_frozen_primary(training_seed, device)
    t, x = _test_bank(driver)
    reference, reference_path, reference_uncertainty = _reference(driver)
    truth = reference.value(t, x @ driver.B)
    common_seed = stable_seed(
        "accuracy_frontier_independent_precision_v1", training_seed, paths, fine_steps
    )
    rows: list[dict[str, Any]] = []
    predictions: dict[int, dict[str, np.ndarray]] = {}
    total_started = time.perf_counter()
    for steps in levels:
        output, details = _evaluate_level(
            driver,
            actor,
            t,
            x,
            paths=paths,
            steps=steps,
            fine_steps=fine_steps,
            seed=common_seed,
            parent_batch=parent_batch,
            path_chunk=path_chunk,
            device=device,
        )
        row = _row(
            output,
            truth,
            details,
            seed=stable_seed("precision_bootstrap", training_seed, steps),
            reference_uncertainty=reference_uncertainty,
            bootstrap_replicates=bootstrap_replicates,
        )
        row.update({"training_seed": training_seed, "actor": "a4"})
        rows.append(row)
        predictions[steps] = output
        print(
            f"seed={training_seed} steps={steps} RMSE={row['value_rmse']:.7g} "
            f"residual={row['residual_non_mc_non_reference_rmse_proxy']:.7g} "
            f"MC={row['mc_rmse_contribution']:.7g} ESS={row['weight_ess_fraction']:.4f}",
            flush=True,
        )

    refinement_rows = []
    for coarse, fine in zip(levels[:-1], levels[1:], strict=True):
        coarse_output, fine_output = predictions[coarse], predictions[fine]
        difference = coarse_output["mcv"] - fine_output["mcv"]
        pair_difference = coarse_output["mcv_pairs"] - fine_output["mcv_pairs"]
        mc_mse = float(np.mean(pair_difference.var(axis=1, ddof=1) / (paths // 2)))
        rms = float(np.sqrt(np.mean(np.square(difference))))
        corrected = math.sqrt(max(rms**2 - mc_mse, 0.0))
        refinement_rows.append(
            {
                "coarse_steps": coarse,
                "fine_steps": fine,
                "coupled_prediction_difference_rms": rms,
                "coupled_difference_mc_rmse": math.sqrt(mc_mse),
                "noise_corrected_weak_difference_rms": corrected,
            }
        )
    for index in range(len(refinement_rows) - 1):
        first = refinement_rows[index]["noise_corrected_weak_difference_rms"]
        second = refinement_rows[index + 1]["noise_corrected_weak_difference_rms"]
        order = None if first <= 0.0 or second <= 0.0 else math.log(first / second, 2.0)
        refinement_rows[index]["observed_order_using_next_difference"] = order
    if refinement_rows:
        refinement_rows[-1]["observed_order_using_next_difference"] = None

    prediction_path = PREDICTION_DIR / (
        f"primary_a4_seed{training_seed}_crn_paths{paths}_fine{fine_steps}.npz"
    )
    arrays: dict[str, np.ndarray] = {"test_t": t, "test_x": x, "truth": truth}
    for steps, output in predictions.items():
        arrays[f"mcv_{steps}"] = output["mcv"]
        arrays[f"plain_{steps}"] = output["plain"]
        arrays[f"mcv_pair_variance_{steps}"] = output["mcv_pair_variance"]
    atomic_npz(prediction_path, **arrays)
    result = {
        "study": "follow-up precision experiment",
        "training_seed": training_seed,
        "actor": "a4",
        "checkpoint": str(checkpoint),
        "checkpoint_selection_reused_unchanged": True,
        "reference_path": str(reference_path),
        "reference_uncertainty_value_rmse": reference_uncertainty,
        "test_bank_sha256": array_sha256(t, x),
        "independent_from_prior_euler_stream": True,
        "common_random_numbers_across_levels": True,
        "rows": rows,
        "refinement_rows": refinement_rows,
        "prediction_path": str(prediction_path),
        "total_wall_seconds": time.perf_counter() - total_started,
        "ground_truth_used_for_training_or_selection": False,
        "configuration": {
            "paths": paths,
            "levels": levels,
            "fine_steps": fine_steps,
            "parent_batch": parent_batch,
            "path_chunk": path_chunk,
            "bootstrap_replicates": bootstrap_replicates,
        },
    }
    atomic_json(
        RAW_DIR / f"precision_primary_a4_seed{training_seed}_paths{paths}.json",
        result,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-seed", type=int, required=True)
    parser.add_argument("--paths", type=int, default=2048)
    parser.add_argument("--levels", type=int, nargs="+", default=list(STEPS))
    parser.add_argument("--fine-steps", type=int, default=512)
    parser.add_argument("--parent-batch", type=int, default=32)
    parser.add_argument("--path-chunk", type=int, default=128)
    parser.add_argument("--bootstrap-replicates", type=int, default=4000)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation is required")
    ensure_dirs()
    atomic_json(
        RAW_DIR / f"environment_precision_seed{args.training_seed}_paths{args.paths}.json",
        environment_payload("precision_evaluation"),
    )
    run(
        args.training_seed,
        paths=args.paths,
        levels=tuple(args.levels),
        fine_steps=args.fine_steps,
        parent_batch=args.parent_batch,
        path_chunk=args.path_chunk,
        bootstrap_replicates=args.bootstrap_replicates,
        device=device,
    )


if __name__ == "__main__":
    main()
