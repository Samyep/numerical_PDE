"""Reference-only effective-dimension and nonlinear-skill characterization."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from typing import Any

import numpy as np
import torch

from .benchmarks import array_sha256, make_benchmark, make_dense_fourier
from .prepare import load_test_bank
from .protocol import (
    BASE_SEED,
    CACHE_DIR,
    CALIBRATION_PATH,
    CHARACTERIZATION_PATH,
    DIMENSIONS,
    PDES,
    REVISION,
    ensure_directories,
    stable_seed,
)


def _geometry(z: np.ndarray) -> tuple[float, float, np.ndarray]:
    centered = z - z.mean(0, keepdims=True)
    covariance = centered.T @ centered / len(centered)
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0)
    trace = float(eigenvalues.sum())
    square = float(np.square(eigenvalues).sum())
    rank = trace * trace / max(square, np.finfo(float).tiny)
    outside = 1.0 - float(eigenvalues[-10:].sum()) / max(trace, np.finfo(float).tiny)
    return rank, outside, eigenvalues


def _dense_heat(
    d: int,
    t: np.ndarray,
    x: np.ndarray,
    *,
    pairs: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    ensure_directories()
    cache = CACHE_DIR / f"dense_fourier_d{d}_heat_antithetic_pairs{pairs}.npz"
    if cache.is_file():
        with np.load(cache) as data:
            return data["mean"], data["standard_error"], json.loads(str(data["metadata"]))
    driver = make_dense_fourier(d)
    omega = torch.as_tensor(driver.omega, dtype=torch.float64, device=device)
    signs = torch.as_tensor(driver.signs, dtype=torch.float64, device=device)
    phases = torch.as_tensor(driver.phases, dtype=torch.float64, device=device)
    state = torch.as_tensor(x, dtype=torch.float64, device=device)
    times = torch.as_tensor(t, dtype=torch.float64, device=device)
    generator = torch.Generator(device=device).manual_seed(
        stable_seed(REVISION, "dense_heat", d, pairs, BASE_SEED)
    )
    normal = torch.randn((pairs, d), dtype=torch.float64, device=device, generator=generator)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
        torch.cuda.synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        noise_projection = normal @ omega.T
        state_phase = state @ omega.T + phases
        scale = torch.sqrt(2.0 * torch.clamp(driver.T - times, min=0.0))
        means: list[torch.Tensor] = []
        errors: list[torch.Tensor] = []
        for begin in range(0, len(t), 2):
            end = min(begin + 2, len(t))
            base = state_phase[begin:end, None, :]
            perturbation = scale[begin:end, None, None] * noise_projection[None, :, :]
            plus = torch.asinh(
                0.5
                + driver.scale
                * torch.sum(signs[None, None, :] * torch.sin(base + perturbation), dim=-1)
            )
            minus = torch.asinh(
                0.5
                + driver.scale
                * torch.sum(signs[None, None, :] * torch.sin(base - perturbation), dim=-1)
            )
            pair_mean = 0.5 * (plus + minus)
            means.append(pair_mean.mean(1))
            errors.append(pair_mean.std(1, unbiased=True) / math.sqrt(pairs))
        mean = torch.cat(means).cpu().numpy()
        standard_error = torch.cat(errors).cpu().numpy()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    metadata = {
        "kind": "antithetic common-random-number heat semigroup Monte Carlo",
        "pairs_per_state": pairs,
        "terminal_evaluations_per_state": 2 * pairs,
        "arithmetic": "float64",
        "seed": stable_seed(REVISION, "dense_heat", d, pairs, BASE_SEED),
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "test_bank_sha256": array_sha256(t, x),
    }
    np.savez_compressed(
        cache,
        mean=mean,
        standard_error=standard_error,
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    return mean, standard_error, metadata


def characterize(device_name: str = "cuda", dense_heat_pairs: int = 4096) -> list[dict[str, Any]]:
    if dense_heat_pairs < 2:
        raise ValueError("dense_heat_pairs must be at least two")
    device = torch.device(device_name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    calibration = json.loads(CALIBRATION_PATH.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for pde in PDES:
        for d in DIMENSIONS:
            benchmark = make_benchmark(pde, d)
            t, x, truth, truth_z, _ = load_test_bank(pde, d)
            rank, outside, eigenvalues = _geometry(truth_z)
            if pde == "b6_fullrank":
                selected = calibration["selections"][str(d)]["selected_row"]
                nonlinear_rmse = float(selected["heat_only_rmse"])
                nonlinear_skill = float(selected["heat_only_skill"])
                heat_uncertainty = float(
                    calibration["selections"][str(d)]["heat_monte_carlo"][
                        "mean_standard_error_by_beta"
                    ][str(float(selected["beta"]))]
                )
                heat_meta = calibration["selections"][str(d)]["heat_monte_carlo"]
                driver_value = -0.5 * np.sum(np.square(truth_z), axis=1)
            else:
                heat, heat_se, heat_meta = _dense_heat(
                    d, t, x, pairs=dense_heat_pairs, device=device
                )
                nonlinear_rmse = float(np.sqrt(np.mean(np.square(truth - heat))))
                nonlinear_skill = nonlinear_rmse / max(float(np.std(truth)), np.finfo(float).tiny)
                heat_uncertainty = float(np.sqrt(np.mean(np.square(heat_se))))
                driver_value = 0.5 * np.tanh(truth) * np.sum(np.square(truth_z), axis=1)
            rows.append(
                {
                    "revision": REVISION,
                    "pde": pde,
                    "dimension": d,
                    "test_states": len(t),
                    "test_bank_sha256": array_sha256(t, x),
                    "value_mean": float(np.mean(truth)),
                    "value_std": float(np.std(truth)),
                    "signal_variance": float(np.var(truth)),
                    "gradient_energy": float(np.mean(np.sum(np.square(truth_z), axis=1))),
                    "gradient_rms": float(np.sqrt(np.mean(np.sum(np.square(truth_z), axis=1)))),
                    "effective_rank": rank,
                    "variance_outside_top10": outside,
                    "nonlinear_heat_rmse": nonlinear_rmse,
                    "S_nonlinear": nonlinear_skill,
                    "heat_reference_mc_se_rms": heat_uncertainty,
                    "driver_mean": float(np.mean(driver_value)),
                    "driver_std": float(np.std(driver_value)),
                    "driver_rms": float(np.sqrt(np.mean(np.square(driver_value)))),
                    "covariance_trace": float(eigenvalues.sum()),
                    "covariance_top_eigenvalue": float(eigenvalues[-1]),
                    "heat_metadata": json.dumps(heat_meta, sort_keys=True),
                    "ground_truth_used_for_training_or_selection": False,
                }
            )
            print(json.dumps(rows[-1], sort_keys=True), flush=True)
    with CHARACTERIZATION_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dense-heat-pairs", type=int, default=4096)
    args = parser.parse_args()
    characterize(args.device, args.dense_heat_pairs)


if __name__ == "__main__":
    main()
