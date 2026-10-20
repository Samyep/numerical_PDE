"""Freeze B6 difficulty parameters before any model training."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import subprocess
import sys
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.gpu_mlp.solver import torch_seed
from .benchmarks import array_sha256, fixed_b6_points, make_b6_equation

from .protocol import (
    BASE_SEED,
    B6_BETAS,
    B6_FIXED_SELECTIONS,
    B6_GAMMAS,
    B6_HEAT_SAMPLES,
    CALIBRATION_PATH,
    DIMENSIONS,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    ensure_directories,
    frozen_protocol,
    write_json_atomic,
)


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=SOURCE_ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()


def _environment(device: torch.device) -> dict[str, Any]:
    result: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "device": str(device),
    }
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(device)
        result.update(
            {
                "gpu_name": properties.name,
                "gpu_total_memory_bytes": properties.total_memory,
                "gpu_compute_capability": f"{properties.major}.{properties.minor}",
            }
        )
    return result


def _heat_only_all_beta(
    d: int,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
    *,
    samples: int = B6_HEAT_SAMPLES,
    state_batch: int = 8,
) -> tuple[dict[float, np.ndarray], dict[str, Any]]:
    base = make_b6_equation(d, 1.0, 0.0)
    directions = torch.as_tensor(
        base.U[:, None] * base.e, dtype=torch.float64, device=device
    )
    state = torch.as_tensor(x, dtype=torch.float64, device=device)
    times = torch.as_tensor(t, dtype=torch.float64, device=device)
    generator = torch.Generator(device=device)
    generator.manual_seed(
        torch_seed(np.random.SeedSequence([BASE_SEED, d, 66, samples, 0x48454154]))
    )
    normal = torch.randn(
        (samples, d), dtype=torch.float64, device=device, generator=generator
    )
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with torch.no_grad():
        state_projection = state @ directions.T
        noise_projection = normal @ directions.T
        scale = torch.sqrt(2.0 * torch.clamp(0.25 - times, min=0.0))
        sums = {beta: torch.zeros(len(t), dtype=torch.float64, device=device) for beta in B6_BETAS}
        squares = {beta: torch.zeros(len(t), dtype=torch.float64, device=device) for beta in B6_BETAS}
        for begin in range(0, len(t), state_batch):
            end = min(begin + state_batch, len(t))
            projected = (
                state_projection[begin:end, None, :]
                + scale[begin:end, None, None] * noise_projection[None, :, :]
            )
            for beta in B6_BETAS:
                values = -(torch.logsumexp(float(beta) * projected, dim=-1) - math.log(base.K))
                sums[beta][begin:end] = values.sum(1)
                squares[beta][begin:end] = values.square().sum(1)
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        estimates = {beta: (sums[beta] / samples).cpu().numpy() for beta in B6_BETAS}
        standard_errors: dict[str, float] = {}
        for beta in B6_BETAS:
            mean = sums[beta] / samples
            variance = torch.clamp(
                (squares[beta] - samples * mean.square()) / (samples - 1), min=0.0
            )
            standard_errors[str(beta)] = float(
                torch.sqrt(variance / samples).mean().cpu()
            )
    return estimates, {
        "samples_per_state": samples,
        "common_random_numbers_across_states_and_configurations": True,
        "arithmetic": "float64",
        "mean_standard_error_by_beta": standard_errors,
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": (
            int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None
        ),
    }


def _global_geometry(z: np.ndarray) -> tuple[float, float]:
    centered = z - z.mean(0, keepdims=True)
    covariance = centered.T @ centered / len(centered)
    eigenvalues = np.maximum(np.linalg.eigvalsh(covariance), 0.0)
    trace = float(eigenvalues.sum())
    square = float(np.square(eigenvalues).sum())
    rank = trace * trace / max(square, np.finfo(float).tiny)
    outside = 1.0 - float(eigenvalues[-10:].sum()) / max(trace, np.finfo(float).tiny)
    return rank, outside


def _local_ranks(
    probabilities: np.ndarray, vectors: np.ndarray, device: torch.device
) -> np.ndarray:
    weights = torch.as_tensor(probabilities, dtype=torch.float64, device=device)
    b = torch.as_tensor(vectors, dtype=torch.float64, device=device)
    gram = b @ b.T
    gram_square = gram.square()
    ranks: list[torch.Tensor] = []
    with torch.no_grad():
        for begin in range(0, len(weights), 64):
            pi = weights[begin : begin + 64]
            pi_gram = pi @ gram
            first = ((pi @ gram_square) * pi).sum(1)
            second = (pi * pi_gram.square()).sum(1)
            mean_square = (pi * pi_gram).sum(1)
            trace = (pi * torch.diag(gram)).sum(1) - mean_square
            frobenius_square = torch.clamp(
                first - 2.0 * second + mean_square.square(), min=0.0
            )
            ranks.append(
                trace.square()
                / torch.clamp(frobenius_square, min=torch.finfo(torch.float64).tiny)
            )
    return torch.cat(ranks).cpu().numpy()


def _violation(row: dict[str, Any]) -> float:
    return float(
        max((5.0 - row["ESS_med"]) / 5.0, 0.0)
        + max((0.5 - row["z_rms"]) / 0.5, 0.0)
        + max((row["z_rms"] - 4.0) / 4.0, 0.0)
        + max((0.15 - row["heat_only_skill"]) / 0.15, 0.0)
        + max((0.05 - row["f_cv"]) / 0.05, 0.0)
    )


def calibrate_dimension(
    d: int, device: torch.device
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    t_all, x_all, is_validation = fixed_b6_points(int(d))
    t, x = t_all[~is_validation], x_all[~is_validation]
    if len(t) != 960:
        raise AssertionError(f"expected 960 calibration/test states, found {len(t)}")
    heat_by_beta, heat_metadata = _heat_only_all_beta(d, t, x, device)
    rows: list[dict[str, Any]] = []
    for beta in B6_BETAS:
        for gamma in B6_GAMMAS:
            started = time.perf_counter()
            equation = make_b6_equation(d, beta, gamma)
            exact = equation.exact_state(t, x)
            value, z = exact[:, 0], exact[:, 1:]
            probabilities = equation.probabilities(t, x)
            entropy = -np.sum(
                probabilities * np.log(np.maximum(probabilities, np.finfo(float).tiny)),
                axis=1,
            )
            generator = -0.5 * np.sum(z * z, axis=1)
            heat = heat_by_beta[beta] - gamma * (x @ equation.v)
            rank, outside = _global_geometry(z)
            local = _local_ranks(probabilities, equation.b, device)
            row: dict[str, Any] = {
                "pde": "b6_fullrank",
                "dimension": d,
                "beta": beta,
                "gamma": gamma,
                "K": equation.K,
                "test_states": len(t),
                "full_bank_sha256": array_sha256(t_all, x_all, is_validation),
                "test_bank_sha256": array_sha256(t, x),
                "ESS_med": float(np.median(np.exp(entropy))),
                "z_rms": float(np.sqrt(np.mean(np.sum(z * z, axis=1)))),
                "f_mean": float(generator.mean()),
                "f_std": float(generator.std()),
                "f_cv": float(generator.std() / max(abs(float(generator.mean())), np.finfo(float).tiny)),
                "heat_only_rmse": float(np.sqrt(np.mean(np.square(heat - value)))),
                "u_std": float(value.std()),
                "heat_only_skill": float(
                    np.sqrt(np.mean(np.square(heat - value)))
                    / max(float(value.std()), np.finfo(float).tiny)
                ),
                "r_glob": rank,
                "z_variance_outside_top10": outside,
                "r_loc": float(np.median(local)),
                "r_loc_min": float(local.min()),
                "r_loc_max": float(local.max()),
                "unequal_norm_range": float(equation.n2.max() - equation.n2.min()),
                "configuration_seconds_excluding_shared_heat": time.perf_counter() - started,
            }
            row["qualifies"] = bool(
                row["ESS_med"] >= 5.0
                and 0.5 <= row["z_rms"] <= 4.0
                and row["heat_only_skill"] >= 0.15
                and row["f_cv"] >= 0.05
            )
            row["normalized_constraint_violation"] = _violation(row)
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    eligible = [row for row in rows if row["qualifies"]]
    if eligible:
        selected = sorted(eligible, key=lambda row: (-row["r_glob"], row["beta"], row["gamma"]))[0]
        status = "qualified"
    else:
        selected = sorted(
            rows,
            key=lambda row: (row["normalized_constraint_violation"], row["beta"], row["gamma"]),
        )[0]
        status = "calibration_limited"
    if d in B6_FIXED_SELECTIONS:
        expected = B6_FIXED_SELECTIONS[d]
        if selected["beta"] != expected["beta"] or selected["gamma"] != expected["gamma"]:
            raise RuntimeError(
                f"B6 d={d} calibration drifted: selected={selected}, expected={expected}"
            )
    return rows, {
        "dimension": d,
        "status": status,
        "selected_beta": selected["beta"],
        "selected_gamma": selected["gamma"],
        "selected_r_glob": selected["r_glob"],
        "selected_row": selected,
        "heat_monte_carlo": heat_metadata,
        "test_bank_sha256": selected["test_bank_sha256"],
        "existing_selection_lock_verified": d in B6_FIXED_SELECTIONS,
    }


def calibrate(device_name: str = "cuda") -> dict[str, Any]:
    ensure_directories()
    device = torch.device(device_name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    if CALIBRATION_PATH.is_file():
        existing = json.loads(CALIBRATION_PATH.read_text(encoding="utf-8"))
        if existing.get("protocol_sha256") != frozen_protocol()["protocol_sha256"]:
            raise RuntimeError("existing B6 calibration belongs to another protocol")
        if set(existing.get("selections", {})) == {str(d) for d in DIMENSIONS}:
            return existing
        raise RuntimeError("partial B6 calibration exists; refusing an ambiguous overwrite")
    rows: list[dict[str, Any]] = []
    selections: dict[str, Any] = {}
    started = time.perf_counter()
    for d in DIMENSIONS:
        current, selection = calibrate_dimension(d, device)
        rows.extend(current)
        selections[str(d)] = selection
    table = RESULTS_DIR / "b6_calibration.csv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "revision": REVISION,
        "protocol_sha256": frozen_protocol()["protocol_sha256"],
        "status": "FROZEN_BEFORE_TRAINING",
        "code_commit_at_calibration": _git("rev-parse", "HEAD"),
        "selection_rule": frozen_protocol()["b6_calibration"]["selection_rule"],
        "selections": selections,
        "grid_rows": len(rows),
        "elapsed_seconds": time.perf_counter() - started,
        "environment": _environment(device),
        "production_solver_touched_b6": False,
        "source_calibration_commit": "281977883e0a5ce4f05895878680b582620461be",
    }
    write_json_atomic(CALIBRATION_PATH, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    print(json.dumps(calibrate(args.device), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
