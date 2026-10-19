"""Reference-only B6 calibration, executed before any production solver run."""

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

from .equations import array_sha256, fixed_points, make_b6_equation
from .protocol import (
    BASE_SEED,
    B6_BETAS,
    B6_GAMMAS,
    B6_HEAT_SAMPLES,
    CALIBRATION_FREEZE_PATH,
    DIMENSIONS,
    RESULTS_ROOT,
    REVISION,
    ensure_directories,
    write_json_atomic,
)


def _git(*args: str) -> str:
    root = RESULTS_ROOT.parents[1]
    return subprocess.run(
        ["git", *args], cwd=root, check=True, capture_output=True, text=True
    ).stdout.strip()


def _environment(device: torch.device) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "device": str(device),
    }
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(device)
        payload.update(
            {
                "gpu_name": properties.name,
                "gpu_total_memory_bytes": properties.total_memory,
                "gpu_compute_capability": f"{properties.major}.{properties.minor}",
            }
        )
    return payload


def _heat_only_all_beta(
    d: int,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
    *,
    samples: int = B6_HEAT_SAMPLES,
    state_batch: int = 8,
) -> tuple[dict[float, np.ndarray], dict[str, Any]]:
    """Common-random-number MC for E[g(X_T)] for every registered beta.

    The same 4,096 Gaussian draws are used for all states/configurations.  This
    is still 4,096 terminal evaluations per state and reduces calibration
    noise in comparisons without changing any marginal estimator.
    """

    base = make_b6_equation(d, 1.0, 0.0)
    directions = torch.tensor(
        base.U[:, None] * base.e, dtype=torch.float64, device=device
    )
    state = torch.tensor(x, dtype=torch.float64, device=device)
    time_tensor = torch.tensor(t, dtype=torch.float64, device=device)
    generator = torch.Generator(device=device)
    generator.manual_seed(
        torch_seed(np.random.SeedSequence([BASE_SEED, d, 66, samples, 0x48454154]))
    )
    normal = torch.randn(
        (samples, d), dtype=torch.float64, device=device, generator=generator
    )
    with torch.no_grad():
        state_projection = torch.matmul(state, directions.T)
        noise_projection = torch.matmul(normal, directions.T)
        scale = torch.sqrt(2.0 * torch.clamp(0.25 - time_tensor, min=0.0))
        sums = {
            beta: torch.zeros(len(t), dtype=torch.float64, device=device)
            for beta in B6_BETAS
        }
        square_sums = {
            beta: torch.zeros(len(t), dtype=torch.float64, device=device)
            for beta in B6_BETAS
        }
        started = time.perf_counter()
        for begin in range(0, len(t), state_batch):
            end = min(begin + state_batch, len(t))
            projected = (
                state_projection[begin:end, None, :]
                + scale[begin:end, None, None] * noise_projection[None, :, :]
            )
            for beta in B6_BETAS:
                values = -(
                    torch.logsumexp(float(beta) * projected, dim=-1)
                    - math.log(base.K)
                )
                sums[beta][begin:end] = torch.sum(values, dim=1)
                square_sums[beta][begin:end] = torch.sum(values.square(), dim=1)
        if device.type == "cuda":
            torch.cuda.synchronize(device)
        elapsed = time.perf_counter() - started
        estimates = {
            beta: (sums[beta] / samples).cpu().numpy() for beta in B6_BETAS
        }
        standard_errors = {}
        for beta in B6_BETAS:
            mean = sums[beta] / samples
            variance = torch.clamp(
                (square_sums[beta] - samples * mean.square()) / (samples - 1),
                min=0.0,
            )
            standard_errors[str(beta)] = float(
                torch.sqrt(variance / samples).mean().cpu().item()
            )
    return estimates, {
        "samples_per_state": samples,
        "common_random_numbers_across_states_and_configurations": True,
        "arithmetic": "float64",
        "mean_standard_error_by_beta": standard_errors,
        "wall_seconds": elapsed,
    }


def _global_geometry(z: np.ndarray) -> tuple[float, float]:
    centered = z - np.mean(z, axis=0, keepdims=True)
    covariance = centered.T @ centered / len(centered)
    eigenvalues = np.linalg.eigvalsh(covariance)
    eigenvalues = np.maximum(eigenvalues, 0.0)
    trace = float(np.sum(eigenvalues))
    square = float(np.sum(eigenvalues**2))
    rank = trace**2 / max(square, np.finfo(np.float64).tiny)
    outside = 1.0 - float(np.sum(eigenvalues[-10:])) / max(
        trace, np.finfo(np.float64).tiny
    )
    return rank, outside


def _local_ranks(
    probabilities: np.ndarray, vectors: np.ndarray, device: torch.device
) -> np.ndarray:
    """Effective rank of Cov_pi(b), evaluated without materializing Hessians."""

    weights = torch.tensor(probabilities, dtype=torch.float64, device=device)
    b = torch.tensor(vectors, dtype=torch.float64, device=device)
    gram = torch.matmul(b, b.T)
    gram_square = gram.square()
    ranks: list[torch.Tensor] = []
    with torch.no_grad():
        for begin in range(0, len(weights), 64):
            pi = weights[begin : begin + 64]
            pi_gram = torch.matmul(pi, gram)
            first = torch.sum(torch.matmul(pi, gram_square) * pi, dim=1)
            second = torch.sum(pi * pi_gram.square(), dim=1)
            mean_square = torch.sum(pi * pi_gram, dim=1)
            trace = torch.sum(pi * torch.diag(gram), dim=1) - mean_square
            frobenius_square = torch.clamp(
                first - 2.0 * second + mean_square.square(), min=0.0
            )
            ranks.append(
                trace.square()
                / torch.clamp(frobenius_square, min=torch.finfo(torch.float64).tiny)
            )
    return torch.cat(ranks).cpu().numpy()


def _violation(row: dict[str, Any]) -> float:
    values = (
        max((5.0 - row["ESS_med"]) / 5.0, 0.0),
        max((0.5 - row["z_rms"]) / 0.5, 0.0),
        max((row["z_rms"] - 4.0) / 4.0, 0.0),
        max((0.15 - row["heat_only_skill"]) / 0.15, 0.0),
        max((0.05 - row["f_cv"]) / 0.05, 0.0),
    )
    return float(sum(values))


def calibrate_dimension(d: int, device: torch.device) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    t_all, x_all, is_validation = fixed_points("B6", d)
    t, x = t_all[~is_validation], x_all[~is_validation]
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
            global_rank, outside_top10 = _global_geometry(z)
            local_rank = _local_ranks(probabilities, equation.b, device)
            row = {
                "pde_id": "B6",
                "dimension": d,
                "beta": beta,
                "gamma": gamma,
                "K": equation.K,
                "test_states": len(t),
                "bank_sha256": array_sha256(t_all, x_all, is_validation),
                "ESS_med": float(np.median(np.exp(entropy))),
                "z_rms": float(np.sqrt(np.mean(np.sum(z * z, axis=1)))),
                "f_mean": float(np.mean(generator)),
                "f_std": float(np.std(generator)),
                "f_cv": float(
                    np.std(generator)
                    / max(abs(float(np.mean(generator))), np.finfo(float).tiny)
                ),
                "heat_only_rmse": float(np.sqrt(np.mean((heat - value) ** 2))),
                "u_std": float(np.std(value)),
                "heat_only_skill": float(
                    np.sqrt(np.mean((heat - value) ** 2))
                    / max(float(np.std(value)), np.finfo(float).tiny)
                ),
                "r_glob": global_rank,
                "z_variance_outside_top10": outside_top10,
                "r_loc": float(np.median(local_rank)),
                "r_loc_min": float(np.min(local_rank)),
                "r_loc_max": float(np.max(local_rank)),
                "unequal_norm_range": float(np.max(equation.n2) - np.min(equation.n2)),
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
        # Primary key is largest r_glob; exact ties use smaller beta and then
        # smaller gamma solely to make the frozen rule deterministic.
        selected = sorted(
            eligible, key=lambda row: (-row["r_glob"], row["beta"], row["gamma"])
        )[0]
        status = "qualified"
    else:
        selected = sorted(
            rows,
            key=lambda row: (
                row["normalized_constraint_violation"], row["beta"], row["gamma"]
            ),
        )[0]
        status = "calibration-limited"
    return rows, {
        "dimension": d,
        "status": status,
        "selected_beta": selected["beta"],
        "selected_gamma": selected["gamma"],
        "selected_r_glob": selected["r_glob"],
        "selected_row": selected,
        "heat_monte_carlo": heat_metadata,
        "test_bank_sha256": selected["bank_sha256"],
    }


def calibrate(device_name: str = "cuda") -> dict[str, Any]:
    ensure_directories()
    device = torch.device(device_name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    all_rows: list[dict[str, Any]] = []
    selections: dict[str, Any] = {}
    started = time.perf_counter()
    for d in DIMENSIONS:
        rows, selection = calibrate_dimension(d, device)
        all_rows.extend(rows)
        selections[str(d)] = selection
    fields = list(all_rows[0])
    with (RESULTS_ROOT / "b6_calibration.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)
    payload = {
        "revision": REVISION,
        "code_commit_at_calibration": _git("rev-parse", "HEAD"),
        "selection_rule": (
            "qualify iff ESS_med>=5, 0.5<=z_rms<=4, heat_only_skill>=0.15, "
            "f_cv>=0.05; maximize r_glob, tie-break smaller beta"
        ),
        "fallback_violation_normalization": {
            "ESS_med": "max((5-x)/5,0)",
            "z_rms_lower": "max((0.5-x)/0.5,0)",
            "z_rms_upper": "max((x-4)/4,0)",
            "heat_only_skill": "max((0.15-x)/0.15,0)",
            "f_cv": "max((0.05-x)/0.05,0)",
        },
        "selections": selections,
        "grid_rows": len(all_rows),
        "elapsed_seconds": time.perf_counter() - started,
        "environment": _environment(device),
        "production_solver_touched_B6": False,
    }
    write_json_atomic(CALIBRATION_FREEZE_PATH, payload)
    write_json_atomic(RESULTS_ROOT / "b6_calibration_summary.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    calibrate(args.device)


if __name__ == "__main__":
    main()
