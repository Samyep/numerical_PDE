"""Stage 0 reference-only calibration and implementation checks."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .common import (
    RESULTS_DIR,
    array_sha256,
    atomic_json,
    ensure_result_dirs,
    fixed_reference_points,
    load_config,
    stable_seed,
)
from .equation import B1Equation, effective_rank
from .mlp_adapter import cpu_gpu_tape_check


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def calibration_row(
    d: int,
    k: int,
    scale: float,
    *,
    point_count: int,
    heat_samples: int,
) -> dict[str, Any]:
    equation = B1Equation(d=d, k=k, coefficient_scale=scale)
    t, x = fixed_reference_points(d, k, point_count)
    value = equation.exact_u(t, x)
    z = equation.exact_z(t, x)
    driver = equation.exact_generator(t, x)
    probabilities = equation.weights(t, x)
    heat = equation.heat_only_mc(
        t,
        x,
        n_samples=heat_samples,
        seed=stable_seed("B1", "heat", d, k, scale, point_count, heat_samples),
        chunk_size=32,
    )
    active_z = z @ equation.B
    covariance = np.atleast_2d(np.cov(active_z, rowvar=False))
    covariance_eigenvalues = np.linalg.eigvalsh(covariance)
    gradient_rank = effective_rank(np.atleast_1d(covariance_eigenvalues))
    hessian_ranks: list[float] = []
    minimum_rank = k
    for hessian in equation.exact_active_hessian(t[: min(128, len(t))], x[: min(128, len(x))]):
        eigenvalues = np.maximum(-np.linalg.eigvalsh(hessian), 0.0)
        hessian_ranks.append(effective_rank(eigenvalues))
        minimum_rank = min(minimum_rank, int(np.linalg.matrix_rank(hessian, tol=1e-11)))
    entropy = -np.sum(probabilities * np.log(np.maximum(probabilities, 1e-300)), axis=1)
    residual = equation.analytic_residual(t, x)
    return {
        "d": d,
        "k": k,
        "K": equation.K,
        "coefficient_scale": None if k == 1 else scale,
        "point_count": point_count,
        "heat_mc_samples": heat_samples,
        "reference_point_sha256": array_sha256(t, x),
        "value_mean": float(np.mean(value)),
        "value_spatial_std": float(np.std(value)),
        "driver_mean": float(np.mean(driver)),
        "driver_spatial_std": float(np.std(driver)),
        "heat_only_rmse": float(np.sqrt(np.mean((heat - value) ** 2))),
        "gradient_covariance_effective_rank": gradient_rank,
        "gradient_rank_ratio": gradient_rank / k,
        "hessian_effective_rank_mean": float(np.mean(hessian_ranks)),
        "hessian_rank_ratio": float(np.mean(hessian_ranks)) / k,
        "hessian_min_algebraic_rank": minimum_rank,
        "softmax_max_probability_mean": float(np.mean(np.max(probabilities, axis=1))),
        "softmax_normalized_entropy_mean": float(np.mean(entropy / math.log(equation.K))),
        "pde_residual_linf": float(np.max(np.abs(residual))),
        "nonlinearity_flag": bool(
            np.std(driver) < 0.05
            or np.sqrt(np.mean((heat - value) ** 2)) < 0.02
            or gradient_rank / k < 0.2
            or np.mean(hessian_ranks) / k < 0.2
        ),
        "basis_sha256": equation.provenance()["basis_sha256"],
        "coefficient_sha256": equation.provenance()["coefficient_sha256"],
    }


def run_reference_calibration() -> tuple[list[dict[str, Any]], list[dict[str, Any]], float]:
    config = load_config()
    candidates = [float(value) for value in config["pde"]["coefficient_scale_candidates"]]
    candidate_rows: list[dict[str, Any]] = []
    for scale in candidates:
        for k in (5, 10, 20, 40):
            candidate_rows.append(
                calibration_row(100, k, scale, point_count=512, heat_samples=1024)
            )
    qualified: list[float] = []
    for scale in candidates:
        rows = [row for row in candidate_rows if row["coefficient_scale"] == scale]
        if rows and all(
            row["softmax_normalized_entropy_mean"] >= 0.55
            and row["hessian_rank_ratio"] >= 0.30
            for row in rows
        ):
            qualified.append(scale)
    selected = max(qualified) if qualified else min(candidates)
    frozen = float(config["pde"]["frozen_coefficient_scale"])
    if selected != frozen:
        raise RuntimeError(
            f"reference-only rule selected scale {selected}, but config freezes {frozen}"
        )

    final_rows: list[dict[str, Any]] = []
    for d in config["pde"]["dimensions"]:
        for k in config["pde"]["effective_dimensions"]:
            final_rows.append(
                calibration_row(d, k, frozen, point_count=512, heat_samples=2048)
            )
    _write_csv(RESULTS_DIR / "scale_calibration_candidates.csv", candidate_rows)
    _write_csv(RESULTS_DIR / "pde_calibration.csv", final_rows)
    provenance = {
        "selection_uses_solver_results": False,
        "selection_rule": config["pde"]["scale_selection_rule"],
        "candidate_scales": candidates,
        "qualified_scales": qualified,
        "selected_scale": selected,
        "basis_seed": config["pde"]["basis_seed"],
        "coefficient_seed": config["pde"]["coefficient_seed"],
        "equations": [
            B1Equation(d, k, selected).provenance()
            for d in config["pde"]["dimensions"]
            for k in config["pde"]["effective_dimensions"]
        ],
    }
    atomic_json(RESULTS_DIR / "frozen_pde_config.json", provenance)
    return candidate_rows, final_rows, selected


def environment() -> dict[str, Any]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    payload: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "cuda_available": torch.cuda.is_available(),
        "device": str(device),
    }
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(device)
        payload.update(
            {
                "gpu_name": properties.name,
                "gpu_total_memory_bytes": properties.total_memory,
                "compute_capability": f"{properties.major}.{properties.minor}",
            }
        )
    try:
        payload["commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parents[3],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        payload["commit"] = None
    return payload


def run_implementation_checks(scale: float) -> list[dict[str, Any]]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    equation = B1Equation(20, 5, scale)
    checks = [
        cpu_gpu_tape_check(equation, method=method, n=2, M=3, device=device)
        for method in ("vanilla_mlp", "double_mlp", "actor_label")
    ]
    _write_csv(RESULTS_DIR / "stage0_mlp_tape_checks.csv", checks)
    if not all(row["passed"] for row in checks):
        raise RuntimeError("CPU/GPU draw-tape equivalence failed")
    return checks


def write_report(
    candidate_rows: list[dict[str, Any]],
    final_rows: list[dict[str, Any]],
    selected: float,
    checks: list[dict[str, Any]] | None,
) -> None:
    flagged = [(row["d"], row["k"]) for row in final_rows if row["nonlinearity_flag"]]
    check_text = (
        "deferred (reference-only invocation)"
        if checks is None
        else f"PASS ({len(checks)}/{len(checks)} exact-draw checks)"
    )
    lines = [
        "# B1 Stage 0 correctness report",
        "",
        "**Verdict:** reference construction PASS; MLP draw replay " + check_text + ".",
        "",
        f"The pre-method analytic calibration selected coefficient scale **{selected:g}**. "
        "The rule and all rejected candidate metrics are saved before solver evaluation.",
        "",
        "## Checks",
        "",
        "- Stable float64 terminal and exact-solution log-sum-exp.",
        "- Analytic PDE residual and terminal identity.",
        "- Closed-form gradient checked independently by PyTorch autodiff in the unit suite.",
        "- Dense coefficient affine rank and active Hessian rank.",
        "- `k=1` exactly reproduces the established P1 direction and coefficients.",
        "- Diffusion and gradient-state convention: `sigma=sqrt(2)`, `z=sqrt(2) grad u`.",
        f"- CPU/Torch recursion on identical draws: {check_text}.",
        "",
        "## Calibration",
        "",
        f"Candidate rows: {len(candidate_rows)}. Frozen-grid rows: {len(final_rows)}.",
        f"Configurations flagged as weak by the preregistered diagnostic thresholds: {flagged or 'none'}.",
        "",
        "Artifacts: `scale_calibration_candidates.csv`, `pde_calibration.csv`, "
        "`frozen_pde_config.json`, and `stage0_mlp_tape_checks.csv` when run with a GPU.",
    ]
    (RESULTS_DIR / "stage0_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-only", action="store_true")
    args = parser.parse_args()
    ensure_result_dirs()
    candidate_rows, final_rows, selected = run_reference_calibration()
    checks = None if args.reference_only else run_implementation_checks(selected)
    atomic_json(RESULTS_DIR / "stage0_environment.json", environment())
    summary = {
        "reference_checks_passed": bool(
            all(row["pde_residual_linf"] <= 1e-12 for row in final_rows)
            and all(row["hessian_min_algebraic_rank"] == row["k"] for row in final_rows)
        ),
        "selected_coefficient_scale": selected,
        "flagged_cells": [
            [row["d"], row["k"]] for row in final_rows if row["nonlinearity_flag"]
        ],
        "mlp_tape_checks_passed": None if checks is None else all(row["passed"] for row in checks),
    }
    if not summary["reference_checks_passed"]:
        raise RuntimeError("reference checks failed")
    atomic_json(RESULTS_DIR / "stage0_summary.json", summary)
    write_report(candidate_rows, final_rows, selected, checks)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
