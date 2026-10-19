"""Frozen calibration table for the spectral HJB family and P1 anchor."""

from __future__ import annotations

import csv
import json
import math
import time
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.double_estimator.equations import (
    fixed_points as anchor_fixed_points,
    make_equation as make_anchor_equation,
)

from .equations import FAMILIES, fixed_points, make_equation, profile_second
from .protocol import RESULTS_ROOT, write_json_atomic
from .solver import current_commit


def _variance_outside(z: np.ndarray, rank: int) -> float:
    centred = np.asarray(z, dtype=np.float64) - np.mean(z, axis=0, keepdims=True)
    singular = np.linalg.svd(centred, compute_uv=False)
    energy = singular * singular
    total = float(np.sum(energy))
    return float(np.sum(energy[rank:]) / total) if total > 0.0 else 0.0


def _row(
    pde_id: str,
    d: int,
    u: np.ndarray,
    z: np.ndarray,
    heat: np.ndarray,
    hessian_eigenvalues: np.ndarray,
) -> dict[str, Any]:
    generator = -0.5 * np.sum(z * z, axis=1)
    value_std = float(np.std(u, ddof=0))
    heat_rmse = float(np.sqrt(np.mean((heat - u) ** 2)))
    squared = hessian_eigenvalues * hessian_eigenvalues
    fourth = squared * squared
    effective_rank = np.sum(squared, axis=1) ** 2 / np.maximum(
        np.sum(fourth, axis=1), np.finfo(np.float64).tiny
    )
    return {
        "pde_id": pde_id,
        "dimension": d,
        "test_points": len(u),
        "std_exact_u": value_std,
        "mean_exact_generator": float(np.mean(generator)),
        "std_exact_generator": float(np.std(generator, ddof=0)),
        "rms_z_norm": float(np.sqrt(np.mean(np.sum(z * z, axis=1)))),
        "mean_z_norm": float(np.linalg.norm(np.mean(z, axis=0))),
        "heat_iterate_skill": heat_rmse / value_std,
        "centred_z_variance_outside_top1": _variance_outside(z, 1),
        "centred_z_variance_outside_top10": _variance_outside(z, 10),
        "mean_terminal_hessian_frobenius": float(
            np.mean(np.sqrt(np.sum(squared, axis=1)))
        ),
        "mean_terminal_hessian_effective_rank": float(np.mean(effective_rank)),
    }


def calibrate() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    started = time.perf_counter()
    for pde_id in FAMILIES:
        for d in (100, 400):
            equation = make_equation(pde_id, d)
            t, x, is_validation = fixed_points(pde_id, d)
            test = ~is_validation
            state = equation.exact_state(t[test], x[test])
            heat = equation.heat_iterate(t[test], x[test])
            hessian = equation.terminal_hessian_eigenvalues(x[test])
            records.append(
                _row(pde_id, d, state[:, 0], state[:, 1:], heat, hessian)
            )

    for d in (100, 400):
        equation = make_anchor_equation("P1", d)
        t, x, is_validation = anchor_fixed_points("P1", d)
        test = ~is_validation
        state = equation.exact_state(t[test], x[test])
        heat = equation.fzero_reference(t[test], x[test])
        projection = x[test] @ equation.w
        hessian_scalar = profile_second(projection)
        hessian = np.zeros((len(projection), d), dtype=np.float64)
        hessian[:, 0] = hessian_scalar
        records.append(_row("P1", d, state[:, 0], state[:, 1:], heat, hessian))

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    fields = list(records[0])
    with (RESULTS_ROOT / "calibration.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    write_json_atomic(
        RESULTS_ROOT / "calibration_summary.json",
        {
            "code_commit": current_commit(),
            "rows": len(records),
            "expected_rows": 8,
            "elapsed_seconds": time.perf_counter() - started,
            "reported_no_verdict": True,
        },
    )
    print(json.dumps(records, indent=2), flush=True)
    return records


if __name__ == "__main__":
    calibrate()
