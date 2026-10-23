"""Freeze query banks, parameters, and independently audited references."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import shutil
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import (
    dense_parameters,
    make_b6 as make_source_b6,
    make_dense_fourier as make_source_dense,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    AllenRadialConfig,
    AllenRadialReference,
    hjb_log_quad_scalar,
    hjb_log_reference,
)

from .benchmarks import benchmark_metadata, exact_values, make_benchmark
from .protocol import (
    B6_SELECTIONS,
    DIMENSIONS,
    FROZEN_PROTOCOL_JSON,
    PARAMETER_DIR,
    PARAMETER_MANIFEST,
    PDES,
    QUERY_COUNT,
    QUERY_SEEDS,
    REFERENCE_VALIDATION,
    RESULTS_DIR,
    REVISION,
    SOURCE_EXPANSION_RESULTS,
    array_sha256,
    ensure_directories,
    file_sha256,
    frozen_protocol,
    query_path,
    write_json_atomic,
    write_npz_atomic,
)


def _save_npy_immutable(path: Path, value: np.ndarray) -> None:
    value = np.asarray(value)
    if path.is_file():
        existing = np.load(path, allow_pickle=False)
        if existing.dtype != value.dtype or existing.shape != value.shape or not np.array_equal(existing, value):
            raise RuntimeError(f"frozen array mismatch: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npy")
    np.save(temporary, value, allow_pickle=False)
    temporary.replace(path)


def _query_bank(pde: str, d: int) -> tuple[np.ndarray, np.ndarray]:
    seed = int(QUERY_SEEDS[pde][d])
    rng = np.random.Generator(np.random.PCG64(seed))
    t = np.zeros(QUERY_COUNT, dtype=np.float64)
    if pde in {"hjb_log", "allen_cahn"}:
        random_x = rng.standard_normal((QUERY_COUNT - 1, d), dtype=np.float64)
        distribution = "N(0,I_d)"
    else:
        random_x = rng.uniform(-1.0, 1.0, (QUERY_COUNT - 1, d)).astype(np.float64)
        distribution = "U([-1,1]^d)"
    x = np.concatenate((np.zeros((1, d), dtype=np.float64), random_x), axis=0)
    if not np.array_equal(x[0], np.zeros(d, dtype=np.float64)):
        raise AssertionError("query row zero is not the exact origin")
    if len(np.unique(x.view(np.dtype((np.void, x.dtype.itemsize * d))))) != QUERY_COUNT:
        raise RuntimeError(f"duplicate query state in {pde} d={d}")
    return t, x


def _hjb_reference(t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    order128 = np.asarray(hjb_log_reference(t, x, order=128), dtype=np.float64)
    order256 = np.asarray(hjb_log_reference(t, x, order=256), dtype=np.float64)
    radius2 = np.square(x).sum(axis=1)
    adaptive = np.asarray(
        [hjb_log_quad_scalar(float(tt), float(rr), d=x.shape[1]) for tt, rr in zip(t, radius2)],
        dtype=np.float64,
    )
    discrepancy = np.maximum(np.abs(order256 - order128), np.abs(order256 - adaptive))
    maximum = float(np.max(discrepancy))
    if maximum > 5e-8:
        raise RuntimeError(f"HJB reference convergence failed: {maximum:.3e}")
    return order256, discrepancy, {
        "kind": "Hopf-Cole scaled Gauss-Laguerre plus independent adaptive quadrature",
        "orders": [128, 256],
        "maximum_order_difference": float(np.max(np.abs(order256 - order128))),
        "maximum_adaptive_difference": float(np.max(np.abs(order256 - adaptive))),
        "maximum_uncertainty": maximum,
        "threshold": 5e-8,
        "passed": True,
    }


def _allen_solve(config: AllenRadialConfig, t: np.ndarray, radius2: np.ndarray) -> np.ndarray:
    return AllenRadialReference(config).evaluate(t, radius2)


def _allen_reference(d: int, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    radius2 = np.square(x).sum(axis=1)
    # The radial solver coordinate is s=||x||^2.  The fixed 20-unit buffer is
    # stricter than the protocol's eight-unit/four-diffusion-SD minimum.
    base_s_max = max(40.0, math.ceil((float(radius2.max()) + 20.0) / 10.0) * 10.0)
    # First converge the far boundary at a fixed fine discretization.  More
    # than one 25% expansion is allowed and was explicitly preregistered.
    domain_levels: list[dict[str, Any]] = []
    previous_s_max = base_s_max
    previous_cfg = AllenRadialConfig(
        s_max=previous_s_max, ds=0.05, dt=5e-5, dimension=d, horizon=0.3,
    )
    previous_value = _allen_solve(previous_cfg, t, radius2)
    domain_levels.append({"name": "domain_0", **asdict(previous_cfg)})
    domain_difference = np.full_like(previous_value, np.inf)
    accepted_s_max = previous_s_max
    accepted_previous_s_max = previous_s_max
    for domain_index in range(1, 7):
        next_s_max = math.ceil(1.25 * previous_s_max / 10.0) * 10.0
        next_cfg = AllenRadialConfig(
            s_max=next_s_max, ds=0.05, dt=5e-5, dimension=d, horizon=0.3,
        )
        next_value = _allen_solve(next_cfg, t, radius2)
        domain_difference = np.abs(next_value - previous_value)
        domain_levels.append({
            "name": f"domain_{domain_index}",
            **asdict(next_cfg),
            "maximum_change_from_previous": float(domain_difference.max()),
        })
        accepted_previous_s_max = previous_s_max
        accepted_s_max = next_s_max
        if float(domain_difference.max()) <= 2e-5:
            break
        previous_s_max, previous_value = next_s_max, next_value
    else:
        raise RuntimeError(f"Allen-Cahn domain convergence failed d={d}")

    # Cancel the first-order radial drift discretization with two independent
    # Richardson estimates.  If their difference is still too large, add one
    # more halved level before applying the frozen threshold.
    coarse_cfg = AllenRadialConfig(
        s_max=accepted_s_max, ds=0.1, dt=1e-4, dimension=d, horizon=0.3,
    )
    fine_cfg = AllenRadialConfig(
        s_max=accepted_s_max, ds=0.05, dt=5e-5, dimension=d, horizon=0.3,
    )
    superfine_cfg = AllenRadialConfig(
        s_max=accepted_s_max, ds=0.025, dt=2.5e-5, dimension=d, horizon=0.3,
    )
    coarse = _allen_solve(coarse_cfg, t, radius2)
    fine = _allen_solve(fine_cfg, t, radius2)
    superfine = _allen_solve(superfine_cfg, t, radius2)
    extrapolated_coarse = 2.0 * fine - coarse
    extrapolated_fine = 2.0 * superfine - fine
    discretization = np.abs(extrapolated_fine - extrapolated_coarse)
    refinement_levels: list[dict[str, Any]] = [
        {"name": "coarse", **asdict(coarse_cfg)},
        {"name": "fine", **asdict(fine_cfg)},
        {"name": "superfine", **asdict(superfine_cfg)},
    ]
    final = extrapolated_fine
    if float(discretization.max()) > 2e-5:
        ultra_cfg = AllenRadialConfig(
            s_max=accepted_s_max, ds=0.0125, dt=1.25e-5,
            dimension=d, horizon=0.3,
        )
        ultra = _allen_solve(ultra_cfg, t, radius2)
        refinement_levels.append({"name": "ultra", **asdict(ultra_cfg)})
        extrapolated_ultra = 2.0 * ultra - superfine
        discretization = np.abs(extrapolated_ultra - extrapolated_fine)
        final = extrapolated_ultra

    previous_domain_superfine_cfg = AllenRadialConfig(
        s_max=accepted_previous_s_max, ds=0.025, dt=2.5e-5,
        dimension=d, horizon=0.3,
    )
    previous_domain_superfine = _allen_solve(
        previous_domain_superfine_cfg, t, radius2
    )
    boundary = np.abs(superfine - previous_domain_superfine)
    uncertainty = np.maximum(discretization, boundary)
    maximum = float(np.max(uncertainty))
    if maximum > 2e-5:
        raise RuntimeError(f"Allen-Cahn reference convergence failed d={d}: {maximum:.3e}")
    return final, uncertainty, {
        "kind": "radial Strang-reaction / Crank-Nicolson convergence audit",
        "coordinate": "s=||x||^2",
        "maximum_query_s": float(radius2.max()),
        "domain_levels": domain_levels,
        "refinement_levels": refinement_levels,
        "richardson_order": 1,
        "maximum_coarse_fine_difference": float(np.abs(fine - coarse).max()),
        "maximum_richardson_difference": float(discretization.max()),
        "maximum_domain_difference": float(boundary.max()),
        "maximum_uncertainty": maximum,
        "threshold": 2e-5,
        "passed": True,
    }


def _analytic_reference(pde: str, d: int, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    primary = exact_values(pde, d, t, x)
    if pde == "b6_fullrank":
        equation = make_source_b6(d).equation
        # Independent stabilized sum with reversed term order.
        logits = x @ equation.b[::-1].T + (equation.T - t)[:, None] * equation.n2[::-1]
        maximum = np.max(logits, axis=1)
        secondary = -(maximum + np.log(np.exp(logits - maximum[:, None]).sum(axis=1)) - math.log(equation.K))
    else:
        driver = make_source_dense(d)
        angle = x @ driver.omega[::-1].T + driver.phases[::-1]
        attenuation = np.exp(-(driver.T - t)[:, None] * driver.omega_norm_squared[::-1])
        terms = attenuation * driver.signs[::-1] * np.sin(angle)
        secondary = np.arcsinh(0.5 + driver.scale * np.sum(terms, axis=1))
    uncertainty = np.abs(primary - secondary)
    return primary, uncertainty, {
        "kind": "analytic formula with independent reversed-order stabilized recomputation",
        "maximum_uncertainty": float(np.max(uncertainty)),
        "passed": bool(np.all(np.isfinite(primary)) and np.max(uncertainty) <= 2e-12),
        "threshold": 2e-12,
    }


def _reference(pde: str, d: int, t: np.ndarray, x: np.ndarray):
    if pde == "hjb_log":
        return _hjb_reference(t, x)
    if pde == "allen_cahn":
        return _allen_reference(d, t, x)
    return _analytic_reference(pde, d, t, x)


def _parameter_arrays(pde: str, d: int) -> dict[str, np.ndarray]:
    if pde == "b6_fullrank":
        benchmark = make_source_b6(d)
        equation = benchmark.equation
        expected = B6_SELECTIONS[d]
        if (equation.beta, equation.gamma) != (expected["beta"], expected["gamma"]):
            raise RuntimeError(f"B6 frozen calibration drift for d={d}")
        return {
            "directions": equation.e,
            "multipliers": equation.U,
            "common_vector": equation.v,
            "vectors": equation.b,
            "squared_norms": equation.n2,
            "beta": np.asarray(equation.beta),
            "gamma": np.asarray(equation.gamma),
            "seed": np.asarray([20261201, d, 66], dtype=np.int64),
        }
    if pde == "dense_fourier":
        omega, signs, phases = dense_parameters(d)
        return {
            "omega": omega,
            "signs": signs,
            "phases": phases,
            "seed": np.asarray([20261201, d, 41], dtype=np.int64),
        }
    return {}


def _source_parameter_rows() -> dict[tuple[str, int], dict[str, Any]]:
    source = SOURCE_EXPANSION_RESULTS / "parameter_and_test_bank_manifest.json"
    payload = json.loads(source.read_text(encoding="utf-8"))
    return {(row["pde"], int(row["dimension"])): row for row in payload["cells"]}


def prepare() -> dict[str, Any]:
    ensure_directories()
    write_json_atomic(FROZEN_PROTOCOL_JSON, frozen_protocol())
    source_rows = _source_parameter_rows()
    parameter_rows: list[dict[str, Any]] = []
    reference_rows: list[dict[str, Any]] = []
    for pde in PDES:
        for d in DIMENSIONS:
            root = query_path(pde, d)
            t, x = _query_bank(pde, d)
            reference, uncertainty, audit = _reference(pde, d, t, x)
            if reference.shape != (QUERY_COUNT,) or uncertainty.shape != (QUERY_COUNT,):
                raise AssertionError((pde, d, reference.shape, uncertainty.shape))
            if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(uncertainty)):
                raise FloatingPointError(f"nonfinite reference for {pde} d={d}")
            _save_npy_immutable(root / "t.npy", t)
            _save_npy_immutable(root / "x.npy", x)
            _save_npy_immutable(root / "reference.npy", reference)
            _save_npy_immutable(root / "reference_uncertainty.npy", uncertainty)
            arrays = _parameter_arrays(pde, d)
            parameter_record: dict[str, Any] = {
                "pde": pde,
                "dimension": d,
                "benchmark": benchmark_metadata(pde, d),
            }
            if arrays:
                parameter_file = PARAMETER_DIR / f"{pde}_d{d}.npz"
                if not parameter_file.is_file():
                    write_npz_atomic(parameter_file, **arrays)
                generated_hash = array_sha256(*arrays.values())
                expected_hash = source_rows[(pde, d)]["parameter_array_sha256"]
                if generated_hash != expected_hash:
                    raise RuntimeError(
                        f"frozen parameter mismatch {pde} d={d}: {generated_hash} != {expected_hash}"
                    )
                parameter_record.update({
                    "file": parameter_file.relative_to(RESULTS_DIR).as_posix(),
                    "file_sha256": file_sha256(parameter_file),
                    "array_sha256": generated_hash,
                    "source_array_sha256": expected_hash,
                    "source_manifest_sha256": file_sha256(
                        SOURCE_EXPANSION_RESULTS / "parameter_and_test_bank_manifest.json"
                    ),
                    "byte_identical_arrays_to_completed_ice_expansion": True,
                })
            parameter_rows.append(parameter_record)
            bank_record = {
                "pde": pde,
                "dimension": d,
                "seed": int(QUERY_SEEDS[pde][d]),
                "distribution": "N(0,I_d)" if pde in {"hjb_log", "allen_cahn"} else "U([-1,1]^d)",
                "origin_row": 0,
                "shape": list(x.shape),
                "dtype": str(x.dtype),
                "t_sha256": file_sha256(root / "t.npy"),
                "x_sha256": file_sha256(root / "x.npy"),
                "reference_sha256": file_sha256(root / "reference.npy"),
                "reference_uncertainty_sha256": file_sha256(root / "reference_uncertainty.npy"),
                "canonical_input_sha256": array_sha256(t, x),
                "canonical_reference_sha256": array_sha256(reference),
                "reference": audit,
                "maximum_reference_uncertainty": float(np.max(uncertainty)),
            }
            write_json_atomic(root / "manifest.json", bank_record)
            reference_rows.append(bank_record)
    parameter_payload = {
        "revision": REVISION,
        "status": "FROZEN_BEFORE_PILOT_OR_FORMAL_TRAINING",
        "source_manifest": (
            SOURCE_EXPANSION_RESULTS / "parameter_and_test_bank_manifest.json"
        ).relative_to(RESULTS_DIR.parents[1]).as_posix(),
        "rows": parameter_rows,
    }
    reference_payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(row["reference"]["passed"] for row in reference_rows) else "FAILED",
        "overall_pass": all(row["reference"]["passed"] for row in reference_rows),
        "banks": reference_rows,
    }
    write_json_atomic(PARAMETER_MANIFEST, parameter_payload)
    write_json_atomic(REFERENCE_VALIDATION, reference_payload)
    if not reference_payload["overall_pass"]:
        raise RuntimeError("reference validation failed")
    return {"parameters": parameter_payload, "references": reference_payload}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()
    result = prepare()
    if args.summary:
        print(json.dumps({
            "parameter_cells": len(result["parameters"]["rows"]),
            "reference_cells": len(result["references"]["banks"]),
            "status": result["references"]["status"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()

