"""Run mathematical compatibility and independent-reference gates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .common import REFERENCE_DIR, RESULTS_DIR, atomic_json, ensure_dirs
from .drivers import SineQuadraticDriver
from .reference2d import save_reference, validate_reference


def _derivative_audit(driver: SineQuadraticDriver) -> dict[str, Any]:
    generator = torch.Generator().manual_seed(7091 + driver.d)
    r = torch.rand(64, generator=generator, dtype=torch.float64)
    q = torch.randn(64, driver.d, generator=generator, dtype=torch.float64) * 0.2
    analytic = driver.tangent_coefficients(r, q)
    autodiff = driver.tangent_coefficients_autodiff(r, q)
    errors = [
        float(torch.max(torch.abs(left - right)).item())
        for left, right in zip(analytic, autodiff, strict=True)
    ]
    point = torch.cat(
        (torch.tensor([0.63], dtype=torch.float64), q[0].detach())
    )

    def scalar(argument: torch.Tensor) -> torch.Tensor:
        return driver.driver(argument[:1], argument[1:].reshape(1, -1))[0]

    actual_hessian = torch.autograd.functional.hessian(scalar, point)
    expected_hessian = torch.zeros_like(actual_hessian)
    expected_hessian[0, 0] = driver.driver_uu(point[:1])[0]
    expected_hessian[1:, 1:] = -torch.eye(driver.d, dtype=torch.float64)
    expected_hessian[0, 1:] = expected_hessian[1:, 0] = -driver.eta * torch.as_tensor(
        driver.e, dtype=torch.float64
    )
    return {
        "tangent_c_max_abs_error": errors[0],
        "tangent_b_max_abs_error": errors[1],
        "tangent_h_max_abs_error": errors[2],
        "joint_hessian_max_abs_error": float(
            torch.max(torch.abs(actual_hessian - expected_hessian)).item()
        ),
        "formal_curvature_at_u_0p25": float(
            driver.driver_uu(torch.tensor([0.25], dtype=torch.float64))[0]
        ),
        "formal_curvature_at_u_0p75": float(
            driver.driver_uu(torch.tensor([0.75], dtype=torch.float64))[0]
        ),
        "mixed_derivative_norm": float(driver.eta),
        "passed": bool(max(errors) < 1e-10 and torch.max(torch.abs(actual_hessian - expected_hessian)) < 1e-10),
    }


def _write_compatibility(audits: dict[str, dict[str, Any]]) -> None:
    lines = [
        "# Mathematical Compatibility Audit",
        "",
        "The implementation uses the existing generic tangent-coefficient and weighted Feynman–Kac code. No PDE-specific linear solver is introduced.",
        "",
        "| Driver | max tangent error | max Hessian error | f_uu>0 fraction | f_uu<0 fraction | invariant [0,1] | gate |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for name, payload in audits.items():
        derivative = payload["derivatives"]
        curvature = payload["reference"]["curvature_coverage"]
        barrier = payload["barriers"]
        passed = payload["passed"]
        lines.append(
            f"| {name} | {max(derivative['tangent_c_max_abs_error'], derivative['tangent_b_max_abs_error'], derivative['tangent_h_max_abs_error']):.3e} | "
            f"{derivative['joint_hessian_max_abs_error']:.3e} | {curvature['f_uu_positive_fraction']:.3f} | "
            f"{curvature['f_uu_negative_fraction']:.3f} | {barrier['invariant_interval_verified_analytically']} | {'PASS' if passed else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            "The interval check uses the forward-time form `w_tau = Delta w + f(w,sqrt(2) grad w)`: 0 is a subsolution and 1 is a supersolution for both drivers. The numerical reference independently remains inside this interval.",
            "",
            "The local Taylor remainder follows from the bounded joint Hessian on the invariant interval and bounded Actor region. This is a local consistency statement; it is not used as a global convergence theorem.",
        ]
    )
    (RESULTS_DIR / "MATHEMATICAL_COMPATIBILITY_AUDIT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def _write_reference_report(audits: dict[str, dict[str, Any]]) -> None:
    lines = [
        "# Reference Convergence Report",
        "",
        "References use a periodic Fourier pseudospectral discretization in the two active coordinates. The production solution uses Strang diffusion splitting and an explicit midpoint nonlinear step. A classical method-of-lines RK4 solve supplies an integration-scheme cross-check.",
        "",
        "| Driver | production (n,dt) | fine→production Value RMSE | Gradient-vector RMSE | PDE residual RMSE | aliasing sensitivity | RK4 scheme difference | range |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for name, payload in audits.items():
        audit = payload["reference"]
        production = audit["production"]
        fine = audit["fine_vs_production"]
        aliasing = audit["aliasing_n96_dealiased_vs_raw"]
        rk4 = audit["independent_mol_rk4_n48_vs_split_n48"]
        lines.append(
            f"| {name} | ({production['n']},{production['dt']:.3g}) | {fine['value_rmse']:.3e} | "
            f"{fine['gradient_vector_rmse']:.3e} | {audit['pde_residual_rmse']:.3e} | "
            f"{aliasing['value_rmse']:.3e} | {rk4['value_rmse']:.3e} | "
            f"[{audit['solution_grid_min']:.6f}, {audit['solution_grid_max']:.6f}] |"
        )
    lines.extend(
        [
            "",
            "Reference arrays are retained as reproducible local artifacts and are never imported by Actor training or checkpoint selection. Only post-freeze evaluation loads them.",
        ]
    )
    (RESULTS_DIR / "REFERENCE_CONVERGENCE_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run(which: str) -> dict[str, Any]:
    ensure_dirs()
    specifications = {
        "primary": SineQuadraticDriver(d=100, lam=0.2, eta=0.0),
        "coupled": SineQuadraticDriver(d=100, lam=0.2, eta=1.0),
    }
    selected = specifications if which == "both" else {which: specifications[which]}
    output_path = RESULTS_DIR / "preflight_audit.json"
    if output_path.exists():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        audits = dict(existing.get("drivers", {}))
    else:
        audits = {}
    for name, driver in selected.items():
        reference, reference_audit = validate_reference(driver)
        path = REFERENCE_DIR / f"{driver.name}.npz"
        save_reference(path, reference)
        derivative_audit = _derivative_audit(driver)
        curvature = reference_audit["curvature_coverage"]
        passed = bool(
            derivative_audit["passed"]
            and driver.barrier_audit()["invariant_interval_verified_analytically"]
            and reference_audit["inside_invariant_interval"]
            and curvature["meaningful_sign_change_gate"]
            and reference_audit["terminal_value_rmse"] < 2e-3
            and reference_audit["fine_vs_production"]["value_rmse"] < 5e-4
        )
        audits[name] = {
            "driver": driver.metadata(),
            "derivatives": derivative_audit,
            "barriers": driver.barrier_audit(),
            "reference": reference_audit,
            "reference_path": str(path),
            "ground_truth_visible_to_training": False,
            "passed": passed,
        }
        print(
            f"{name}: gate={'PASS' if passed else 'FAIL'} "
            f"ref={reference_audit['fine_vs_production']['value_rmse']:.3e} "
            f"curvature=({curvature['f_uu_positive_fraction']:.3f},"
            f"{curvature['f_uu_negative_fraction']:.3f})",
            flush=True,
        )
    payload = {
        "stage": "A mathematical and CPU correctness",
        "drivers": audits,
        "all_requested_passed": all(item["passed"] for item in audits.values()),
        "tests": "pytest test_correctness.py",
    }
    atomic_json(output_path, payload)
    _write_compatibility(audits)
    _write_reference_report(audits)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", choices=("primary", "coupled", "both"), default="primary")
    args = parser.parse_args()
    result = run(args.driver)
    if not result["all_requested_passed"]:
        raise SystemExit("preflight gate failed")


if __name__ == "__main__":
    main()
