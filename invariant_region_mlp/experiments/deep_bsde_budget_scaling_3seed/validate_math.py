"""Independent PDE, terminal, parameter, and source-fidelity audit."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import math

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import (
    exact_z,
    make_benchmark,
    metadata,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    OfficialDeepBSDE,
)

from .protocol import DIMENSION, RESULTS_DIR, TRAJECTORY_BY_ID, write_json_atomic
from .solver import VariableStartDeepBSDE, deep_config


EXPECTED_PARAMETER_HASHES = {
    "b6_fullrank": "8d661b71b57f0f37c388a0c1ad81d97897f5838e315eabda8a499bf4abe7927a",
    "dense_fourier": "b96804e2f7c20cbd4e72c0f961a4ebfc21bf46cae120d198b6afe6ce2e2684b6",
}


def _exact_torch(pde: str, benchmark, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    if pde == "dense_fourier":
        return benchmark.exact_u_torch(t, x)
    equation = benchmark.equation
    vectors = torch.as_tensor(equation.b, dtype=x.dtype, device=x.device)
    norms = torch.as_tensor(equation.n2, dtype=x.dtype, device=x.device)
    logits = x @ vectors.T + (benchmark.T - t.reshape(-1, 1)) * norms
    return -(torch.logsumexp(logits, dim=1) - math.log(equation.K))


def _cell(pde: str) -> dict:
    benchmark = make_benchmark(pde, DIMENSION)
    generator = torch.Generator().manual_seed(4817 if pde == "b6_fullrank" else 4818)
    x = (2.0 * torch.rand((1, DIMENSION), generator=generator, dtype=torch.float64) - 1.0).requires_grad_(True)
    t_value = 0.37 * benchmark.T
    t = torch.tensor([t_value], dtype=torch.float64, requires_grad=True)
    u = _exact_torch(pde, benchmark, t, x)
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    diagonal = []
    for coordinate in range(DIMENSION):
        second = torch.autograd.grad(
            grad[0, coordinate], x, retain_graph=True
        )[0][0, coordinate]
        diagonal.append(second)
    laplacian = torch.stack(diagonal).sum()
    z = math.sqrt(2.0) * grad
    driver = benchmark.driver(u, z)
    residual = u_t + laplacian + driver

    terminal_t = torch.tensor([benchmark.T], dtype=torch.float64)
    terminal_u = _exact_torch(pde, benchmark, terminal_t, x.detach())
    terminal = benchmark.terminal(x.detach())
    z_numpy = exact_z(benchmark, np.asarray([t_value]), x.detach().numpy())
    z_error = float(np.max(np.abs(z.detach().numpy() - z_numpy)))
    info = metadata(benchmark)
    return {
        "pde": pde,
        "T": benchmark.T,
        "sigma": benchmark.sigma,
        "driver_formula": (
            "-0.5*sum(z**2)" if pde == "b6_fullrank"
            else "0.5*tanh(y)*sum(z**2)"
        ),
        "parameter_sha256": info["parameter_sha256"],
        "expected_parameter_sha256": EXPECTED_PARAMETER_HASHES[pde],
        "parameter_hash_match": info["parameter_sha256"] == EXPECTED_PARAMETER_HASHES[pde],
        "pde_residual_abs": abs(float(residual.detach())),
        "terminal_max_abs_error": float(torch.max(torch.abs(terminal_u - terminal))),
        "z_sqrt2_grad_max_abs_error": z_error,
    }


def validate() -> dict:
    cells = [_cell(pde) for pde in ("b6_fullrank", "dense_fourier")]
    benchmark = make_benchmark("b6_fullrank", DIMENSION)
    config = deep_config(benchmark, 0.1, TRAJECTORY_BY_ID["c2_steps10"], updates=1)
    model = VariableStartDeepBSDE(
        benchmark, config, torch.zeros(DIMENSION, dtype=torch.float64), 0.1
    )
    source_checks = {
        "inherits_official_architecture": isinstance(model, OfficialDeepBSDE),
        "scalar_initial_value_shape": list(model.y_init.shape) == [1],
        "learned_initial_z_shape": list(model.z_init.shape) == [1, DIMENSION],
        "time_slice_subnet_count": len(model.subnets) == config.intervals - 1,
        "float64_profile": config.dtype == "float64",
        "remaining_horizon": math.isclose(config.total_time, benchmark.T - 0.1, abs_tol=2e-15),
        "euler_dt": math.isclose(model.dt, (benchmark.T - 0.1) / 10, abs_tol=2e-15),
    }
    failures = []
    for cell in cells:
        if not cell["parameter_hash_match"]:
            failures.append(f"{cell['pde']}: parameter hash")
        if cell["pde_residual_abs"] > 2e-10:
            failures.append(f"{cell['pde']}: PDE residual")
        if cell["terminal_max_abs_error"] > 2e-12:
            failures.append(f"{cell['pde']}: terminal identity")
        if cell["z_sqrt2_grad_max_abs_error"] > 2e-11:
            failures.append(f"{cell['pde']}: z convention")
    failures.extend(name for name, passed in source_checks.items() if not passed)
    payload = {
        "status": "PASS" if not failures else "FAIL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "cells": cells,
        "source_fidelity_checks": source_checks,
        "failures": failures,
    }
    write_json_atomic(RESULTS_DIR / "mathematical_adaptation_audit.json", payload)
    return payload


if __name__ == "__main__":
    result = validate()
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)

