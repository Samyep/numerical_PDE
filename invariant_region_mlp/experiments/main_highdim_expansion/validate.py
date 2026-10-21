"""Independent mathematical and stochastic pre-training validation gate."""

from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    evaluate_policy as evaluate_gradient_policy,
    generate_labels as generate_gradient_labels,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy as evaluate_joint_policy,
)
from invariant_region_mlp.experiments.tangent_actor_mcv_training.labels import (
    mcv_tangent_labels,
)

from .benchmarks import (
    B6Equation,
    array_sha256,
    dense_parameters,
    make_benchmark,
    make_dense_fourier,
    sample_cube_states,
)
from .joint_actor import JointTerminalActor
from .prepare import file_sha256, prepare
from .protocol import (
    DIMENSIONS,
    GRADIENT_ACTOR_CONFIG,
    JOINT_ACTOR_CONFIG,
    MANIFEST_PATH,
    REVISION,
    VALIDATION_PATH,
    frozen_protocol,
    stable_seed,
    write_json_atomic,
)


def _check(name: str, value: float, threshold: float, relation: str = "le") -> dict[str, Any]:
    passed = value <= threshold if relation == "le" else value >= threshold
    return {
        "name": name,
        "value": float(value),
        "threshold": float(threshold),
        "relation": relation,
        "passed": bool(passed),
    }


def _laplacian(value: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True)[0]
    result = torch.zeros_like(value)
    for index in range(x.shape[1]):
        result = result + torch.autograd.grad(
            gradient[:, index].sum(), x, create_graph=True
        )[0][:, index]
    return result


def _math_checks() -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    rng = np.random.default_rng(20261020)

    b6 = B6Equation(6, beta=1.5, gamma=1.0)
    t_np = rng.uniform(0.0, 0.24, 4)
    x_np = rng.uniform(-1.0, 1.0, (4, 6))
    state = b6.exact_state(t_np, x_np)
    terminal_error = np.max(
        np.abs(b6.exact_u(np.full(4, b6.T), x_np) - b6.terminal(x_np))
    )
    checks.append(_check("b6_terminal_identity", terminal_error, 2e-13))
    h = 1e-5
    finite = np.empty_like(x_np)
    for j in range(b6.d):
        plus, minus = x_np.copy(), x_np.copy()
        plus[:, j] += h
        minus[:, j] -= h
        finite[:, j] = (b6.exact_u(t_np, plus) - b6.exact_u(t_np, minus)) / (2 * h)
    checks.append(
        _check("b6_analytic_z_finite_difference", np.max(np.abs(math.sqrt(2.0) * finite - state[:, 1:])), 3e-9)
    )
    b = torch.tensor(b6.b.copy(), dtype=torch.float64)
    n2 = torch.tensor(b6.n2.copy(), dtype=torch.float64)
    t = torch.tensor(t_np, dtype=torch.float64, requires_grad=True)
    x = torch.tensor(x_np, dtype=torch.float64, requires_grad=True)
    u = -(torch.logsumexp(x @ b.T + (b6.T - t)[:, None] * n2, 1) - math.log(b6.K))
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    residual = u_t + _laplacian(u, x) - grad.square().sum(1)
    checks.append(_check("b6_autodiff_pde_residual", float(residual.detach().abs().max()), 2e-11))
    hopf_cole = torch.exp(-u)
    hopf_t = torch.autograd.grad(hopf_cole.sum(), t, create_graph=True)[0]
    hopf_residual = hopf_t + _laplacian(hopf_cole, x)
    checks.append(
        _check(
            "b6_hopf_cole_heat_identity",
            float(hopf_residual.detach().abs().max()),
            2e-11,
        )
    )

    dense = make_dense_fourier(8)
    t_np = rng.uniform(0.0, 0.49, 4)
    x_np = rng.uniform(-1.0, 1.0, (4, 8))
    state = dense.exact_state(t_np, x_np)
    checks.append(
        _check(
            "dense_terminal_identity",
            np.max(np.abs(dense.exact_u(np.full(4, dense.T), x_np) - dense.terminal_numpy(x_np))),
            2e-13,
        )
    )
    finite = np.empty_like(x_np)
    for j in range(dense.d):
        plus, minus = x_np.copy(), x_np.copy()
        plus[:, j] += h
        minus[:, j] -= h
        finite[:, j] = (dense.exact_u(t_np, plus) - dense.exact_u(t_np, minus)) / (2 * h)
    checks.append(
        _check("dense_analytic_z_finite_difference", np.max(np.abs(math.sqrt(2.0) * finite - state[:, 1:])), 3e-9)
    )
    t = torch.tensor(t_np, dtype=torch.float64, requires_grad=True)
    x = torch.tensor(x_np, dtype=torch.float64, requires_grad=True)
    u = dense.exact_u_torch(t, x)
    u_t = torch.autograd.grad(u.sum(), t, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), x, create_graph=True)[0]
    z = math.sqrt(2.0) * grad
    residual = u_t + _laplacian(u, x) + 0.5 * torch.tanh(u) * z.square().sum(1)
    checks.append(_check("dense_autodiff_pde_residual", float(residual.detach().abs().max()), 2e-11))
    omega = torch.as_tensor(dense.omega, dtype=torch.float64)
    signs = torch.as_tensor(dense.signs, dtype=torch.float64)
    phases = torch.as_tensor(dense.phases, dtype=torch.float64)
    attenuation = torch.exp(
        -(dense.T - t)[:, None]
        * torch.as_tensor(dense.omega_norm_squared, dtype=torch.float64)
    )
    transformed = 0.5 + dense.scale * torch.sum(
        attenuation * signs * torch.sin(x @ omega.T + phases), dim=1
    )
    checks.append(
        _check(
            "dense_sinh_transform_value_identity",
            float((torch.sinh(u) - transformed).detach().abs().max()),
            2e-13,
        )
    )
    transformed_t = torch.autograd.grad(
        transformed.sum(), t, create_graph=True
    )[0]
    transformed_residual = transformed_t + _laplacian(transformed, x)
    checks.append(
        _check(
            "dense_sinh_transform_heat_identity",
            float(transformed_residual.detach().abs().max()),
            2e-11,
        )
    )

    generator = torch.Generator().manual_seed(991)
    r = torch.randn(5, dtype=torch.float64, generator=generator, requires_grad=True)
    q = torch.randn(5, dense.d, dtype=torch.float64, generator=generator, requires_grad=True)
    f = dense.driver(r, q)
    expected_c = torch.autograd.grad(f.sum(), r, create_graph=True)[0]
    expected_b = torch.autograd.grad(f.sum(), q, create_graph=True)[0]
    c, drift, source = dense.tangent_coefficients(r, q)
    checks.append(_check("dense_tangent_c_autodiff", float((c - expected_c).detach().abs().max()), 2e-13))
    checks.append(_check("dense_tangent_b_autodiff", float((drift - expected_b).detach().abs().max()), 2e-13))
    reconstructed = c * r + (drift * q).sum(1) + source
    checks.append(_check("dense_tangent_reconstruction", float((reconstructed - f).detach().abs().max()), 2e-13))
    return checks


def _artifact_checks(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    seen: set[tuple[str, int]] = set()
    for cell in manifest["cells"]:
        key = (cell["pde"], int(cell["dimension"]))
        seen.add(key)
        checks.append(_check(f"{key}_parameter_nonfinite", cell["nonfinite_parameters"], 0.0))
        checks.append(_check(f"{key}_bank_nonfinite", cell["nonfinite_bank"], 0.0))
        checks.append(_check(f"{key}_test_states", cell["test_states"], 960.0, "ge"))
        checks.append(
            {
                "name": f"{key}_parameter_file_hash",
                "value": file_sha256(cell_path := Path(cell["parameter_file"])),
                "threshold": cell["parameter_file_sha256"],
                "relation": "eq",
                "passed": file_sha256(cell_path) == cell["parameter_file_sha256"],
            }
        )
        bank_path = Path(cell["test_bank_file"])
        checks.append(
            {
                "name": f"{key}_bank_file_hash",
                "value": file_sha256(bank_path),
                "threshold": cell["test_bank_file_sha256"],
                "relation": "eq",
                "passed": file_sha256(bank_path) == cell["test_bank_file_sha256"],
            }
        )
    expected = {(pde, d) for pde in ("b6_fullrank", "dense_fourier") for d in DIMENSIONS}
    checks.append(
        {
            "name": "complete_parameter_and_bank_grid",
            "value": sorted(seen),
            "threshold": sorted(expected),
            "relation": "eq",
            "passed": seen == expected,
        }
    )
    return checks


def _stochastic_checks(device: torch.device) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {}
    # Small independent points and budgets exercise the production FK/MCV
    # paths without consuming any formal training/test states.
    b6 = make_benchmark("b6_fullrank", 20)
    seed_every = stable_seed(REVISION, "validation", "b6_actor")
    seed_everything(seed_every)
    b6_actor = TerminalAwareActor(
        b6,
        width=16,
        depth=2,
        activation=GRADIENT_ACTOR_CONFIG.activation,
    ).to(device)
    t, x = b6.sample_states(8, seed=stable_seed(REVISION, "validation", "b6_points"))
    q_label, label_meta = generate_gradient_labels(
        b6,
        b6_actor,
        t,
        x,
        paths=16,
        time_steps=4,
        seed=stable_seed(REVISION, "validation", "b6_labels"),
        device=device,
        parent_batch=4,
        path_chunk=8,
        estimator="mcv",
    )
    output, eval_meta = evaluate_gradient_policy(
        b6,
        b6_actor,
        t,
        x,
        paths=32,
        time_steps=4,
        seed=stable_seed(REVISION, "validation", "b6_eval"),
        device=device,
        parent_batch=4,
        path_chunk=8,
    )
    checks.append(_check("b6_gradient_labels_nonfinite", np.size(q_label) - np.count_nonzero(np.isfinite(q_label)), 0.0))
    checks.append(_check("b6_mcv_evaluation_nonfinite", eval_meta["nonfinite_count"], 0.0))
    checks.append(_check("b6_reaction_weight_ess", eval_meta["weight_ess_fraction"], 0.999999, "ge"))
    metadata["b6"] = {"labels": label_meta, "evaluation": eval_meta}

    dense = make_dense_fourier(20)
    seed_everything(stable_seed(REVISION, "validation", "dense_actor"))
    dense_actor = JointTerminalActor(dense, width=16, depth=2).to(device)
    t, x = sample_cube_states(
        "dense_fourier",
        20,
        8,
        seed=stable_seed(REVISION, "validation", "dense_points"),
    )
    r_label, q_label, label_meta = mcv_tangent_labels(
        dense,
        dense_actor,
        t,
        x,
        paths=16,
        time_steps=4,
        seed=stable_seed(REVISION, "validation", "dense_labels"),
        device=device,
        parent_batch=4,
        path_chunk=8,
    )
    if q_label is None:
        raise RuntimeError("joint validation labels omitted q")
    output, eval_meta = evaluate_joint_policy(
        dense,
        dense_actor,
        t,
        x,
        paths=32,
        time_steps=4,
        seed=stable_seed(REVISION, "validation", "dense_eval"),
        device=device,
        parent_batch=4,
        path_chunk=8,
    )
    nonfinite = (
        np.size(r_label) + np.size(q_label)
        - np.count_nonzero(np.isfinite(r_label))
        - np.count_nonzero(np.isfinite(q_label))
    )
    checks.append(_check("dense_joint_labels_nonfinite", nonfinite, 0.0))
    checks.append(_check("dense_mcv_evaluation_nonfinite", eval_meta["nonfinite_count"], 0.0))
    checks.append(_check("dense_reaction_weight_positive", label_meta["weight_min"], 0.999999, "ge"))
    checks.append(_check("dense_reaction_weight_ess", eval_meta["weight_ess_fraction"], 0.8, "ge"))
    # Plain and MCV estimate the same discretized policy value.  Use their
    # paired statewise uncertainty rather than requiring either to be closer
    # to the analytic nonlinear solution.
    difference = output["plain"] - output["mcv"]
    se = np.sqrt(
        output["plain_pair_variance"] / 16.0
        + output["mcv_pair_variance"] / 16.0
    )
    standardized = float(np.sqrt(np.mean(np.square(difference / np.maximum(se, 1e-12)))))
    checks.append(_check("dense_plain_mcv_consistency", standardized, 5.0))
    metadata["dense_fourier"] = {"labels": label_meta, "evaluation": eval_meta}
    return checks, metadata


def validate(device_name: str = "cpu") -> dict[str, Any]:
    started = time.perf_counter()
    device = torch.device(device_name)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    manifest = prepare()
    if manifest["protocol"] != frozen_protocol():
        raise RuntimeError("frozen input protocol mismatch")
    checks = _math_checks() + _artifact_checks(manifest)
    stochastic, stochastic_metadata = _stochastic_checks(device)
    checks.extend(stochastic)
    overall = all(item["passed"] for item in checks)
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if overall else "FAILED",
        "overall_pass": overall,
        "checks": checks,
        "failures": [item for item in checks if not item["passed"]],
        "stochastic_metadata": stochastic_metadata,
        "manifest_path": str(MANIFEST_PATH),
        "manifest_sha256": file_sha256(MANIFEST_PATH),
        "ground_truth_used_for_training_or_selection": False,
        "elapsed_seconds": time.perf_counter() - started,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "device": str(device),
            "cuda_runtime": torch.version.cuda,
        },
    }
    write_json_atomic(VALIDATION_PATH, payload)
    if not overall:
        raise RuntimeError(json.dumps(payload["failures"], indent=2))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    result = validate(args.device)
    print(json.dumps({"status": result["status"], "checks": len(result["checks"])}, sort_keys=True))


if __name__ == "__main__":
    main()
