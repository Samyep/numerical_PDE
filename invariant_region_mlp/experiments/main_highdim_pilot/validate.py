"""Mathematical validation and pre-training difficulty calibration."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from typing import Any, Callable

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    evaluate_policy as evaluate_gradient_policy,
)
from invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.solver import (
    evaluate_policy as evaluate_joint_policy,
)

from .actors import JointTerminalActor
from .common import (
    environment_record,
    load_bank,
    load_dense,
    load_lqg,
    synchronize,
)
from .equations import CoupledLQGEquation, DenseFourierDriver, LQGActorBenchmark
from .prepare import prepare
from .protocol import (
    CALIBRATION_PATH,
    GRADIENT_CONFIG,
    HEAT_MC_PATHS,
    HEAT_MC_REPLICATES,
    JOINT_CONFIG,
    REFERENCE_AUDIT_ORDERS,
    REFERENCE_QUADRATURE_ORDER,
    REVISION,
    VALIDATION_PATH,
    ensure_directories,
    write_json_atomic,
)


def _check(value: float, threshold: float, *, relation: str = "le") -> dict[str, Any]:
    passed = value <= threshold if relation == "le" else value >= threshold
    return {
        "value": float(value),
        "threshold": float(threshold),
        "relation": relation,
        "passed": bool(passed),
    }


def _autodiff_residual(
    exact: Callable[[torch.Tensor, torch.Tensor], torch.Tensor],
    driver: Callable[[torch.Tensor, torch.Tensor], torch.Tensor],
    t: np.ndarray,
    x: np.ndarray,
) -> tuple[float, float]:
    tt = torch.tensor(t, dtype=torch.float64, requires_grad=True)
    xx = torch.tensor(x, dtype=torch.float64, requires_grad=True)
    u = exact(tt, xx)
    ut = torch.autograd.grad(u.sum(), tt, create_graph=True)[0]
    grad = torch.autograd.grad(u.sum(), xx, create_graph=True)[0]
    laplacian = torch.zeros_like(u)
    for coordinate in range(xx.shape[1]):
        second = torch.autograd.grad(
            grad[:, coordinate].sum(), xx, create_graph=True, retain_graph=True
        )[0][:, coordinate]
        laplacian = laplacian + second
    z = math.sqrt(2.0) * grad
    residual = ut + laplacian + driver(u, z)
    absolute = float(torch.max(torch.abs(residual)).item())
    scale = float(
        torch.max(torch.abs(ut) + torch.abs(laplacian) + torch.abs(driver(u, z))).item()
    )
    return absolute, absolute / max(scale, 1e-14)


def _effective_dimension(z: np.ndarray) -> tuple[float, float]:
    centered = np.asarray(z, dtype=np.float64) - np.mean(z, axis=0, keepdims=True)
    covariance = centered.T @ centered / max(len(centered) - 1, 1)
    eigenvalues = np.linalg.eigvalsh(covariance)[::-1]
    eigenvalues = np.maximum(eigenvalues, 0.0)
    total = float(np.sum(eigenvalues))
    rank = total**2 / max(float(np.sum(eigenvalues**2)), np.finfo(float).tiny)
    outside = float(np.sum(eigenvalues[10:]) / max(total, np.finfo(float).tiny))
    return rank, outside


def _heat_mc(
    terminal: Callable[[torch.Tensor], torch.Tensor],
    t: np.ndarray,
    x: np.ndarray,
    *,
    T: float,
    paths: int,
    seed: int,
    device: torch.device,
    parent_batch: int = 16,
    path_chunk: int = 1024,
) -> tuple[np.ndarray, np.ndarray, float, int | None]:
    if paths % path_chunk or path_chunk % 2:
        raise ValueError("heat MC paths must be divisible by an even chunk")
    generator = torch.Generator(device=device).manual_seed(seed)
    means = np.empty(len(t), dtype=np.float64)
    pair_variance = np.empty(len(t), dtype=np.float64)
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with torch.no_grad():
        for begin in range(0, len(t), parent_batch):
            end = min(begin + parent_batch, len(t))
            bx = torch.as_tensor(x[begin:end], dtype=torch.float32, device=device)
            bt = torch.as_tensor(t[begin:end], dtype=torch.float32, device=device)
            tau_scale = torch.sqrt(torch.clamp(2.0 * (T - bt), min=0.0))
            pairs: list[torch.Tensor] = []
            for _ in range(0, paths, path_chunk):
                half = path_chunk // 2
                normal = torch.randn(
                    (end - begin, half, x.shape[1]),
                    generator=generator,
                    dtype=torch.float32,
                    device=device,
                )
                positive = bx[:, None, :] + tau_scale[:, None, None] * normal
                negative = bx[:, None, :] - tau_scale[:, None, None] * normal
                pair = 0.5 * (
                    terminal(positive.reshape(-1, x.shape[1])).reshape(end - begin, half)
                    + terminal(negative.reshape(-1, x.shape[1])).reshape(end - begin, half)
                )
                pairs.append(pair.double())
            local = torch.cat(pairs, dim=1)
            means[begin:end] = local.mean(dim=1).cpu().numpy()
            pair_variance[begin:end] = local.var(dim=1, unbiased=True).cpu().numpy()
    synchronize(device)
    elapsed = time.perf_counter() - started
    peak = (
        int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None
    )
    return means, pair_variance, elapsed, peak


def _lqg_hopf_cole_mc(
    equation: CoupledLQGEquation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    samples: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray]:
    generator = torch.Generator(device=device).manual_seed(seed)
    bx = torch.as_tensor(x, dtype=torch.float64, device=device)
    bt = torch.as_tensor(t, dtype=torch.float64, device=device)
    tau = torch.sqrt(torch.clamp(2.0 * (equation.T - bt), min=0.0))
    sum_value = torch.zeros(len(t), dtype=torch.float64, device=device)
    sum_square = torch.zeros_like(sum_value)
    count = 0
    chunk = 4096
    c1 = torch.as_tensor(equation.c1, dtype=torch.float64, device=device)
    c2 = torch.as_tensor(equation.c2, dtype=torch.float64, device=device)
    with torch.no_grad():
        for _ in range(0, samples, chunk):
            local = min(chunk, samples - count)
            normal = torch.randn(
                (len(t), local, equation.d),
                generator=generator,
                dtype=torch.float64,
                device=device,
            )
            state = bx[:, None, :] + tau[:, None, None] * normal
            q = torch.sum(
                c1 * (state[..., :-1] - state[..., 1:]).square()
                + c2 * state[..., 1:].square(),
                dim=-1,
            )
            value = 2.0 / (1.0 + q)
            sum_value += value.sum(dim=1)
            sum_square += value.square().sum(dim=1)
            count += local
    mean = sum_value / count
    variance = torch.clamp(sum_square / count - mean.square(), min=0.0)
    se = torch.sqrt(variance / count)
    return mean.cpu().numpy(), se.cpu().numpy()


def _difficulty_row(
    pde: str,
    truth: np.ndarray,
    driver_values: np.ndarray,
    heat: np.ndarray,
    heat_pair_variance: np.ndarray,
    heat_pairs_total: int,
    heat_seconds: float,
    reference_uncertainty: float,
    benchmark_label: str,
) -> dict[str, Any]:
    u, z = truth[:, 0], truth[:, 1:]
    rank, outside = _effective_dimension(z)
    u_std = float(np.std(u))
    heat_error = heat - u
    return {
        "pde": pde,
        "benchmark_label": benchmark_label,
        "dimension": z.shape[1],
        "test_states": len(u),
        "u_mean": float(np.mean(u)),
        "u_std": u_std,
        "u_min": float(np.min(u)),
        "u_max": float(np.max(u)),
        "z_rms": float(np.sqrt(np.mean(np.sum(z * z, axis=1)))),
        "effective_rank": rank,
        "z_variance_outside_top10": outside,
        "driver_rms": float(np.sqrt(np.mean(driver_values**2))),
        "driver_mean": float(np.mean(driver_values)),
        "driver_std": float(np.std(driver_values)),
        "nonlinear_heat_rmse": float(np.sqrt(np.mean(heat_error**2))),
        "nonlinear_skill": float(
            np.sqrt(np.mean(heat_error**2)) / max(u_std, np.finfo(float).tiny)
        ),
        "heat_mc_pointwise_se_rms": float(
            np.sqrt(np.mean(heat_pair_variance / heat_pairs_total))
        ),
        "heat_mc_paths_total": int(2 * heat_pairs_total),
        "heat_mc_seconds": float(heat_seconds),
        "reference_uncertainty_rmse": float(reference_uncertainty),
    }


def validate(device: torch.device) -> dict[str, Any]:
    ensure_directories()
    prepare()
    dense, lqg = load_dense(), load_lqg()
    dense_bank, lqg_bank = load_bank("dense_fourier"), load_bank("coupled_lqg")
    checks: dict[str, dict[str, Any]] = {}

    rng = np.random.default_rng(90210)
    dense_t = rng.uniform(0.03, 0.47, 3)
    dense_x = rng.uniform(-0.8, 0.8, (3, dense.d))
    dense_terminal_error = float(
        np.max(
            np.abs(
                dense.exact_state(np.full(3, dense.T), dense_x)[:, 0]
                - dense.terminal_numpy(dense_x)
            )
        )
    )
    checks["dense_terminal_identity"] = _check(dense_terminal_error, 2e-13)

    tx = torch.tensor(dense_x, dtype=torch.float64, requires_grad=True)
    tt = torch.tensor(dense_t, dtype=torch.float64, requires_grad=True)
    tu = dense.exact_u_torch(tt, tx)
    tz = math.sqrt(2.0) * torch.autograd.grad(tu.sum(), tx)[0]
    dense_gradient_error = float(
        np.max(np.abs(tz.detach().numpy() - dense.exact_z(dense_t, dense_x)))
    )
    checks["dense_exact_gradient_autodiff"] = _check(dense_gradient_error, 2e-11)
    residual, relative = _autodiff_residual(
        dense.exact_u_torch, dense.driver, dense_t[:2], dense_x[:2]
    )
    checks["dense_pde_residual_absolute"] = _check(residual, 2e-10)
    checks["dense_pde_residual_relative"] = _check(relative, 2e-10)

    r = torch.tensor(rng.normal(size=5), dtype=torch.float64)
    q = torch.tensor(rng.normal(size=(5, dense.d)), dtype=torch.float64)
    c, b, h = dense.tangent_coefficients(r, q)
    ac, ab, ah = dense.tangent_coefficients_autodiff(r, q)
    tangent_derivative_error = max(
        float(torch.max(torch.abs(c - ac)).item()),
        float(torch.max(torch.abs(b - ab)).item()),
        float(torch.max(torch.abs(h - ah)).item()),
    )
    reconstruction = c * r + torch.sum(b * q, dim=-1) + h
    reconstruction_error = float(
        torch.max(torch.abs(reconstruction - dense.driver(r, q))).item()
    )
    checks["dense_tangent_derivatives"] = _check(tangent_derivative_error, 2e-12)
    checks["dense_tangent_reconstruction"] = _check(reconstruction_error, 2e-12)

    dense_actor = JointTerminalActor(
        dense, width=JOINT_CONFIG.width, depth=JOINT_CONFIG.depth
    ).to(device)
    dense_actor.eval()
    terminal_x = torch.as_tensor(dense_x, dtype=torch.float32, device=device)
    terminal_t = torch.full((len(dense_x),), dense.T, device=device)
    with torch.no_grad():
        ar, aq = dense_actor.predict(terminal_t, terminal_x, dense.T)
    actor_terminal_error = max(
        float(np.max(np.abs(ar.cpu().numpy() - dense.terminal_numpy(dense_x)))),
        float(
            np.max(
                np.abs(
                    aq.cpu().numpy()
                    - dense.sigma * dense.terminal_gradient_numpy(dense_x)
                )
            )
        ),
    )
    checks["dense_a0_terminal_anchor"] = _check(actor_terminal_error, 2e-5)

    # LQG matrix construction is independently compared to the published sum.
    lqg_x = rng.normal(size=(8, lqg.d))
    direct_q = np.sum(
        lqg.c1 * np.square(lqg_x[:, :-1] - lqg_x[:, 1:])
        + lqg.c2 * np.square(lqg_x[:, 1:]),
        axis=1,
    )
    matrix_error = float(np.max(np.abs(direct_q - lqg.quadratic(lqg_x))))
    checks["lqg_matrix_construction"] = _check(matrix_error, 5e-12)
    checks["lqg_matrix_spd"] = _check(float(lqg.eigenvalues[0]), 0.0, relation="ge")

    audit_t = np.array([0.03, 0.17, 0.31, 0.46])
    audit_x = lqg_bank["test_x"][:4]
    references = {
        order: lqg.exact_state(audit_t, audit_x, order=order)
        for order in REFERENCE_AUDIT_ORDERS
    }
    refine_96_144 = float(np.max(np.abs(references[96] - references[144])))
    refine_144_192 = float(np.max(np.abs(references[144] - references[192])))
    refine_192_256 = float(np.max(np.abs(references[192] - references[256])))
    checks["lqg_quadrature_finest_difference"] = _check(refine_192_256, 2e-10)
    checks["lqg_quadrature_refines"] = {
        "value": refine_192_256,
        "previous": refine_144_192,
        "coarsest": refine_96_144,
        "criterion": (
            "monotone decrease before the 2e-10 convergence floor, or all "
            "successive differences already below that floor"
        ),
        "passed": bool(
            (refine_192_256 <= refine_144_192 <= refine_96_144)
            or max(refine_96_144, refine_144_192, refine_192_256) <= 2e-10
        ),
    }

    terminal_state = lqg.exact_state(
        np.full(len(audit_x), lqg.T), audit_x, order=REFERENCE_QUADRATURE_ORDER
    )
    lqg_terminal_error = max(
        float(np.max(np.abs(terminal_state[:, 0] - lqg.terminal(audit_x)))),
        float(
            np.max(
                np.abs(
                    terminal_state[:, 1:]
                    - lqg.sigma * lqg.terminal_gradient(audit_x)
                )
            )
        ),
    )
    checks["lqg_terminal_and_gradient"] = _check(lqg_terminal_error, 2e-12)

    # Finite differences audit the independently derived integral gradient.
    fd_errors: list[float] = []
    epsilon = 2e-5
    for row in range(2):
        for coordinate in (0, 1, 17, 63, 99):
            plus, minus = audit_x[row].copy(), audit_x[row].copy()
            plus[coordinate] += epsilon
            minus[coordinate] -= epsilon
            fd = (
                lqg.exact_u(audit_t[row], plus[None, :])[0]
                - lqg.exact_u(audit_t[row], minus[None, :])[0]
            ) / (2.0 * epsilon)
            analytic = lqg.exact_z(audit_t[row], audit_x[row : row + 1])[0, coordinate] / lqg.sigma
            fd_errors.append(abs(fd - analytic))
    checks["lqg_integral_gradient_finite_difference"] = _check(max(fd_errors), 2e-7)

    mc_w, mc_se = _lqg_hopf_cole_mc(
        lqg,
        audit_t,
        audit_x,
        samples=262144,
        seed=991827,
        device=device,
    )
    quadrature_w = np.exp(-references[REFERENCE_QUADRATURE_ORDER][:, 0])
    standardized = np.abs(mc_w - quadrature_w) / np.maximum(mc_se, 1e-15)
    checks["lqg_hopf_cole_independent_mc"] = _check(float(np.max(standardized)), 4.5)

    residual, relative = _autodiff_residual(
        lambda t, x: lqg.exact_u_torch(t, x, order=REFERENCE_QUADRATURE_ORDER),
        lambda u, z: -0.5 * torch.sum(z.square(), dim=-1),
        audit_t[:1],
        audit_x[:1],
    )
    checks["lqg_pde_residual_absolute"] = _check(residual, 2e-8)
    checks["lqg_pde_residual_relative"] = _check(relative, 2e-8)

    # Cheap stochastic checks exercise the actual imported FK/MCV paths.
    dense_st = dense_bank["selector_t"][:8]
    dense_sx = dense_bank["selector_x"][:8]
    dense_out, dense_meta = evaluate_joint_policy(
        dense,
        dense_actor,
        dense_st,
        dense_sx,
        paths=256,
        time_steps=16,
        seed=77123,
        device=device,
        parent_batch=4,
        path_chunk=32,
    )
    dense_difference = dense_out["plain"] - dense_out["mcv"]
    dense_difference_se = np.sqrt(
        (
            dense_out["plain_pair_variance"]
            + dense_out["mcv_pair_variance"]
        )
        / 128.0
    )
    dense_standardized_mean = abs(float(np.mean(dense_difference))) / max(
        float(np.sqrt(np.sum(dense_difference_se**2)) / len(dense_difference)),
        1e-14,
    )
    checks["dense_plain_mcv_consistency"] = _check(dense_standardized_mean, 5.0)
    checks["dense_mcv_variance_reduction"] = _check(
        float(dense_meta["variance_reduction"]), 1.0, relation="ge"
    )
    checks["dense_discount_weight_sign"] = _check(
        float(dense_meta["weight_min"]), 1.0 - 2e-6, relation="ge"
    )

    lqg_benchmark = LQGActorBenchmark(lqg)
    lqg_actor = TerminalAwareActor(
        lqg_benchmark,
        width=GRADIENT_CONFIG.width,
        depth=GRADIENT_CONFIG.depth,
        activation=GRADIENT_CONFIG.activation,
    ).to(device)
    lqg_actor.eval()
    lqg_out, lqg_meta = evaluate_gradient_policy(
        lqg_benchmark,
        lqg_actor,
        lqg_bank["selector_t"][:8],
        lqg_bank["selector_x"][:8],
        paths=256,
        time_steps=16,
        seed=77124,
        device=device,
        parent_batch=4,
        path_chunk=32,
    )
    lqg_difference = lqg_out["plain"] - lqg_out["mcv"]
    lqg_se = np.sqrt(
        (lqg_out["plain_pair_variance"] + lqg_out["mcv_pair_variance"]) / 128.0
    )
    lqg_standardized_mean = abs(float(np.mean(lqg_difference))) / max(
        float(np.sqrt(np.sum(lqg_se**2)) / len(lqg_se)), 1e-14
    )
    checks["lqg_plain_mcv_consistency"] = _check(lqg_standardized_mean, 5.0)
    checks["lqg_mcv_variance_reduction"] = _check(
        float(lqg_meta["variance_reduction"]), 1.0, relation="ge"
    )

    calibration_rows: list[dict[str, Any]] = []
    calibration_start = time.perf_counter()
    for pde, equation, bank, terminal, benchmark_label in (
        (
            "Dense Fourier joint",
            dense,
            dense_bank,
            dense.terminal,
            "ANALYTICALLY_INTEGRABLE_HIGH_DIMENSIONAL_STRESS_TEST",
        ),
        (
            "Published coupled LQG",
            lqg,
            lqg_bank,
            lqg.terminal_torch,
            "RECREATED_PUBLISHED_BENCHMARK",
        ),
    ):
        test_t, test_x = bank["test_t"], bank["test_x"]
        if isinstance(equation, DenseFourierDriver):
            truth = equation.exact_state(test_t, test_x)
            driver_values = 0.5 * np.tanh(truth[:, 0]) * np.sum(
                truth[:, 1:] ** 2, axis=1
            )
            reference_uncertainty = 2e-13
        else:
            coarse = equation.exact_state(test_t, test_x, order=192)
            truth = equation.exact_state(test_t, test_x, order=256)
            driver_values = -0.5 * np.sum(truth[:, 1:] ** 2, axis=1)
            reference_uncertainty = float(
                np.sqrt(np.mean(np.square(coarse[:, 0] - truth[:, 0])))
            )
        heat_predictions, variances, seconds, _peak = [], [], 0.0, None
        for replicate in range(HEAT_MC_REPLICATES):
            prediction, variance, elapsed, peak = _heat_mc(
                terminal,
                test_t,
                test_x,
                T=equation.T,
                paths=HEAT_MC_PATHS,
                seed=882100 + 97 * replicate + (0 if pde.startswith("Dense") else 1),
                device=device,
            )
            heat_predictions.append(prediction)
            variances.append(variance)
            seconds += elapsed
            _peak = max(_peak or 0, peak or 0)
        heat = np.mean(heat_predictions, axis=0)
        mean_pair_variance = np.mean(variances, axis=0)
        total_pairs = HEAT_MC_REPLICATES * (HEAT_MC_PATHS // 2)
        row = _difficulty_row(
            pde,
            truth,
            driver_values,
            heat,
            mean_pair_variance,
            total_pairs,
            seconds,
            reference_uncertainty,
            benchmark_label,
        )
        row["heat_mc_peak_gpu_bytes"] = _peak
        calibration_rows.append(row)

    with CALIBRATION_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(calibration_rows[0]))
        writer.writeheader()
        writer.writerows(calibration_rows)

    required = [name for name in checks if not name.endswith("variance_reduction")]
    overall_pass = all(bool(checks[name]["passed"]) for name in required)
    payload = {
        "study_revision": REVISION,
        "status": "VALIDATED" if overall_pass else "REFERENCE_NOT_VALIDATED",
        "overall_pass": overall_pass,
        "checks": checks,
        "dense_fourier": dense.metadata(),
        "coupled_lqg": lqg.metadata(),
        "lqg_quadrature_refinement": {
            "96_to_144_max_state_difference": refine_96_144,
            "144_to_192_max_state_difference": refine_144_192,
            "192_to_256_max_state_difference": refine_192_256,
            "independent_mc_max_standard_errors": float(np.max(standardized)),
        },
        "stochastic_audits": {
            "dense": dense_meta,
            "lqg": lqg_meta,
        },
        "calibration_rows": calibration_rows,
        "calibration_total_seconds": time.perf_counter() - calibration_start,
        "environment": environment_record(),
        "ground_truth_used_for_training_or_selection": False,
    }
    write_json_atomic(VALIDATION_PATH, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    result = validate(device)
    print(json.dumps({"status": result["status"], "checks": result["checks"]}, indent=2))
    if not result["overall_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()


__all__ = ["validate"]

