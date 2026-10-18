"""Generate frozen transfer inputs and run mathematical/numerical preflight."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    JointConcaveDriver, JointConvexOriginalDriver, SignTransformedDriver,
    dense_orthonormal_basis,
)
from invariant_region_mlp.experiments.nonconvex_tangent_actor.drivers import (
    SineQuadraticDriver,
)

from .benchmarks import (
    AllenCahnBenchmark, HJBLogBenchmark, MLPFailureBenchmark, OUGMMBenchmark,
)
from .prepare_suite import sha256_file, write_json
from .protocol import DIMENSIONS, MLP_FAILURE_DIMENSIONS, RESULTS_DIR, stable_seed


INPUT_DIR = RESULTS_DIR / "inputs"


def _array_hash(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str(value.dtype).encode())
        digest.update(str(value.shape).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _generate_ou_inputs(source_root: Path) -> dict[str, Any]:
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    official_dir = source_root / "invariant_region_mlp/experiments/original_pde_benchmark_comparison/official_inputs"
    records: dict[str, Any] = {}
    for d in DIMENSIONS:
        output = INPUT_DIR / f"ou_gmm_d{d}.npz"
        if d == 100:
            mean = torch.load(
                official_dir / "mean_100d_ms=1.0_vs=2.0_5.pt",
                map_location="cpu", weights_only=False,
            ).double().numpy()
            pi = torch.load(
                official_dir / "pi_100d_ms=1.0_vs=2.0_5.pt",
                map_location="cpu", weights_only=False,
            ).double().numpy()
            source = "official tracked tensors; covariance restored to 2I"
            seed = None
        else:
            seed = stable_seed("ou_gmm_transfer_constructor", d)
            generator = torch.Generator(device="cpu").manual_seed(seed)
            mean = (2.0 * torch.rand((5, d), generator=generator) - 1.0).double().numpy()
            raw_pi = torch.rand(5, generator=generator).double()
            pi = (raw_pi / raw_pi.sum()).numpy()
            source = "official constructor distributions: means U[-1,1], covariance 2I, normalized U[0,1] weights"
        np.savez_compressed(
            output, dimension=np.asarray(d), mean=mean, pi=pi,
            variance_diagonal=np.full((5, d), 2.0, dtype=np.float64),
            constructor_seed=np.asarray(-1 if seed is None else seed),
        )
        records[str(d)] = {
            "path": str(output.relative_to(source_root)),
            "sha256": sha256_file(output),
            "parameter_array_sha256": _array_hash(mean, pi),
            "source": source,
            "seed": seed,
            "mean_min": float(mean.min()), "mean_max": float(mean.max()),
            "pi": pi.tolist(), "pi_sum": float(pi.sum()),
        }
    return records


def _test_bank(
    pde: str, d: int, count: int, *, T: float,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(stable_seed("icml_final_test_bank", pde, d))
    t = rng.uniform(0.0, np.nextafter(T, 0.0), count)
    x = rng.standard_normal((count, d))
    if pde == "ou_gmm":
        x *= np.sqrt(4.0 + t[:, None])
    elif pde == "allen_cahn":
        x *= np.sqrt(2.0 * t[:, None])
    return t.astype(np.float64), x.astype(np.float64)


def _generate_banks(source_root: Path) -> dict[str, Any]:
    historical = {
        ("hjb_log", 100): source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison/evaluation_bank_hjb_log_d100.npz",
        ("ou_gmm", 100): source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison/evaluation_bank_dpi_ou_gmm_d100.npz",
        ("allen_cahn", 100): source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison/evaluation_bank_allen_cahn_d100.npz",
        ("nonconvex_coupled", 100): source_root / "invariant_region_mlp/results/tangent_actor_accuracy_frontier/predictions/post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz",
    }
    settings = {
        "hjb_log": (DIMENSIONS, 1.0), "ou_gmm": (DIMENSIONS, 1.0),
        "joint_concave": (DIMENSIONS, 0.5), "joint_convex": (DIMENSIONS, 0.5),
        "nonconvex_coupled": (DIMENSIONS, 0.5),
        "mlp_failure": (MLP_FAILURE_DIMENSIONS, 1.0),
        "allen_cahn": (DIMENSIONS, 0.3),
    }
    records: dict[str, Any] = {}
    for pde, (dims, T) in settings.items():
        for d in dims:
            key = f"{pde}_d{d}"
            if (pde, d) in historical:
                path = historical[(pde, d)]
                records[key] = {
                    "path": str(path.relative_to(source_root)), "sha256": sha256_file(path),
                    "historical_exact_reuse": True,
                }
                continue
            t, x = _test_bank(pde, d, 960, T=T)
            path = INPUT_DIR / f"test_bank_{key}.npz"
            np.savez_compressed(path, t=t, x=x)
            records[key] = {
                "path": str(path.relative_to(source_root)), "sha256": sha256_file(path),
                "array_sha256": _array_hash(t, x), "states": len(t),
                "historical_exact_reuse": False,
            }
    return records


def _smooth_tangent_check(benchmark: Any, *, seed: int) -> dict[str, float]:
    generator = torch.Generator().manual_seed(seed)
    count = 64
    r = torch.rand(count, generator=generator, dtype=torch.float64, requires_grad=True)
    q = torch.randn((count, benchmark.d), generator=generator, dtype=torch.float64, requires_grad=True)
    value = benchmark.driver(r, q)
    derivative_r, derivative_q = torch.autograd.grad(
        value.sum(), (r, q), allow_unused=True
    )
    if derivative_r is None:
        derivative_r = torch.zeros_like(r)
    if derivative_q is None:
        derivative_q = torch.zeros_like(q)
    c, b, h = benchmark.tangent_coefficients(r.detach(), q.detach())
    reconstructed = c * r.detach() + (b * q.detach()).sum(dim=-1) + h
    return {
        "c_max_abs": float((c - derivative_r).abs().max()),
        "b_max_abs": float((b - derivative_q).abs().max()),
        "tangent_reconstruction_max_abs": float((reconstructed - value.detach()).abs().max()),
    }


def _joint_tangent_check(driver: Any, *, seed: int) -> dict[str, float]:
    generator = torch.Generator().manual_seed(seed)
    r = driver.r_low + (driver.r_high - driver.r_low) * torch.rand(
        64, generator=generator, dtype=torch.float64
    )
    q = torch.randn((64, driver.d), generator=generator, dtype=torch.float64)
    c, b, h = driver.tangent_coefficients(r, q)
    ca, ba, ha = driver.tangent_coefficients_autodiff(r, q)
    return {
        "c_max_abs": float((c - ca).abs().max()),
        "b_max_abs": float((b - ba).abs().max()),
        "h_max_abs": float((h - ha).abs().max()),
    }


def _math_audit(source_root: Path) -> dict[str, Any]:
    ou = OUGMMBenchmark.from_npz(INPUT_DIR / "ou_gmm_d100.npz")
    smooth = {
        "hjb_log": _smooth_tangent_check(HJBLogBenchmark(10), seed=11),
        "ou_gmm": _smooth_tangent_check(ou, seed=12),
        "allen_cahn": _smooth_tangent_check(AllenCahnBenchmark(10), seed=13),
    }
    joint_concave = JointConcaveDriver(d=10, k=2)
    joint_convex_original = JointConvexOriginalDriver(d=10, k=2)
    joint_convex_working = SignTransformedDriver(joint_convex_original)
    nonconvex = SineQuadraticDriver(d=10, lam=0.2, eta=1.0)
    joint = {
        "joint_concave": _joint_tangent_check(joint_concave, seed=21),
        "joint_convex_working": _joint_tangent_check(joint_convex_working, seed=22),
        "nonconvex_coupled": _joint_tangent_check(nonconvex, seed=23),
    }
    y = torch.linspace(0.05, 0.95, 64, dtype=torch.float64)
    z = torch.randn((64, 10), generator=torch.Generator().manual_seed(24), dtype=torch.float64)
    sign_identity = float(
        (joint_convex_working.driver(-y, -z) + joint_convex_original.driver(y, z)).abs().max()
    )
    failure = MLPFailureBenchmark(20)
    q_zero = torch.zeros((8, 20), dtype=torch.float64)
    c0, b0, h0 = failure.tangent_coefficients(None, q_zero)
    q_nonzero = torch.randn((32, 20), generator=torch.Generator().manual_seed(25), dtype=torch.float64)
    _, b1, h1 = failure.tangent_coefficients(None, q_nonzero)
    failure_check = {
        "zero_subgradient_norm": float(torch.linalg.vector_norm(b0, dim=-1).max()),
        "zero_h_max_abs": float(h0.abs().max()),
        "nonzero_inactive_direction_norm_error": float(
            (torch.linalg.vector_norm(b1[:, 1:], dim=-1) - 1.0).abs().max()
        ),
        "nonzero_tangent_remainder_max_abs": float(h1.abs().max()),
        "selection_rule": "minimum-norm subgradient b(0)=0",
    }
    terminal_checks: dict[str, Any] = {}
    for d in DIMENSIONS:
        basis = dense_orthonormal_basis(d, 2, seed=24017)
        terminal_checks[str(d)] = {
            "basis_orthogonality_max_abs": float(np.max(np.abs(basis.T @ basis - np.eye(2)))),
            "basis_max_abs_entry": float(np.max(np.abs(basis))),
            "dense_nonzero_fraction": float(np.mean(basis != 0.0)),
        }
    existing_references = {
        "original": source_root / "invariant_region_mlp/results/original_pde_benchmark_comparison/reference_accuracy_audit.json",
        "joint": source_root / "invariant_region_mlp/results/convex_concave_semilinear_extended/numerical_reference_validation.csv",
        "nonconvex": source_root / "invariant_region_mlp/results/tangent_actor_accuracy_frontier/references/coupled_refinement_audit.json",
        "mlp_failure": source_root / "invariant_region_mlp/results/counterexample_rescue/independent_validation.json",
    }
    reference_hashes = {key: sha256_file(path) for key, path in existing_references.items()}
    all_small = all(max(check.values()) < 1e-10 for check in smooth.values())
    all_joint_small = all(max(check.values()) < 1e-10 for check in joint.values())
    basis_ok = all(row["basis_orthogonality_max_abs"] < 1e-12 for row in terminal_checks.values())
    return {
        "smooth_tangent_checks": smooth,
        "joint_tangent_checks": joint,
        "convex_sign_transform_identity_max_abs": sign_identity,
        "mlp_failure_subgradient": failure_check,
        "dense_basis": terminal_checks,
        "existing_reference_artifact_sha256": reference_hashes,
        "passed": bool(
            all_small and all_joint_small and sign_identity < 1e-12 and basis_ok
            and failure_check["zero_subgradient_norm"] == 0.0
            and failure_check["nonzero_tangent_remainder_max_abs"] < 1e-10
        ),
    }


def main() -> None:
    source_root = Path(__file__).resolve().parents[3]
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    ou = _generate_ou_inputs(source_root)
    banks = _generate_banks(source_root)
    math_audit = _math_audit(source_root)
    payload = {
        "passed": math_audit["passed"],
        "ou_gmm_transfer_parameters": ou,
        "test_banks": banks,
        "mathematical_audit": math_audit,
        "ground_truth_used_for_training_or_selection": False,
    }
    write_json(RESULTS_DIR / "PREFLIGHT_AUDIT.json", payload)
    print(json.dumps({
        "passed": payload["passed"], "ou_instances": len(ou), "test_banks": len(banks),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
