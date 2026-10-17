"""Create immutable test banks, official inputs, and numerical reference audit."""

from __future__ import annotations

import json
import math
from pathlib import Path
import shutil

import numpy as np
import torch

from .benchmarks import BENCHMARK_NAMES, DPIOUBenchmark, make_benchmark
from .common import (
    OFFICIAL_WORK_DIR,
    PACKAGE_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    sha256_file,
    stable_seed,
)
from .references import (
    AllenRadialConfig,
    AllenRadialReference,
    allen_reference_extrapolated,
    dpi_exact_reference,
    hjb_log_quad_scalar,
    hjb_log_reference,
)


OFFICIAL_ROOT = Path("/work/nvme/bgvf/ywang151/original_pde_source_audit")


def _hjb_audit() -> dict[str, object]:
    benchmark = make_benchmark("hjb_log_d100")
    t, x = benchmark.sample_states(
        24, seed=stable_seed("hjb_reference_audit"), purpose="function"
    )
    order64 = hjb_log_reference(t, x, order=64)
    order128 = hjb_log_reference(t, x, order=128)
    order256 = hjb_log_reference(t, x, order=256)
    adaptive = np.asarray(
        [
            hjb_log_quad_scalar(float(ti), float(np.dot(xi, xi)))
            for ti, xi in zip(t[:8], x[:8])
        ]
    )
    native = float(
        hjb_log_reference(np.asarray([0.0]), np.zeros((1, 100)), order=128)[0]
    )
    return {
        "method": "Hopf-Cole plus scaled Gauss-Laguerre quadrature",
        "native_u_0_0": native,
        "order64_vs_128_max_abs": float(np.max(np.abs(order64 - order128))),
        "order128_vs_256_max_abs": float(np.max(np.abs(order128 - order256))),
        "order128_vs_adaptive_quad_max_abs_first8": float(
            np.max(np.abs(order128[:8] - adaptive))
        ),
        "states": len(t),
        "passed": bool(
            np.max(np.abs(order64 - order128)) < 1e-10
            and np.max(np.abs(order128[:8] - adaptive)) < 1e-10
        ),
    }


def _dpi_residual(benchmark: DPIOUBenchmark) -> dict[str, float | bool]:
    rng = np.random.default_rng(stable_seed("dpi_exact_residual"))
    t = torch.as_tensor(
        rng.uniform(0.05, 0.95, 4), dtype=torch.float64
    ).requires_grad_(True)
    x = torch.as_tensor(rng.normal(size=(4, 100)), dtype=torch.float64).requires_grad_(True)
    value = benchmark.exact_value(t, x)
    u_t = torch.autograd.grad(value.sum(), t, create_graph=True)[0]
    gradient = torch.autograd.grad(value.sum(), x, create_graph=True)[0]
    laplacian = torch.zeros_like(value)
    for coordinate in range(benchmark.d):
        laplacian = laplacian + torch.autograd.grad(
            gradient[:, coordinate].sum(), x, create_graph=False,
            retain_graph=coordinate + 1 < benchmark.d,
        )[0][:, coordinate]
    F = benchmark.theta * (benchmark.mu - x)
    residual = (
        u_t
        + 0.5 * benchmark.alpha * laplacian
        - (F * gradient).sum(dim=1)
        - 0.5 * benchmark.alpha * gradient.square().sum(dim=1)
        - benchmark.d * benchmark.theta
    )
    return {
        "maximum_exact_pde_residual": float(residual.detach().abs().max().item()),
        "rms_exact_pde_residual": float(
            torch.sqrt(residual.detach().square().mean()).item()
        ),
        "passed": bool(residual.detach().abs().max().item() < 2e-9),
    }


def _dpi_audit() -> dict[str, object]:
    benchmark = DPIOUBenchmark()
    t, x = benchmark.sample_states(
        32, seed=stable_seed("dpi_reference_audit"), purpose="official"
    )
    numpy_value, numpy_z = dpi_exact_reference(
        benchmark, t, x, return_gradient=True
    )
    tt = torch.as_tensor(t, dtype=torch.float64)
    xx = torch.as_tensor(x, dtype=torch.float64)
    torch_value = benchmark.exact_value(tt, xx).detach().numpy()
    torch_z = (benchmark.sigma * benchmark.exact_gradient(tt, xx)).detach().numpy()
    return {
        "method": "closed-form evolving diagonal Gaussian mixture",
        "numpy_vs_torch_value_max_abs": float(
            np.max(np.abs(numpy_value - torch_value))
        ),
        "numpy_vs_torch_z_max_abs": float(np.max(np.abs(numpy_z - torch_z))),
        "terminal_value_identity_max_abs": float(
            np.max(
                np.abs(
                    dpi_exact_reference(
                        benchmark, np.ones(len(x)), x
                    )
                    - benchmark.terminal(xx).detach().numpy()
                )
            )
        ),
        **_dpi_residual(benchmark),
    }


def _allen_audit() -> dict[str, object]:
    configs = [
        AllenRadialConfig(s_max=180.0, ds=0.1, dt=5e-5),
        AllenRadialConfig(s_max=200.0, ds=0.05, dt=2.5e-5),
        AllenRadialConfig(s_max=200.0, ds=0.025, dt=1.25e-5),
        AllenRadialConfig(s_max=200.0, ds=0.0125, dt=6.25e-6),
    ]
    raw: list[dict[str, object]] = []
    for config in configs:
        value = AllenRadialReference(config).evaluate(
            np.asarray([0.0]), np.asarray([0.0])
        )[0]
        raw.append({"config": config.__dict__, "u_0_0": float(value)})
    dummy_x = np.zeros((1, 100), dtype=np.float64)
    extrapolated, metadata = allen_reference_extrapolated(
        np.asarray([0.0]), dummy_x
    )
    successive = [
        abs(float(raw[index]["u_0_0"]) - float(raw[index - 1]["u_0_0"]))
        for index in range(1, len(raw))
    ]
    return {
        "method": "squared-radius PDE with monotone CN/Strang solver and Richardson extrapolation",
        "refinements": raw,
        "successive_native_differences": successive,
        "native_u_0_0_extrapolated": float(extrapolated[0]),
        "extrapolation": metadata,
        "reference_uncertainty_proxy": float(successive[-1]),
        "invariant_region_verified_numerically": bool(
            0.0 < extrapolated[0] < 1.0
        ),
        "passed": bool(successive[-1] < 1.1e-5),
    }


def _write_bank(name: str, count: int) -> dict[str, object]:
    benchmark = make_benchmark(name)
    purpose = "function" if name == "hjb_log_d100" else (
        "official" if name == "dpi_ou_gmm_d100" else "native"
    )
    t, x = benchmark.sample_states(
        count, seed=stable_seed("original_benchmark_final_bank", name), purpose=purpose
    )
    if name == "hjb_log_d100":
        value, z = hjb_log_reference(t, x, order=128, return_gradient=True)
        reference_meta: dict[str, object] = {"kind": "Hopf-Cole quadrature", "order": 128}
    elif name == "dpi_ou_gmm_d100":
        value, z = dpi_exact_reference(benchmark, t, x, return_gradient=True)
        reference_meta = {"kind": "analytic evolving GMM"}
    else:
        value, reference_meta = allen_reference_extrapolated(t, x)
        z = np.empty((0, benchmark.d), dtype=np.float64)
    rng = np.random.default_rng(stable_seed("original_benchmark_final_splits", name))
    permutation = rng.permutation(count)
    validation_count = count // 5
    validation = np.sort(permutation[:validation_count])
    test = np.sort(permutation[validation_count:])
    path = RESULTS_DIR / f"evaluation_bank_{name}.npz"
    np.savez_compressed(
        path,
        t=t,
        x=x,
        value=value,
        z=z,
        validation_indices=validation,
        test_indices=test,
    )
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "count": count,
        "validation_count": len(validation),
        "test_count": len(test),
        "sampling_distribution": purpose,
        "reference": reference_meta,
    }


def _prepare_sdgd_bank(count: int = 20_000) -> dict[str, object]:
    benchmark = make_benchmark("hjb_log_d100")
    with np.load(RESULTS_DIR / "evaluation_bank_hjb_log_d100.npz") as common:
        common_t = common["t"]
        common_x = common["x"]
        common_value = common["value"]
    remaining = count - len(common_t)
    if remaining < 0:
        raise ValueError("SDGD bank is smaller than the common bank")
    extra_t, extra_x = benchmark.sample_states(
        remaining, seed=stable_seed("official_sdgd_fixed_test_bank"), purpose="sdgd"
    )
    extra_value = hjb_log_reference(extra_t, extra_x, order=128)
    t = np.concatenate((common_t, extra_t))
    x = np.concatenate((common_x, extra_x))
    value = np.concatenate((common_value, extra_value))
    payload = {
        "X": torch.as_tensor(np.column_stack((x, t)), dtype=torch.float32),
        "Y": torch.as_tensor(value, dtype=torch.float32),
    }
    directory = OFFICIAL_WORK_DIR / "sdgd"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "HJB_Log_100_grad.pkl"
    torch.save(payload, path)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "count": count,
        "common_bank_prefix_count": len(common_t),
        "common_bank_is_exact_prefix": True,
        "labels_used_only_for_official_periodic_L2_reporting": True,
        "training_loss_uses_labels": False,
        "reference": "order-128 independent Hopf-Cole quadrature",
    }


def _prepare_dpi_inputs() -> dict[str, object]:
    source = PACKAGE_DIR / "official_inputs"
    destination = OFFICIAL_WORK_DIR / "dpi_inputs"
    destination.mkdir(parents=True, exist_ok=True)
    records: dict[str, object] = {}
    for name in (
        "mean_100d_ms=1.0_vs=2.0_5.pt",
        "pi_100d_ms=1.0_vs=2.0_5.pt",
    ):
        shutil.copy2(source / name, destination / name)
        records[name] = {
            "sha256": sha256_file(destination / name),
            "origin": "tracked official DPI repository tensor",
        }
    variance = torch.stack(
        [2.0 * torch.eye(100, dtype=torch.float64) for _ in range(5)]
    )
    variance_path = destination / "var_100d_ms=1.0_vs=2.0_5.pt"
    torch.save(variance, variance_path)
    records[variance_path.name] = {
        "sha256": sha256_file(variance_path),
        "origin": "deterministically restored from official constructor var_scale=2",
        "reason": "file referenced by source but absent from official git tree",
    }
    return {"directory": str(destination), "files": records}


def main() -> None:
    ensure_dirs()
    definitions = {
        "source_commits": {
            "DeepBSDE": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
            "DeepPicardIteration": "e8616c2f531c9f95299b246a13c1180b8f02822a",
            "SDGD_PINN": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
            "HTE_PINN": "8b9ca6310cc9d050a54f92528401573413567fe4",
            "SCaSML_FullHistoryMLP": "c2412b942b1720084accf9ac9651c6085e073919",
        },
        "benchmarks": {name: make_benchmark(name).metadata() for name in BENCHMARK_NAMES},
        "excluded": {
            "hjb_rosenbrock": {
                "compatible_driver": True,
                "included": False,
                "reason": (
                    "official repository omits its generated test/reference data; the supplied "
                    "generator uses only 1000 effective MC paths per point and does not quantify "
                    "reference error, so it fails the independent-reference eligibility rule"
                ),
            }
        },
    }
    atomic_json(RESULTS_DIR / "pde_definitions_and_references.json", definitions)
    audits = {
        "hjb_log": _hjb_audit(),
        "dpi_ou_gmm": _dpi_audit(),
        "allen_cahn": _allen_audit(),
    }
    audits["passed"] = all(bool(item["passed"]) for item in audits.values())
    atomic_json(RESULTS_DIR / "reference_accuracy_audit.json", audits)
    banks = {
        "hjb_log_d100": _write_bank("hjb_log_d100", 1200),
        "dpi_ou_gmm_d100": _write_bank("dpi_ou_gmm_d100", 1200),
        "allen_cahn_d100": _write_bank("allen_cahn_d100", 512),
        "sdgd_official_test_bank": _prepare_sdgd_bank(),
        "dpi_official_inputs": _prepare_dpi_inputs(),
    }
    atomic_json(RESULTS_DIR / "prepared_inputs_manifest.json", banks)
    print(json.dumps({"status": "complete", "reference_passed": audits["passed"]}))


if __name__ == "__main__":
    main()

