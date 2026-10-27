"""Mathematical, implementation-fidelity, and isolation gate."""

from __future__ import annotations

import argparse
import copy
import inspect
import json
import math
import platform
import sys
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    _sample_dw,
    clipped_terminal_loss,
)

from .benchmark import frozen_basis, make_benchmark, sample_states
from .deep_bsde import FixedInitialDeepBSDE
from .fh_mlp import AuditedFullHistoryMLP, UnprojectedDriver
from .protocol import (
    DIMENSIONS,
    PACKAGE_DIR,
    PARAMETER_MANIFEST,
    PDE,
    PREFLIGHT_AUDIT,
    QUERY_COUNT,
    RAW_DIR,
    REFERENCE_VALIDATION,
    REVISION,
    file_sha256,
    load_query_bank,
    stable_seed,
    write_json_atomic,
)


def _terminal_and_tangent() -> dict[str, Any]:
    driver = make_benchmark(20)
    rng = np.random.default_rng(stable_seed("preflight", "terminal"))
    x = rng.standard_normal((64, 20))
    direction = rng.standard_normal((64, 20))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    eps = 2e-6
    finite_difference = (
        driver.terminal_numpy(x + eps * direction)
        - driver.terminal_numpy(x - eps * direction)
    ) / (2.0 * eps)
    analytic = np.sum(driver.terminal_gradient_numpy(x) * direction, axis=1)
    terminal_gradient_error = float(np.max(np.abs(finite_difference - analytic)))

    generator = torch.Generator().manual_seed(stable_seed("preflight", "tangent"))
    r = torch.rand(256, dtype=torch.float64, generator=generator)
    q = torch.randn((256, 20), dtype=torch.float64, generator=generator)
    analytic_coefficients = driver.tangent_coefficients(r, q)
    autodiff_coefficients = driver.tangent_coefficients_autodiff(r, q)
    tangent_error = max(
        float(torch.max(torch.abs(left - right)).item())
        for left, right in zip(analytic_coefficients, autodiff_coefficients)
    )
    return {
        "terminal_gradient_directional_max_abs": terminal_gradient_error,
        "terminal_gradient_threshold": 2e-10,
        "tangent_coefficients_max_abs": tangent_error,
        "tangent_coefficients_threshold": 1e-11,
        "passed": terminal_gradient_error <= 2e-10 and tangent_error <= 1e-11,
    }


def _joint_geometry() -> dict[str, Any]:
    rng = np.random.default_rng(stable_seed("preflight", "joint_geometry"))
    y = rng.uniform(0.0, 1.0, 1024)
    # In the span of (value,v), the Hessian is
    # [[-6y-1,-1],[-1,-1]]; all other z eigenvalues equal -1.
    a = -6.0 * y - 1.0
    largest = 0.5 * (a - 1.0 + np.sqrt(np.square(a + 1.0) + 4.0))
    maximum_eigenvalue = float(np.max(largest))
    minimum_eigenvalue = float(np.min(0.5 * (a - 1.0 - np.sqrt(np.square(a + 1.0) + 4.0))))
    mixed_norm = float(np.linalg.norm(make_benchmark(20).v))
    return {
        "states": len(y),
        "maximum_joint_hessian_eigenvalue": maximum_eigenvalue,
        "minimum_joint_hessian_eigenvalue": minimum_eigenvalue,
        "mixed_derivative_norm": mixed_norm,
        "analytic_mixed_derivative": "nabla_z f_y=-v",
        "passed": maximum_eigenvalue <= 1e-12 and mixed_norm >= 0.999999999999,
    }


def _reduction() -> dict[str, Any]:
    driver = make_benchmark(50)
    rng = np.random.default_rng(stable_seed("preflight", "reduction"))
    y = rng.uniform(0.0, 1.0, 256)
    reduced_z = rng.standard_normal((256, 2))
    ambient_z = reduced_z @ driver.B.T
    ambient = driver.driver(
        torch.as_tensor(y, dtype=torch.float64),
        torch.as_tensor(ambient_z, dtype=torch.float64),
    ).numpy()
    reduced = driver.driver_reduced_numpy(y, reduced_z)
    error = float(np.max(np.abs(ambient - reduced)))
    return {"driver_reduction_max_abs": error, "threshold": 1e-12, "passed": error <= 1e-12}


def _query_isolation() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for dimension in DIMENSIONS:
        _, query_x, _, _ = load_query_bank(dimension)
        query_bytes = {np.ascontiguousarray(row).tobytes() for row in query_x}
        overlaps = 0
        checked = 0
        for split, count in (
            ("train_a1", 8192),
            ("train_a2", 8192),
            ("train_a3", 8192),
            ("train_a4", 8192),
            ("fixed_validation", 512),
            ("selector", 128),
        ):
            _, bank = sample_states(
                dimension,
                count,
                seed=stable_seed("ta_states", dimension, split, 0),
            )
            checked += len(bank)
            overlaps += sum(
                np.ascontiguousarray(row).tobytes() in query_bytes for row in bank
            )
        rows.append({
            "dimension": dimension,
            "rows_checked": checked,
            "exact_overlap_count": overlaps,
            "passed": overlaps == 0,
        })
    return {"rows": rows, "passed": all(row["passed"] for row in rows)}


def _fh_replay(device: torch.device) -> dict[str, Any]:
    driver = UnprojectedDriver(make_benchmark(20))
    outputs = []
    stats = []
    for _ in range(2):
        solver = AuditedFullHistoryMLP(
            driver,
            4,
            device,
            seed=stable_seed("preflight", "fh_replay"),
            dtype=torch.float64,
        )
        with torch.no_grad():
            output = solver.solve(
                2,
                torch.zeros(2, dtype=torch.float64, device=device),
                torch.zeros((2, 20), dtype=torch.float64, device=device),
            )
        outputs.append(output.detach().cpu())
        stats.append(vars(solver.stats))
    difference = float(torch.max(torch.abs(outputs[0] - outputs[1])).item())
    return {
        "draw_tape_replay_max_abs": difference,
        "stats_identical": stats[0] == stats[1],
        "finite": bool(torch.all(torch.isfinite(outputs[0])).item()),
        "passed": difference == 0.0 and stats[0] == stats[1],
    }


def _step(model, optimizer, dw) -> None:
    terminal_y, terminal_x = model.rollout(dw)
    loss = clipped_terminal_loss(terminal_y - model.benchmark.terminal(terminal_x))
    optimizer.zero_grad(set_to_none=True)
    loss.backward()
    optimizer.step()


def _deep_resume(device: torch.device) -> dict[str, Any]:
    driver = make_benchmark(4)
    config = DeepBSDEConfig(
        benchmark=f"{PDE}_d4",
        total_time=driver.T,
        intervals=2,
        hidden_widths=(8, 8),
        learning_rate=1e-2,
        iterations=2,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.0, 1.0),
    )
    seed_everything(stable_seed("preflight", "deep_init"))
    template = FixedInitialDeepBSDE(driver, config, torch.zeros(4)).to(
        device=device, dtype=torch.float64
    )
    initial = copy.deepcopy(template.state_dict())

    continuous = FixedInitialDeepBSDE(driver, config, torch.zeros(4)).to(
        device=device, dtype=torch.float64
    )
    continuous.load_state_dict(initial)
    optimizer = torch.optim.Adam(continuous.parameters(), lr=1e-2, eps=1e-8)
    generator = torch.Generator(device=device).manual_seed(
        stable_seed("preflight", "deep_dw")
    )
    first = _sample_dw(config, 4, 4, generator=generator, device=device)
    _step(continuous, optimizer, first)
    checkpoint_model = copy.deepcopy(continuous.state_dict())
    checkpoint_optimizer = copy.deepcopy(optimizer.state_dict())
    checkpoint_generator = generator.get_state()
    second = _sample_dw(config, 4, 4, generator=generator, device=device)
    _step(continuous, optimizer, second)

    resumed = FixedInitialDeepBSDE(driver, config, torch.zeros(4)).to(
        device=device, dtype=torch.float64
    )
    resumed_optimizer = torch.optim.Adam(resumed.parameters(), lr=1e-2, eps=1e-8)
    resumed.load_state_dict(checkpoint_model)
    resumed_optimizer.load_state_dict(checkpoint_optimizer)
    resumed_generator = torch.Generator(device=device)
    resumed_generator.set_state(checkpoint_generator)
    repeated = _sample_dw(config, 4, 4, generator=resumed_generator, device=device)
    _step(resumed, resumed_optimizer, repeated)
    maximum = max(
        float(torch.max(torch.abs(continuous.state_dict()[key] - resumed.state_dict()[key])).item())
        for key in continuous.state_dict()
    )
    return {
        "second_draw_bitwise_equal": bool(torch.equal(second, repeated)),
        "maximum_parameter_difference": maximum,
        "threshold": 1e-12,
        "passed": torch.equal(second, repeated) and maximum <= 1e-12,
    }


def _reference_access_audit() -> dict[str, Any]:
    from invariant_region_mlp.experiments.single_point_32_ice import deep_bsde, fh_mlp
    from invariant_region_mlp.experiments.main_highdim_expansion import run_dense_actor

    deep = inspect.getsource(deep_bsde.execute)
    fh = inspect.getsource(fh_mlp.execute)
    ta_train = inspect.getsource(run_dense_actor.train)
    ta_eval = inspect.getsource(run_dense_actor.evaluate)
    checks = {
        "deep_reference_after_completed_schedule": (
            deep.index("reference =") > deep.index("if completed_updates < config.iterations")
        ),
        "fh_reference_after_solver_loop": fh.index("reference =") > fh.index("for point in range"),
        "ta_training_has_no_test_bank": "load_test_bank" not in ta_train,
        "ta_truth_only_in_evaluation": "load_test_bank" in ta_eval,
    }
    return {"checks": checks, "passed": all(checks.values())}


def run(device: torch.device) -> dict[str, Any]:
    started = time.perf_counter()
    reference = json.loads(REFERENCE_VALIDATION.read_text(encoding="utf-8"))
    parameters = json.loads(PARAMETER_MANIFEST.read_text(encoding="utf-8"))
    formal_root = RAW_DIR / "formal"
    formal_files = list(formal_root.rglob("*")) if formal_root.exists() else []
    checks: dict[str, Any] = {
        "reference": {
            "passed": reference.get("overall_pass") is True,
            "manifest_sha256": file_sha256(REFERENCE_VALIDATION),
            "gates": reference.get("gates"),
        },
        "parameters": {
            "passed": parameters.get("overall_pass") is True,
            "manifest_sha256": file_sha256(PARAMETER_MANIFEST),
        },
        "terminal_and_tangent": _terminal_and_tangent(),
        "joint_geometry": _joint_geometry(),
        "ambient_reduction": _reduction(),
        "query_isolation": _query_isolation(),
        "fh_draw_tape_replay": _fh_replay(device),
        "deep_bsde_resume": _deep_resume(device),
        "reference_access": _reference_access_audit(),
        "no_formal_outputs_before_gate": {
            "count": sum(path.is_file() for path in formal_files),
            "passed": not any(path.is_file() for path in formal_files),
        },
    }
    mathematical_pass = all(item["passed"] for item in checks.values())
    cuda_pass = device.type == "cuda" and torch.cuda.is_available()
    overall = mathematical_pass and cuda_pass
    code_files = (
        "protocol.py", "benchmark.py", "prepare.py", "fh_mlp.py",
        "deep_bsde.py", "tangent_actor.py", "preflight.py",
    )
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if overall else (
            "CPU_VALIDATED_PENDING_CUDA" if mathematical_pass else "FAILED"
        ),
        "overall_pass": overall,
        "mathematical_checks_pass": mathematical_pass,
        "device": str(device),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu_name": torch.cuda.get_device_name(device) if cuda_pass else None,
        },
        "checks": checks,
        "code_sha256": {name: file_sha256(PACKAGE_DIR / name) for name in code_files},
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(PREFLIGHT_AUDIT, payload)
    if not mathematical_pass:
        raise RuntimeError("joint-extension preflight failed")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    result = run(device)
    print(json.dumps({
        "status": result["status"],
        "overall_pass": result["overall_pass"],
        "mathematical_checks_pass": result["mathematical_checks_pass"],
        "elapsed_seconds": result["elapsed_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["run"]
