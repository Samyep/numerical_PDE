"""H100, source-fidelity, equation, and smoke gates for the joint addendum."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.joint_actor import JointTerminalActor
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import _sample_dw, clipped_terminal_loss

from . import deep_bsde, fh_mlp, tangent_actor
from .benchmark import frozen_basis, make_benchmark
from .protocol import (
    DIMENSION,
    JOINT_ACTOR_CONFIG,
    PREFLIGHT_AUDIT,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    SOURCE_VALIDATED_COMMIT,
    array_sha256,
    bank_path,
    file_sha256,
    load_bank,
    require_reference,
    stable_seed,
    write_json_atomic,
)


SOURCE_FILES = (
    "invariant_region_mlp/experiments/single_point_joint_extension/benchmark.py",
    "invariant_region_mlp/experiments/single_point_joint_extension/tangent_actor.py",
    "invariant_region_mlp/experiments/single_point_joint_extension/fh_mlp.py",
    "invariant_region_mlp/experiments/single_point_joint_extension/deep_bsde.py",
    "invariant_region_mlp/experiments/main_highdim_expansion/run_dense_actor.py",
    "invariant_region_mlp/experiments/single_point_32_ice/fh_mlp.py",
    "invariant_region_mlp/experiments/single_point_32_ice/deep_bsde.py",
)


def _git_blob(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_VALIDATED_COMMIT}:{path}"], cwd=SOURCE_ROOT)


def _source_audit() -> dict[str, Any]:
    rows = []
    for relative in SOURCE_FILES:
        # Git's Windows checkout may materialize CRLF while the validated blob
        # is stored with LF.  Compare canonical source bytes, not checkout EOLs.
        current = (SOURCE_ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
        expected = _git_blob(relative)
        rows.append({
            "path": relative,
            "current_sha256": sha256(current).hexdigest(),
            "validated_commit_sha256": sha256(expected).hexdigest(),
            "byte_identical": current == expected,
        })
    return {"overall_pass": all(row["byte_identical"] for row in rows), "rows": rows}


def run(device: torch.device) -> dict[str, Any]:
    reference = require_reference()
    gpu = tangent_actor._require_h100(device)
    calibration = load_bank("calibration", require_truth=True)
    final = load_bank("final", require_truth=True)
    timing = load_bank("timing", require_truth=False)
    driver = make_benchmark(DIMENSION)
    basis = frozen_basis(DIMENSION)

    rng = np.random.default_rng(stable_seed("preflight", "equation"))
    x_np = rng.standard_normal((8, DIMENSION))
    y_np = rng.uniform(0.2, 0.8, 8)
    z_np = rng.standard_normal((8, DIMENSION)) * 0.1
    x = torch.as_tensor(x_np, dtype=torch.float64, device=device)
    y = torch.as_tensor(y_np, dtype=torch.float64, device=device)
    z = torch.as_tensor(z_np, dtype=torch.float64, device=device)
    terminal_torch = driver.terminal(x).detach().cpu().numpy()
    terminal_numpy = driver.terminal_numpy(x_np)
    expected_driver = y_np - y_np**3 - 0.5 * np.square(z_np + y_np[:, None] * driver.v).sum(axis=1)
    actual_driver = driver.driver(y, z).detach().cpu().numpy()
    expected_z = -(z_np + y_np[:, None] * driver.v)
    actual_z = driver.driver_z(y, z).detach().cpu().numpy()

    # Deterministic replay of the exact FH estimator implementation.
    benchmark = fh_mlp.UnprojectedDriver(driver)
    first, first_stats = fh_mlp._source._solve_point(
        benchmark, 0.0, calibration["x"][0], n=2, M=8,
        seed=stable_seed("preflight", "fh_replay"), device=device,
    )
    second, second_stats = fh_mlp._source._solve_point(
        benchmark, 0.0, calibration["x"][0], n=2, M=8,
        seed=stable_seed("preflight", "fh_replay"), device=device,
    )

    # Forward/backward smoke for the exact native Deep-BSDE architecture.
    deep_bsde._configure()
    torch.manual_seed(stable_seed("preflight", "deep") % (2**31))
    model = deep_bsde.FixedInitialDeepBSDE(
        driver, deep_bsde.CONFIG,
        torch.as_tensor(calibration["x"][1], dtype=torch.float64),
    ).to(device=device, dtype=torch.float64)
    generator = torch.Generator(device=device).manual_seed(stable_seed("preflight", "deep_dw"))
    dw = _sample_dw(deep_bsde.CONFIG, DIMENSION, 4, generator=generator, device=device)
    terminal_y, terminal_x = model.rollout(dw)
    loss = clipped_terminal_loss(terminal_y - driver.terminal(terminal_x))
    loss.backward()
    deep_finite = bool(torch.isfinite(loss) and all(
        parameter.grad is None or torch.isfinite(parameter.grad).all()
        for parameter in model.parameters()
    ))

    # Direct and tiny MCV smoke for the exact joint Actor class/evaluator.
    torch.manual_seed(stable_seed("preflight", "actor") % (2**31))
    actor = JointTerminalActor(driver, width=JOINT_ACTOR_CONFIG.width, depth=JOINT_ACTOR_CONFIG.depth).to(device)
    direct, direct_q = tangent_actor._source._predict(
        actor, driver, calibration["t"][:2], calibration["x"][:2], device
    )
    mcv, mcv_meta = tangent_actor._source.evaluate_policy(
        driver, actor, calibration["t"][:2], calibration["x"][:2],
        paths=4, time_steps=2, seed=stable_seed("preflight", "actor_mcv"),
        device=device, parent_batch=2, path_chunk=4,
    )

    source = _source_audit()
    checks = {
        "reference_validated": reference.get("status") == "VALIDATED",
        "gpu_is_h100_80gb": "H100" in gpu.upper() and torch.cuda.get_device_properties(device).total_memory >= 75 * 1024**3,
        "calibration_shape": calibration["x"].shape == (32, DIMENSION),
        "final_shape": final["x"].shape == (64, DIMENSION),
        "timing_shape": timing["x"].shape == (1_048_576, DIMENSION),
        "timing_has_no_truth": "reference" not in timing,
        "basis_orthonormal": float(np.max(np.abs(basis.T @ basis - np.eye(2)))) <= 1e-12,
        "terminal_numpy_torch": float(np.max(np.abs(terminal_torch - terminal_numpy))) <= 1e-12,
        "driver_formula": float(np.max(np.abs(actual_driver - expected_driver))) <= 1e-12,
        "driver_z_formula": float(np.max(np.abs(actual_z - expected_z))) <= 1e-12,
        "mixed_derivative_nonzero": float(np.linalg.norm(driver.v)) >= 0.999999999,
        "fh_exact_replay": first == second and asdict(first_stats) == asdict(second_stats),
        "deep_bsde_forward_backward_finite": deep_finite,
        "joint_actor_direct_finite": bool(np.isfinite(direct).all() and np.isfinite(direct_q).all()),
        "joint_actor_mcv_finite": bool(np.isfinite(mcv["mcv"]).all() and mcv_meta["nonfinite_count"] == 0),
        "validated_sources_byte_identical": source["overall_pass"],
    }
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(checks.values()) else "FAILED",
        "overall_pass": all(checks.values()),
        "gpu": gpu,
        "checks": checks,
        "diagnostics": {
            "terminal_max_abs": float(np.max(np.abs(terminal_torch - terminal_numpy))),
            "driver_max_abs": float(np.max(np.abs(actual_driver - expected_driver))),
            "driver_z_max_abs": float(np.max(np.abs(actual_z - expected_z))),
            "fh_replay_value": first,
            "deep_smoke_loss": float(loss.detach()),
            "actor_direct_sha256": array_sha256(direct, direct_q),
            "actor_mcv_sha256": array_sha256(np.asarray(mcv["mcv"])),
        },
        "source_fidelity": source,
        "bank_files": {name: file_sha256(bank_path(name)) for name in ("calibration", "final", "timing")},
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(PREFLIGHT_AUDIT, payload)
    if not payload["overall_pass"]:
        raise RuntimeError("joint break-even preflight failed")
    return payload


def main() -> None:
    started = time.perf_counter()
    result = run(torch.device("cuda"))
    result["preflight_wall_seconds"] = time.perf_counter() - started
    write_json_atomic(PREFLIGHT_AUDIT, result)
    print(json.dumps({"status": result["status"], "gpu": result["gpu"], "checks": result["checks"]}, sort_keys=True))


if __name__ == "__main__":
    main()
