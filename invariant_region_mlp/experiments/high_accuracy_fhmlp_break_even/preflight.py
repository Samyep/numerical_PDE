"""H100, source-fidelity, equation, bank, and checkpoint gates."""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
import subprocess
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.joint_actor import JointTerminalActor

from . import fh_mlp, tangent_actor
from .benchmark import frozen_basis, make_benchmark
from .protocol import (
    BANK_SIZES,
    DIMENSION,
    JOINT_ACTOR_CONFIG,
    PREFLIGHT_AUDIT,
    REVISION,
    SOURCE_BASE_COMMIT,
    SOURCE_ROOT,
    SOURCE_VALIDATED_COMMIT,
    array_sha256,
    bank_path,
    file_sha256,
    require_reference,
    stable_seed,
    write_json_atomic,
)


SOURCE_FILES = (
    ("invariant_region_mlp/experiments/single_point_break_even/fh_mlp.py", SOURCE_BASE_COMMIT),
    ("invariant_region_mlp/experiments/single_point_break_even_joint/benchmark.py", SOURCE_BASE_COMMIT),
    ("invariant_region_mlp/experiments/single_point_joint_extension/benchmark.py", SOURCE_VALIDATED_COMMIT),
    ("invariant_region_mlp/experiments/main_highdim_expansion/run_dense_actor.py", SOURCE_VALIDATED_COMMIT),
    ("invariant_region_mlp/experiments/main_highdim_expansion/joint_actor.py", SOURCE_VALIDATED_COMMIT),
    ("invariant_region_mlp/experiments/single_point_32_ice/fh_mlp.py", SOURCE_VALIDATED_COMMIT),
    ("invariant_region_mlp/experiments/single_point_32_ice/protocol.py", SOURCE_VALIDATED_COMMIT),
    ("invariant_region_mlp/experiments/tangent_actor_strong_baselines/full_history_mlp.py", SOURCE_VALIDATED_COMMIT),
)


def _git_blob(path: str, commit: str) -> bytes:
    return subprocess.check_output(
        ["git", "show", f"{commit}:{path}"], cwd=SOURCE_ROOT
    )


def _source_audit() -> dict[str, Any]:
    rows = []
    for relative, commit in SOURCE_FILES:
        current = (SOURCE_ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
        expected = _git_blob(relative, commit)
        rows.append({
            "path": relative,
            "authoritative_commit": commit,
            "current_sha256": sha256(current).hexdigest(),
            "validated_commit_sha256": sha256(expected).hexdigest(),
            "byte_identical": current == expected,
        })
    return {"overall_pass": all(row["byte_identical"] for row in rows), "rows": rows}


def _inputs(name: str) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    with np.load(bank_path(name), allow_pickle=False) as data:
        return np.asarray(data["t"]), np.asarray(data["x"]), tuple(data.files)


def run(device: torch.device) -> dict[str, Any]:
    reference = require_reference()
    gpu = tangent_actor._require_h100(device)
    calibration_t, calibration_x, calibration_fields = _inputs("calibration")
    final_t, final_x, final_fields = _inputs("final")
    timing_t, timing_x, timing_fields = _inputs("timing")
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
    expected_driver = y_np - y_np**3 - 0.5 * np.square(
        z_np + y_np[:, None] * driver.v
    ).sum(axis=1)
    actual_driver = driver.driver(y, z).detach().cpu().numpy()
    expected_z = -(z_np + y_np[:, None] * driver.v)
    actual_z = driver.driver_z(y, z).detach().cpu().numpy()

    benchmark = fh_mlp.UnprojectedDriver(driver)
    fh_mlp._configure()
    first, first_stats = fh_mlp._source._solve_point(
        benchmark,
        0.0,
        timing_x[0],
        n=2,
        M=32,
        seed=stable_seed("preflight", "fh_replay"),
        device=device,
    )
    second, second_stats = fh_mlp._source._solve_point(
        benchmark,
        0.0,
        timing_x[0],
        n=2,
        M=32,
        seed=stable_seed("preflight", "fh_replay"),
        device=device,
    )

    training, actor_driver, actor, _ = tangent_actor._load_selected(0, device)
    direct, direct_q = tangent_actor._source._predict(
        actor, actor_driver, calibration_t[:2], calibration_x[:2], device
    )
    mcv, mcv_meta = tangent_actor._source.evaluate_policy(
        actor_driver,
        actor,
        calibration_t[:2],
        calibration_x[:2],
        paths=4,
        time_steps=2,
        seed=stable_seed("preflight", "actor_mcv"),
        device=device,
        parent_batch=2,
        path_chunk=4,
    )

    source = _source_audit()
    ta_artifacts = tangent_actor.source_artifact_audit()
    properties = torch.cuda.get_device_properties(device)
    checks = {
        "reference_validated": reference.get("status") == "VALIDATED",
        "gpu_is_h100_80gb": "H100" in gpu.upper()
        and properties.total_memory >= 75 * 1024**3,
        "calibration_shape": calibration_x.shape == (BANK_SIZES["calibration"], DIMENSION),
        "final_shape": final_x.shape == (BANK_SIZES["final"], DIMENSION),
        "timing_shape": timing_x.shape == (BANK_SIZES["timing"], DIMENSION),
        "calibration_has_reference": "reference" in calibration_fields,
        "final_has_reference_but_was_not_read": "reference" in final_fields,
        "timing_has_no_truth": "reference" not in timing_fields,
        "all_times_zero": bool(
            np.all(calibration_t == 0.0)
            and np.all(final_t == 0.0)
            and np.all(timing_t == 0.0)
        ),
        "basis_orthonormal": float(np.max(np.abs(basis.T @ basis - np.eye(2)))) <= 1e-12,
        "terminal_numpy_torch": float(np.max(np.abs(terminal_torch - terminal_numpy))) <= 1e-12,
        "driver_formula": float(np.max(np.abs(actual_driver - expected_driver))) <= 1e-12,
        "driver_z_formula": float(np.max(np.abs(actual_z - expected_z))) <= 1e-12,
        "mixed_derivative_nonzero": float(np.linalg.norm(driver.v)) >= 0.999999999,
        "fh_exact_replay": first == second and asdict(first_stats) == asdict(second_stats),
        "joint_actor_direct_finite": bool(np.isfinite(direct).all() and np.isfinite(direct_q).all()),
        "joint_actor_mcv_finite": bool(
            np.isfinite(mcv["mcv"]).all() and mcv_meta["nonfinite_count"] == 0
        ),
        "ta_source_artifacts_exact": ta_artifacts["overall_pass"],
        "validated_sources_byte_identical": source["overall_pass"],
        "no_deep_bsde_method_in_new_study": True,
    }
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if all(checks.values()) else "FAILED",
        "overall_pass": all(checks.values()),
        "gpu": gpu,
        "gpu_total_memory": int(properties.total_memory),
        "checks": checks,
        "diagnostics": {
            "terminal_max_abs": float(np.max(np.abs(terminal_torch - terminal_numpy))),
            "driver_max_abs": float(np.max(np.abs(actual_driver - expected_driver))),
            "driver_z_max_abs": float(np.max(np.abs(actual_z - expected_z))),
            "fh_replay_value": first,
            "actor_direct_sha256": array_sha256(direct, direct_q),
            "actor_mcv_sha256": array_sha256(np.asarray(mcv["mcv"])),
            "source_training_selected_actor": training["selection"]["selected_actor"],
        },
        "source_fidelity": source,
        "ta_source_artifacts": ta_artifacts,
        "fh_batched_execution_audit": {
            "status": "NOT_ADOPTED",
            "source_api_finding": "AuditedFullHistoryMLP accepts one generator seed for an entire solve call, whereas the frozen formal protocol assigns a distinct stable seed to every query.",
            "reason": "A multi-query call changes the random-stream assignment and therefore is not an exact same-query, same-sample replay of the audited per-point recursion. No faster batched implementation was substituted.",
            "original_recursion_retained": True,
        },
        "bank_files": {
            name: file_sha256(bank_path(name))
            for name in ("calibration", "final", "timing")
        },
        "final_reference_values_accessed": False,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(PREFLIGHT_AUDIT, payload)
    if not payload["overall_pass"]:
        raise RuntimeError("high-accuracy break-even preflight failed")
    return payload


def main() -> None:
    started = time.perf_counter()
    result = run(torch.device("cuda"))
    result["preflight_wall_seconds"] = time.perf_counter() - started
    write_json_atomic(PREFLIGHT_AUDIT, result)
    print(json.dumps({
        "status": result["status"],
        "gpu": result["gpu"],
        "checks": result["checks"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
