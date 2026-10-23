"""Stage-A mathematical, fidelity, isolation, and runtime gates."""

from __future__ import annotations

import argparse
import copy
from dataclasses import asdict
import inspect
import json
import math
import platform
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
import scipy
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    OfficialDeepBSDE,
    _sample_dw,
    clipped_terminal_loss,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.models import (
    TerminalAwareActor,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.tangent_solver import (
    antithetic_normals,
)

from .benchmarks import make_benchmark
from .deep_bsde import FixedInitialDeepBSDE
from .fh_mlp import AuditedFullHistoryMLP, _UnprojectedDriver
from .prepare import prepare
from .protocol import (
    DIMENSIONS,
    METHOD_FIDELITY_AUDIT,
    PACKAGE_DIR,
    PARAMETER_MANIFEST,
    PDES,
    QUERY_COUNT,
    RAW_DIR,
    REFERENCE_REVALIDATION,
    REFERENCE_VALIDATION,
    REVISION,
    SOURCE_ROOT,
    STAGE_A_AUDIT,
    TRAINING_SEEDS,
    array_sha256,
    file_sha256,
    load_query_bank,
    query_path,
    stable_seed,
    write_json_atomic,
)
from .tangent_actor import (
    _configure_b6,
    _configure_dense,
    _configure_standard,
    _value_mcv_rollout_chunk,
)


def _row_bytes(array: np.ndarray) -> set[bytes]:
    value = np.ascontiguousarray(array, dtype=np.float64)
    return {row.tobytes() for row in value}


def _query_isolation() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for pde in PDES:
        for d in DIMENSIONS:
            query_rows = _row_bytes(np.load(query_path(pde, d) / "x.npy", allow_pickle=False))
            checked = 0
            overlaps = 0
            for seed in TRAINING_SEEDS:
                if pde in {"hjb_log", "allen_cahn"}:
                    impl = _configure_standard(False)
                    benchmark = impl.make_benchmark(pde, d)
                    splits = [(f"train_a{i}", 8192) for i in range(1, 5)]
                    splits += [("fixed_validation", 512), ("fixed_point", 128)]
                    banks = [
                        impl._states(pde, benchmark, count, split=split, seed=seed)[1]
                        for split, count in splits
                    ]
                elif pde == "b6_fullrank":
                    impl = _configure_b6(False)
                    benchmark = impl.make_benchmark(pde, d)
                    splits = [(f"train_a{i}", 8192) for i in range(1, 5)]
                    splits += [("fixed_validation", 512), ("fixed_point", 128)]
                    banks = [
                        impl._states(pde, benchmark, count, split=split, seed=seed)[1]
                        for split, count in splits
                    ]
                else:
                    impl = _configure_dense(False, seed)
                    splits = [(f"train_a{i}", 8192) for i in range(1, 5)]
                    splits += [("fixed_validation", 512), ("selector", 128)]
                    banks = [
                        impl._states(d, count, split=split, seed=seed)[1]
                        for split, count in splits
                    ]
                for bank in banks:
                    checked += len(bank)
                    overlaps += len(query_rows.intersection(_row_bytes(bank)))
            rows.append({
                "pde": pde,
                "dimension": d,
                "training_validation_selection_rows_checked": checked,
                "exact_query_overlap_count": overlaps,
                "passed": overlaps == 0,
            })
    return {"rows": rows, "passed": all(row["passed"] for row in rows)}


def _native_rollout_anchor(pde: str, device: torch.device) -> dict[str, Any]:
    benchmark = make_benchmark(pde, 100)
    allen = pde == "allen_cahn"
    config = DeepBSDEConfig(
        benchmark=f"{pde}_d100",
        total_time=0.3 if allen else 1.0,
        intervals=20,
        hidden_widths=(110, 110),
        learning_rate=5e-4 if allen else 1e-2,
        iterations=1,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.3, 0.5) if allen else (0.0, 1.0),
    )
    seed_everything(stable_seed("native_anchor", pde))
    native = OfficialDeepBSDE(benchmark, config).to(device=device, dtype=torch.float64)
    fixed = FixedInitialDeepBSDE(benchmark, config, torch.zeros(100)).to(
        device=device, dtype=torch.float64
    )
    missing, unexpected = fixed.load_state_dict(native.state_dict(), strict=False)
    if set(missing) != {"fixed_x0"} or unexpected:
        raise RuntimeError((missing, unexpected))
    native.eval()
    fixed.eval()
    generator = torch.Generator(device=device).manual_seed(stable_seed("native_anchor_dw", pde))
    dw = _sample_dw(config, 100, 4, generator=generator, device=device)
    with torch.no_grad():
        native_y, native_x = native.rollout(dw)
        fixed_y, fixed_x = fixed.rollout(dw)
    y_error = float(torch.max(torch.abs(native_y - fixed_y)).item())
    x_error = float(torch.max(torch.abs(native_x - fixed_x)).item())
    return {
        "pde": pde,
        "dimension": 100,
        "x0": "zero",
        "max_y_rollout_difference": y_error,
        "max_x_rollout_difference": x_error,
        "passed": y_error == 0.0 and x_error == 0.0,
    }


def _optimizer_step(model, optimizer, dw) -> float:
    model.train()
    terminal_y, terminal_x = model.rollout(dw)
    loss = clipped_terminal_loss(terminal_y - model.benchmark.terminal(terminal_x))
    optimizer.zero_grad(set_to_none=True)
    loss.backward()
    optimizer.step()
    return float(loss.detach())


def _resume_equivalence(device: torch.device) -> dict[str, Any]:
    benchmark = make_benchmark("hjb_log", 2)
    config = DeepBSDEConfig(
        benchmark="hjb_log_d2",
        total_time=1.0,
        intervals=2,
        hidden_widths=(4, 4),
        learning_rate=1e-2,
        iterations=2,
        batch_size=4,
        validation_size=4,
        logging_frequency=1,
        y_init_range=(0.0, 1.0),
    )
    seed_everything(stable_seed("resume_equivalence_init"))
    template = FixedInitialDeepBSDE(benchmark, config, torch.tensor([0.3, -0.4])).to(
        device=device, dtype=torch.float64
    )
    initial = {key: value.detach().clone() for key, value in template.state_dict().items()}

    continuous = FixedInitialDeepBSDE(benchmark, config, torch.tensor([0.3, -0.4])).to(
        device=device, dtype=torch.float64
    )
    continuous.load_state_dict(initial)
    continuous_optimizer = torch.optim.Adam(continuous.parameters(), lr=1e-2, eps=1e-8)
    continuous_generator = torch.Generator(device=device).manual_seed(stable_seed("resume_dw"))
    first_dw = _sample_dw(config, 2, 4, generator=continuous_generator, device=device)
    _optimizer_step(continuous, continuous_optimizer, first_dw)
    saved_generator_state = continuous_generator.get_state()
    second_dw = _sample_dw(config, 2, 4, generator=continuous_generator, device=device)
    _optimizer_step(continuous, continuous_optimizer, second_dw)

    interrupted = FixedInitialDeepBSDE(benchmark, config, torch.tensor([0.3, -0.4])).to(
        device=device, dtype=torch.float64
    )
    interrupted.load_state_dict(initial)
    interrupted_optimizer = torch.optim.Adam(interrupted.parameters(), lr=1e-2, eps=1e-8)
    interrupted_generator = torch.Generator(device=device).manual_seed(stable_seed("resume_dw"))
    first_repeated = _sample_dw(config, 2, 4, generator=interrupted_generator, device=device)
    if not torch.equal(first_dw, first_repeated):
        raise RuntimeError("resume test first Brownian batch mismatch")
    _optimizer_step(interrupted, interrupted_optimizer, first_repeated)
    saved_model = {key: value.detach().clone() for key, value in interrupted.state_dict().items()}
    saved_optimizer = copy.deepcopy(interrupted_optimizer.state_dict())

    resumed = FixedInitialDeepBSDE(benchmark, config, torch.tensor([0.3, -0.4])).to(
        device=device, dtype=torch.float64
    )
    resumed_optimizer = torch.optim.Adam(resumed.parameters(), lr=1e-2, eps=1e-8)
    resumed.load_state_dict(saved_model)
    resumed_optimizer.load_state_dict(saved_optimizer)
    resumed_generator = torch.Generator(device=device)
    resumed_generator.set_state(saved_generator_state)
    second_repeated = _sample_dw(config, 2, 4, generator=resumed_generator, device=device)
    if not torch.equal(second_dw, second_repeated):
        raise RuntimeError("resume test second Brownian batch mismatch")
    _optimizer_step(resumed, resumed_optimizer, second_repeated)
    maximum = max(
        float(torch.max(torch.abs(continuous.state_dict()[key] - resumed.state_dict()[key])).item())
        for key in continuous.state_dict()
    )
    threshold = 0.0 if device.type == "cpu" else 1e-12
    return {
        "device": str(device),
        "maximum_parameter_difference": maximum,
        "threshold": threshold,
        "generator_continuation_bitwise_equal": True,
        "passed": maximum <= threshold,
    }


def _fh_joint_audit(device: torch.device) -> dict[str, Any]:
    driver = _UnprojectedDriver(make_benchmark("dense_fourier", 20))
    solver = AuditedFullHistoryMLP(
        driver, 4, device, seed=stable_seed("fh_joint_audit"), dtype=torch.float64
    )
    with torch.no_grad():
        state = solver.solve(
            2,
            torch.zeros(1, dtype=torch.float64, device=device),
            torch.zeros((1, 20), dtype=torch.float64, device=device),
        )
    y = torch.tensor([0.2, 0.6], dtype=torch.float64, device=device)
    z = torch.ones((2, 20), dtype=torch.float64, device=device)
    base = driver.driver(y, z)
    changed_value = driver.driver(y + 0.1, z)
    changed_gradient = driver.driver(y, 2.0 * z)
    source = inspect.getsource(AuditedFullHistoryMLP.solve)
    passed = bool(
        torch.all(torch.isfinite(state)).item()
        and not torch.equal(base, changed_value)
        and torch.allclose(changed_gradient, 4.0 * base)
        and "1e-7" not in source
    )
    return {
        "state_finite": bool(torch.all(torch.isfinite(state)).item()),
        "driver_depends_on_value": not torch.equal(base, changed_value),
        "driver_quadratic_gradient_scaling": bool(torch.allclose(changed_gradient, 4.0 * base)),
        "untruncated_bismut_denominator": "1e-7" not in source,
        "stats": asdict(solver.stats),
        "passed": passed,
    }


def _value_mcv_audit(device: torch.device) -> dict[str, Any]:
    benchmark = make_benchmark("allen_cahn", 20).impl
    seed_everything(stable_seed("value_mcv_actor"))
    actor = TerminalAwareActor(benchmark, width=16, depth=2, activation="silu").to(device)
    actor.eval()
    generator = torch.Generator(device=device).manual_seed(stable_seed("value_mcv_paths"))
    normals = antithetic_normals(1, 1024, 8, 20, generator=generator, device=device)
    plain, corrected, _, martingale, diagnostics = _value_mcv_rollout_chunk(
        benchmark,
        actor,
        torch.zeros(1, dtype=torch.float32, device=device),
        torch.zeros((1, 20), dtype=torch.float32, device=device),
        normals,
        use_mcv=True,
    )
    samples = martingale.double().reshape(-1)
    mean = float(samples.mean())
    se = float(samples.std(unbiased=True) / math.sqrt(len(samples)))
    identity_error = float(torch.max(torch.abs(corrected - (plain - martingale))).item())
    return {
        "martingale_mean": mean,
        "martingale_standard_error": se,
        "mean_z_score": mean / max(se, 1e-300),
        "pathwise_identity_max_error": identity_error,
        "nonfinite_count": diagnostics["nonfinite_count"],
        "passed": abs(mean) <= 6.0 * se + 1e-12 and identity_error == 0.0
        and diagnostics["nonfinite_count"] == 0,
    }


def _copy_fidelity_audit() -> str:
    source = PACKAGE_DIR / "METHOD_FIDELITY_AUDIT.md"
    content = source.read_text(encoding="utf-8")
    METHOD_FIDELITY_AUDIT.parent.mkdir(parents=True, exist_ok=True)
    temporary = METHOD_FIDELITY_AUDIT.with_suffix(".md.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(METHOD_FIDELITY_AUDIT)
    return file_sha256(METHOD_FIDELITY_AUDIT)


def run(device: torch.device, *, generate_inputs: bool = True) -> dict[str, Any]:
    started = time.perf_counter()
    if generate_inputs:
        prepare()
    reference = json.loads(REFERENCE_VALIDATION.read_text(encoding="utf-8"))
    revalidation = json.loads(
        REFERENCE_REVALIDATION.read_text(encoding="utf-8")
    )
    parameters = json.loads(PARAMETER_MANIFEST.read_text(encoding="utf-8"))
    formal_root = RAW_DIR / "formal"
    formal_artifacts = list(formal_root.rglob("*")) if formal_root.exists() else []
    checks: dict[str, Any] = {
        "references": {
            "passed": reference.get("overall_pass") is True,
            "cells": len(reference.get("banks", [])),
            "manifest_sha256": file_sha256(REFERENCE_VALIDATION),
        },
        "reference_revalidation": {
            "passed": revalidation.get("overall_pass") is True,
            "cells": len(revalidation.get("banks", [])),
            "manifest_sha256": file_sha256(REFERENCE_REVALIDATION),
            "maximum_frozen_vs_recomputed_reference_difference": max(
                float(row["frozen_vs_recomputed_reference_max_abs"])
                for row in revalidation.get("banks", [])
            ),
            "maximum_frozen_vs_recomputed_uncertainty_difference": max(
                float(row["frozen_vs_recomputed_uncertainty_max_abs"])
                for row in revalidation.get("banks", [])
            ),
        },
        "parameters": {
            "passed": all(
                row.get("byte_identical_arrays_to_completed_ice_expansion", True)
                for row in parameters["rows"]
            ),
            "cells": len(parameters["rows"]),
            "manifest_sha256": file_sha256(PARAMETER_MANIFEST),
        },
        "query_isolation": _query_isolation(),
        "deep_bsde_native_anchors": [
            _native_rollout_anchor("hjb_log", device),
            _native_rollout_anchor("allen_cahn", device),
        ],
        "deep_bsde_resume": _resume_equivalence(device),
        "fh_mlp_joint": _fh_joint_audit(device),
        "allen_value_actor_mcv": _value_mcv_audit(device),
        "no_formal_artifacts_before_gate": {
            "count": len([path for path in formal_artifacts if path.is_file()]),
            "passed": not any(path.is_file() for path in formal_artifacts),
        },
    }
    anchors_pass = all(item["passed"] for item in checks["deep_bsde_native_anchors"])
    overall = bool(
        checks["references"]["passed"]
        and checks["reference_revalidation"]["passed"]
        and checks["parameters"]["passed"]
        and checks["query_isolation"]["passed"]
        and anchors_pass
        and checks["deep_bsde_resume"]["passed"]
        and checks["fh_mlp_joint"]["passed"]
        and checks["allen_value_actor_mcv"]["passed"]
        and checks["no_formal_artifacts_before_gate"]["passed"]
    )
    fidelity_hash = _copy_fidelity_audit()
    code_files = [
        "protocol.py", "benchmarks.py", "prepare.py", "deep_bsde.py",
        "fh_mlp.py", "tangent_actor.py", "preflight.py",
    ]
    payload = {
        "revision": REVISION,
        "status": "VALIDATED" if overall and device.type == "cuda" else (
            "CPU_VALIDATED_PENDING_CUDA" if overall else "FAILED"
        ),
        "overall_pass": overall and device.type == "cuda",
        "mathematical_checks_pass": overall,
        "device": str(device),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu_name": torch.cuda.get_device_name(device) if device.type == "cuda" else None,
        },
        "checks": checks,
        "method_fidelity_audit_sha256": fidelity_hash,
        "code_sha256": {
            name: file_sha256(PACKAGE_DIR / name) for name in code_files
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(STAGE_A_AUDIT, payload)
    if not overall:
        raise RuntimeError("Stage-A mathematical preflight failed")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--skip-prepare", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA preflight requested but CUDA is unavailable")
    result = run(device, generate_inputs=not args.skip_prepare)
    print(json.dumps({
        "status": result["status"],
        "overall_pass": result["overall_pass"],
        "mathematical_checks_pass": result["mathematical_checks_pass"],
        "elapsed_seconds": result["elapsed_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

