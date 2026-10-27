"""Build disjoint query banks and apply the frozen Allen--Cahn truth gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    AllenRadialConfig,
    AllenRadialReference,
)
from invariant_region_mlp.experiments.single_point_32_ice.prepare import (
    _allen_reference as audited_allen_reference,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.benchmarks import (
    AllenCahnBenchmark,
)

from .protocol import (
    BANK_DIR,
    BANK_SIZES,
    DIMENSION,
    HORIZON,
    NAMESPACE_SEED,
    RESULTS_DIR,
    REVISION,
    SOURCE_ROOT,
    array_sha256,
    bank_path,
    ensure_directories,
    file_sha256,
    protocol_payload,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


OLD_BANK = (
    SOURCE_ROOT
    / "invariant_region_mlp"
    / "results"
    / "single_point_32_ice"
    / "frozen_query_banks"
    / "allen_cahn_d100"
)


def _row_bytes(x: np.ndarray) -> set[bytes]:
    value = np.ascontiguousarray(x, dtype=np.float64)
    return {row.tobytes() for row in value}


def _new_states(name: str, count: int) -> tuple[np.ndarray, int]:
    bank_seed = stable_seed("query_bank", name)
    generator = np.random.Generator(np.random.PCG64(bank_seed))
    x = generator.standard_normal((count, DIMENSION), dtype=np.float64)
    if name != "timing":
        x[0] = 0.0
    return x, bank_seed


def _assert_bank_separation(banks: dict[str, np.ndarray], old_x: np.ndarray) -> dict[str, Any]:
    """Audit exact row separation, except for the protocol-mandated origin.

    The frozen prose simultaneously mandates an origin in calibration/final
    and exact disjointness from each other and the old bank, which also has an
    origin.  Those requirements cannot all hold.  We preserve the mandated
    origin and make the unavoidable exception explicit rather than concealing
    it; every stochastic row remains exactly disjoint.
    """
    origin = np.zeros(DIMENSION, dtype=np.float64).tobytes()
    sets = {name: _row_bytes(value) for name, value in banks.items()}
    old = _row_bytes(old_x)
    rows: list[dict[str, Any]] = []
    for left, right, a, b in (
        ("calibration", "final", sets["calibration"], sets["final"]),
        ("calibration", "timing", sets["calibration"], sets["timing"]),
        ("final", "timing", sets["final"], sets["timing"]),
        ("calibration", "previous_32", sets["calibration"], old),
        ("final", "previous_32", sets["final"], old),
        ("timing", "previous_32", sets["timing"], old),
    ):
        overlap = a & b
        permitted = {origin} if left in {"calibration", "final"} and right in {
            "calibration", "final", "previous_32"
        } else set()
        unexpected = overlap - permitted
        rows.append({
            "left": left,
            "right": right,
            "overlap_count": len(overlap),
            "mandated_origin_overlap": bool(origin in overlap),
            "unexpected_overlap_count": len(unexpected),
        })
        if unexpected:
            raise RuntimeError(f"unexpected exact query overlap: {left} / {right}")
    return {
        "status": "VALIDATED_WITH_MANDATED_ORIGIN_EXCEPTION",
        "internal_protocol_conflict": (
            "calibration, final, and the previous bank are each required to "
            "contain x=0, so literal pairwise state disjointness is impossible"
        ),
        "comparisons": rows,
    }


def _residual_at_scale(
    solver: AllenRadialReference,
    t: np.ndarray,
    s: np.ndarray,
    *,
    ht: float,
    hs: float,
) -> tuple[np.ndarray, np.ndarray]:
    count = len(t)
    query_t = np.concatenate((t, t + ht, t - ht, t, t))
    query_s = np.concatenate((s, s, s, s + hs, s - hs))
    values = solver.evaluate(query_t, query_s).reshape(5, count)
    center, later, earlier, right, left = values
    u_t = (later - earlier) / (2.0 * ht)
    u_s = (right - left) / (2.0 * hs)
    u_ss = (right - 2.0 * center + left) / (hs * hs)
    reaction = center - center**3
    residual = u_t + 4.0 * s * u_ss + 2.0 * DIMENSION * u_s + reaction
    scale = np.abs(u_t) + np.abs(4.0 * s * u_ss) + np.abs(2.0 * DIMENSION * u_s) + np.abs(reaction)
    return residual, np.maximum(scale, 1e-12)


def _pde_identity_audit() -> dict[str, Any]:
    seed = stable_seed("reference", "pde_residual_audit")
    generator = np.random.Generator(np.random.PCG64(seed))
    x = generator.standard_normal((200, DIMENSION), dtype=np.float64)
    s = np.square(x).sum(axis=1)
    t = generator.uniform(0.01, HORIZON - 0.01, len(s)).astype(np.float64)
    s_max = math.ceil((float(s.max()) + 20.0) / 10.0) * 10.0
    config = AllenRadialConfig(
        s_max=s_max,
        ds=0.0125,
        dt=1.25e-5,
        dimension=DIMENSION,
        horizon=HORIZON,
    )
    solver = AllenRadialReference(config)
    fine, fine_scale = _residual_at_scale(solver, t, s, ht=0.001, hs=0.2)
    coarse, _ = _residual_at_scale(solver, t, s, ht=0.002, hs=0.4)
    relative = np.abs(fine) / fine_scale

    # At t=T the radial grid is initialized with the exact terminal function.
    terminal_indices = generator.integers(0, len(solver.s), size=200)
    terminal_s = solver.s[terminal_indices]
    terminal = solver.evaluate(np.full(200, HORIZON), terminal_s)
    terminal_exact = 0.5 / (1.0 + 0.2 * terminal_s)
    terminal_solver_max = float(np.max(np.abs(terminal - terminal_exact)))

    benchmark = AllenCahnBenchmark(DIMENSION)
    torch_terminal = benchmark.terminal(torch.as_tensor(x, dtype=torch.float64)).numpy()
    formula_terminal = 0.5 / (1.0 + 0.2 * s)
    terminal_impl_max = float(np.max(np.abs(torch_terminal - formula_terminal)))

    residual_rms = float(np.sqrt(np.mean(fine**2)))
    residual_max = float(np.max(np.abs(fine)))
    relative_rms = float(np.sqrt(np.mean(relative**2)))
    refinement_rms = float(np.sqrt(np.mean((fine - coarse) ** 2)))
    # These numerical-consistency tolerances were fixed in code before the
    # new reference was evaluated.  The direct truth-uncertainty gate below
    # is substantially tighter and remains the main accuracy safeguard.
    passed = bool(
        terminal_solver_max <= 1e-14
        and terminal_impl_max <= 1e-14
        and residual_rms <= 2e-3
        and residual_max <= 2e-2
        and relative_rms <= 2e-3
        and refinement_rms <= 2e-3
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "seed": seed,
        "states": 200,
        "radial_config": asdict(config),
        "fine_difference_steps": {"dt": 0.001, "ds": 0.2},
        "coarse_difference_steps": {"dt": 0.002, "ds": 0.4},
        "residual_rms": residual_rms,
        "residual_max_abs": residual_max,
        "relative_residual_rms": relative_rms,
        "coarse_fine_residual_rms": refinement_rms,
        "thresholds": {
            "residual_rms": 2e-3,
            "residual_max_abs": 2e-2,
            "relative_residual_rms": 2e-3,
            "coarse_fine_residual_rms": 2e-3,
            "terminal_max_abs": 1e-14,
        },
        "terminal_solver_max_abs": terminal_solver_max,
        "terminal_implementation_max_abs": terminal_impl_max,
    }


def prepare(*, force_verify: bool = False) -> dict[str, Any]:
    ensure_directories()
    manifest_path = RESULTS_DIR / "bank_manifest.json"
    audit_path = RESULTS_DIR / "reference_audit.json"
    if manifest_path.is_file() and audit_path.is_file() and not force_verify:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") != "VALIDATED":
            raise RuntimeError("existing reference gate is not validated")
        return {"manifest": manifest, "reference_audit": audit}

    if not (OLD_BANK / "x.npy").is_file():
        raise FileNotFoundError(f"missing frozen previous bank: {OLD_BANK}")
    old_x = np.load(OLD_BANK / "x.npy", allow_pickle=False)
    if old_x.shape != (32, DIMENSION):
        raise RuntimeError(f"unexpected previous bank shape: {old_x.shape}")

    states: dict[str, np.ndarray] = {}
    bank_seeds: dict[str, int] = {}
    for name, count in BANK_SIZES.items():
        states[name], bank_seeds[name] = _new_states(name, count)
        if len(_row_bytes(states[name])) != count:
            raise RuntimeError(f"duplicate states within {name} bank")
    separation = _assert_bank_separation(states, old_x)

    cal_t = np.zeros(BANK_SIZES["calibration"], dtype=np.float64)
    final_t = np.zeros(BANK_SIZES["final"], dtype=np.float64)
    combined_t = np.concatenate((cal_t, final_t))
    combined_x = np.concatenate((states["calibration"], states["final"]), axis=0)
    reference, uncertainty, source_audit = audited_allen_reference(
        DIMENSION, combined_t, combined_x
    )
    reference = np.asarray(reference, dtype=np.float64)
    uncertainty = np.asarray(uncertainty, dtype=np.float64)
    uncertainty_rms = float(np.sqrt(np.mean(uncertainty**2)))
    uncertainty_max = float(np.max(uncertainty))
    identities = _pde_identity_audit()
    finite = bool(np.all(np.isfinite(reference)) and np.all(np.isfinite(uncertainty)))
    gate_passed = bool(
        finite
        and uncertainty_rms <= 5e-6
        and uncertainty_max <= 1e-5
        and identities["passed"]
    )

    cal_count = BANK_SIZES["calibration"]
    cal_reference, final_reference = reference[:cal_count], reference[cal_count:]
    cal_uncertainty, final_uncertainty = uncertainty[:cal_count], uncertainty[cal_count:]
    if gate_passed:
        write_npz_atomic(
            bank_path("calibration"),
            t=cal_t,
            x=states["calibration"],
            reference=cal_reference,
            reference_uncertainty=cal_uncertainty,
        )
        write_npz_atomic(
            bank_path("final"),
            t=final_t,
            x=states["final"],
            reference=final_reference,
            reference_uncertainty=final_uncertainty,
        )
        # Deliberately no reference field in the timing artifact.
        write_npz_atomic(
            bank_path("timing"),
            t=np.zeros(BANK_SIZES["timing"], dtype=np.float64),
            x=states["timing"],
        )

    audit = {
        "revision": REVISION,
        "status": "VALIDATED" if gate_passed else "BLOCKED_REFERENCE_GATE",
        "overall_pass": gate_passed,
        "finite": finite,
        "uncertainty_rms": uncertainty_rms,
        "uncertainty_max": uncertainty_max,
        "thresholds": {"uncertainty_rms": 5e-6, "uncertainty_max": 1e-5},
        "source_reference_audit": source_audit,
        "pde_and_terminal_audit": identities,
        "combined_input_sha256": array_sha256(combined_t, combined_x),
        "combined_reference_sha256": array_sha256(reference),
        "combined_uncertainty_sha256": array_sha256(uncertainty),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(audit_path, audit)
    if not gate_passed:
        raise RuntimeError(
            "BLOCKED_REFERENCE_GATE: "
            f"uncertainty_rms={uncertainty_rms:.3e}, "
            f"uncertainty_max={uncertainty_max:.3e}, "
            f"identity={identities['status']}"
        )

    rows = {}
    for name in BANK_SIZES:
        path = bank_path(name)
        bank_arrays = {"t": np.zeros(BANK_SIZES[name]), "x": states[name]}
        rows[name] = {
            "size": BANK_SIZES[name],
            "seed": bank_seeds[name],
            "law": "t=0; x~N(0,I_100)" + (
                "; row zero replaced by exact origin" if name != "timing" else ""
            ),
            "path": str(path.relative_to(SOURCE_ROOT)),
            "file_sha256": file_sha256(path),
            "input_array_sha256": array_sha256(bank_arrays["t"], bank_arrays["x"]),
            "contains_reference": name != "timing",
        }
    manifest = {
        "revision": REVISION,
        "namespace_seed": NAMESPACE_SEED,
        "protocol": protocol_payload(),
        "banks": rows,
        "separation_audit": separation,
        "previous_bank": {
            "path": str(OLD_BANK.relative_to(SOURCE_ROOT)),
            "input_array_sha256": array_sha256(
                np.load(OLD_BANK / "t.npy", allow_pickle=False), old_x
            ),
        },
        "truth_policy": "calibration/final only; timing bank contains no truth",
    }
    write_json_atomic(manifest_path, manifest)
    return {"manifest": manifest, "reference_audit": audit}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-verify", action="store_true")
    args = parser.parse_args()
    result = prepare(force_verify=args.force_verify)
    print(json.dumps({
        "status": result["reference_audit"]["status"],
        "uncertainty_rms": result["reference_audit"]["uncertainty_rms"],
        "uncertainty_max": result["reference_audit"]["uncertainty_max"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

