"""Resumable runner for the nonlinear-generator mechanism benchmark suite.

Every repetition is an independent, atomic NPZ artifact.  Existing artifacts
are validated and skipped, so a long run can be resumed without overwriting
completed work.  The E0 gate file controls admission to E1--E4.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable

import numpy as np
import scipy

from .equations import (
    BASE_SEED,
    BurgersFisher,
    MultiDirectionLSEHJB,
    RidgeLSEHJB,
    VBa,
    make_points,
    published_vb,
)
from .gates import run_analytic_gates
from .mechanism_mlp import (
    BALL,
    BATCH_BOX,
    BOX,
    CENTRE,
    F_ZERO,
    ORACLE_STATE,
    ORACLE_U,
    ORACLE_Z,
    RAW,
    SEGMENT,
    SIGN_ONLY,
    SPAN_ONLY,
    U_DOSE_SCALE,
    Z_ZERO,
    dose_method,
    dose_u_method,
    load_repetition,
    method_from_dict,
    run_single_repetition,
    save_repetition,
    tuning_methods,
)
from .norm_hjb import NormDriverHJB, build_norm_reference, load_norm_reference


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
FULL_HISTORY_SOURCE = HERE.parent / "active_vb_high_budget" / "vb_mlp_methods.py"
RESULTS_ROOT = PACKAGE_ROOT / "results" / "mechanism_suite"
RAW_ROOT = RESULTS_ROOT / "raw"
GATES_PATH = RESULTS_ROOT / "gates.json"
REFERENCE_ROOT = RESULTS_ROOT / "reference_cache"
P4_REFERENCE = REFERENCE_ROOT / "norm_hjb_beta2_lambda1_T0p5.npz"
P4_REJECTED_REFERENCE = REFERENCE_ROOT / "norm_hjb_beta2_lambda1_T0p25.npz"
POINT_COUNT = 1200
P4_IMPLEMENTATION_REVISION = "p4_cached_reference_derivative_v1"
HEADLINE_CONFIGS = ((3, 6), (4, 6))
FULL_CONFIGS = ((3, 6), (3, 10), (4, 3), (4, 6))
DOSE_CONFIGS = ((4, 3), (3, 6))
DOSE_LEVELS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)


def _git(command: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *command.split()],
            cwd=PACKAGE_ROOT.parent,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _json_safe(value: Any) -> Any:
    """Convert NumPy and non-finite values to strict-JSON representations."""

    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if math.isnan(number):
            return "NaN"
        if math.isinf(number):
            return "Infinity" if number > 0 else "-Infinity"
        return number
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def ensure_p4_references() -> None:
    """Create only missing pre-registered P4 reference artifacts."""

    REFERENCE_ROOT.mkdir(parents=True, exist_ok=True)
    if not P4_REFERENCE.exists():
        print(f"building P4 reference: {P4_REFERENCE}", flush=True)
        build_norm_reference(P4_REFERENCE, beta=2.0, lambda_f=1.0, horizon=0.5)
    if not P4_REJECTED_REFERENCE.exists():
        print(f"building P4 screen reference: {P4_REJECTED_REFERENCE}", flush=True)
        build_norm_reference(
            P4_REJECTED_REFERENCE, beta=2.0, lambda_f=1.0, horizon=0.25
        )
    for path, label in (
        (P4_REFERENCE, "T=0.5"),
        (P4_REJECTED_REFERENCE, "T=0.25"),
    ):
        metadata = load_norm_reference(str(path.resolve())).metadata
        if not metadata.get("refinement_passed", False):
            raise RuntimeError(f"P4 {label} reference did not pass its refinement gate")


def equation_specs() -> list[dict[str, Any]]:
    """All formal E0 candidates at their pre-registered dimensions."""

    specs: list[dict[str, Any]] = []
    specs.extend({"pde_id": "P1", "d": d, "kind": "p1"} for d in (20, 50, 100, 200))
    specs.extend(
        {"pde_id": "P2_a4", "d": d, "kind": "p2", "a": 4.0}
        for d in (20, 50, 100)
    )
    specs.extend(
        {"pde_id": "P2_a8", "d": d, "kind": "p2", "a": 8.0}
        for d in (20, 50, 100)
    )
    specs.extend(
        {"pde_id": "P3_rho1", "d": d, "kind": "p3", "rho": 1.0}
        for d in (20, 50)
    )
    specs.extend(
        {"pde_id": "P3_rho2", "d": d, "kind": "p3", "rho": 2.0}
        for d in (20, 50)
    )
    # P4 did not list dimensions.  Use the common 20/50/100 grid, with the
    # two smallest dimensions used in the compute-limited E3 control panel.
    specs.extend({"pde_id": "P4", "d": d, "kind": "p4"} for d in (20, 50, 100))
    specs.append({"pde_id": "N3", "d": 100, "kind": "n3"})
    specs.extend({"pde_id": "N4", "d": d, "kind": "n4"} for d in (20, 50, 100))
    return specs


def make_equation(spec: dict[str, Any]) -> Any:
    kind = spec["kind"]
    d = int(spec["d"])
    if kind == "p1":
        return RidgeLSEHJB(d=d)
    if kind == "p2":
        a = float(spec["a"])
        return VBa(d=d, a=a, name=f"P2_vba_a{a:g}")
    if kind == "p3":
        rho = float(spec["rho"])
        return BurgersFisher(d=d, rho=rho, name=f"P3_burgers_fisher_rho{rho:g}")
    if kind == "p4":
        return NormDriverHJB(d=d, reference_path=str(P4_REFERENCE.resolve()), T=0.5)
    if kind == "n3":
        return MultiDirectionLSEHJB(d=d, n_directions=5, strength=0.5)
    if kind == "n4":
        return published_vb(d)
    raise ValueError(f"unknown equation kind {kind!r}")


def chunk_size_for_dimension(d: int, stage: str = "e0") -> int:
    # E0 was launched with the conservative smoke-screen setting.  Once its
    # memory audit passed, the full studies use 16 at d=20 to amortize the
    # fixed cost of numerical-reference spline calls.  Pairing is preserved
    # because every method within a study receives the same stage-specific
    # chunk size.  High-dimensional cells retain the pre-registered guide of 4.
    if stage == "e2" and 50 <= d < 200:
        return 8
    if d >= 50:
        return 4
    return 8 if stage == "e0" else 16


def _safe_name(value: str) -> str:
    return value.replace("-", "m").replace(".", "p").replace("+", "")


def task_path(task: dict[str, Any]) -> Path:
    method = _safe_name(str(task["method"]["name"]))
    return (
        RAW_ROOT
        / str(task["stage"])
        / str(task["spec"]["pde_id"])
        / f"d{int(task['spec']['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / method
        / f"rep{int(task['repetition']):02d}.npz"
    )


def _artifact_matches(path: Path, task: dict[str, Any]) -> bool:
    if not path.exists():
        return False
    loaded = load_repetition(path)
    metadata = loaded["metadata"]
    expected = (
        metadata["pde_id"] == task["spec"]["pde_id"]
        and int(metadata["dimension"]) == int(task["spec"]["d"])
        and int(metadata["n"]) == int(task["n"])
        and int(metadata["M"]) == int(task["M"])
        and int(metadata["repetition"]) == int(task["repetition"])
        and metadata["method"]["name"] == task["method"]["name"]
        and metadata["method"]["transform"] == task["method"]["transform"]
        and float(metadata["method"]["factor"]) == float(task["method"]["factor"])
        and int(metadata["chunk_size"]) == int(task["chunk_size"])
        and int(metadata["base_seed"]) == BASE_SEED
        and metadata["dtype"] == "float64"
    )
    if not expected:
        raise RuntimeError(f"existing artifact does not match requested task: {path}")
    if (
        loaded["prediction_u"].shape != (POINT_COUNT,)
        or loaded["truth_u"].shape != (POINT_COUNT,)
        or loaded["is_validation"].shape != (POINT_COUNT,)
        or loaded["prediction_u"].dtype != np.float64
        or loaded["truth_u"].dtype != np.float64
        or loaded["is_validation"].dtype != np.bool_
    ):
        raise RuntimeError(f"existing artifact has wrong shape or dtype: {path}")
    if (
        task["spec"]["kind"] == "p4"
        and metadata.get("implementation_revision") != P4_IMPLEMENTATION_REVISION
    ):
        return False
    return True


def _execute_task(task: dict[str, Any]) -> dict[str, Any]:
    """Worker entry point; arguments and return value are spawn-pickle safe."""

    path = task_path(task)
    if _artifact_matches(path, task):
        return {"status": "skipped", "path": str(path)}
    if path.exists():
        backup = RAW_ROOT / "superseded" / path.relative_to(RAW_ROOT)
        backup.parent.mkdir(parents=True, exist_ok=True)
        if not backup.exists():
            shutil.copy2(path, backup)
    equation = make_equation(task["spec"])
    points = make_points(equation, n_points=POINT_COUNT, seed=BASE_SEED)
    result = run_single_repetition(
        pde_id=str(task["spec"]["pde_id"]),
        equation=equation,
        method=method_from_dict(task["method"]),
        n=int(task["n"]),
        M=int(task["M"]),
        repetition=int(task["repetition"]),
        t=points["t"],
        x=points["x"],
        is_validation=points["is_validation"],
        base_seed=BASE_SEED,
        chunk_size=int(task["chunk_size"]),
    )
    save_repetition(path, result)
    return {
        "status": "completed",
        "path": str(path),
        "seconds": result["metadata"]["wall_clock_seconds"],
        "skill": result["metadata"]["metrics"]["test"]["skill"],
    }


def reuse_source_task(task: dict[str, Any]) -> dict[str, Any] | None:
    """Return an exactly identical earlier computation, when one exists."""

    stage = task["stage"]
    method = task["method"]["name"]
    source_stage: str | None = None
    if stage == "e3" and method in {
        "raw",
        "box",
        "segment",
        "oracle_state",
        "oracle_z",
        "oracle_u",
    }:
        source_stage = "e1"
    elif stage == "e6" and method in {"raw", "box", "segment", "oracle_state"}:
        source_stage = "e1"
    elif stage == "e6" and method == "centre":
        dimension = int(task["spec"]["d"])
        if dimension in {20, 50}:
            source_stage = "e3"
        elif int(task["repetition"]) < 3:
            # E0 already evaluated this exact (4,6) centre task at d>=50,
            # where its conservative chunk size is the same as E6's.
            source_stage = "e0"
    elif (
        stage == "e5"
        and int(task["repetition"]) < 3
        and int(task["spec"]["d"]) >= 50
        and method in {"raw", "box", "oracle_state", "centre"}
    ):
        source_stage = "e0"
    if source_stage is None:
        return None
    source = dict(task)
    source["stage"] = source_stage
    source["chunk_size"] = chunk_size_for_dimension(
        int(source["spec"]["d"]), source_stage
    )
    if int(source["chunk_size"]) != int(task["chunk_size"]):
        return None
    return source


def materialize_exact_reuses(tasks: list[dict[str, Any]]) -> int:
    """Hard-link identical earlier artifacts, falling back to a byte copy."""

    reused = 0
    for task in tasks:
        destination = task_path(task)
        if destination.exists():
            continue
        source_task = reuse_source_task(task)
        if source_task is None:
            continue
        source = task_path(source_task)
        if not source.exists() or not _artifact_matches(source, source_task):
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(source, destination)
        except OSError:
            shutil.copy2(source, destination)
        reused += 1
    return reused


def _task(
    stage: str,
    spec: dict[str, Any],
    n: int,
    M: int,
    method: Any,
    repetition: int,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "spec": dict(spec),
        "n": int(n),
        "M": int(M),
        "method": method.to_dict(),
        "repetition": int(repetition),
        "chunk_size": chunk_size_for_dimension(int(spec["d"]), stage),
    }


def _cross_tasks(
    stage: str,
    specs: Iterable[dict[str, Any]],
    configs: Iterable[tuple[int, int]],
    method_fn: Any,
    repetitions: int,
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for spec in specs:
        for n, M in configs:
            for method in method_fn(spec):
                for repetition in range(repetitions):
                    tasks.append(_task(stage, spec, n, M, method, repetition))
    return tasks


def e0_tasks() -> list[dict[str, Any]]:
    return _cross_tasks(
        "e0",
        equation_specs(),
        HEADLINE_CONFIGS,
        lambda spec: [RAW, BOX, ORACLE_STATE, CENTRE],
        3,
    )


def _load_gates() -> dict[str, Any]:
    if not GATES_PATH.exists():
        raise RuntimeError("gates.json is missing; run --stage e0 first")
    return json.loads(GATES_PATH.read_text(encoding="utf-8"))


def admitted_specs() -> list[dict[str, Any]]:
    gates = _load_gates()
    admitted = {
        key for key, value in gates["candidate_summary"].items() if value["admitted"]
    }
    return [spec for spec in equation_specs() if spec["pde_id"] in admitted]


def e1_methods(spec: dict[str, Any]) -> list[Any]:
    common = [RAW, BOX]
    if spec["kind"] in {"p1", "p4"}:
        common.append(SEGMENT)
    return common + [ORACLE_STATE, ORACLE_Z, ORACLE_U]


def e1_tasks() -> list[dict[str, Any]]:
    return _cross_tasks("e1", admitted_specs(), FULL_CONFIGS, e1_methods, 10)


def e2_methods(spec: dict[str, Any]) -> list[Any]:
    methods: list[Any] = []
    for level in DOSE_LEVELS:
        methods.extend([dose_method(level), dose_method(level, projected="box")])
        if spec["kind"] in {"p1", "p4"}:
            methods.append(dose_method(level, projected="segment"))
        if spec["kind"] == "p3":
            methods.extend(
                [dose_u_method(level), dose_u_method(level, clipped=True)]
            )
    return methods


def e2_tasks() -> list[dict[str, Any]]:
    return _cross_tasks("e2", admitted_specs(), DOSE_CONFIGS, e2_methods, 10)


def _two_smallest_dimensions(specs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_pde: dict[str, list[dict[str, Any]]] = {}
    for spec in specs:
        by_pde.setdefault(str(spec["pde_id"]), []).append(spec)
    selected: list[dict[str, Any]] = []
    for group in by_pde.values():
        selected.extend(sorted(group, key=lambda item: int(item["d"]))[:2])
    return selected


def e3_methods(spec: dict[str, Any]) -> list[Any]:
    methods: list[Any] = [RAW, BOX]
    if spec["kind"] in {"p1", "p4"}:
        methods.extend([SEGMENT, BALL, SPAN_ONLY, BATCH_BOX])
    else:
        methods.extend([SIGN_ONLY, BALL, BATCH_BOX])
    methods.extend([ORACLE_STATE, ORACLE_Z, ORACLE_U, Z_ZERO, F_ZERO, CENTRE])
    methods.extend(tuning_methods())
    return methods


def e3_tasks() -> list[dict[str, Any]]:
    specs = _two_smallest_dimensions(admitted_specs())
    return _cross_tasks("e3", specs, FULL_CONFIGS, e3_methods, 10)


def e5_tasks() -> list[dict[str, Any]]:
    specs = [spec for spec in equation_specs() if spec["pde_id"] in {"N3", "N4"}]
    return _cross_tasks(
        "e5", specs, HEADLINE_CONFIGS, lambda spec: [RAW, BOX, ORACLE_STATE, CENTRE], 10
    )


def e6_tasks() -> list[dict[str, Any]]:
    specs = [spec for spec in equation_specs() if spec["pde_id"] == "P1"]
    return _cross_tasks(
        "e6",
        specs,
        ((4, 6),),
        lambda spec: [RAW, BOX, SEGMENT, ORACLE_STATE, CENTRE],
        10,
    )


def run_tasks(tasks: list[dict[str, Any]], *, workers: int, dry_run: bool) -> None:
    reused = 0 if dry_run else materialize_exact_reuses(tasks)
    missing: list[dict[str, Any]] = []
    for task in tasks:
        path = task_path(task)
        if not _artifact_matches(path, task):
            missing.append(task)
    print(
        f"tasks total={len(tasks)} complete={len(tasks)-len(missing)} "
        f"reused_now={reused} missing={len(missing)}",
        flush=True,
    )
    if dry_run or not missing:
        return
    started = time.perf_counter()
    completed = 0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_execute_task, task): task for task in missing}
        for future in as_completed(futures):
            task = futures[future]
            result = future.result()
            completed += 1
            if completed == 1 or completed % 10 == 0 or completed == len(missing):
                elapsed = time.perf_counter() - started
                rate = completed / max(elapsed, 1e-9)
                remaining = (len(missing) - completed) / max(rate, 1e-9)
                print(
                    f"[{completed}/{len(missing)}] {task['stage']} "
                    f"{task['spec']['pde_id']} d={task['spec']['d']} "
                    f"n={task['n']} M={task['M']} {task['method']['name']} "
                    f"rep={task['repetition']} eta={remaining/60:.1f}m",
                    flush=True,
                )


def _metric(stage: str, spec: dict[str, Any], n: int, M: int, method: str, rep: int) -> float:
    task = {
        "stage": stage,
        "spec": spec,
        "n": n,
        "M": M,
        "method": {"name": method, "transform": "unused", "factor": 1.0},
        "repetition": rep,
        "chunk_size": chunk_size_for_dimension(int(spec["d"]), stage),
    }
    path = task_path(task)
    return float(load_repetition(path)["metadata"]["metrics"]["test"]["skill"])


def analytic_gates() -> dict[str, dict[str, Any]]:
    ensure_p4_references()
    results: dict[str, dict[str, Any]] = {}
    for index, spec in enumerate(equation_specs(), start=1):
        equation = make_equation(spec)
        print(
            f"analytic gate [{index}/{len(equation_specs())}] {spec['pde_id']} d={spec['d']}",
            flush=True,
        )
        results.setdefault(str(spec["pde_id"]), {})[str(spec["d"])] = run_analytic_gates(
            equation
        )
    return results


def _p4_screen_record() -> dict[str, Any]:
    record: dict[str, Any] = {
        "selected": {"beta": 2.0, "lambda_f": 1.0, "T": 0.5},
        "selection_rule": "beta=2 and the stronger suggested lambda_f=1 were fixed before MLP inspection; T=0.25 and T=0.5 were compared by E0, retaining the setting that passed G2",
        "not_evaluated": {
            "parameters": {"lambda_f": 0.5, "T": [0.25, 0.5]},
            "reason": "the protocol described these as suggested values rather than requiring an exhaustive parameter sweep; no MLP result was used to choose lambda_f=1",
        },
    }
    if P4_REJECTED_REFERENCE.exists():
        equation = NormDriverHJB(
            d=20,
            reference_path=str(P4_REJECTED_REFERENCE.resolve()),
            T=0.25,
        )
        rejected = run_analytic_gates(equation)
        record["screened_out"] = {
            "parameters": {"beta": 2.0, "lambda_f": 1.0, "T": 0.25},
            "reason": "G2 below 0.15",
            "gates": rejected,
        }
    return record


def finalize_gates(analytic: dict[str, dict[str, Any]]) -> dict[str, Any]:
    expected_failure = {"N3": "G2", "N4": "G3"}
    summaries: dict[str, Any] = {}
    for pde_id, dimensions in analytic.items():
        instance_summaries: dict[str, Any] = {}
        for dimension_text, analytic_result in dimensions.items():
            d = int(dimension_text)
            spec = next(
                item
                for item in equation_specs()
                if item["pde_id"] == pde_id and int(item["d"]) == d
            )
            cells: list[dict[str, Any]] = []
            for n, M in HEADLINE_CONFIGS:
                values: dict[str, list[float]] = {}
                for method in ("raw", "box", "oracle_state", "centre"):
                    values[method] = [
                        _metric("e0", spec, n, M, method, repetition)
                        for repetition in range(3)
                    ]
                means = {key: float(np.mean(value)) for key, value in values.items()}
                medians = {key: float(np.median(value)) for key, value in values.items()}
                ratio = means["raw"] / means["oracle_state"]
                cells.append(
                    {
                        "n": n,
                        "M": M,
                        "mean_test_skill": means,
                        "median_test_skill": medians,
                        "raw_over_oracle": ratio,
                        "G4_cell_passed": bool(ratio >= 2.0),
                    }
                )
            g3 = any(
                cell["mean_test_skill"]["oracle_state"] <= 0.10 for cell in cells
            )
            g4 = all(cell["G4_cell_passed"] for cell in cells)
            instance_summaries[dimension_text] = {
                "analytic": analytic_result,
                "G3": {
                    "passed": g3,
                    "threshold": "mean oracle_state test skill <=0.10 in at least one headline configuration",
                },
                "G4": {
                    "passed": g4,
                    "threshold": "mean raw/oracle_state skill >=2 in both headline cells",
                },
                "headline_cells": cells,
            }
        gate_status = {
            gate: all(
                bool(instance["analytic"][gate]["passed"])
                if gate in {"G1", "G2", "G5", "G6"}
                else bool(instance[gate]["passed"])
                for instance in instance_summaries.values()
            )
            for gate in ("G1", "G2", "G3", "G4", "G5", "G6")
        }
        failures = [gate for gate, passed in gate_status.items() if not passed]
        predicted = expected_failure.get(pde_id)
        negative_eligible = len(failures) == 1 and failures[0] == predicted
        summaries[pde_id] = {
            "gate_status": gate_status,
            "failed_gates": failures,
            "predicted_failure": predicted,
            "negative_control_eligible": negative_eligible,
            "admitted": pde_id.startswith("P") and not failures,
            "instances": instance_summaries,
        }
    payload = {
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "protocol": {
            "point_count": POINT_COUNT,
            "test_distributions": {
                "P1_N3_P4": "t uniform on [0,T), x uniform on [-1,1]^d",
                "P2_P3_N4": "inherited VB make_test_points: 1000 interior + 200 face-boundary points at n=1200",
            },
            "analytic_G2_point_count": 2000,
            "containment_point_count": 100000,
            "headline_configurations": [list(value) for value in HEADLINE_CONFIGS],
            "e0_repetitions": 3,
            "base_seed": BASE_SEED,
            "dtype": "float64",
            "p3_u_dose_scale": U_DOSE_SCALE,
            "p4_dimensions": {
                "E0_E1_E2": [20, 50, 100],
                "E3": [20, 50],
                "note": "the supplied P4 specification did not list dimensions",
            },
            "chunk_sizes": {
                "E0": {"d20": 8, "d_ge_50": 4},
                "E1_E3_E5_E6": {"d20": 16, "d_ge_50": 4},
                "E2": {"d20": 16, "d50_d100": 8, "d200": 4},
            },
            "G4_interpretation": "both pre-registered headline cells must pass",
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_count": os.cpu_count(),
            "git_branch": _git("branch --show-current"),
            "git_commit": _git("rev-parse HEAD"),
        },
        "unchanged_recursion_source": {
            "path": "invariant_region_mlp/experiments/active_vb_high_budget/vb_mlp_methods.py",
            "sha256": hashlib.sha256(FULL_HISTORY_SOURCE.read_bytes()).hexdigest(),
            "class": "FullHistoryMLP",
        },
        "reference_bundle": {
            "reference_code_present": False,
            "note": "The supplied attachment contained only Pasted text.txt; no reference_code directory was present. N3 was reconstructed from the mathematical specification.",
        },
        "p4_parameter_screen": _p4_screen_record(),
        "candidate_summary": summaries,
        "historical_negative_controls": {
            "N1": {
                "problem": "Rosenbrock HJB",
                "failed_gate": "G2",
                "prediction": "f=0 is the best method",
                "rerun": False,
                "evidence": [
                    "invariant_region_mlp/docs/HJB_LIFE_OR_DEATH_ABLATION.md",
                    "invariant_region_mlp/results/hjb_life_or_death_summary.json",
                ],
            },
            "N2": {
                "problem": "100D funding",
                "failed_gate": "G2",
                "prediction": "z=0 beats all corrections",
                "rerun": False,
                "evidence": [
                    "invariant_region_mlp/docs/FUNDING_LIFE_OR_DEATH_ABLATION.md",
                    "invariant_region_mlp/results/funding_life_or_death_summary.json",
                ],
            },
            "N5": {
                "problem": "linear driver f=0",
                "failed_gate": "G4",
                "prediction": "all corrections identical in value",
                "rerun": False,
                "evidence": [
                    "invariant_region_mlp/docs/batchir_negative_controls.md"
                ],
            },
        },
    }
    _write_json_atomic(GATES_PATH, payload)
    return payload


def run_e0(*, workers: int, dry_run: bool) -> None:
    ensure_p4_references()
    run_tasks(e0_tasks(), workers=workers, dry_run=dry_run)
    if dry_run:
        return
    analytic = analytic_gates()
    payload = finalize_gates(analytic)
    print("E0 gate verdicts:", flush=True)
    for pde_id, summary in payload["candidate_summary"].items():
        verdict = "ADMIT" if summary["admitted"] else "EXCLUDE"
        print(f"  {pde_id}: {verdict}; failures={summary['failed_gates']}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("e0", "e1", "e2", "e3", "e5", "e6", "all"),
        default="e0",
    )
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    stages = ("e0", "e1", "e2", "e3", "e5", "e6") if args.stage == "all" else (args.stage,)
    for stage in stages:
        print(f"=== {stage.upper()} ===", flush=True)
        if stage == "e0":
            run_e0(workers=args.workers, dry_run=args.dry_run)
        else:
            builder = {
                "e1": e1_tasks,
                "e2": e2_tasks,
                "e3": e3_tasks,
                "e5": e5_tasks,
                "e6": e6_tasks,
            }[stage]
            run_tasks(builder(), workers=args.workers, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
