"""Resumable runner for the confirmatory round-2 mechanism suite."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time
from typing import Any, Iterable

import numpy as np
import scipy

from .equations import make_points
from .mechanism_mlp import (
    BOX,
    CENTRE,
    F_ZERO,
    ORACLE_STATE,
    ORACLE_Z,
    RAW,
    SEGMENT,
    SHRINK_CENTRE_FACTORS,
    SHRINK_FACTORS,
    Z_ZERO,
    dose_method,
    illegal_method,
    load_repetition,
    method_from_dict,
    shrink_centre_method,
    shrink_method,
)
from .round2_certificate import (
    BOUND_CACHE_PATH,
    ROUND2_BASE_SEED,
    ensure_p4_bound_cache,
    run_p4_containment_audit,
)
from .round2_mlp import (
    ROUND2_IMPLEMENTATION_REVISION,
    ROUND2_PREREGISTRATION_COMMIT,
    TIGHT_CENTRE,
    TIGHT_SEGMENT,
    run_single_repetition_r2,
    save_round2_repetition,
    tight_illegal_method,
    tightness_method,
)
from .run_suite import P4_REFERENCE, chunk_size_for_dimension, make_equation


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parent
FULL_HISTORY_SOURCE = HERE.parent / "active_vb_high_budget" / "vb_mlp_methods.py"
RESULTS_ROOT = PACKAGE_ROOT / "results" / "mechanism_suite_r2"
RAW_ROOT = RESULTS_ROOT / "raw"
CONTAINMENT_PATH = RESULTS_ROOT / "r3_containment.json"
CONTAINMENT_INITIAL_AUDIT_PATH = (
    RESULTS_ROOT / "audit_history" / "r3_containment_initial_four_grid.json"
)
POINT_COUNT = 1200
ROUND2_CONFIGS = ((3, 6), (3, 10), (4, 3), (4, 6))
R1_DEEP_CONFIGS = ((4, 3), (4, 4), (4, 6), (5, 2))
R1_SHALLOW_CONFIGS = ((3, 6), (3, 10), (3, 16), (2, 32))
R1_CONFIGS = R1_DEEP_CONFIGS + R1_SHALLOW_CONFIGS
R3_CONFIGS = ((3, 6), (4, 6))
R4_CONFIGS = ((4, 3), (3, 6))
R4_DOSES = (0.5, 1.0, 2.0, 4.0)
R6_CONFIGS = ((4, 6),)
ILLEGAL_R2_FACTORS = (0.25, 0.5, 0.75)
TIGHTNESS_LEVELS = (0.0, 0.25, 0.5, 0.75, 1.0)


def _git(*args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args],
            cwd=REPOSITORY_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if math.isnan(number):
            return "NaN"
        if math.isinf(number):
            return "Infinity" if number > 0 else "-Infinity"
        return number
    return value


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(
            _json_safe(payload), indent=2, sort_keys=True, allow_nan=False
        )
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _spec(pde_id: str, d: int) -> dict[str, Any]:
    if pde_id == "P1":
        return {"pde_id": pde_id, "d": d, "kind": "p1"}
    if pde_id == "P2_a4":
        return {"pde_id": pde_id, "d": d, "kind": "p2", "a": 4.0}
    if pde_id == "P2_a8":
        return {"pde_id": pde_id, "d": d, "kind": "p2", "a": 8.0}
    if pde_id == "P3_rho1":
        return {"pde_id": pde_id, "d": d, "kind": "p3", "rho": 1.0}
    if pde_id == "P3_rho2":
        return {"pde_id": pde_id, "d": d, "kind": "p3", "rho": 2.0}
    if pde_id == "P4":
        return {"pde_id": pde_id, "d": d, "kind": "p4"}
    raise ValueError(f"unknown round-2 PDE {pde_id!r}")


def r1_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    specs.extend(_spec("P1", d) for d in (20, 50, 100))
    for pde in ("P2_a4", "P2_a8", "P3_rho1", "P3_rho2", "P4"):
        specs.extend(_spec(pde, d) for d in (20, 50))
    return specs


def r2_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for pde in (
        "P1",
        "P2_a4",
        "P2_a8",
        "P3_rho1",
        "P3_rho2",
        "P4",
    ):
        specs.extend(_spec(pde, d) for d in (20, 50))
    return specs


def primary_method(spec: dict[str, Any]) -> Any:
    return TIGHT_SEGMENT if spec["kind"] == "p4" else BOX


def r1_methods(spec: dict[str, Any]) -> list[Any]:
    if spec["kind"] == "p1":
        return [RAW, BOX, ORACLE_STATE, ORACLE_Z, SEGMENT]
    if spec["kind"] in {"p2", "p3"}:
        return [RAW, BOX, ORACLE_STATE, ORACLE_Z]
    return [
        RAW,
        TIGHT_SEGMENT,
        ORACLE_STATE,
        ORACLE_Z,
        SEGMENT,
        BOX,
    ]


def r2_methods(spec: dict[str, Any]) -> list[Any]:
    centre = TIGHT_CENTRE if spec["kind"] == "p4" else CENTRE
    methods = [RAW, primary_method(spec), centre]
    methods.extend(shrink_method(value) for value in SHRINK_FACTORS)
    methods.extend((Z_ZERO, F_ZERO))
    if spec["kind"] == "p4":
        methods.extend(tight_illegal_method(value) for value in ILLEGAL_R2_FACTORS)
        methods.extend((SEGMENT, BOX))
    else:
        methods.extend(illegal_method(value) for value in ILLEGAL_R2_FACTORS)
        if spec["kind"] == "p1":
            methods.append(SEGMENT)
    by_name: dict[str, Any] = {}
    for method in methods:
        previous = by_name.get(method.name)
        if previous is not None and previous.to_dict() != method.to_dict():
            raise RuntimeError(f"conflicting method definitions for {method.name}")
        by_name[method.name] = method
    return list(by_name.values())


def _safe_name(value: str) -> str:
    return value.replace("-", "m").replace(".", "p").replace("+", "")


def task_path(task: dict[str, Any]) -> Path:
    return (
        RAW_ROOT
        / str(task["stage"])
        / str(task["spec"]["pde_id"])
        / f"d{int(task['spec']['d'])}"
        / f"n{int(task['n'])}_M{int(task['M'])}"
        / _safe_name(str(task["method"]["name"]))
        / f"rep{int(task['repetition']):02d}.npz"
    )


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
    methods: Any,
    repetitions: int = 10,
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for spec in specs:
        for n, M in configs:
            for method in methods(spec):
                for repetition in range(repetitions):
                    tasks.append(
                        _task(stage, spec, n, M, method, repetition)
                    )
    return tasks


def r1_tasks() -> list[dict[str, Any]]:
    return _cross_tasks("r1", r1_specs(), R1_CONFIGS, r1_methods)


def r2_tasks() -> list[dict[str, Any]]:
    return _cross_tasks("r2", r2_specs(), ROUND2_CONFIGS, r2_methods)


def r3_tasks() -> list[dict[str, Any]]:
    specs = [_spec("P4", d) for d in (20, 50)]
    return _cross_tasks(
        "r3",
        specs,
        R3_CONFIGS,
        lambda spec: [tightness_method(theta) for theta in TIGHTNESS_LEVELS],
    )


def r4_tasks() -> list[dict[str, Any]]:
    specs = [_spec("P4", d) for d in (20, 50)]
    return _cross_tasks(
        "r4",
        specs,
        R4_CONFIGS,
        lambda spec: [dose_method(level) for level in R4_DOSES],
    )


def r5_tasks() -> list[dict[str, Any]]:
    specs = [
        spec
        for spec in r1_specs()
        if spec["kind"] in {"p2", "p3"}
    ]
    return _cross_tasks(
        "r5",
        specs,
        R1_CONFIGS,
        lambda spec: [
            shrink_centre_method(value) for value in SHRINK_CENTRE_FACTORS
        ],
    )


def r6_tasks() -> list[dict[str, Any]]:
    specs = [_spec("P1", d) for d in (20, 50, 100, 200, 400)]
    return _cross_tasks(
        "r6",
        specs,
        R6_CONFIGS,
        lambda spec: [RAW, BOX, SEGMENT, ORACLE_STATE, CENTRE],
    )


def _artifact_matches(path: Path, task: dict[str, Any]) -> bool:
    if not path.exists():
        return False
    try:
        loaded = load_repetition(path)
    except Exception as error:
        raise RuntimeError(f"cannot validate existing artifact {path}") from error
    metadata = loaded["metadata"]
    expected = (
        metadata.get("round") == 2
        and metadata.get("implementation_revision")
        == ROUND2_IMPLEMENTATION_REVISION
        and metadata.get("preregistration_commit")
        == ROUND2_PREREGISTRATION_COMMIT
        and metadata["pde_id"] == task["spec"]["pde_id"]
        and int(metadata["dimension"]) == int(task["spec"]["d"])
        and int(metadata["n"]) == int(task["n"])
        and int(metadata["M"]) == int(task["M"])
        and int(metadata["repetition"]) == int(task["repetition"])
        and metadata["method"]["name"] == task["method"]["name"]
        and metadata["method"]["transform"]
        == task["method"]["transform"]
        and float(metadata["method"]["factor"])
        == float(task["method"]["factor"])
        and int(metadata["chunk_size"]) == int(task["chunk_size"])
        and int(metadata["base_seed"]) == ROUND2_BASE_SEED
        and metadata["dtype"] == "float64"
    )
    if not expected:
        raise RuntimeError(
            f"existing round-2 artifact does not match its task: {path}"
        )
    if (
        loaded["prediction_u"].shape != (POINT_COUNT,)
        or loaded["truth_u"].shape != (POINT_COUNT,)
        or loaded["is_validation"].shape != (POINT_COUNT,)
        or loaded["prediction_u"].dtype != np.float64
        or loaded["truth_u"].dtype != np.float64
        or loaded["is_validation"].dtype != np.bool_
    ):
        raise RuntimeError(f"artifact has wrong shape or dtype: {path}")
    return True


def _execute_task(task: dict[str, Any]) -> dict[str, Any]:
    path = task_path(task)
    if _artifact_matches(path, task):
        return {"status": "skipped", "path": str(path)}
    if path.exists():
        raise RuntimeError(
            f"refusing to overwrite a non-matching round-2 artifact: {path}"
        )
    equation = make_equation(task["spec"])
    points = make_points(
        equation, n_points=POINT_COUNT, seed=ROUND2_BASE_SEED
    )
    result = run_single_repetition_r2(
        pde_id=str(task["spec"]["pde_id"]),
        equation=equation,
        method=method_from_dict(task["method"]),
        n=int(task["n"]),
        M=int(task["M"]),
        repetition=int(task["repetition"]),
        t=points["t"],
        x=points["x"],
        is_validation=points["is_validation"],
        base_seed=ROUND2_BASE_SEED,
        chunk_size=int(task["chunk_size"]),
        study=str(task["stage"]),
    )
    result["metadata"]["point_set"] = {
        "seed": ROUND2_BASE_SEED,
        "count": POINT_COUNT,
        "validation_fraction": 0.2,
    }
    result["metadata"]["code_commit_at_execution"] = _git("rev-parse", "HEAD")
    save_round2_repetition(path, result)
    return {
        "status": "completed",
        "path": str(path),
        "seconds": result["metadata"]["wall_clock_seconds"],
        "skill": result["metadata"]["metrics"]["test"]["skill"],
    }


def _source_candidates(task: dict[str, Any]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    stage = str(task["stage"])
    pde = str(task["spec"]["pde_id"])
    d = int(task["spec"]["d"])
    config = (int(task["n"]), int(task["M"]))
    name = str(task["method"]["name"])

    def source(source_stage: str) -> dict[str, Any]:
        item = dict(task)
        item["stage"] = source_stage
        item["chunk_size"] = chunk_size_for_dimension(d, source_stage)
        return item

    if (
        pde == "P4"
        and d in {20, 50}
        and config in R3_CONFIGS
        and name in {"segment", "tight_segment"}
        and stage in {"r1", "r2"}
    ):
        candidates.append(source("r3"))
    if stage == "r1" and d in {20, 50} and config in ROUND2_CONFIGS:
        overlap = {
            "P1": {"raw", "box", "segment"},
            "P2_a4": {"raw", "box"},
            "P2_a8": {"raw", "box"},
            "P3_rho1": {"raw", "box"},
            "P3_rho2": {"raw", "box"},
            "P4": {"raw", "box", "segment", "tight_segment"},
        }
        if name in overlap[pde]:
            candidates.append(source("r2"))
    if (
        stage == "r6"
        and d in {20, 50, 100}
        and name in {"raw", "box", "segment", "oracle_state"}
    ):
        candidates.append(source("r1"))
    return candidates


def materialize_exact_reuses(tasks: list[dict[str, Any]]) -> int:
    reused = 0
    for task in tasks:
        destination = task_path(task)
        if destination.exists():
            continue
        for source_task in _source_candidates(task):
            if int(source_task["chunk_size"]) != int(task["chunk_size"]):
                continue
            source_path = task_path(source_task)
            if not source_path.exists() or not _artifact_matches(
                source_path, source_task
            ):
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(source_path, destination)
            except OSError:
                shutil.copy2(source_path, destination)
            reused += 1
            break
    return reused


def run_tasks(
    tasks: list[dict[str, Any]], *, workers: int, dry_run: bool
) -> None:
    reused = 0 if dry_run else materialize_exact_reuses(tasks)
    missing: list[dict[str, Any]] = []
    for task in tasks:
        if not _artifact_matches(task_path(task), task):
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
        futures = {
            executor.submit(_execute_task, task): task for task in missing
        }
        for future in as_completed(futures):
            task = futures[future]
            result = future.result()
            completed += 1
            if (
                completed == 1
                or completed % 10 == 0
                or completed == len(missing)
            ):
                elapsed = time.perf_counter() - started
                rate = completed / max(elapsed, 1e-9)
                remaining = (len(missing) - completed) / max(rate, 1e-9)
                print(
                    f"[{completed}/{len(missing)}] {task['stage']} "
                    f"{task['spec']['pde_id']} d={task['spec']['d']} "
                    f"n={task['n']} M={task['M']} "
                    f"{task['method']['name']} rep={task['repetition']} "
                    f"eta={remaining/60:.1f}m",
                    flush=True,
                )
    summed_seconds = sum(
        float(future.result().get("seconds", 0.0)) for future in futures
    )
    print(f"worker-summed seconds={summed_seconds:.1f}")


def _containment_is_current(payload: dict[str, Any]) -> bool:
    return (
        int(payload.get("base_seed", -1)) == ROUND2_BASE_SEED
        and payload.get("implementation_revision")
        == ROUND2_IMPLEMENTATION_REVISION
        and payload.get("preregistration_commit")
        == ROUND2_PREREGISTRATION_COMMIT
    )


def run_r3_gate(*, dry_run: bool) -> bool:
    if CONTAINMENT_PATH.exists():
        payload = json.loads(CONTAINMENT_PATH.read_text(encoding="utf-8"))
        if not _containment_is_current(payload):
            raise RuntimeError(
                "existing R3 containment artifact belongs to another revision"
            )
        print(
            f"R3 containment existing proceed={payload['proceed']}",
            flush=True,
        )
        return bool(payload["proceed"])
    if dry_run:
        print("R3 containment missing", flush=True)
        return False
    if not P4_REFERENCE.exists():
        raise RuntimeError(
            "round-1 P4 reference is missing; refusing to modify round-1 results"
        )
    print(f"building/auditing P4 bound cache: {BOUND_CACHE_PATH}", flush=True)
    cache_metadata = ensure_p4_bound_cache()
    print("running R3 four-grid sharp-bound containment audit", flush=True)
    payload = run_p4_containment_audit(P4_REFERENCE)
    payload.update(
        {
            "implementation_revision": ROUND2_IMPLEMENTATION_REVISION,
            "preregistration_commit": ROUND2_PREREGISTRATION_COMMIT,
            "code_commit": _git("rev-parse", "HEAD"),
            "git_branch": _git("branch", "--show-current"),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_count": os.cpu_count(),
            "bound_cache": cache_metadata,
            "prior_failed_audit": (
                CONTAINMENT_INITIAL_AUDIT_PATH.relative_to(RESULTS_ROOT).as_posix()
                if CONTAINMENT_INITIAL_AUDIT_PATH.exists()
                else None
            ),
        }
    )
    _write_json_atomic(CONTAINMENT_PATH, payload)
    print(
        f"R3 containment proceed={payload['proceed']} "
        f"finest={payload['successive_grid_results'][-1]['maximum_violation']:.3e}",
        flush=True,
    )
    return bool(payload["proceed"])


def _p4_only(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [task for task in tasks if task["spec"]["pde_id"] == "P4"]


def run_stage(stage: str, *, workers: int, dry_run: bool) -> None:
    print(f"=== {stage.upper()} ===", flush=True)
    if stage == "r3":
        proceed = run_r3_gate(dry_run=dry_run)
        if not proceed and not dry_run:
            raise RuntimeError(
                "R3 G5 refinement rule failed; stopping before P4 MLP runs"
            )
        run_tasks(r3_tasks(), workers=workers, dry_run=dry_run)
        print("R3 prerequisite: P4 R2 cells", flush=True)
        run_tasks(_p4_only(r2_tasks()), workers=workers, dry_run=dry_run)
        print("R3 prerequisite: P4 R1 cells", flush=True)
        run_tasks(_p4_only(r1_tasks()), workers=workers, dry_run=dry_run)
        return
    mapping = {
        "r2": r2_tasks,
        "r1": r1_tasks,
        "r4": r4_tasks,
        "r6": r6_tasks,
        "r5": r5_tasks,
    }
    if stage == "r7":
        print("R7 is derived from R1; auditing R1 completeness", flush=True)
        run_tasks(r1_tasks(), workers=workers, dry_run=True)
        return
    if stage not in mapping:
        raise ValueError(stage)
    if stage in {"r2", "r1", "r4"}:
        gate_passed = run_r3_gate(dry_run=True)
        if not gate_passed and not dry_run:
            raise RuntimeError("R3 containment must pass before later stages")
        if not gate_passed and dry_run:
            print(
                "R3 containment is not yet available; continuing task-count "
                "audit only",
                flush=True,
            )
    run_tasks(mapping[stage](), workers=workers, dry_run=dry_run)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("r3", "r2", "r1", "r4", "r6", "r5", "r7", "all"),
        default="all",
    )
    parser.add_argument(
        "--workers", type=int, default=min(8, os.cpu_count() or 1)
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.workers < 1:
        raise ValueError("workers must be positive")
    stages = (
        ("r3", "r2", "r1", "r4", "r6", "r5", "r7")
        if arguments.stage == "all"
        else (arguments.stage,)
    )
    for stage in stages:
        run_stage(stage, workers=arguments.workers, dry_run=arguments.dry_run)


if __name__ == "__main__":
    main()
