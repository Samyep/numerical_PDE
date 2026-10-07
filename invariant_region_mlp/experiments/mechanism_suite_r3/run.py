"""Resumable runner for the pre-registered round-3 mechanism suite."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any, Iterable

import numpy as np
import scipy

from invariant_region_mlp.experiments.mechanism_suite.mechanism_mlp import (
    method_from_dict,
)

from .equations import (
    BASE_SEED,
    C3_REFERENCE,
    POINT_COUNT,
    RESULTS_ROOT,
    make_equation,
    make_points,
)
from .gates import analytic_gates
from .mlp import (
    IMPLEMENTATION_REVISION,
    PREREGISTRATION_COMMIT,
    load_repetition,
    run_single_repetition,
    save_repetition,
)
from .protocol import (
    D400_PILOT_PATH,
    GATES_PATH,
    RAW_ROOT,
    S1_CONFIGS,
    make_task,
    primary_method,
    regular_methods,
    s0_tasks,
    s1_tasks,
    s2_tasks,
    s3_tasks,
    s4_tasks,
    task_identity,
    task_path,
)
from .reference_c3 import build_and_audit_reference, load_reference


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parent
FULL_HISTORY_SOURCE = (
    HERE.parent / "active_vb_high_budget" / "vb_mlp_methods.py"
)
ROUND1_GATES = PACKAGE_ROOT / "results" / "mechanism_suite" / "gates.json"
INITIAL_GATES_AUDIT = (
    RESULTS_ROOT / "audit_history" / "gates_initial_c3_terminal_spline.json"
)


def _git(*args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=REPOSITORY_ROOT, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


EXECUTION_COMMIT = _git("rev-parse", "HEAD")


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
            return "Infinity" if number > 0.0 else "-Infinity"
        return number
    return value


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def ensure_c3_reference(*, dry_run: bool = False) -> dict[str, Any] | None:
    if C3_REFERENCE.exists():
        _, _, _, metadata = load_reference(C3_REFERENCE)
        parameters = metadata.get("parameters", {})
        expected = {"kappa": 0.5, "A": 2.0, "beta": 2.0, "T": 0.25}
        if any(
            not math.isclose(float(parameters.get(key, float("nan"))), value)
            for key, value in expected.items()
        ):
            raise RuntimeError("existing C3 reference has incompatible parameters")
        print(
            f"C3 reference exists: passed={metadata.get('passed')} "
            f"finest={metadata.get('finest_richardson_difference_L8')}",
            flush=True,
        )
        return metadata
    if dry_run:
        print(f"C3 reference missing: {C3_REFERENCE}", flush=True)
        return None
    print("building C3 four-grid, two-domain production reference", flush=True)
    metadata = build_and_audit_reference(C3_REFERENCE, seed=BASE_SEED)
    print(
        "C3 reference gate "
        f"passed={metadata['passed']} "
        f"finest={metadata['finest_richardson_difference_L8']:.3e} "
        f"domain={metadata['domain_L6_L8_max_difference']:.3e} "
        f"fd={metadata['independent_fd_residual']['max_abs_residual']:.3e}",
        flush=True,
    )
    return metadata


def _artifact_matches(path: Path, task: dict[str, Any]) -> bool:
    if not path.exists():
        return False
    try:
        loaded = load_repetition(path)
    except Exception as error:
        raise RuntimeError(f"cannot validate existing artifact {path}") from error
    metadata = loaded["metadata"]
    matches = (
        metadata.get("round") == 3
        and metadata.get("implementation_revision") == IMPLEMENTATION_REVISION
        and metadata.get("preregistration_commit") == PREREGISTRATION_COMMIT
        and metadata.get("pde_id") == task["pde_id"]
        and int(metadata.get("dimension", -1)) == int(task["d"])
        and int(metadata.get("n", -1)) == int(task["n"])
        and int(metadata.get("M", -1)) == int(task["M"])
        and int(metadata.get("repetition", -1)) == int(task["repetition"])
        and metadata.get("method") == task["method"]
        and int(metadata.get("base_seed", -1)) == BASE_SEED
        and int(metadata.get("chunk_size", -1)) == int(task["chunk_size"])
        and metadata.get("dtype") == "float64"
    )
    if not matches:
        raise RuntimeError(f"refusing to overwrite non-matching artifact {path}")
    if (
        loaded["prediction_u"].shape != (POINT_COUNT,)
        or loaded["truth_u"].shape != (POINT_COUNT,)
        or loaded["is_validation"].shape != (POINT_COUNT,)
        or loaded["prediction_u"].dtype != np.float64
        or loaded["truth_u"].dtype != np.float64
        or loaded["is_validation"].dtype != np.bool_
    ):
        raise RuntimeError(f"artifact shape/dtype mismatch: {path}")
    return True


def _execute_task(task: dict[str, Any]) -> dict[str, Any]:
    path = task_path(task)
    if _artifact_matches(path, task):
        return {"status": "skipped", "path": str(path), "seconds": 0.0}
    if path.exists():
        raise RuntimeError(f"refusing to overwrite {path}")
    equation = make_equation(str(task["pde_id"]), int(task["d"]))
    points = make_points(
        equation, pde_id=str(task["pde_id"]), n_points=POINT_COUNT,
        seed=BASE_SEED,
    )
    result = run_single_repetition(
        pde_id=str(task["pde_id"]),
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
        study=str(task["study"]),
        code_commit=EXECUTION_COMMIT,
    )
    result["metadata"]["point_set"] = {
        "seed": BASE_SEED,
        "count": POINT_COUNT,
        "validation_fraction": 0.2,
        "pde_id_seed_component": str(task["pde_id"]),
    }
    save_repetition(path, result)
    return {
        "status": "completed",
        "path": str(path),
        "seconds": result["metadata"]["wall_clock_seconds"],
        "skill": result["metadata"]["metrics"]["test"]["skill"],
    }


def _deduplicate(tasks: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[tuple[Any, ...], dict[str, Any]] = {}
    for task in tasks:
        key = task_identity(task)
        if key in unique:
            continue
        unique[key] = task
    return list(unique.values())


def _run_pool(
    tasks: list[dict[str, Any]], *, workers: int, label: str
) -> None:
    missing = [task for task in tasks if not _artifact_matches(task_path(task), task)]
    print(
        f"{label}: total={len(tasks)} complete={len(tasks)-len(missing)} "
        f"missing={len(missing)} workers={workers}",
        flush=True,
    )
    if not missing:
        return
    # Start expensive cells first so short tails do not strand worker slots.
    missing.sort(
        key=lambda task: (
            int(task["M"]) ** int(task["n"]) * int(task["d"]),
            int(task["n"]),
        ),
        reverse=True,
    )
    started = time.perf_counter()
    worker_seconds = 0.0
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_execute_task, task): task for task in missing}
        for completed, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            result = future.result()
            worker_seconds += float(result.get("seconds", 0.0))
            if completed == 1 or completed % 10 == 0 or completed == len(missing):
                elapsed = time.perf_counter() - started
                rate = completed / max(elapsed, 1.0e-12)
                eta = (len(missing) - completed) / max(rate, 1.0e-12)
                print(
                    f"[{completed}/{len(missing)}] {task['pde_id']} "
                    f"d={task['d']} n={task['n']} M={task['M']} "
                    f"{task['method']['name']} rep={task['repetition']} "
                    f"elapsed={elapsed/60:.1f}m naive_eta={eta/60:.1f}m",
                    flush=True,
                )
    print(
        f"{label}: worker-summed={worker_seconds/3600:.2f}h "
        f"wall={(time.perf_counter()-started)/3600:.2f}h",
        flush=True,
    )


def _materialize_n1_aliases(tasks: list[dict[str, Any]]) -> int:
    aliases = 0
    for task in tasks:
        if int(task["n"]) != 1 or task["method"]["name"] == "raw":
            continue
        destination = task_path(task)
        if _artifact_matches(destination, task):
            continue
        source_task = dict(task)
        source_task["method"] = {"name": "raw", "transform": "raw", "factor": 1.0}
        source = task_path(source_task)
        if not _artifact_matches(source, source_task):
            continue
        loaded = load_repetition(source)
        loaded["metadata"] = dict(loaded["metadata"])
        loaded["metadata"]["method"] = dict(task["method"])
        loaded["metadata"]["study_origin"] = str(task["study"])
        loaded["metadata"]["exact_reuse"] = {
            "source": source.relative_to(RESULTS_ROOT).as_posix(),
            "reason": "n=1 has no generator evaluation, so every transform is exactly inactive",
        }
        save_repetition(destination, loaded)
        aliases += 1
    return aliases


def run_tasks(
    tasks: Iterable[dict[str, Any]], *, workers: int, dry_run: bool, label: str
) -> None:
    tasks = _deduplicate(tasks)
    if dry_run:
        complete = sum(_artifact_matches(task_path(t), t) for t in tasks)
        print(f"{label}: total={len(tasks)} complete={complete} missing={len(tasks)-complete}")
        return
    n1_raw: dict[tuple[Any, ...], dict[str, Any]] = {}
    for task in tasks:
        if int(task["n"]) == 1:
            raw = dict(task)
            raw["method"] = {"name": "raw", "transform": "raw", "factor": 1.0}
            n1_raw[task_identity(raw)] = raw
    if n1_raw:
        _run_pool(list(n1_raw.values()), workers=workers, label=f"{label}/n1-raw")
        aliases = _materialize_n1_aliases(tasks)
        print(f"{label}: materialized {aliases} exact n=1 aliases", flush=True)
    _run_pool(tasks, workers=workers, label=label)


def _mean_test_skill(task: dict[str, Any]) -> float:
    return float(load_repetition(task_path(task))["metadata"]["metrics"]["test"]["skill"])


def _carry_anchor(pde_id: str, round1: dict[str, Any]) -> dict[str, Any]:
    source = round1["candidate_summary"][pde_id]
    return {
        "source": "round 1 gates carried without rerun",
        "source_path": "invariant_region_mlp/results/mechanism_suite/gates.json",
        "gate_status": source["gate_status"],
        "failed_gates": source["failed_gates"],
        "admitted": source["admitted"],
        "instances": source["instances"],
    }


def finalize_gates(analytic: dict[str, dict[str, Any]]) -> dict[str, Any]:
    summaries: dict[str, Any] = {}
    reference_passed = bool(load_reference(C3_REFERENCE)[3]["passed"])
    for pde_id, dimensions in analytic.items():
        instances: dict[str, Any] = {}
        for dimension_text, result in dimensions.items():
            d = int(dimension_text)
            cells = []
            if pde_id != "C3" or reference_passed:
                for n, M in ((3, 6), (4, 6)):
                    methods = {
                        name: next(
                            method for method in (
                                [
                                    {"name": "raw", "transform": "raw", "factor": 1.0},
                                    primary_method(pde_id).to_dict(),
                                    {"name": "oracle_state", "transform": "oracle_state", "factor": 1.0},
                                ]
                            ) if method["name"] == name
                        )
                        for name in ("raw", "oracle_state")
                    }
                    values: dict[str, list[float]] = {}
                    for name, method in methods.items():
                        values[name] = []
                        for repetition in range(3):
                            task = make_task("s0", pde_id, d, n, M, method_from_dict(method), repetition)
                            if not _artifact_matches(task_path(task), task):
                                raise RuntimeError(f"missing S0 gate artifact {task_path(task)}")
                            values[name].append(_mean_test_skill(task))
                    means = {key: float(np.mean(value)) for key, value in values.items()}
                    ratio = means["raw"] / means["oracle_state"]
                    cells.append({
                        "n": n,
                        "M": M,
                        "mean_test_skill": means,
                        "raw_over_oracle": ratio,
                        "G4_cell_passed": bool(ratio >= 2.0),
                    })
            g3 = bool(cells and any(
                cell["mean_test_skill"]["oracle_state"] <= 0.10 for cell in cells
            ))
            g4 = bool(cells and all(cell["G4_cell_passed"] for cell in cells))
            instances[dimension_text] = {
                "analytic": result,
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
                for instance in instances.values()
            )
            for gate in ("G1", "G2", "G3", "G4", "G5", "G6")
        }
        failures = [gate for gate, passed in gate_status.items() if not passed]
        summaries[pde_id] = {
            "gate_status": gate_status,
            "failed_gates": failures,
            "admitted": not failures,
            "instances": instances,
        }

    round1 = json.loads(ROUND1_GATES.read_text(encoding="utf-8"))
    summaries = {
        "P1": _carry_anchor("P1", round1),
        "P4": _carry_anchor("P4", round1),
        **summaries,
    }
    payload = {
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "protocol": {
            "base_seed": BASE_SEED,
            "point_count": POINT_COUNT,
            "validation_fraction": 0.2,
            "headline_configurations": [[3, 6], [4, 6]],
            "repetitions": 3,
            "dtype": "float64",
            "tree_seed": "SeedSequence([20261207,d,n,M,rep,chunk_index])",
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_count": os.cpu_count(),
            "git_branch": _git("branch", "--show-current"),
            "git_commit": _git("rev-parse", "HEAD"),
        },
        "preregistration_commit": PREREGISTRATION_COMMIT,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "unchanged_recursion": {
            "path": FULL_HISTORY_SOURCE.relative_to(REPOSITORY_ROOT).as_posix(),
            "sha256": hashlib.sha256(FULL_HISTORY_SOURCE.read_bytes()).hexdigest(),
        },
        "candidates": summaries,
        "LQG_G_prediction": "may fail G2 because the d-scaled time-only term dominates std(u)",
    }
    write_json_atomic(GATES_PATH, payload)
    return payload


def run_s0(*, workers: int, dry_run: bool) -> None:
    reference = ensure_c3_reference(dry_run=dry_run)
    if dry_run:
        run_tasks(s0_tasks(), workers=workers, dry_run=True, label="S0")
        return
    analytic: dict[str, dict[str, Any]] = {}
    for pde_id in ("C1", "C2-convex", "C2-cancel", "C2-flip", "C3", "LQG"):
        for d in (20, 100):
            print(f"analytic gates {pde_id} d={d}", flush=True)
            analytic.setdefault(pde_id, {})[str(d)] = analytic_gates(pde_id, d)
    tasks = s0_tasks()
    if not reference or not bool(reference["passed"]):
        print("C3 production reference failed; stopping C3 before MLP S0", flush=True)
        tasks = [task for task in tasks if task["pde_id"] != "C3"]
    run_tasks(tasks, workers=workers, dry_run=False, label="S0")
    if GATES_PATH.exists() and not INITIAL_GATES_AUDIT.exists():
        INITIAL_GATES_AUDIT.parent.mkdir(parents=True, exist_ok=True)
        INITIAL_GATES_AUDIT.write_bytes(GATES_PATH.read_bytes())
    payload = finalize_gates(analytic)
    for pde_id, summary in payload["candidates"].items():
        print(
            f"{pde_id}: {'ADMIT' if summary['admitted'] else 'EXCLUDE'} "
            f"failed={summary['failed_gates']}",
            flush=True,
        )


def _load_gate_admission() -> dict[str, bool]:
    if not GATES_PATH.exists():
        raise RuntimeError("S0 gates must be completed first")
    payload = json.loads(GATES_PATH.read_text(encoding="utf-8"))
    return {
        pde_id: bool(summary["admitted"])
        for pde_id, summary in payload["candidates"].items()
    }


def _filter_admitted(
    tasks: Iterable[dict[str, Any]], admission: dict[str, bool]
) -> list[dict[str, Any]]:
    return [task for task in tasks if admission.get(task["pde_id"], False)]


def _run_d400_pilots(*, dry_run: bool) -> dict[str, Any]:
    if D400_PILOT_PATH.exists():
        payload = json.loads(D400_PILOT_PATH.read_text(encoding="utf-8"))
        if payload.get("implementation_revision") != IMPLEMENTATION_REVISION:
            raise RuntimeError("existing d=400 pilot file has another revision")
        return payload
    pilots = []
    if dry_run:
        return {"dropped_configs": [], "pilots": [], "dry_run": True}
    for n, M in S1_CONFIGS:
        task = make_task(
            "s1_pilot", "P1", 400, n, M,
            method_from_dict({"name": "raw", "transform": "raw", "factor": 1.0}),
            0,
        )
        print(f"d=400 cutoff pilot n={n} M={M}", flush=True)
        _execute_task(task)
        metadata = load_repetition(task_path(task))["metadata"]
        seconds = float(metadata["wall_clock_seconds"])
        pilots.append({"n": n, "M": M, "seconds": seconds, "drop": seconds > 1800.0})
        print(f"d=400 cutoff pilot n={n} M={M}: {seconds:.1f}s", flush=True)
    dropped = [[item["n"], item["M"]] for item in pilots if item["drop"]]
    payload = {
        "schema_version": 1,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "code_commit": _git("rev-parse", "HEAD"),
        "representative_pde": "P1",
        "representative_method": "raw",
        "dimension": 400,
        "threshold_seconds": 1800.0,
        "pilots": pilots,
        "dropped_configs": dropped,
    }
    write_json_atomic(D400_PILOT_PATH, payload)
    return payload


def run_stage(stage: str, *, workers: int, dry_run: bool) -> None:
    if stage == "reference":
        ensure_c3_reference(dry_run=dry_run)
        return
    if stage == "s0":
        run_s0(workers=workers, dry_run=dry_run)
        return
    admission = _load_gate_admission()
    if stage == "s2":
        tasks = _filter_admitted(s2_tasks(), admission)
        run_tasks(tasks, workers=workers, dry_run=dry_run, label="S2")
        if not dry_run:
            from .onestep import OUTPUT_PATH, run_all

            if OUTPUT_PATH.exists():
                print(f"S2 one-step check exists: {OUTPUT_PATH}", flush=True)
            else:
                run_all(OUTPUT_PATH)
        return
    if stage == "s1":
        include_lqg = admission.get("LQG", False)
        tasks = _filter_admitted(s1_tasks(include_lqg=include_lqg), admission)
        lower = [task for task in tasks if int(task["d"]) in {20, 100}]
        run_tasks(lower, workers=workers, dry_run=dry_run, label="S1 d=20,100")
        pilots = _run_d400_pilots(dry_run=dry_run)
        dropped = {tuple(value) for value in pilots.get("dropped_configs", [])}
        high = [
            task for task in tasks
            if int(task["d"]) == 400
            and (int(task["n"]), int(task["M"])) not in dropped
        ]
        run_tasks(high, workers=max(1, min(workers, 4)), dry_run=dry_run, label="S1 d=400")
        return
    if stage == "s3":
        tasks = _filter_admitted(s3_tasks(), admission)
        run_tasks(tasks, workers=workers, dry_run=dry_run, label="S3")
        return
    if stage == "s4":
        tasks = _filter_admitted(s4_tasks(), admission)
        run_tasks(tasks, workers=workers, dry_run=dry_run, label="S4")
        return
    raise ValueError(stage)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("reference", "s0", "s2", "s1", "s3", "s4", "all"),
        default="all",
    )
    parser.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    stages = (
        ("s0", "s2", "s1", "s3", "s4")
        if args.stage == "all" else (args.stage,)
    )
    for stage in stages:
        print(f"=== {stage.upper()} ===", flush=True)
        run_stage(stage, workers=args.workers, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
