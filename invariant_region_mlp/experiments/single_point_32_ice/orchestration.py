"""Deterministic four-GPU dispatch and stage gates for the ICE experiment."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import socket
import time
import traceback
from typing import Any, Iterator

import torch

from .protocol import (
    DIMENSIONS,
    FH_CONFIGURATIONS,
    FH_MAIN_CONFIGURATION,
    FH_REPLICATES,
    PDES,
    QUERY_COUNT,
    QUEUE_DIR,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    STAGE_B_AUDIT,
    STATUS_DIR,
    TA_AUDIT_EULER_STEPS,
    TA_INFERENCE_REPLICATES,
    TA_PRIMARY_EULER_STEPS,
    TRAINING_SEEDS,
    ensure_directories,
    require_stage_a,
    require_stage_b,
    write_json_atomic,
)


SLOTS = 4
DIMENSION_AUDIT_DIR = RESULTS_DIR / "dimension_audits"
PIPELINE_JOBS = RESULTS_DIR / "pipeline_jobs.json"
_ALGORITHM_FILES = (
    "protocol.py",
    "benchmarks.py",
    "deep_bsde.py",
    "fh_mlp.py",
    "tangent_actor.py",
    "orchestration.py",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256(encoded).hexdigest()


def algorithm_hashes() -> dict[str, str]:
    root = Path(__file__).resolve().parent
    return {
        name: sha256(
            (root / name).read_text(encoding="utf-8").replace("\r\n", "\n").encode("utf-8")
        ).hexdigest()
        for name in _ALGORITHM_FILES
    }


def _task_id(task: dict[str, Any]) -> str:
    stage = task["stage"].lower()
    method = task["method"]
    pde = task["pde"]
    d = task["dimension"]
    if method == "tangent_actor":
        suffix = f"s{task['seed']}"
    elif method == "deep_bsde":
        suffix = f"q{task['point']:02d}-s{task['seed']}"
    else:
        suffix = f"n{task['n']}-M{task['M']}-r{task['replicate']:02d}"
    return f"{stage}-{method}-{pde}-d{d}-{suffix}"


def _ta_weight(pde: str, d: int, pilot: bool) -> float:
    if pilot:
        return 1.0
    family = {
        "hjb_log": 1.0,
        "allen_cahn": 1.15,
        "b6_fullrank": 0.6,
        "dense_fourier": 1.2,
    }[pde]
    return family * (1.0 + d / 800.0) * 15.0


def _deep_weight(pde: str, d: int, pilot: bool) -> float:
    if pilot:
        return 0.15
    iterations = 2.0 if pde == "allen_cahn" else 1.0
    return iterations * (1.0 + d / 250.0)


def _fh_weight(pde: str, d: int, n: int, M: int, pilot: bool) -> float:
    family = 1.5 if pde in {"hjb_log", "b6_fullrank", "dense_fourier"} else 1.0
    point_factor = 4 if pilot else QUERY_COUNT
    return family * point_factor * (M ** n) * (1.0 + d / 800.0) / 5000.0


def stage_b_tasks() -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for pde in PDES:
        tasks.append({
            "stage": "B",
            "method": "tangent_actor",
            "pde": pde,
            "dimension": 100,
            "seed": 0,
            "pilot": True,
            "weight": _ta_weight(pde, 100, True),
        })
        for point in range(4):
            tasks.append({
                "stage": "B",
                "method": "deep_bsde",
                "pde": pde,
                "dimension": 100,
                "point": point,
                "seed": 0,
                "pilot": True,
                "weight": _deep_weight(pde, 100, True),
            })
        for n, M in FH_CONFIGURATIONS:
            tasks.append({
                "stage": "B",
                "method": "fh_mlp_standard",
                "pde": pde,
                "dimension": 100,
                "n": n,
                "M": M,
                "replicate": 0,
                "pilot": True,
                "weight": _fh_weight(pde, 100, n, M, True),
            })
    for task in tasks:
        task["task_id"] = _task_id(task)
    return tasks


def formal_tasks(dimension: int) -> list[dict[str, Any]]:
    if dimension not in DIMENSIONS:
        raise ValueError(dimension)
    tasks: list[dict[str, Any]] = []
    for pde in PDES:
        for seed in TRAINING_SEEDS:
            tasks.append({
                "stage": "C",
                "method": "tangent_actor",
                "pde": pde,
                "dimension": dimension,
                "seed": seed,
                "pilot": False,
                "weight": _ta_weight(pde, dimension, False),
            })
        for point in range(QUERY_COUNT):
            for seed in TRAINING_SEEDS:
                tasks.append({
                    "stage": "C",
                    "method": "deep_bsde",
                    "pde": pde,
                    "dimension": dimension,
                    "point": point,
                    "seed": seed,
                    "pilot": False,
                    "weight": _deep_weight(pde, dimension, False),
                })
        for n, M in FH_CONFIGURATIONS:
            for replicate in range(FH_REPLICATES):
                tasks.append({
                    "stage": "C",
                    "method": "fh_mlp_standard",
                    "pde": pde,
                    "dimension": dimension,
                    "n": n,
                    "M": M,
                    "replicate": replicate,
                    "pilot": False,
                    "weight": _fh_weight(pde, dimension, n, M, False),
                })
    for task in tasks:
        task["task_id"] = _task_id(task)
    return tasks


def assign_slots(tasks: list[dict[str, Any]], slots: int = SLOTS) -> list[dict[str, Any]]:
    """Accuracy-blind longest-processing-time assignment with stable ties."""
    if slots <= 0:
        raise ValueError(slots)
    loads = [0.0] * slots
    sequences = [0] * slots
    assigned: list[dict[str, Any]] = []
    for source in sorted(tasks, key=lambda row: (-float(row["weight"]), row["task_id"])):
        slot = min(range(slots), key=lambda index: (loads[index], index))
        task = dict(source)
        task["slot"] = slot
        task["sequence"] = sequences[slot]
        sequences[slot] += 1
        loads[slot] += float(task["weight"])
        assigned.append(task)
    return sorted(assigned, key=lambda row: (row["slot"], row["sequence"]))


def build_manifest(stage: str, dimension: int | None, slots: int = SLOTS) -> dict[str, Any]:
    normalized = stage.upper()
    if normalized == "B":
        if dimension not in {None, 100}:
            raise ValueError("Stage B is fixed at d=100")
        tasks = stage_b_tasks()
        dimension = 100
    elif normalized == "C":
        if dimension is None:
            raise ValueError("Stage C requires a dimension")
        tasks = formal_tasks(dimension)
    else:
        raise KeyError(stage)
    assigned = assign_slots(tasks, slots)
    payload: dict[str, Any] = {
        "revision": REVISION,
        "stage": normalized,
        "dimension": dimension,
        "slots": slots,
        "task_count": len(assigned),
        "slot_task_counts": [
            sum(row["slot"] == slot for row in assigned) for slot in range(slots)
        ],
        "slot_weights": [
            sum(
                float(row["weight"])
                for row in assigned
                if row["slot"] == slot
            )
            for slot in range(slots)
        ],
        "ordering": "accuracy-blind deterministic LPT with task-id tie break",
        "algorithm_hashes": algorithm_hashes(),
        "tasks": assigned,
    }
    payload["manifest_sha256"] = _digest(payload)
    return payload


def manifest_path(stage: str, dimension: int | None = None) -> Path:
    if stage.upper() == "B":
        return QUEUE_DIR / "stage_b_d100.json"
    if dimension is None:
        raise ValueError("dimension required")
    return QUEUE_DIR / f"stage_c_d{dimension}.json"


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    if path.is_file():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != payload:
            raise RuntimeError(f"refusing to rewrite frozen queue: {path}")
        return
    write_json_atomic(path, payload)


def load_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    claimed = payload.get("manifest_sha256")
    unsigned = dict(payload)
    unsigned.pop("manifest_sha256", None)
    if claimed != _digest(unsigned):
        raise RuntimeError(f"manifest hash mismatch: {path}")
    if payload.get("revision") != REVISION:
        raise RuntimeError(f"manifest revision mismatch: {path}")
    if payload.get("algorithm_hashes") != algorithm_hashes():
        raise RuntimeError(f"manifest algorithm hash mismatch: {path}")
    return payload


def _ta_records(task: dict[str, Any]) -> list[Path]:
    pilot = bool(task["pilot"])
    tag = "pilot" if pilot else "formal"
    pde, d, seed = task["pde"], task["dimension"], task["seed"]
    pairs = [(0, 4)] if pilot else [
        *((replicate, TA_PRIMARY_EULER_STEPS)
          for replicate in range(TA_INFERENCE_REPLICATES)),
        *((0, steps) for steps in TA_AUDIT_EULER_STEPS),
    ]
    root = RAW_DIR / tag / "tangent_actor_evaluations"
    return [
        root
        / f"ta__{pde}__d{d}__s{seed}__rep{replicate:02d}__steps{steps}__{tag}.json"
        for replicate, steps in pairs
    ]


def expected_records(task: dict[str, Any]) -> list[Path]:
    tag = "pilot" if task["pilot"] else "formal"
    pde, d = task["pde"], task["dimension"]
    if task["method"] == "tangent_actor":
        return _ta_records(task)
    if task["method"] == "deep_bsde":
        run_id = (
            f"deep_bsde__{pde}__d{d}__q{task['point']:02d}"
            f"__s{task['seed']}__{tag}"
        )
        return [RAW_DIR / tag / "deep_bsde" / f"{run_id}.json"]
    run_id = (
        f"fh_mlp_standard__{pde}__d{d}__n{task['n']}_M{task['M']}"
        f"__r{task['replicate']:02d}__{tag}"
    )
    return [RAW_DIR / tag / "fh_mlp" / f"{run_id}.json"]


def _record_terminal(
    path: Path, task: dict[str, Any]
) -> tuple[bool, dict[str, Any] | None]:
    if not path.is_file():
        return False, None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False, None
    expected = "pilot_complete" if task["pilot"] else "complete"
    accepted = {expected}
    if task["method"] == "fh_mlp_standard":
        accepted.add("numerical_failure")
    if record.get("status") not in accepted:
        return False, record
    prediction = Path(str(record.get("prediction_file", "")))
    if not prediction.is_file():
        return False, record
    return True, record


def task_complete(task: dict[str, Any]) -> bool:
    return all(_record_terminal(path, task)[0] for path in expected_records(task))


@contextmanager
def task_lock(task_id: str) -> Iterator[None]:
    root = STATUS_DIR / "locks"
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{task_id}.lock"
    handle = path.open("a+b")
    try:
        if os.name == "nt":
            import msvcrt
            if path.stat().st_size == 0:
                handle.write(b"0")
                handle.flush()
            handle.seek(0)
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as error:
                raise RuntimeError(f"duplicate writer for {task_id}") from error
        else:
            import fcntl
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError(f"duplicate writer for {task_id}") from error
        yield
    finally:
        if os.name == "nt":
            import msvcrt
            handle.seek(0)
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            except OSError:
                pass
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def _task_status_path(task: dict[str, Any]) -> Path:
    return STATUS_DIR / task["stage"].lower() / f"{task['task_id']}.json"


def _environment() -> dict[str, Any]:
    result: dict[str, Any] = {
        "hostname": socket.gethostname(),
        "pid": os.getpid(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "cuda_available": torch.cuda.is_available(),
    }
    if torch.cuda.is_available():
        index = torch.cuda.current_device()
        properties = torch.cuda.get_device_properties(index)
        result.update({
            "cuda_device_index": index,
            "cuda_device_name": properties.name,
            "cuda_total_memory": properties.total_memory,
        })
    return result


def _run_solver(task: dict[str, Any], device: torch.device) -> list[dict[str, Any]]:
    pde = task["pde"]
    d = int(task["dimension"])
    pilot = bool(task["pilot"])
    if task["method"] == "tangent_actor":
        from .tangent_actor import evaluate, train
        records = [train(pde, d, int(task["seed"]), device, pilot=pilot)]
        if pilot:
            records.append(
                evaluate(pde, d, int(task["seed"]), 0, device, pilot=True)
            )
        else:
            for replicate in range(TA_INFERENCE_REPLICATES):
                records.append(evaluate(
                    pde,
                    d,
                    int(task["seed"]),
                    replicate,
                    device,
                    pilot=False,
                    time_steps=TA_PRIMARY_EULER_STEPS,
                ))
            for steps in TA_AUDIT_EULER_STEPS:
                records.append(evaluate(
                    pde,
                    d,
                    int(task["seed"]),
                    0,
                    device,
                    pilot=False,
                    time_steps=steps,
                ))
        return records
    if task["method"] == "deep_bsde":
        from .deep_bsde import execute
        return [execute(
            pde,
            d,
            int(task["point"]),
            int(task["seed"]),
            device,
            pilot=pilot,
        )]
    from .fh_mlp import execute
    return [execute(
        pde,
        d,
        int(task["n"]),
        int(task["M"]),
        int(task["replicate"]),
        device,
        pilot=pilot,
    )]


def execute_task(task: dict[str, Any], device: torch.device) -> dict[str, Any]:
    require_stage_a()
    if not task["pilot"]:
        require_stage_b()
    ensure_directories()
    status_path = _task_status_path(task)
    with task_lock(task["task_id"]):
        if task_complete(task):
            if status_path.is_file():
                existing = json.loads(status_path.read_text(encoding="utf-8"))
                if existing.get("status") == "COMPLETE":
                    return existing
            payload = {
                "revision": REVISION,
                "task": task,
                "status": "COMPLETE",
                "skipped_existing": True,
                "checked_utc": utc_now(),
                "environment": _environment(),
            }
            write_json_atomic(status_path, payload)
            return payload
        started_utc = utc_now()
        started = time.perf_counter()
        try:
            records = _run_solver(task, device)
            if not task_complete(task):
                raise RuntimeError("solver returned without all terminal artifacts")
            payload = {
                "revision": REVISION,
                "task": task,
                "status": "COMPLETE",
                "started_utc": started_utc,
                "finished_utc": utc_now(),
                "wall_seconds": time.perf_counter() - started,
                "environment": _environment(),
                "record_summaries": [
                    {
                        key: record.get(key)
                        for key in (
                            "run_id",
                            "status",
                            "nonfinite_count",
                            "training_seconds",
                            "inference_seconds",
                            "peak_gpu_bytes",
                        )
                    }
                    for record in records
                ],
            }
            write_json_atomic(status_path, payload)
            return payload
        except BaseException as error:
            payload = {
                "revision": REVISION,
                "task": task,
                "status": "FAILED",
                "started_utc": started_utc,
                "finished_utc": utc_now(),
                "wall_seconds": time.perf_counter() - started,
                "environment": _environment(),
                "exception_type": type(error).__name__,
                "exception": str(error),
                "traceback": traceback.format_exc(),
            }
            write_json_atomic(status_path, payload)
            raise


def run_slot(manifest: dict[str, Any], slot: int, device: torch.device) -> None:
    if not 0 <= slot < int(manifest["slots"]):
        raise ValueError(slot)
    tasks = [row for row in manifest["tasks"] if int(row["slot"]) == slot]
    tasks.sort(key=lambda row: int(row["sequence"]))
    print(json.dumps({
        "event": "slot_start",
        "slot": slot,
        "tasks": len(tasks),
        "manifest_sha256": manifest["manifest_sha256"],
        "utc": utc_now(),
    }, sort_keys=True), flush=True)
    for index, task in enumerate(tasks):
        print(json.dumps({
            "event": "task_start",
            "slot": slot,
            "index": index,
            "task_id": task["task_id"],
            "utc": utc_now(),
        }, sort_keys=True), flush=True)
        result = execute_task(task, device)
        print(json.dumps({
            "event": "task_end",
            "slot": slot,
            "index": index,
            "task_id": task["task_id"],
            "status": result["status"],
            "utc": utc_now(),
        }, sort_keys=True), flush=True)


def run_index(manifest: dict[str, Any], index: int, device: torch.device) -> None:
    """Run one immutable task; Slurm's array throttle supplies four GPUs."""
    tasks = sorted(
        manifest["tasks"],
        key=lambda row: (int(row["sequence"]), int(row["slot"])),
    )
    if not 0 <= index < len(tasks):
        raise ValueError(index)
    task = tasks[index]
    print(json.dumps({
        "event": "task_start",
        "array_index": index,
        "task_id": task["task_id"],
        "frozen_slot": task["slot"],
        "frozen_sequence": task["sequence"],
        "manifest_sha256": manifest["manifest_sha256"],
        "utc": utc_now(),
    }, sort_keys=True), flush=True)
    result = execute_task(task, device)
    print(json.dumps({
        "event": "task_end",
        "array_index": index,
        "task_id": task["task_id"],
        "status": result["status"],
        "utc": utc_now(),
    }, sort_keys=True), flush=True)


def _audit_task(task: dict[str, Any]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for path in expected_records(task):
        terminal, record = _record_terminal(path, task)
        rows.append({
            "file": path.relative_to(RESULTS_DIR).as_posix(),
            "terminal": terminal,
            "status": None if record is None else record.get("status"),
            "nonfinite_count": (
                None if record is None else record.get("nonfinite_count")
            ),
        })
    return {
        "task_id": task["task_id"],
        "method": task["method"],
        "pde": task["pde"],
        "terminal": all(row["terminal"] for row in rows),
        "records": rows,
    }


def audit_manifest(manifest: dict[str, Any]) -> dict[str, Any]:
    audits = [_audit_task(task) for task in manifest["tasks"]]
    missing = [row["task_id"] for row in audits if not row["terminal"]]
    checks: dict[str, Any] = {
        "all_tasks_terminal": {"passed": not missing, "missing": missing},
        "task_count": {
            "expected": int(manifest["task_count"]),
            "observed": len(audits),
            "passed": len(audits) == int(manifest["task_count"]),
        },
    }
    if manifest["stage"] == "B":
        ta_bad: list[str] = []
        deep_bad: list[str] = []
        fh_main_bad: list[str] = []
        for task, audit in zip(manifest["tasks"], audits, strict=True):
            if not audit["terminal"]:
                continue
            if task["method"] == "tangent_actor":
                if any(
                    row["status"] != "pilot_complete"
                    or row["nonfinite_count"] not in {0, None}
                    for row in audit["records"]
                ):
                    ta_bad.append(task["task_id"])
            elif task["method"] == "deep_bsde":
                if any(
                    row["status"] != "pilot_complete"
                    or row["nonfinite_count"] not in {0, None}
                    for row in audit["records"]
                ):
                    deep_bad.append(task["task_id"])
            elif (task["n"], task["M"]) == FH_MAIN_CONFIGURATION:
                if any(
                    row["status"] != "pilot_complete"
                    for row in audit["records"]
                ):
                    fh_main_bad.append(task["task_id"])
        checks["tangent_actor_finite"] = {
            "passed": not ta_bad,
            "failed": ta_bad,
        }
        checks["deep_bsde_finite"] = {
            "passed": not deep_bad,
            "failed": deep_bad,
        }
        checks["fh_main_configuration_finite"] = {
            "passed": not fh_main_bad,
            "failed": fh_main_bad,
        }
    else:
        non_fh_bad = [
            task["task_id"]
            for task, audit in zip(manifest["tasks"], audits, strict=True)
            if task["method"] != "fh_mlp_standard" and not audit["terminal"]
        ]
        checks["trained_methods_complete"] = {
            "passed": not non_fh_bad,
            "failed": non_fh_bad,
        }
    overall = all(bool(check["passed"]) for check in checks.values())
    payload = {
        "revision": REVISION,
        "stage": manifest["stage"],
        "dimension": manifest["dimension"],
        "manifest_sha256": manifest["manifest_sha256"],
        "algorithm_hashes": algorithm_hashes(),
        "checked_utc": utc_now(),
        "checks": checks,
        "method_terminal_counts": {
            method: sum(
                row["terminal"] and row["method"] == method for row in audits
            )
            for method in (
                "tangent_actor",
                "deep_bsde",
                "fh_mlp_standard",
            )
        },
        "fh_numerical_failure_records": sum(
            record["status"] == "numerical_failure"
            for audit in audits
            if audit["method"] == "fh_mlp_standard"
            for record in audit["records"]
        ),
        "overall_pass": overall,
        "status": "VALIDATED" if overall else "FAILED",
    }
    output = (
        STAGE_B_AUDIT
        if manifest["stage"] == "B"
        else DIMENSION_AUDIT_DIR / f"d{manifest['dimension']}.json"
    )
    write_json_atomic(output, payload)
    if not overall:
        raise RuntimeError(f"{manifest['stage']} audit failed: {output}")
    return payload


def record_pipeline(entries: list[str]) -> None:
    jobs: dict[str, str] = {}
    for entry in entries:
        name, separator, job_id = entry.partition("=")
        if not separator or not name or not job_id:
            raise ValueError(entry)
        jobs[name] = job_id
    write_json_atomic(PIPELINE_JOBS, {
        "revision": REVISION,
        "recorded_utc": utc_now(),
        "jobs": jobs,
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    freeze = subparsers.add_parser("freeze-manifest")
    freeze.add_argument("--stage", choices=("B", "C"), required=True)
    freeze.add_argument("--dimension", type=int)
    freeze.add_argument("--slots", type=int, default=SLOTS)
    freeze.add_argument("--output", type=Path)
    worker = subparsers.add_parser("run-slot")
    worker.add_argument("--manifest", type=Path, required=True)
    worker.add_argument("--slot", type=int, required=True)
    worker.add_argument("--device", default="cuda")
    task_parser = subparsers.add_parser("run-task")
    task_parser.add_argument("--manifest", type=Path, required=True)
    task_parser.add_argument("--index", type=int, required=True)
    task_parser.add_argument("--device", default="cuda")
    audit = subparsers.add_parser("audit")
    audit.add_argument("--manifest", type=Path, required=True)
    jobs = subparsers.add_parser("record-pipeline")
    jobs.add_argument("entries", nargs="+")
    args = parser.parse_args()

    if args.action == "freeze-manifest":
        payload = build_manifest(args.stage, args.dimension, args.slots)
        output = args.output or manifest_path(args.stage, args.dimension)
        write_manifest(output, payload)
        print(json.dumps({
            "path": str(output),
            "manifest_sha256": payload["manifest_sha256"],
            "task_count": payload["task_count"],
            "slot_task_counts": payload["slot_task_counts"],
            "slot_weights": payload["slot_weights"],
        }, indent=2, sort_keys=True))
    elif args.action in {"run-slot", "run-task"}:
        device = torch.device(args.device)
        if device.type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA device requested but unavailable")
        manifest = load_manifest(args.manifest)
        if args.action == "run-slot":
            run_slot(manifest, args.slot, device)
        else:
            run_index(manifest, args.index, device)
    elif args.action == "audit":
        print(json.dumps(
            audit_manifest(load_manifest(args.manifest)),
            indent=2,
            sort_keys=True,
        ))
    else:
        record_pipeline(args.entries)


if __name__ == "__main__":
    main()
