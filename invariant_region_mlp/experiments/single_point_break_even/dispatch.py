"""Small manifest-based runner for H100 calibration, final, and timing stages."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import traceback
from typing import Any

import torch

from . import deep_bsde, fh_mlp, tangent_actor
from .protocol import (
    DB_TIMING_COUNTS,
    FH_CONFIGS,
    FINAL_SEEDS,
    RAW_DIR,
    RESULTS_DIR,
    REVISION,
    TA_PATHS,
    TARGETS,
    require_frozen_selections,
    write_json_atomic,
)


MANIFEST_DIR = RESULTS_DIR / "manifests"
TASK_STATUS_DIR = RAW_DIR / "task_status"


def _manifest(stage: str, tasks: list[dict[str, Any]]) -> dict[str, Any]:
    payload = {
        "revision": REVISION,
        "stage": stage,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tasks": [{"index": index, **task} for index, task in enumerate(tasks)],
    }
    write_json_atomic(MANIFEST_DIR / f"{stage}.json", payload)
    return payload


def make_calibration_manifest() -> dict[str, Any]:
    tasks: list[dict[str, Any]] = []
    tasks += [{"kind": "ta_calibration_bundle", "seed": seed} for seed in (0, 1, 2)]
    tasks += [
        {"kind": "fh_calibration_bundle", "n": n, "M": M}
        for n, M in FH_CONFIGS
    ]
    tasks += [{"kind": "db_calibration_bundle", "point": point} for point in range(32)]
    return _manifest("calibration", tasks)


def _selected_candidates(frozen: dict[str, Any], method: str) -> list[str]:
    values = []
    for target in TARGETS:
        item = frozen["targets"][target][method]
        candidate = item.get("selected_candidate")
        if item.get("status") == "SELECTED" and candidate and candidate not in values:
            values.append(candidate)
    return values


def make_final_manifest() -> dict[str, Any]:
    frozen = require_frozen_selections()
    tasks: list[dict[str, Any]] = []
    ta_candidates = _selected_candidates(frozen, "tangent_actor")
    for seed in FINAL_SEEDS:
        tasks.append({
            "kind": "ta_final_bundle", "seed": seed, "candidates": ta_candidates
        })
    for candidate in _selected_candidates(frozen, "fh_mlp_standard"):
        n_text, M_text = candidate.split("_")
        tasks.append({
            "kind": "fh_final_bundle",
            "n": int(n_text[1:]),
            "M": int(M_text[1:]),
        })
    db_updates = sorted({
        int(candidate.split("_")[1])
        for candidate in _selected_candidates(frozen, "deep_bsde")
    })
    if db_updates:
        for point in range(64):
            tasks.append({
                "kind": "db_final_bundle",
                "point": point,
                "capture_updates": db_updates,
                "max_updates": max(db_updates),
            })
    return _manifest("final_accuracy", tasks)


def make_timing_manifest() -> dict[str, Any]:
    summary_path = RESULTS_DIR / "final_accuracy_summary.json"
    if not summary_path.is_file():
        raise RuntimeError("final accuracy must be aggregated before timing")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    tasks: list[dict[str, Any]] = []
    eligible: dict[str, set[str]] = {
        "tangent_actor": set(), "fh_mlp_standard": set(), "deep_bsde": set()
    }
    for target in TARGETS:
        for method in eligible:
            item = summary["targets"][target][method]
            if item["status"] == "FINAL_ACCURACY_PASSED":
                eligible[method].add(item["selected_candidate"])
    for candidate in sorted(eligible["tangent_actor"]):
        for seed in FINAL_SEEDS:
            tasks.append({"kind": "ta_timing", "seed": seed, "candidate": candidate})
    for candidate in sorted(eligible["fh_mlp_standard"]):
        n_text, M_text = candidate.split("_")
        for replicate in range(3):
            tasks.append({
                "kind": "fh_timing",
                "n": int(n_text[1:]),
                "M": int(M_text[1:]),
                "replicate": replicate,
            })
    for candidate in sorted(eligible["deep_bsde"]):
        updates = int(candidate.split("_")[1])
        for replicate in range(3):
            tasks.append({
                "kind": "db_timing", "updates": updates, "replicate": replicate
            })
    return _manifest("timing", tasks)


def _run(task: dict[str, Any], device: torch.device) -> list[dict[str, Any]]:
    kind = task["kind"]
    records = []
    if kind == "ta_calibration_bundle":
        seed = int(task["seed"])
        records.append(tangent_actor.train(seed, device))
        records.append(tangent_actor.evaluate("calibration", seed, "direct", 0, device))
        for paths in TA_PATHS:
            for replicate in range(2):
                records.append(tangent_actor.evaluate(
                    "calibration", seed, f"mcv_{paths}", replicate, device
                ))
    elif kind == "fh_calibration_bundle":
        for replicate in range(5):
            records.append(fh_mlp.evaluate(
                "calibration", int(task["n"]), int(task["M"]), replicate, device
            ))
    elif kind == "db_calibration_bundle":
        for seed in (0, 1, 2):
            records.append(deep_bsde.evaluate(
                "calibration",
                int(task["point"]),
                seed,
                device,
                max_updates=8000,
                capture_updates=(2000, 4000, 8000),
            ))
    elif kind == "ta_final_bundle":
        seed = int(task["seed"])
        records.append(tangent_actor.train(seed, device))
        for candidate in task["candidates"]:
            replicates = (0,) if candidate == "direct" else (0, 1)
            for replicate in replicates:
                records.append(tangent_actor.evaluate(
                    "final", seed, candidate, replicate, device
                ))
    elif kind == "fh_final_bundle":
        for replicate in range(5):
            records.append(fh_mlp.evaluate(
                "final", int(task["n"]), int(task["M"]), replicate, device
            ))
    elif kind == "db_final_bundle":
        capture = tuple(int(value) for value in task["capture_updates"])
        for seed in FINAL_SEEDS:
            records.append(deep_bsde.evaluate(
                "final",
                int(task["point"]),
                seed,
                device,
                max_updates=int(task["max_updates"]),
                capture_updates=capture,
            ))
    elif kind == "ta_timing":
        records.append(tangent_actor.timing(
            int(task["seed"]), task["candidate"], device
        ))
    elif kind == "fh_timing":
        records.append(fh_mlp.timing(
            int(task["n"]), int(task["M"]), int(task["replicate"]), device
        ))
    elif kind == "db_timing":
        records.append(deep_bsde.timing(
            int(task["updates"]), int(task["replicate"]), device
        ))
    else:
        raise KeyError(kind)
    return records


def run_task(manifest_path: Path, index: int, device: torch.device) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["revision"] != REVISION:
        raise RuntimeError("task manifest revision mismatch")
    task = manifest["tasks"][index]
    if int(task["index"]) != index:
        raise RuntimeError("task manifest index mismatch")
    task_id = f"{manifest['stage']}__{index:04d}__{task['kind']}"
    status_path = TASK_STATUS_DIR / f"{task_id}.json"
    started = time.perf_counter()
    payload: dict[str, Any]
    try:
        records = _run(task, device)
        payload = {
            "revision": REVISION,
            "task_id": task_id,
            "stage": manifest["stage"],
            "task": task,
            "status": "complete",
            "record_ids": [row.get("run_id", row.get("training_run_id")) for row in records],
            "task_wall_seconds": time.perf_counter() - started,
            "slurm_job_id": __import__("os").environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": __import__("os").environ.get("SLURM_ARRAY_TASK_ID"),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
    except Exception as error:
        payload = {
            "revision": REVISION,
            "task_id": task_id,
            "stage": manifest["stage"],
            "task": task,
            "status": "failed",
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
            "task_wall_seconds": time.perf_counter() - started,
            "slurm_job_id": __import__("os").environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": __import__("os").environ.get("SLURM_ARRAY_TASK_ID"),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json_atomic(status_path, payload)
        raise
    write_json_atomic(status_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    for name in ("make-calibration", "make-final", "make-timing"):
        subparsers.add_parser(name)
    run = subparsers.add_parser("run-task")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--index", type=int, required=True)
    run.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.action == "make-calibration":
        result = make_calibration_manifest()
    elif args.action == "make-final":
        result = make_final_manifest()
    elif args.action == "make-timing":
        result = make_timing_manifest()
    else:
        result = run_task(args.manifest, args.index, torch.device(args.device))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

