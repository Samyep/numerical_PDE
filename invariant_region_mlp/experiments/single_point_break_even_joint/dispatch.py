"""Manifest-based H100 runner for calibration, final accuracy, and timing."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
import traceback
from typing import Any

import torch

from . import deep_bsde, fh_mlp, tangent_actor
from .protocol import (
    DB_CANDIDATE_UPDATES,
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
        "revision": REVISION, "stage": stage,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tasks": [{"index": index, **task} for index, task in enumerate(tasks)],
    }
    write_json_atomic(MANIFEST_DIR / f"{stage}.json", payload)
    return payload


def make_calibration_manifest() -> dict[str, Any]:
    tasks: list[dict[str, Any]] = []
    tasks += [{"kind": "ta_calibration_bundle", "seed": seed} for seed in (0, 1, 2)]
    tasks += [{"kind": "fh_calibration_bundle", "n": n, "M": M} for n, M in FH_CONFIGS]
    tasks += [{"kind": "db_calibration_bundle", "point": point} for point in range(32)]
    return _manifest("calibration", tasks)


def _selected(frozen: dict[str, Any], method: str) -> list[str]:
    values: list[str] = []
    for target in TARGETS:
        item = frozen["targets"][target][method]
        candidate = item.get("selected_candidate")
        if item.get("status") == "SELECTED" and candidate and candidate not in values:
            values.append(candidate)
    return values


def make_final_manifest() -> dict[str, Any]:
    frozen = require_frozen_selections()
    tasks: list[dict[str, Any]] = []
    ta_candidates = _selected(frozen, "tangent_actor")
    tasks += [{"kind": "ta_final_bundle", "seed": seed, "candidates": ta_candidates} for seed in FINAL_SEEDS]
    for candidate in _selected(frozen, "fh_mlp_standard"):
        n_text, m_text = candidate.split("_")
        tasks.append({"kind": "fh_final_bundle", "n": int(n_text[1:]), "M": int(m_text[1:])})
    db_updates = sorted({int(value.split("_")[1]) for value in _selected(frozen, "deep_bsde")})
    if db_updates:
        tasks += [{"kind": "db_final_bundle", "point": point, "capture_updates": db_updates, "max_updates": max(db_updates)} for point in range(64)]
    return _manifest("final_accuracy", tasks)


def make_timing_manifest() -> dict[str, Any]:
    path = RESULTS_DIR / "final_accuracy_summary.json"
    if not path.is_file():
        raise RuntimeError("final accuracy must be aggregated before timing")
    final = json.loads(path.read_text(encoding="utf-8"))
    tasks: list[dict[str, Any]] = []
    coverage: dict[str, Any] = {"revision": REVISION, "policy": "primary target first; secondary runs only when it reuses the primary configuration", "targets": {}}
    primary = final["targets"]["epsilon_1"]
    primary_candidates: dict[str, str] = {}
    for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
        item = primary[method]
        if item["status"] == "FINAL_ACCURACY_PASSED":
            primary_candidates[method] = item["selected_candidate"]
    for target in TARGETS:
        coverage["targets"][target] = {}
        for method in ("tangent_actor", "fh_mlp_standard", "deep_bsde"):
            item = final["targets"][target][method]
            if item["status"] != "FINAL_ACCURACY_PASSED":
                coverage["targets"][target][method] = item["status"]
            elif primary_candidates.get(method) == item.get("selected_candidate"):
                coverage["targets"][target][method] = "RUN"
            elif target == "epsilon_1":
                coverage["targets"][target][method] = "RUN"
            else:
                coverage["targets"][target][method] = "NOT_EVALUATED_COMPUTE_CAP"
    if "tangent_actor" in primary_candidates:
        tasks += [{"kind": "ta_timing", "seed": seed, "candidate": primary_candidates["tangent_actor"]} for seed in FINAL_SEEDS]
    if "fh_mlp_standard" in primary_candidates:
        candidate = primary_candidates["fh_mlp_standard"]
        n_text, m_text = candidate.split("_")
        tasks += [{"kind": "fh_timing", "n": int(n_text[1:]), "M": int(m_text[1:]), "replicate": replicate} for replicate in range(3)]
    if "deep_bsde" in primary_candidates:
        updates = int(primary_candidates["deep_bsde"].split("_")[1])
        tasks += [{"kind": "db_timing", "updates": updates, "replicate": replicate} for replicate in range(3)]
    write_json_atomic(RESULTS_DIR / "timing_coverage.json", coverage)
    return _manifest("timing", tasks)


def _run(task: dict[str, Any], device: torch.device) -> list[dict[str, Any]]:
    kind = task["kind"]
    rows = []
    if kind == "ta_calibration_bundle":
        seed = int(task["seed"]); rows.append(tangent_actor.train(seed, device)); rows.append(tangent_actor.evaluate("calibration", seed, "direct", 0, device))
        for paths in TA_PATHS:
            for replicate in range(2): rows.append(tangent_actor.evaluate("calibration", seed, f"mcv_{paths}", replicate, device))
    elif kind == "fh_calibration_bundle":
        for replicate in range(5): rows.append(fh_mlp.evaluate("calibration", int(task["n"]), int(task["M"]), replicate, device))
    elif kind == "db_calibration_bundle":
        for seed in (0, 1, 2): rows.append(deep_bsde.evaluate("calibration", int(task["point"]), seed, device, max_updates=4000, capture_updates=DB_CANDIDATE_UPDATES))
    elif kind == "ta_final_bundle":
        seed = int(task["seed"]); rows.append(tangent_actor.train(seed, device))
        for candidate in task["candidates"]:
            for replicate in ((0,) if candidate == "direct" else (0, 1)): rows.append(tangent_actor.evaluate("final", seed, candidate, replicate, device))
    elif kind == "fh_final_bundle":
        for replicate in range(5): rows.append(fh_mlp.evaluate("final", int(task["n"]), int(task["M"]), replicate, device))
    elif kind == "db_final_bundle":
        capture = tuple(int(value) for value in task["capture_updates"])
        for seed in FINAL_SEEDS: rows.append(deep_bsde.evaluate("final", int(task["point"]), seed, device, max_updates=int(task["max_updates"]), capture_updates=capture))
    elif kind == "ta_timing": rows.append(tangent_actor.timing(int(task["seed"]), task["candidate"], device))
    elif kind == "fh_timing": rows.append(fh_mlp.timing(int(task["n"]), int(task["M"]), int(task["replicate"]), device))
    elif kind == "db_timing": rows.append(deep_bsde.timing(int(task["updates"]), int(task["replicate"]), device))
    else: raise KeyError(kind)
    return rows


def run_task(manifest_path: Path, index: int, device: torch.device) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["revision"] != REVISION: raise RuntimeError("task manifest revision mismatch")
    task = manifest["tasks"][index]
    task_id = f"{manifest['stage']}__{index:04d}__{task['kind']}"
    attempt = f"{os.environ.get('SLURM_JOB_ID', 'local')}_{os.environ.get('SLURM_ARRAY_TASK_ID', index)}"
    status_path = TASK_STATUS_DIR / f"{task_id}__{attempt}.json"
    started = time.perf_counter()
    try:
        rows = _run(task, device)
        payload = {"revision": REVISION, "task_id": task_id, "stage": manifest["stage"], "task": task, "status": "complete", "record_ids": [row.get("run_id", row.get("training_run_id")) for row in rows], "task_wall_seconds": time.perf_counter() - started, "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"), "gpu_count": 1, "finished_utc": datetime.now(timezone.utc).isoformat()}
    except Exception as error:
        payload = {"revision": REVISION, "task_id": task_id, "stage": manifest["stage"], "task": task, "status": "failed", "error_type": type(error).__name__, "error": str(error), "traceback": traceback.format_exc(), "task_wall_seconds": time.perf_counter() - started, "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"), "gpu_count": 1, "finished_utc": datetime.now(timezone.utc).isoformat()}
        write_json_atomic(status_path, payload); raise
    write_json_atomic(status_path, payload); return payload


def main() -> None:
    parser = argparse.ArgumentParser(); sub = parser.add_subparsers(dest="action", required=True)
    for name in ("make-calibration", "make-final", "make-timing"): sub.add_parser(name)
    p = sub.add_parser("run-task"); p.add_argument("--manifest", type=Path, required=True); p.add_argument("--index", type=int, required=True); p.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.action == "make-calibration": result = make_calibration_manifest()
    elif args.action == "make-final": result = make_final_manifest()
    elif args.action == "make-timing": result = make_timing_manifest()
    else: result = run_task(args.manifest, args.index, torch.device(args.device))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
