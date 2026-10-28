"""Manifest-based H100 runner for profiling, accuracy, and actual-prefix timing."""

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

from . import fh_mlp, tangent_actor
from .protocol import (
    FH_ACCURACY_REPLICATES,
    FH_CONFIGS,
    FH_TIMING_REPLICATES,
    FORMAL_TIMING_TASK_CAP_SECONDS,
    MANIFEST_DIR,
    MCV_TIMING_TASK_CAP_SECONDS,
    PROFILE_REPLICATES,
    RESULTS_DIR,
    REVISION,
    TA_MCV_ACCURACY_REPLICATES,
    TA_PATHS,
    TA_TIMING_REPLICATES,
    TA_TRAINING_SEEDS,
    TARGETS,
    TASK_STATUS_DIR,
    require_cost_gate,
    require_frozen_selections,
    require_frozen_timing_grid,
    write_json_atomic,
)


def _manifest(stage: str, tasks: list[dict[str, Any]]) -> dict[str, Any]:
    payload = {
        "revision": REVISION,
        "stage": stage,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tasks": [{"index": index, **task} for index, task in enumerate(tasks)],
    }
    write_json_atomic(MANIFEST_DIR / f"{stage}.json", payload)
    return payload


def make_profile_manifest() -> dict[str, Any]:
    tasks = [
        {"kind": "fh_profile", "n": n, "M": M, "replicate": replicate}
        for n, M in FH_CONFIGS
        for replicate in range(PROFILE_REPLICATES)
    ]
    return _manifest("profile", tasks)


def make_calibration_manifest() -> dict[str, Any]:
    gate = require_cost_gate()
    tasks: list[dict[str, Any]] = [
        {"kind": "ta_calibration_bundle", "seed": seed}
        for seed in TA_TRAINING_SEEDS
    ]
    for key, row in gate["eligible_fh_configs"].items():
        tasks.append({"kind": "fh_calibration_bundle", "candidate": key, "n": row["n"], "M": row["M"]})
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
    ta_candidates = ["direct"] if any(
        frozen["targets"][target]["ta_direct"].get("status") == "SELECTED"
        for target in TARGETS
    ) else []
    for candidate in _selected(frozen, "ta_fk_mcv"):
        if candidate not in ta_candidates:
            ta_candidates.append(candidate)
    tasks += [
        {"kind": "ta_final_bundle", "seed": seed, "candidates": ta_candidates}
        for seed in TA_TRAINING_SEEDS
    ]
    for candidate in _selected(frozen, "fh_mlp_standard"):
        n_text, m_text = candidate.split("_")
        tasks.append({
            "kind": "fh_final_bundle",
            "candidate": candidate,
            "n": int(n_text[1:]),
            "M": int(m_text[1:]),
        })
    return _manifest("final_accuracy", tasks)


def make_timing_manifest() -> dict[str, Any]:
    require_frozen_timing_grid()
    final = json.loads(
        (RESULTS_DIR / "final_accuracy_summary.json").read_text(encoding="utf-8")
    )
    tasks: list[dict[str, Any]] = []
    coverage: dict[str, Any] = {"revision": REVISION, "targets": {}}

    direct_needed = False
    fh_candidates: list[str] = []
    mcv_candidates: list[str] = []
    for target in TARGETS:
        local = final["targets"][target]
        coverage["targets"][target] = {}
        direct_pass = local["ta_direct"]["status"] == "FINAL_ACCURACY_PASSED"
        fh_pass = local["fh_mlp_standard"]["status"] == "FINAL_ACCURACY_PASSED"
        mcv_pass = local["ta_fk_mcv"]["status"] == "FINAL_ACCURACY_PASSED"
        coverage["targets"][target]["ta_direct"] = "RUN" if direct_pass else local["ta_direct"]["status"]
        coverage["targets"][target]["fh_mlp_standard"] = "RUN" if fh_pass else local["fh_mlp_standard"]["status"]
        coverage["targets"][target]["ta_fk_mcv"] = "RUN" if mcv_pass else local["ta_fk_mcv"]["status"]
        if direct_pass:
            direct_needed = True
        if fh_pass:
            candidate = local["fh_mlp_standard"]["selected_candidate"]
            if candidate not in fh_candidates:
                fh_candidates.append(candidate)
        if mcv_pass:
            candidate = local["ta_fk_mcv"]["selected_candidate"]
            if candidate not in mcv_candidates:
                mcv_candidates.append(candidate)

    if direct_needed:
        tasks += [
            {"kind": "ta_timing", "candidate": "direct", "replicate": replicate}
            for replicate in range(TA_TIMING_REPLICATES)
        ]
    for candidate in mcv_candidates:
        tasks += [
            {"kind": "ta_timing", "candidate": candidate, "replicate": replicate}
            for replicate in range(TA_TIMING_REPLICATES)
        ]
    for candidate in fh_candidates:
        n_text, m_text = candidate.split("_")
        tasks += [
            {
                "kind": "fh_timing",
                "candidate": candidate,
                "n": int(n_text[1:]),
                "M": int(m_text[1:]),
                "replicate": replicate,
            }
            for replicate in range(FH_TIMING_REPLICATES)
        ]
    write_json_atomic(RESULTS_DIR / "timing_coverage.json", coverage)
    return _manifest("timing", tasks)


def _run(task: dict[str, Any], device: torch.device) -> list[dict[str, Any]]:
    kind = task["kind"]
    rows: list[dict[str, Any]] = []
    if kind == "fh_profile":
        rows.append(fh_mlp.profile(task["n"], task["M"], task["replicate"], device))
    elif kind == "ta_calibration_bundle":
        seed = int(task["seed"])
        rows.append(tangent_actor.evaluate("calibration", seed, "direct", 0, device))
        for paths in TA_PATHS:
            for replicate in range(TA_MCV_ACCURACY_REPLICATES):
                rows.append(tangent_actor.evaluate(
                    "calibration", seed, f"mcv_{paths}", replicate, device
                ))
    elif kind == "fh_calibration_bundle":
        for replicate in range(FH_ACCURACY_REPLICATES):
            rows.append(fh_mlp.evaluate(
                "calibration", task["n"], task["M"], replicate, device
            ))
    elif kind == "ta_final_bundle":
        seed = int(task["seed"])
        for candidate in task["candidates"]:
            replicates = (0,) if candidate == "direct" else range(TA_MCV_ACCURACY_REPLICATES)
            for replicate in replicates:
                rows.append(tangent_actor.evaluate(
                    "final", seed, candidate, replicate, device
                ))
    elif kind == "fh_final_bundle":
        for replicate in range(FH_ACCURACY_REPLICATES):
            rows.append(fh_mlp.evaluate(
                "final", task["n"], task["M"], replicate, device
            ))
    elif kind == "ta_timing":
        cap = (
            FORMAL_TIMING_TASK_CAP_SECONDS
            if task["candidate"] == "direct"
            else MCV_TIMING_TASK_CAP_SECONDS
        )
        rows.append(tangent_actor.timing(
            task["replicate"], task["candidate"], device,
            max_active_seconds=cap,
        ))
    elif kind == "fh_timing":
        rows.append(fh_mlp.timing(
            task["n"], task["M"], task["replicate"], device,
            max_active_seconds=FORMAL_TIMING_TASK_CAP_SECONDS,
        ))
    else:
        raise KeyError(kind)
    return rows


def run_task(manifest_path: Path, index: int, device: torch.device) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["revision"] != REVISION:
        raise RuntimeError("task manifest revision mismatch")
    task = manifest["tasks"][index]
    task_id = f"{manifest['stage']}__{index:04d}__{task['kind']}"
    attempt = f"{os.environ.get('SLURM_JOB_ID', 'local')}_{os.environ.get('SLURM_ARRAY_TASK_ID', index)}"
    status_path = TASK_STATUS_DIR / f"{task_id}__{attempt}.json"
    started = time.perf_counter()
    try:
        rows = _run(task, device)
        payload = {
            "revision": REVISION,
            "task_id": task_id,
            "stage": manifest["stage"],
            "task": task,
            "status": "complete",
            "record_ids": [row.get("run_id") for row in rows],
            "task_wall_seconds": time.perf_counter() - started,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "gpu_count": 1,
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
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "gpu_count": 1,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_json_atomic(status_path, payload)
        raise
    write_json_atomic(status_path, payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    for name in ("make-profile", "make-calibration", "make-final", "make-timing"):
        sub.add_parser(name)
    item = sub.add_parser("run-task")
    item.add_argument("--manifest", type=Path, required=True)
    item.add_argument("--index", type=int, required=True)
    item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.action == "make-profile":
        result = make_profile_manifest()
    elif args.action == "make-calibration":
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
