"""Seed-0 Tangent-Actor policy and cell-bundled Stage-C dispatch.

The v1 scientific task manifests remain immutable.  This module selects the
amended primary task set and groups those tasks into natural GPU allocations
without changing any solver call, seed, or task-level output.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import torch

from .orchestration import (
    DIMENSION_AUDIT_DIR,
    _audit_task,
    algorithm_hashes,
    execute_task,
    load_manifest,
    manifest_path,
    task_complete,
    utc_now,
)
from .protocol import (
    DIMENSIONS,
    FH_CONFIGURATIONS,
    FH_REPLICATES,
    PDES,
    QUERY_COUNT,
    RESULTS_DIR,
    TRAINING_SEEDS,
    write_json_atomic,
)


AMENDMENT_REVISION = "single_point_32_ice_v2_ta_seed0_bundled"
PRIMARY_TA_TRAINING_SEEDS = (0,)
AMENDMENT_PATH = Path(__file__).with_name("PROTOCOL_AMENDMENT_2026-10-23.md")
BUNDLE_QUEUE_DIR = RESULTS_DIR / "dispatch_bundles"
AMENDED_PIPELINE_JOBS = RESULTS_DIR / "pipeline_jobs_amended.json"


def _digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def _portable_text_sha256(path: Path) -> str:
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
    return sha256(text.encode("utf-8")).hexdigest()


def _selected_tasks(source: dict[str, Any]) -> list[dict[str, Any]]:
    if source.get("stage") != "C":
        raise ValueError("amended dispatch only accepts Stage-C manifests")
    return [
        task
        for task in source["tasks"]
        if task["method"] != "tangent_actor"
        or int(task["seed"]) in PRIMARY_TA_TRAINING_SEEDS
    ]


def _bundle_key(task: dict[str, Any]) -> tuple[Any, ...]:
    common = (task["method"], task["pde"], int(task["dimension"]))
    if task["method"] == "tangent_actor":
        return (*common, int(task["seed"]))
    if task["method"] == "deep_bsde":
        return (*common, int(task["seed"]))
    return (*common, int(task["n"]), int(task["M"]))


def _bundle_id(key: tuple[Any, ...]) -> str:
    method, pde, dimension, *suffix = key
    if method == "tangent_actor":
        tail = f"s{suffix[0]}"
    elif method == "deep_bsde":
        tail = f"s{suffix[0]}"
    else:
        tail = f"n{suffix[0]}-M{suffix[1]}"
    return f"bundle-{method}-{pde}-d{dimension}-{tail}"


def _task_order(task: dict[str, Any]) -> tuple[int, ...]:
    if task["method"] == "deep_bsde":
        return (int(task["point"]),)
    if task["method"] == "fh_mlp_standard":
        return (int(task["replicate"]),)
    return (int(task["seed"]),)


def build_bundle_manifest(source: dict[str, Any]) -> dict[str, Any]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for task in _selected_tasks(source):
        grouped[_bundle_key(task)].append(task)

    bundles: list[dict[str, Any]] = []
    for key, tasks in grouped.items():
        ordered = sorted(tasks, key=_task_order)
        bundles.append({
            "bundle_id": _bundle_id(key),
            "method": ordered[0]["method"],
            "pde": ordered[0]["pde"],
            "dimension": int(ordered[0]["dimension"]),
            "weight": round(sum(float(task["weight"]) for task in ordered), 9),
            "task_ids": [task["task_id"] for task in ordered],
        })

    # Descending weight gives Slurm's normally ascending array dispatch an LPT
    # order while leaving every scientific task and seed untouched.
    bundles.sort(key=lambda row: (-float(row["weight"]), row["bundle_id"]))
    for index, bundle in enumerate(bundles):
        bundle["index"] = index

    payload: dict[str, Any] = {
        "amendment_revision": AMENDMENT_REVISION,
        "amendment_sha256": _portable_text_sha256(AMENDMENT_PATH),
        "source_revision": source["revision"],
        "source_manifest_sha256": source["manifest_sha256"],
        "stage": "C",
        "dimension": int(source["dimension"]),
        "primary_tangent_actor_training_seeds": list(
            PRIMARY_TA_TRAINING_SEEDS
        ),
        "deep_bsde_training_seeds": list(TRAINING_SEEDS),
        "dispatch_policy": (
            "one H100 allocation per TA training, Deep-BSDE PDE/seed 32-point "
            "cell, or FH-MLP PDE/configuration 10-replicate cell"
        ),
        "selected_task_count": sum(len(row["task_ids"]) for row in bundles),
        "bundle_count": len(bundles),
        "bundles": bundles,
    }
    payload["bundle_manifest_sha256"] = _digest(payload)
    return payload


def bundle_manifest_path(dimension: int) -> Path:
    return BUNDLE_QUEUE_DIR / f"stage_c_d{dimension}.json"


def write_bundle_manifest(path: Path, payload: dict[str, Any]) -> None:
    if path.is_file():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != payload:
            raise RuntimeError(f"refusing to rewrite amended queue: {path}")
        return
    write_json_atomic(path, payload)


def load_bundle_manifest(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    claimed = payload.get("bundle_manifest_sha256")
    unsigned = dict(payload)
    unsigned.pop("bundle_manifest_sha256", None)
    if claimed != _digest(unsigned):
        raise RuntimeError(f"bundle manifest hash mismatch: {path}")
    if payload.get("amendment_revision") != AMENDMENT_REVISION:
        raise RuntimeError(f"amendment revision mismatch: {path}")
    dimension = int(payload["dimension"])
    source = load_manifest(manifest_path("C", dimension))
    rebuilt = build_bundle_manifest(source)
    if payload != rebuilt:
        raise RuntimeError(f"bundle manifest is not reproducible: {path}")
    return payload, source


def _task_lookup(source: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {task["task_id"]: task for task in source["tasks"]}


def bundle_complete(
    bundle: dict[str, Any], source: dict[str, Any]
) -> bool:
    lookup = _task_lookup(source)
    return all(task_complete(lookup[task_id]) for task_id in bundle["task_ids"])


def pending_indices(
    payload: dict[str, Any], source: dict[str, Any]
) -> list[int]:
    return [
        int(bundle["index"])
        for bundle in payload["bundles"]
        if not bundle_complete(bundle, source)
    ]


def run_bundle(
    payload: dict[str, Any],
    source: dict[str, Any],
    index: int,
    device: torch.device,
) -> None:
    bundles = payload["bundles"]
    if not 0 <= index < len(bundles):
        raise ValueError(index)
    bundle = bundles[index]
    lookup = _task_lookup(source)
    print(json.dumps({
        "event": "bundle_start",
        "bundle_index": index,
        "bundle_id": bundle["bundle_id"],
        "tasks": len(bundle["task_ids"]),
        "utc": utc_now(),
    }, sort_keys=True), flush=True)
    for position, task_id in enumerate(bundle["task_ids"]):
        task = lookup[task_id]
        print(json.dumps({
            "event": "bundle_task_start",
            "bundle_id": bundle["bundle_id"],
            "position": position,
            "task_id": task_id,
            "utc": utc_now(),
        }, sort_keys=True), flush=True)
        result = execute_task(task, device)
        print(json.dumps({
            "event": "bundle_task_end",
            "bundle_id": bundle["bundle_id"],
            "position": position,
            "task_id": task_id,
            "status": result["status"],
            "utc": utc_now(),
        }, sort_keys=True), flush=True)
    print(json.dumps({
        "event": "bundle_end",
        "bundle_index": index,
        "bundle_id": bundle["bundle_id"],
        "utc": utc_now(),
    }, sort_keys=True), flush=True)


def audit_bundle_manifest(
    payload: dict[str, Any], source: dict[str, Any]
) -> dict[str, Any]:
    lookup = _task_lookup(source)
    selected = [
        lookup[task_id]
        for bundle in payload["bundles"]
        for task_id in bundle["task_ids"]
    ]
    audits = [_audit_task(task) for task in selected]
    missing = [row["task_id"] for row in audits if not row["terminal"]]
    trained_nonfinite = [
        audit["task_id"]
        for task, audit in zip(selected, audits, strict=True)
        if task["method"] != "fh_mlp_standard"
        and any(
            record["status"] != "complete"
            or record["nonfinite_count"] not in {0, None}
            for record in audit["records"]
        )
    ]
    expected_counts = {
        "tangent_actor": len(PDES) * len(PRIMARY_TA_TRAINING_SEEDS),
        "deep_bsde": len(PDES) * QUERY_COUNT * len(TRAINING_SEEDS),
        "fh_mlp_standard": (
            len(PDES) * len(FH_CONFIGURATIONS) * FH_REPLICATES
        ),
    }
    observed_counts = {
        method: sum(
            audit["terminal"] and audit["method"] == method
            for audit in audits
        )
        for method in expected_counts
    }
    checks = {
        "all_selected_tasks_terminal": {
            "passed": not missing,
            "missing": missing,
        },
        "selected_task_count": {
            "passed": len(selected) == int(payload["selected_task_count"]),
            "expected": int(payload["selected_task_count"]),
            "observed": len(selected),
        },
        "method_terminal_counts": {
            "passed": observed_counts == expected_counts,
            "expected": expected_counts,
            "observed": observed_counts,
        },
        "trained_methods_finite": {
            "passed": not trained_nonfinite,
            "failed": trained_nonfinite,
        },
    }
    overall = all(bool(check["passed"]) for check in checks.values())
    result = {
        "revision": AMENDMENT_REVISION,
        "source_revision": source["revision"],
        "stage": "C",
        "dimension": int(payload["dimension"]),
        "source_manifest_sha256": source["manifest_sha256"],
        "bundle_manifest_sha256": payload["bundle_manifest_sha256"],
        "algorithm_hashes": algorithm_hashes(),
        "checked_utc": utc_now(),
        "checks": checks,
        "method_terminal_counts": observed_counts,
        "fh_numerical_failure_records": sum(
            record["status"] == "numerical_failure"
            for audit in audits
            if audit["method"] == "fh_mlp_standard"
            for record in audit["records"]
        ),
        "overall_pass": overall,
        "status": "VALIDATED" if overall else "FAILED",
    }
    output = DIMENSION_AUDIT_DIR / f"d{payload['dimension']}.json"
    write_json_atomic(output, result)
    if not overall:
        raise RuntimeError(f"amended Stage-C audit failed: {output}")
    return result


def freeze_all() -> None:
    BUNDLE_QUEUE_DIR.mkdir(parents=True, exist_ok=True)
    for dimension in DIMENSIONS:
        source = load_manifest(manifest_path("C", dimension))
        payload = build_bundle_manifest(source)
        write_bundle_manifest(bundle_manifest_path(dimension), payload)
        print(json.dumps({
            "dimension": dimension,
            "selected_tasks": payload["selected_task_count"],
            "bundles": payload["bundle_count"],
            "path": str(bundle_manifest_path(dimension)),
        }, sort_keys=True))


def record_pipeline(job_id: str, dimension: int) -> None:
    write_json_atomic(AMENDED_PIPELINE_JOBS, {
        "revision": AMENDMENT_REVISION,
        "recorded_utc": utc_now(),
        "initial_dimension": dimension,
        "initial_launcher_job_id": job_id,
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("freeze-all")
    pending = subparsers.add_parser("pending-spec")
    pending.add_argument("--manifest", type=Path, required=True)
    run = subparsers.add_parser("run-bundle")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--index", type=int, required=True)
    run.add_argument("--device", default="cuda")
    audit = subparsers.add_parser("audit")
    audit.add_argument("--manifest", type=Path, required=True)
    record = subparsers.add_parser("record-pipeline")
    record.add_argument("--job-id", required=True)
    record.add_argument("--dimension", type=int, required=True)
    args = parser.parse_args()

    if args.action == "freeze-all":
        freeze_all()
        return
    if args.action == "record-pipeline":
        record_pipeline(args.job_id, args.dimension)
        return
    payload, source = load_bundle_manifest(args.manifest)
    if args.action == "pending-spec":
        print(",".join(str(value) for value in pending_indices(payload, source)))
    elif args.action == "run-bundle":
        run_bundle(payload, source, args.index, torch.device(args.device))
    elif args.action == "audit":
        print(json.dumps(audit_bundle_manifest(payload, source), sort_keys=True))


if __name__ == "__main__":
    main()
