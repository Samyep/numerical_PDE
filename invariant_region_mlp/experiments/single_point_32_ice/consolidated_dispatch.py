"""V3 Stage-C dispatch with one FH-MLP allocation per PDE cell."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
from typing import Any

import torch

from .amended_dispatch import (
    AMENDMENT_PATH,
    AMENDMENT_REVISION,
    PRIMARY_TA_TRAINING_SEEDS,
    _digest,
    _portable_text_sha256,
    _selected_tasks,
    audit_bundle_manifest,
    pending_indices,
    run_bundle,
)
from .orchestration import load_manifest, manifest_path
from .protocol import (
    DIMENSIONS,
    RESULTS_DIR,
    TRAINING_SEEDS,
    write_json_atomic,
)


DISPATCH_REVISION = "single_point_32_ice_dispatch_v3_fh_by_pde"
SCHEDULING_AMENDMENT_PATH = Path(__file__).with_name(
    "SCHEDULING_AMENDMENT_2026-10-23_V3.md"
)
CONSOLIDATED_QUEUE_DIR = RESULTS_DIR / "dispatch_bundles_v3"


def _bundle_key(task: dict[str, Any]) -> tuple[Any, ...]:
    common = (task["method"], task["pde"], int(task["dimension"]))
    if task["method"] in {"tangent_actor", "deep_bsde"}:
        return (*common, int(task["seed"]))
    return common


def _bundle_id(key: tuple[Any, ...]) -> str:
    method, pde, dimension, *suffix = key
    tail = f"-s{suffix[0]}" if suffix else "-all-configs-replicates"
    return f"bundle-{method}-{pde}-d{dimension}{tail}"


def _task_order(task: dict[str, Any]) -> tuple[Any, ...]:
    if task["method"] == "deep_bsde":
        return (int(task["point"]),)
    if task["method"] == "fh_mlp_standard":
        # Complete expensive configurations first so a timeout leaves only
        # cheaper restart-safe work behind.
        return (
            -float(task["weight"]),
            -int(task["n"]),
            -int(task["M"]),
            int(task["replicate"]),
        )
    return (int(task["seed"]),)


def build_consolidated_manifest(source: dict[str, Any]) -> dict[str, Any]:
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
    bundles.sort(key=lambda row: (-float(row["weight"]), row["bundle_id"]))
    for index, bundle in enumerate(bundles):
        bundle["index"] = index

    payload: dict[str, Any] = {
        "amendment_revision": AMENDMENT_REVISION,
        "dispatch_revision": DISPATCH_REVISION,
        "amendment_sha256": _portable_text_sha256(AMENDMENT_PATH),
        "scheduling_amendment_sha256": _portable_text_sha256(
            SCHEDULING_AMENDMENT_PATH
        ),
        "source_revision": source["revision"],
        "source_manifest_sha256": source["manifest_sha256"],
        "stage": "C",
        "dimension": int(source["dimension"]),
        "primary_tangent_actor_training_seeds": list(
            PRIMARY_TA_TRAINING_SEEDS
        ),
        "deep_bsde_training_seeds": list(TRAINING_SEEDS),
        "dispatch_policy": (
            "one H100 allocation per TA PDE, Deep-BSDE PDE/seed 32-point "
            "cell, or FH-MLP PDE all-configuration/all-replicate cell"
        ),
        "selected_task_count": sum(len(row["task_ids"]) for row in bundles),
        "bundle_count": len(bundles),
        "bundles": bundles,
    }
    payload["bundle_manifest_sha256"] = _digest(payload)
    return payload


def consolidated_manifest_path(dimension: int) -> Path:
    return CONSOLIDATED_QUEUE_DIR / f"stage_c_d{dimension}.json"


def write_consolidated_manifest(path: Path, payload: dict[str, Any]) -> None:
    if path.is_file():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != payload:
            raise RuntimeError(f"refusing to rewrite consolidated queue: {path}")
        return
    write_json_atomic(path, payload)


def load_consolidated_manifest(
    path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    claimed = payload.get("bundle_manifest_sha256")
    unsigned = dict(payload)
    unsigned.pop("bundle_manifest_sha256", None)
    if claimed != _digest(unsigned):
        raise RuntimeError(f"consolidated manifest hash mismatch: {path}")
    if payload.get("amendment_revision") != AMENDMENT_REVISION:
        raise RuntimeError(f"amendment revision mismatch: {path}")
    if payload.get("dispatch_revision") != DISPATCH_REVISION:
        raise RuntimeError(f"dispatch revision mismatch: {path}")
    dimension = int(payload["dimension"])
    source = load_manifest(manifest_path("C", dimension))
    if payload != build_consolidated_manifest(source):
        raise RuntimeError(f"consolidated manifest is not reproducible: {path}")
    return payload, source


def freeze_all() -> None:
    CONSOLIDATED_QUEUE_DIR.mkdir(parents=True, exist_ok=True)
    for dimension in DIMENSIONS:
        source = load_manifest(manifest_path("C", dimension))
        payload = build_consolidated_manifest(source)
        write_consolidated_manifest(consolidated_manifest_path(dimension), payload)
        print(json.dumps({
            "dimension": dimension,
            "selected_tasks": payload["selected_task_count"],
            "bundles": payload["bundle_count"],
            "path": str(consolidated_manifest_path(dimension)),
        }, sort_keys=True))


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
    args = parser.parse_args()

    if args.action == "freeze-all":
        freeze_all()
        return
    payload, source = load_consolidated_manifest(args.manifest)
    if args.action == "pending-spec":
        print(",".join(str(value) for value in pending_indices(payload, source)))
    elif args.action == "run-bundle":
        run_bundle(payload, source, args.index, torch.device(args.device))
    elif args.action == "audit":
        print(json.dumps(audit_bundle_manifest(payload, source), sort_keys=True))


if __name__ == "__main__":
    main()
