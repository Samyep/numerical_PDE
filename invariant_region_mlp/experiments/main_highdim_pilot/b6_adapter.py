"""Read-only adapter for the frozen B6 benchmark and mature results.

The pilot must reuse B6 rather than train it again.  This module provides a
small, auditable interface to the frozen calibration, iteration diagnostics,
timing records, and (when the original worktree is present) selected
checkpoints.  It never writes to the mature experiment directory.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = PROJECT_ROOT / "results" / "highdim_nosym"
SOURCE_COMMIT = "281977883e0a5ce4f05895878680b582620461be"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _checkpoint_from_training(training: dict[str, Any], actor: str) -> Path | None:
    if actor == "a0":
        return None
    iteration = int(actor[1:])
    phase = next(item for item in training["phases"] if item["iteration"] == iteration)
    return Path(phase["checkpoint"])


def load_b6_reuse(dimension: int) -> dict[str, Any]:
    """Load and cross-check the frozen seed-0, 64-step B6 records."""

    dimension = int(dimension)
    if dimension not in (100, 400):
        raise ValueError("the frozen B6 study contains dimensions 100 and 400")
    required = [
        SOURCE_ROOT / "B6_CALIBRATION_FREEZE.json",
        SOURCE_ROOT / "b6_calibration.csv",
        SOURCE_ROOT / "iteration_diagnostics.csv",
        SOURCE_ROOT / "resource_accounting.csv",
        SOURCE_ROOT / "raw" / "actor" / f"training_B6_d{dimension}_mature_s0.json",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing mature B6 artifacts: {missing}")

    freeze = json.loads(required[0].read_text(encoding="utf-8"))
    selection = freeze["selections"][str(dimension)]
    calibration = selection["selected_row"]
    if calibration["bank_sha256"] != selection["test_bank_sha256"]:
        raise RuntimeError("B6 calibration/test-bank hash mismatch")

    diagnostics = [
        row
        for row in _read_csv(required[2])
        if row["pde_id"] == "B6"
        and int(row["dimension"]) == dimension
        and int(row["training_seed"]) == 0
        and int(row["steps"]) == 64
        and row["actor_request"] in {"a0", "a1", "a2", "a3", "a4"}
    ]
    if {row["actor_request"] for row in diagnostics} != {
        "a0",
        "a1",
        "a2",
        "a3",
        "a4",
    }:
        raise RuntimeError("incomplete B6 iteration diagnostics")

    evaluation_hashes: dict[str, list[str]] = {}
    for actor in ("a0", "a1", "a2", "a3", "a4"):
        local = []
        for replicate in range(4):
            path = (
                SOURCE_ROOT
                / "raw"
                / "actor_evaluations"
                / f"B6_d{dimension}_mature_s0"
                / actor
                / f"steps64_rep{replicate:02d}.json"
            )
            record = json.loads(path.read_text(encoding="utf-8"))
            if record["test_bank_sha256"] != selection["test_bank_sha256"]:
                raise RuntimeError(f"B6 {actor} evaluation bank hash mismatch")
            local.append(_sha256(path))
        evaluation_hashes[actor] = local

    training = json.loads(required[4].read_text(encoding="utf-8"))
    if training["status"] != "complete" or training["selection"]["selected_actor"] != "a4":
        raise RuntimeError("unexpected mature B6 training/selection state")
    selected_checkpoint = _checkpoint_from_training(training, "a4")
    checkpoint_audit: dict[str, Any] = {
        "recorded_path": str(selected_checkpoint),
        "available": bool(selected_checkpoint and selected_checkpoint.is_file()),
        "sha256": None,
    }
    if selected_checkpoint and selected_checkpoint.is_file():
        checkpoint_audit["sha256"] = _sha256(selected_checkpoint)

    resources = [
        row
        for row in _read_csv(required[3])
        if row["pde_id"] == "B6"
        and int(row["dimension"]) == dimension
        and row["seed_or_repetitions"] == "0"
        and row["method"] == "TA-training"
    ]
    if len(resources) != 1:
        raise RuntimeError("missing or duplicate B6 seed-0 training resource row")

    return {
        "label": "REUSED_MATURE_RESULT",
        "source_commit": SOURCE_COMMIT,
        "dimension": dimension,
        "calibration": calibration,
        "iterations": sorted(diagnostics, key=lambda row: int(row["actor_request"][1:])),
        "training": training,
        "training_resources": resources[0],
        "selected_checkpoint": checkpoint_audit,
        "evaluation_record_sha256": evaluation_hashes,
        "artifact_sha256": {str(path.relative_to(PROJECT_ROOT)): _sha256(path) for path in required},
        "audit_status": "PASS",
    }


__all__ = ["SOURCE_COMMIT", "SOURCE_ROOT", "load_b6_reuse"]
