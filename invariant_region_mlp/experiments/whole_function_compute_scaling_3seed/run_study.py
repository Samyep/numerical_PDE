"""Small CLI for calibration, formal training, and post-training evaluation."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess

import torch

from invariant_region_mlp.experiments.whole_function_compute_frontier.common import (
    file_sha256,
    write_json_atomic,
)
from invariant_region_mlp.experiments.whole_function_compute_frontier.prepare import (
    historical_results_root,
)

from .calibration import build_manifest
from .engine import evaluate_spec, train_spec, trajectory_path
from .protocol import (
    BASE_REVISION,
    CALIBRATION_MANIFEST,
    EXPECTED_INPUTS,
    PDES,
    RESULTS_DIR,
    RUN_MANIFEST,
    SOURCE_ROOT,
    STATUS_DIR,
    calibration_specs,
    ensure_directories,
    formal_specs,
)


def _device(name: str) -> torch.device:
    if name == "auto":
        name = "cuda" if torch.cuda.is_available() else "cpu"
    result = torch.device(name)
    if result.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    return result


def prepare() -> dict:
    ensure_directories()
    source = historical_results_root()
    inputs = []
    for pde in PDES:
        bank = source / "inputs" / f"{pde}_d100_test_bank.npz"
        parameters = source / "inputs" / f"{pde}_d100_parameters.npz"
        if file_sha256(bank) != EXPECTED_INPUTS[pde]["test_bank_file_sha256"]:
            raise RuntimeError(f"frozen test-bank file drift: {pde}")
        if file_sha256(parameters) != EXPECTED_INPUTS[pde]["parameter_file_sha256"]:
            raise RuntimeError(f"frozen parameter file drift: {pde}")
        inputs.append({
            "pde": pde, "test_bank": str(bank),
            "test_bank_file_sha256": file_sha256(bank),
            "parameter_file": str(parameters),
            "parameter_file_sha256": file_sha256(parameters),
        })
    payload = {
        "status": "PREPARED_WITHOUT_READING_TEST_TRUTH_ARRAY",
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "base_revision": BASE_REVISION,
        "checked_out_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=SOURCE_ROOT, text=True
        ).strip(),
        "historical_results_root": str(source),
        "inputs": inputs,
    }
    write_json_atomic(RESULTS_DIR / "preflight_manifest.json", payload)
    return payload


def _supplemental_specs():
    path = RESULTS_DIR / "calibration_supplemental_specs.json"
    payloads = json.loads(path.read_text(encoding="utf-8"))
    from .protocol import RunSpec
    result = []
    for payload in payloads:
        payload.pop("run_id", None)
        payload["checkpoint_values"] = tuple(payload["checkpoint_values"])
        result.append(RunSpec(**payload))
    return result


def training_complete() -> dict:
    specs = formal_specs()
    incomplete = []
    for spec in specs:
        path = trajectory_path(spec)
        if not path.is_file():
            incomplete.append({"run_id": spec.run_id, "reason": "missing"})
            continue
        status = json.loads(path.read_text(encoding="utf-8")).get("status")
        if status not in {"TRAINED_NO_TRUTH", "COMPLETE"}:
            incomplete.append({"run_id": spec.run_id, "reason": status})
    if incomplete:
        raise RuntimeError(f"formal training incomplete: {incomplete}")
    payload = {
        "status": "ALL_FORMAL_TRAINING_COMPLETE_TEST_EVALUATION_ALLOWED",
        "opened_utc": datetime.now(timezone.utc).isoformat(),
        "formal_run_count": len(specs),
        "runs": [
            {"run_id": spec.run_id, "trajectory_sha256": file_sha256(trajectory_path(spec))}
            for spec in specs
        ],
    }
    write_json_atomic(STATUS_DIR / "formal_training_complete.json", payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare")
    list_parser = sub.add_parser("list")
    list_parser.add_argument("--phase", choices=("calibration", "supplemental", "formal"), required=True)
    train = sub.add_parser("train")
    train.add_argument("--phase", choices=("calibration", "supplemental", "formal"), required=True)
    train.add_argument("--index", type=int, required=True)
    train.add_argument("--device", default="auto")
    sub.add_parser("build-calibration")
    sub.add_parser("mark-training-complete")
    evaluate = sub.add_parser("evaluate")
    evaluate.add_argument("--index", type=int, required=True)
    evaluate.add_argument("--device", default="auto")
    args = parser.parse_args()

    if args.command == "prepare":
        result = prepare()
    elif args.command == "build-calibration":
        result = build_manifest()
    elif args.command == "mark-training-complete":
        result = training_complete()
    elif args.command == "list":
        specs = (
            calibration_specs() if args.phase == "calibration"
            else _supplemental_specs() if args.phase == "supplemental"
            else formal_specs()
        )
        result = [spec.payload() for spec in specs]
    elif args.command == "train":
        specs = (
            calibration_specs() if args.phase == "calibration"
            else _supplemental_specs() if args.phase == "supplemental"
            else formal_specs()
        )
        result = train_spec(specs[args.index], _device(args.device))
    else:
        if not (STATUS_DIR / "formal_training_complete.json").is_file():
            raise RuntimeError("formal training has not completed")
        specs = formal_specs()
        result = evaluate_spec(specs[args.index], _device(args.device))
    print(json.dumps(result, indent=2, sort_keys=True, default=str), flush=True)


if __name__ == "__main__":
    main()
