"""Command-line entry point for staged, globally truth-gated execution."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import torch

from .common import file_sha256, write_json_atomic
from .prepare import freeze
from .protocol import (
    CONFIGURATIONS,
    METHODS,
    PDES,
    RAW_DIR,
    REVISION,
    STATUS_DIR,
    ensure_directories,
    protocol_payload,
)


def _record_path(method: str, pde: str, configuration: str) -> Path:
    short = "ta" if method == "tangent_actor" else method
    return (
        RAW_DIR / "trajectories"
        / f"{short}_{pde}_d100_{configuration.lower()}_s0.json"
    )


def training_inventory() -> list[dict[str, str]]:
    inventory: list[dict[str, str]] = []
    for pde in PDES:
        for method in METHODS:
            for configuration in CONFIGURATIONS:
                path = _record_path(method, pde, configuration)
                status = "MISSING"
                if path.is_file():
                    status = str(json.loads(path.read_text(encoding="utf-8")).get("status"))
                inventory.append({
                    "pde": pde, "method": method,
                    "configuration_id": configuration,
                    "path": str(path), "status": status,
                })
    return inventory


def mark_all_training_complete() -> dict:
    ensure_directories()
    inventory = training_inventory()
    incomplete = [
        row for row in inventory
        if row["status"] not in {"TRAINED_NO_TRUTH", "COMPLETE"}
    ]
    if incomplete:
        raise RuntimeError(f"truth gate remains closed; incomplete trajectories: {incomplete}")
    protocol = protocol_payload()
    payload = {
        "revision": REVISION,
        "status": "ALL_18_TRAINING_TRAJECTORIES_STOPPED_TRUTH_GATE_OPEN",
        "opened_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": protocol["protocol_sha256"],
        "trajectory_count": len(inventory),
        "trajectories": [
            {**row, "sha256": file_sha256(Path(row["path"]))}
            for row in inventory
        ],
    }
    write_json_atomic(STATUS_DIR / "all_training_complete.json", payload)
    return payload


def _device(name: str) -> torch.device:
    if name == "auto":
        name = "cuda" if torch.cuda.is_available() else "cpu"
    result = torch.device(name)
    if result.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("freeze")
    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--method", choices=METHODS, required=True)
    train_parser.add_argument("--pde", choices=PDES, required=True)
    train_parser.add_argument("--configuration", choices=CONFIGURATIONS, required=True)
    train_parser.add_argument("--device", default="auto")
    subparsers.add_parser("mark-ready")
    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--method", choices=METHODS, required=True)
    evaluate_parser.add_argument("--pde", choices=PDES, required=True)
    evaluate_parser.add_argument("--configuration", choices=CONFIGURATIONS, required=True)
    evaluate_parser.add_argument("--device", default="auto")
    subparsers.add_parser("inventory")
    args = parser.parse_args()

    if args.command == "freeze":
        result = freeze()
    elif args.command == "inventory":
        result = training_inventory()
    elif args.command == "mark-ready":
        result = mark_all_training_complete()
    elif args.command == "train":
        device = _device(args.device)
        if args.method == "tangent_actor":
            from .ta_frontier import train
            result = train(args.pde, args.configuration, device)
        else:
            from .baseline_frontier import train
            result = train(args.method, args.pde, args.configuration, device)
    else:
        device = _device(args.device)
        if args.method == "tangent_actor":
            from .ta_frontier import evaluate
            result = evaluate(args.pde, args.configuration, device)
        else:
            from .baseline_frontier import evaluate
            result = evaluate(args.method, args.pde, args.configuration, device)
    print(json.dumps(result, indent=2, sort_keys=True, default=str), flush=True)


if __name__ == "__main__":
    main()
