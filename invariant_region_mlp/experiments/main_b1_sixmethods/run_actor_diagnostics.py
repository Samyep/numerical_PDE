"""Run the paired oracle-actor and pathwise-baseline diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import traceback
from typing import Any

import numpy as np
import torch

from .actor_diagnostics import (
    DIAGNOSTIC_METHODS,
    PAIRED_SEED_NAMESPACE,
    evaluate_actor_diagnostic,
)
from .common import RESULTS_DIR, array_sha256, atomic_json, fixed_test_points, load_config
from .equation import B1Equation
from .metrics import prediction_metrics


DIAGNOSTIC_RESULTS_DIR = RESULTS_DIR / "actor_bottleneck_diagnostics"
RAW_DIR = DIAGNOSTIC_RESULTS_DIR / "raw"


def _paths(method: str, d: int, k: int) -> tuple[Path, Path]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"{method}_d{d}_k{k}"
    return RAW_DIR / f"{stem}.json", RAW_DIR / f"{stem}_predictions.npz"


def _save_predictions(path: Path, predictions: list[dict[str, Any]]) -> None:
    payload: dict[str, np.ndarray] = {}
    for item in predictions:
        repetition = int(item["repetition"])
        payload[f"value_seed{repetition}"] = item["value"]
        payload[f"z_seed{repetition}"] = item["z"]
    np.savez_compressed(path, **payload)


def _learned_actor_record(d: int, k: int, test_hash: str) -> tuple[Path, str]:
    path = RESULTS_DIR / "raw" / "stage_b" / f"actor_mlp_d{d}_k{k}.json"
    if not path.is_file():
        raise FileNotFoundError(f"missing completed learned Actor--MLP record: {path}")
    raw = path.read_bytes()
    payload = json.loads(raw)
    if payload.get("status") != "complete" or len(payload.get("rows", [])) != 10:
        raise RuntimeError(f"learned Actor--MLP record is incomplete: {path}")
    hashes = {row["test_point_sha256"] for row in payload["rows"]}
    if hashes != {test_hash}:
        raise RuntimeError(
            f"held-out point hash differs from learned Actor--MLP for d={d}, k={k}"
        )
    return path, hashlib.sha256(raw).hexdigest()


def run_diagnostic(
    method: str,
    d: int,
    k: int,
    *,
    force: bool = False,
) -> dict[str, Any]:
    if method not in DIAGNOSTIC_METHODS:
        raise ValueError(method)
    output_path, prediction_path = _paths(method, d, k)
    if output_path.is_file() and not force:
        return json.loads(output_path.read_text())

    config = load_config()
    stage = config["stage_b"]
    equation = B1Equation(
        d=d,
        k=k,
        coefficient_scale=float(config["pde"]["frozen_coefficient_scale"]),
    )
    points = fixed_test_points(d, k, int(stage["test_points"]))
    test_hash = array_sha256(points)
    learned_path, learned_sha256 = _learned_actor_record(d, k, test_hash)
    t = np.zeros(len(points), dtype=np.float64)
    reference_u = equation.exact_u(t, points)
    reference_z = equation.exact_z(t, points)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("actor diagnostics require an allocated GPU")

    rows: list[dict[str, Any]] = []
    predictions: list[dict[str, Any]] = []
    started = time.perf_counter()
    for repetition in range(int(stage["mlp_seeds"])):
        state, cost = evaluate_actor_diagnostic(
            equation,
            t,
            points,
            method=method,
            n=int(stage["mlp_n"]),
            M=int(stage["mlp_M"]),
            repetition=repetition,
            device=device,
        )
        row: dict[str, Any] = {
            "stage": "actor_bottleneck_diagnostic",
            "method": method,
            "d": d,
            "k": k,
            "K": equation.K,
            "coefficient_scale": None if k == 1 else equation.coefficient_scale,
            "repetition": repetition,
            "test_points": len(points),
            "test_point_sha256": test_hash,
            "test_time": 0.0,
            "test_distribution": "x~N(0,I_d)",
            "training_seconds": 0.0,
            "inference_seconds": cost["wall_seconds"],
            "parameter_count": 0,
        }
        row.update(
            prediction_metrics(
                state[:, 0], state[:, 1:], reference_u, reference_z
            )
        )
        row.update({f"cost_{key}": value for key, value in cost.items()})
        rows.append(row)
        predictions.append(
            {
                "repetition": repetition,
                "value": state[:, 0],
                "z": state[:, 1:],
            }
        )

    payload = {
        "schema_version": 1,
        "status": "complete",
        "stage": "actor_bottleneck_diagnostic",
        "method": method,
        "d": d,
        "k": k,
        "device": str(device),
        "elapsed_seconds": time.perf_counter() - started,
        "oracle_reference_used_inside_solver": method == "oracle_actor",
        "paired_with_learned_actor": True,
        "paired_seed_namespace": PAIRED_SEED_NAMESPACE,
        "learned_actor_record": str(learned_path.relative_to(RESULTS_DIR)),
        "learned_actor_record_sha256": learned_sha256,
        "rows": rows,
    }
    _save_predictions(prediction_path, predictions)
    atomic_json(output_path, payload)
    print(
        json.dumps(
            {
                key: payload[key]
                for key in ("status", "method", "d", "k", "elapsed_seconds")
            },
            indent=2,
        )
    )
    return payload


def _run_all(force: bool) -> None:
    config = load_config()
    dimensions = [int(value) for value in config["pde"]["dimensions"]]
    effective_dimensions = [
        int(value) for value in config["pde"]["effective_dimensions"]
    ]
    # Audit the originally preregistered ablation cell first, then resume the
    # full grid. Completed JSON cells are immutable checkpoints.
    ordered_cells = [(100, 10)] + [
        (d, k)
        for d in dimensions
        for k in effective_dimensions
        if (d, k) != (100, 10)
    ]
    for d, k in ordered_cells:
        for method in DIAGNOSTIC_METHODS:
            run_diagnostic(method, d, k, force=force)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=DIAGNOSTIC_METHODS)
    parser.add_argument("--d", type=int)
    parser.add_argument("--k", type=int)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        if args.aggregate_only:
            from .actor_diagnostics_report import aggregate_actor_diagnostics

            print(json.dumps(aggregate_actor_diagnostics(), indent=2))
        elif args.pilot:
            for method in DIAGNOSTIC_METHODS:
                run_diagnostic(method, 100, 10, force=args.force)
        elif args.all:
            _run_all(args.force)
            from .actor_diagnostics_report import aggregate_actor_diagnostics

            print(json.dumps(aggregate_actor_diagnostics(), indent=2))
        elif args.method is not None and args.d is not None and args.k is not None:
            run_diagnostic(args.method, args.d, args.k, force=args.force)
        else:
            parser.error("use --all, --pilot, or provide --method, --d, and --k")
    except Exception as error:
        DIAGNOSTIC_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        atomic_json(
            DIAGNOSTIC_RESULTS_DIR / "failure.json",
            {
                "status": "failed",
                "error": repr(error),
                "traceback": traceback.format_exc(),
            },
        )
        raise


if __name__ == "__main__":
    main()

