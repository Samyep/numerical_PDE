"""Native scalar-Y0 Deep BSDE adapter for the joint nonseparable PDE."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import DeepBSDEConfig
from invariant_region_mlp.experiments.single_point_break_even import deep_bsde as _source

from .benchmark import make_benchmark
from .protocol import (
    ACTIVE_CAP_SECONDS,
    CHECKPOINT_DIR,
    DB_CANDIDATE_UPDATES,
    DB_TIMING_COUNTS,
    DIMENSION,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    array_sha256,
    bank_path,
    ensure_directories,
    file_sha256,
    load_bank,
    require_preflight,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


CONFIG = DeepBSDEConfig(
    benchmark="joint_nonseparable_d100",
    total_time=0.5,
    intervals=20,
    hidden_widths=(110, 110),
    learning_rate=1e-2,
    iterations=4000,
    batch_size=64,
    validation_size=256,
    logging_frequency=500,
    y_init_range=(0.0, 1.0),
    dtype="float64",
)
FixedInitialDeepBSDE = _source.FixedInitialDeepBSDE


def _configure() -> None:
    _source.CONFIG = CONFIG
    _source.CHECKPOINT_DIR = CHECKPOINT_DIR
    _source.DB_CANDIDATE_UPDATES = DB_CANDIDATE_UPDATES
    _source.DB_TIMING_COUNTS = DB_TIMING_COUNTS
    _source.DIMENSION = DIMENSION
    _source.PREDICTION_DIR = PREDICTION_DIR
    _source.RAW_DIR = RAW_DIR
    _source.REVISION = REVISION
    _source.array_sha256 = array_sha256
    _source.ensure_directories = ensure_directories
    _source.file_sha256 = file_sha256
    _source.load_bank = load_bank
    _source.make_benchmark = lambda _pde, dimension: make_benchmark(dimension)
    _source.stable_seed = stable_seed
    _source.write_json_atomic = write_json_atomic
    _source.write_npz_atomic = write_npz_atomic


def evaluate(
    stage: str,
    point: int,
    seed: int,
    device: torch.device,
    *,
    max_updates: int,
    capture_updates: tuple[int, ...],
) -> dict[str, Any]:
    if stage not in {"calibration", "final"}:
        raise KeyError(stage)
    if stage == "calibration" and (
        max_updates != DB_CANDIDATE_UPDATES[-1]
        or tuple(capture_updates) != DB_CANDIDATE_UPDATES
    ):
        raise ValueError("calibration must run the frozen 500/1000/2000/4000 trajectory")
    require_preflight()
    _configure()
    ensure_directories()
    gpu = _source._require_h100(device)
    with np.load(bank_path(stage), allow_pickle=False) as data:
        count = len(data["x"])
        if not 0 <= point < count:
            raise ValueError(point)
        x0 = np.asarray(data["x"][point], dtype=np.float64)
        t0 = np.asarray(data["t"][point:point + 1], dtype=np.float64)
    capture_tag = "-".join(map(str, capture_updates))
    result_id = f"db__{stage}__q{point:03d}__s{seed}__u{max_updates}__c{capture_tag}"
    output_path = RAW_DIR / stage / "deep_bsde" / f"{result_id}.json"
    if output_path.is_file():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        if existing.get("status") in {"complete", "numerical_failure"}:
            return existing
    resume_path = CHECKPOINT_DIR / stage / "deep_bsde" / f"{result_id}_resume.pt"
    checkpoint_path = CHECKPOINT_DIR / stage / "deep_bsde" / f"{result_id}.pt"
    result = _source._train_one(
        x0,
        full_seed=stable_seed("deep_bsde", stage, point, seed),
        capture_updates=capture_updates,
        max_updates=max_updates,
        device=device,
        resume_path=resume_path,
        checkpoint_path=checkpoint_path,
    )
    # Truth is opened only after the complete, predeclared training trajectory.
    with np.load(bank_path(stage), allow_pickle=False) as data:
        reference = float(data["reference"][point])
    candidates = []
    for row in result["candidates"]:
        local = dict(row)
        local.update({
            "reference": reference,
            "error": float(local["y0"]) - reference,
            "absolute_error": abs(float(local["y0"]) - reference),
        })
        candidates.append(local)
    prediction_path = PREDICTION_DIR / stage / "deep_bsde" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        updates=np.asarray([row["updates"] for row in candidates], dtype=np.int64),
        prediction=np.asarray([row["y0"] for row in candidates], dtype=np.float64),
        reference=np.asarray(reference, dtype=np.float64),
    )
    record = {
        "revision": REVISION, "run_id": result_id, "status": result["status"],
        "stage": stage, "method": "deep_bsde_native_scalar_y0", "point": point,
        "seed": seed, "max_updates": max_updates, "capture_updates": list(capture_updates),
        "configuration": asdict(CONFIG), "candidates": candidates,
        "active_seconds": result["active_seconds"], "completed_updates": result["completed_updates"],
        "nonfinite_count": result["nonfinite_count"], "peak_gpu_bytes": result["peak_gpu_bytes"],
        "parameter_count": result["parameter_count"], "validation_dw_sha256": result["validation_dw_sha256"],
        "checkpoint": str(checkpoint_path), "checkpoint_sha256": file_sha256(checkpoint_path),
        "resume_checkpoint": str(resume_path), "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(np.asarray([row["y0"] for row in candidates], dtype=np.float64)),
        "query_input_sha256": array_sha256(t0, x0[None, :]), "gpu": gpu,
        "native_scalar_y0": True, "amortized_decoder": False, "fresh_model_per_query": True,
        "truth_loaded_only_after_frozen_training_schedule": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(output_path, record)
    return record


def timing(updates: int, replicate: int, device: torch.device, *, max_active_seconds: float = ACTIVE_CAP_SECONDS):
    require_preflight()
    _configure()
    return _source.timing(int(updates), int(replicate), device, max_active_seconds=max_active_seconds)


def _parse_updates(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("evaluate"); p.add_argument("--stage", choices=("calibration", "final"), required=True); p.add_argument("--point", type=int, required=True); p.add_argument("--seed", type=int, required=True); p.add_argument("--max-updates", type=int, required=True); p.add_argument("--capture-updates", type=_parse_updates, required=True)
    p = sub.add_parser("timing"); p.add_argument("--updates", type=int, required=True); p.add_argument("--replicate", type=int, required=True); p.add_argument("--max-active-seconds", type=float, default=ACTIVE_CAP_SECONDS)
    for item in sub.choices.values(): item.add_argument("--device", default="cuda")
    args = parser.parse_args(); device = torch.device(args.device)
    if args.action == "evaluate": result = evaluate(args.stage, args.point, args.seed, device, max_updates=args.max_updates, capture_updates=args.capture_updates)
    else: result = timing(args.updates, args.replicate, device, max_active_seconds=args.max_active_seconds)
    print(json.dumps({key: result.get(key) for key in ("run_id", "status", "active_seconds", "completed_queries")}, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["CONFIG", "FixedInitialDeepBSDE", "evaluate", "timing"]
