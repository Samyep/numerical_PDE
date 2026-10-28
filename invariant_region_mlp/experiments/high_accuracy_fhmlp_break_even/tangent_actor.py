"""Frozen-checkpoint Tangent-Actor accuracy and actual-prefix timing."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion import run_dense_actor as _source
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import synchronize

from .benchmark import make_benchmark
from .protocol import (
    DIMENSION,
    EXPECTED_TA_ARTIFACTS,
    JOINT_ACTOR_CONFIG,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    TA_MCV_ACCURACY_REPLICATES,
    TA_PATHS,
    TA_STEPS,
    TA_TRAINING_SEEDS,
    array_sha256,
    ensure_directories,
    file_sha256,
    load_bank,
    require_frozen_timing_grid,
    require_preflight,
    source_results_dir,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


def _require_h100(device: torch.device) -> str:
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("formal break-even computation requires CUDA")
    name = torch.cuda.get_device_name(device)
    if "H100" not in name.upper():
        raise RuntimeError(f"formal timing requires NVIDIA H100, found {name!r}")
    return name


def _run_id(seed: int) -> str:
    return f"joint_break_even_d{DIMENSION}_s{int(seed)}"


def _source_training_path(seed: int):
    return source_results_dir() / "raw" / "actor" / f"training_{_run_id(seed)}.json"


def _source_checkpoint_path(seed: int):
    selected = EXPECTED_TA_ARTIFACTS[int(seed)]["selected_actor"]
    return source_results_dir() / "checkpoints" / f"{_run_id(seed)}_{selected}.pt"


def source_artifact_audit() -> dict[str, Any]:
    rows = []
    for seed in TA_TRAINING_SEEDS:
        expected = EXPECTED_TA_ARTIFACTS[int(seed)]
        training_path = _source_training_path(seed)
        checkpoint_path = _source_checkpoint_path(seed)
        training = json.loads(training_path.read_text(encoding="utf-8"))
        checks = {
            "training_exists": training_path.is_file(),
            "checkpoint_exists": checkpoint_path.is_file(),
            "training_sha256": file_sha256(training_path) == expected["training_sha256"],
            "checkpoint_sha256": file_sha256(checkpoint_path) == expected["checkpoint_sha256"],
            "training_complete": training.get("status") == "complete",
            "selected_actor": training.get("selection", {}).get("selected_actor")
            == expected["selected_actor"],
            "truth_free_selection": training.get(
                "ground_truth_used_for_training_or_selection"
            ) is False,
            "charged_seconds_match": abs(
                float(training["method_active_training_seconds"])
                - float(expected["charged_training_seconds"])
            ) <= 1.0e-9,
        }
        rows.append({
            "seed": int(seed),
            "training_path": str(training_path),
            "checkpoint_path": str(checkpoint_path),
            "checks": checks,
            "overall_pass": all(checks.values()),
        })
    return {"rows": rows, "overall_pass": all(row["overall_pass"] for row in rows)}


def _configure(seed: int) -> None:
    _source.make_dense_fourier = make_benchmark
    _source._run_id = lambda _dimension, local_seed: _run_id(local_seed)
    _source.CHECKPOINT_DIR = source_results_dir() / "checkpoints"
    _source.CACHE_DIR = source_results_dir() / "cache"
    _source.JOINT_ACTOR_CONFIG = JOINT_ACTOR_CONFIG
    _source.PREDICTION_DIR = PREDICTION_DIR
    _source.PRIMARY_EULER_STEPS = TA_STEPS
    _source.RAW_DIR = RAW_DIR
    _source.REVISION = REVISION
    _source.TRAINING_SEED = int(seed)
    _source.ensure_directories = ensure_directories
    _source.stable_seed = stable_seed
    _source.write_json_atomic = write_json_atomic


def _fixed_training_seconds(training: dict[str, Any]) -> float:
    return float(training["method_active_training_seconds"])


def _load_selected(seed: int, device: torch.device):
    _configure(seed)
    training_path = _source_training_path(seed)
    training = json.loads(training_path.read_text(encoding="utf-8"))
    expected = EXPECTED_TA_ARTIFACTS[int(seed)]
    if training.get("status") != "complete":
        raise RuntimeError("source training record is incomplete")
    if file_sha256(training_path) != expected["training_sha256"]:
        raise RuntimeError("source Tangent-Actor training hash mismatch")
    checkpoint_path = _source_checkpoint_path(seed)
    if file_sha256(checkpoint_path) != expected["checkpoint_sha256"]:
        raise RuntimeError("source Tangent-Actor checkpoint hash mismatch")
    driver = make_benchmark(DIMENSION)
    synchronize(device)
    started = time.perf_counter()
    actor = _source._actors_for_evaluation(
        _run_id(seed), driver, training, ("selected",), device
    )["selected"]
    synchronize(device)
    return training, driver, actor, time.perf_counter() - started


def _direct_predict(actor, driver, t, x, device):
    synchronize(device)
    torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    prediction, _ = _source._predict(actor, driver, t, x, device)
    synchronize(device)
    value = np.asarray(prediction, dtype=np.float64).reshape(-1)
    return value, {
        "wall_seconds": time.perf_counter() - started,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "nonfinite_count": int((~np.isfinite(value)).sum()),
        "estimator": "selected_frozen_joint_value_actor_head",
        "paths": 0,
        "time_steps": 0,
    }


def _predict_candidate(actor, driver, t, x, candidate: str, seed: int, device):
    if candidate == "direct":
        return _direct_predict(actor, driver, t, x, device)
    paths = int(candidate.split("_")[1])
    output, metadata = _source.evaluate_policy(
        driver,
        actor,
        t,
        x,
        paths=paths,
        time_steps=TA_STEPS,
        seed=seed,
        device=device,
        parent_batch=8,
        path_chunk=min(64, paths),
    )
    return np.asarray(output["mcv"], dtype=np.float64).reshape(-1), metadata


def evaluate(
    stage: str,
    seed: int,
    candidate: str,
    replicate: int,
    device: torch.device,
) -> dict[str, Any]:
    if stage not in {"calibration", "final"}:
        raise KeyError(stage)
    if seed not in TA_TRAINING_SEEDS:
        raise ValueError(seed)
    if candidate != "direct" and candidate not in {f"mcv_{p}" for p in TA_PATHS}:
        raise ValueError(candidate)
    if candidate == "direct" and replicate != 0:
        raise ValueError("direct evaluation is deterministic conditional on training")
    if candidate != "direct" and not 0 <= replicate < TA_MCV_ACCURACY_REPLICATES:
        raise ValueError(replicate)
    require_preflight()
    gpu = _require_h100(device)
    bank = load_bank(stage, require_truth=True)
    training, driver, actor, load_seconds = _load_selected(seed, device)
    prediction, metadata = _predict_candidate(
        actor,
        driver,
        bank["t"],
        bank["x"],
        candidate,
        stable_seed("ta", stage, seed, candidate, replicate),
        device,
    )
    error = prediction - bank["reference"]
    run_id = f"ta__{stage}__s{seed}__{candidate}__r{replicate}"
    prediction_path = PREDICTION_DIR / stage / "tangent_actor" / f"{run_id}.npz"
    write_npz_atomic(
        prediction_path,
        t=bank["t"],
        x=bank["x"],
        prediction=prediction,
        reference=bank["reference"],
    )
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "status": "complete" if np.isfinite(prediction).all() else "numerical_failure",
        "stage": stage,
        "method": "tangent_actor",
        "training_seed": int(seed),
        "candidate": candidate,
        "inference_replicate": int(replicate),
        "selected_actor": training["selection"]["selected_actor"],
        "source_training_artifact": str(_source_training_path(seed)),
        "source_training_artifact_sha256": file_sha256(_source_training_path(seed)),
        "source_checkpoint": str(_source_checkpoint_path(seed)),
        "source_checkpoint_sha256": file_sha256(_source_checkpoint_path(seed)),
        "charged_training_seconds": _fixed_training_seconds(training),
        "points": len(prediction),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
        "inference_seconds": float(metadata["wall_seconds"]),
        "checkpoint_load_seconds": load_seconds,
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "peak_gpu_bytes": int(metadata.get("peak_gpu_bytes") or 0),
        "metadata": metadata,
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"], bank["x"]),
        "gpu": gpu,
        "truth_used_for_training_or_actor_selection": False,
        "inference_mode": "direct_value_head" if candidate == "direct" else "fk_mcv",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RAW_DIR / stage / "tangent_actor" / f"{run_id}.json", record)
    return record


def _next_prefix_fits(
    fixed: float, inference: float, previous: int, next_count: int, cap: float
) -> bool:
    if fixed + inference >= cap:
        return False
    if previous <= 0 or inference <= 0:
        return True
    projected = inference / previous * (next_count - previous)
    return fixed + inference + projected <= cap


def timing(
    replicate: int,
    candidate: str,
    device: torch.device,
    *,
    max_active_seconds: float,
) -> dict[str, Any]:
    if not 0 <= replicate < 3:
        raise ValueError(replicate)
    if candidate != "direct" and candidate not in {f"mcv_{p}" for p in TA_PATHS}:
        raise ValueError(candidate)
    require_preflight()
    frozen = require_frozen_timing_grid()
    counts = tuple(int(value) for value in frozen["timing_counts"])
    seed = TA_TRAINING_SEEDS[replicate]
    gpu = _require_h100(device)
    bank = load_bank("timing", require_truth=False)
    if "reference" in bank:
        raise RuntimeError("truth leaked into timing bank")
    training, driver, actor, load_seconds = _load_selected(seed, device)
    training_seconds = _fixed_training_seconds(training)
    predictions: list[np.ndarray] = []
    prefix_rows: list[dict[str, Any]] = []
    inference_seconds = 0.0
    previous = 0
    peak = 0
    for count in counts:
        fixed = training_seconds + load_seconds
        if not _next_prefix_fits(
            fixed, inference_seconds, previous, count, max_active_seconds
        ):
            break
        local, metadata = _predict_candidate(
            actor,
            driver,
            bank["t"][previous:count],
            bank["x"][previous:count],
            candidate,
            stable_seed("ta", "timing", replicate, candidate, previous, count),
            device,
        )
        predictions.append(local)
        inference_seconds += float(metadata["wall_seconds"])
        peak = max(peak, int(metadata.get("peak_gpu_bytes") or 0))
        prefix_rows.append({
            "queries": count,
            "increment_queries": count - previous,
            "inference_seconds": inference_seconds,
            "checkpoint_load_seconds": load_seconds,
            "charged_training_and_selection_seconds": training_seconds,
            "total_seconds": fixed + inference_seconds,
        })
        previous = count
    prediction = np.concatenate(predictions) if predictions else np.empty(0, dtype=np.float64)
    status = "complete" if previous == counts[-1] else "TIME_CAP"
    if not np.isfinite(prediction).all():
        status = "numerical_failure"
    run_id = f"ta__timing__{candidate}__r{replicate}"
    prediction_path = PREDICTION_DIR / "timing" / "tangent_actor" / f"{run_id}.npz"
    write_npz_atomic(
        prediction_path,
        prediction=prediction,
        prefix_counts=np.asarray([row["queries"] for row in prefix_rows], dtype=np.int64),
    )
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "status": status,
        "method": "tangent_actor",
        "timing_replicate": int(replicate),
        "training_seed": int(seed),
        "candidate": candidate,
        "selected_actor": training["selection"]["selected_actor"],
        "source_training_artifact": str(_source_training_path(seed)),
        "source_checkpoint": str(_source_checkpoint_path(seed)),
        "prefix_rows": prefix_rows,
        "completed_queries": previous,
        "charged_training_seconds": training_seconds,
        "checkpoint_load_seconds": load_seconds,
        "inference_seconds": inference_seconds,
        "active_seconds": training_seconds + load_seconds + inference_seconds,
        "actual_new_experiment_seconds": load_seconds + inference_seconds,
        "peak_gpu_bytes": peak,
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"][:previous], bank["x"][:previous]),
        "gpu": gpu,
        "truth_loaded": False,
        "actual_prefix_execution": True,
        "inference_mode": "direct_value_head" if candidate == "direct" else "fk_mcv",
    }
    write_json_atomic(RAW_DIR / "timing" / "tangent_actor" / f"{run_id}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    item = sub.add_parser("evaluate")
    item.add_argument("--stage", choices=("calibration", "final"), required=True)
    item.add_argument("--seed", type=int, required=True)
    item.add_argument("--candidate", required=True)
    item.add_argument("--replicate", type=int, required=True)
    item = sub.add_parser("timing")
    item.add_argument("--replicate", type=int, required=True)
    item.add_argument("--candidate", required=True)
    item.add_argument("--max-active-seconds", type=float, required=True)
    for item in sub.choices.values():
        item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if args.action == "evaluate":
        result = evaluate(args.stage, args.seed, args.candidate, args.replicate, device)
    else:
        result = timing(
            args.replicate,
            args.candidate,
            device,
            max_active_seconds=args.max_active_seconds,
        )
    print(json.dumps({
        key: result.get(key)
        for key in ("run_id", "status", "rmse", "active_seconds", "completed_queries")
    }, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["evaluate", "source_artifact_audit", "timing"]
