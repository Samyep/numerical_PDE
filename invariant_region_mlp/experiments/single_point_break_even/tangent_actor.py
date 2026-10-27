"""Isolated adapter around the frozen mature Allen--Cahn Tangent-Actor."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    synchronize,
)
from invariant_region_mlp.experiments.single_point_32_ice import (
    tangent_actor as source_ta,
)

from .protocol import (
    CACHE_DIR,
    CHECKPOINT_DIR,
    DIMENSION,
    FINAL_SEEDS,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    TA_INFERENCE_REPLICATES,
    TA_PATHS,
    TA_STEPS,
    TIMING_COUNTS,
    array_sha256,
    ensure_directories,
    file_sha256,
    load_bank,
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


def _next_prefix_fits(
    *,
    fixed_seconds: float,
    inference_seconds: float,
    previous_queries: int,
    next_queries: int,
    max_active_seconds: float,
) -> bool:
    """Truth-free guard for the frozen active-compute timing cap."""
    if fixed_seconds + inference_seconds >= max_active_seconds:
        return False
    if previous_queries <= 0 or inference_seconds <= 0.0:
        return True
    seconds_per_query = inference_seconds / previous_queries
    projected_increment = seconds_per_query * (next_queries - previous_queries)
    return fixed_seconds + inference_seconds + projected_increment <= max_active_seconds


def _configure():
    """Redirect every source artifact into this isolated result namespace."""
    source_ta.RAW_DIR = RAW_DIR
    source_ta.CACHE_DIR = CACHE_DIR
    source_ta.CHECKPOINT_DIR = CHECKPOINT_DIR
    source_ta.PREDICTION_DIR = PREDICTION_DIR
    source_ta.REVISION = REVISION
    source_ta.stable_seed = stable_seed
    source_ta.file_sha256 = file_sha256
    source_ta.write_json_atomic = write_json_atomic
    return source_ta._configure_standard(False)


def train(seed: int, device: torch.device) -> dict[str, Any]:
    """Run the frozen four-update, reference-free mature configuration."""
    ensure_directories()
    gpu = _require_h100(device)
    implementation = _configure()
    record = source_ta._train_allen_value_mcv(
        implementation, DIMENSION, int(seed), device, False
    )
    result = dict(record)
    result.update({
        "study_revision": REVISION,
        "gpu": gpu,
        "training_seed": int(seed),
        "formal_protocol": True,
    })
    # Preserve the source training record and add a study-level immutable view.
    path = RAW_DIR / "tangent_actor" / f"training_seed{seed}.json"
    write_json_atomic(path, result)
    return result


def _load_selected(seed: int, device: torch.device):
    implementation = _configure()
    started = time.perf_counter()
    data = source_ta._selected_standard(
        implementation, "allen_cahn", DIMENSION, int(seed), device, False
    )
    synchronize(device)
    load_seconds = time.perf_counter() - started
    return (*data, load_seconds)


def _direct_predict(
    implementation,
    actor,
    benchmark,
    t: np.ndarray,
    x: np.ndarray,
    device: torch.device,
) -> tuple[np.ndarray, dict[str, Any]]:
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    prediction = np.asarray(
        implementation._predict(actor, benchmark, t, x, device), dtype=np.float64
    ).reshape(-1)
    synchronize(device)
    elapsed = time.perf_counter() - started
    return prediction, {
        "wall_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "estimator": "selected_frozen_value_actor_head",
        "paths": 0,
        "time_steps": 0,
    }


def evaluate(
    stage: str,
    seed: int,
    candidate: str,
    replicate: int,
    device: torch.device,
) -> dict[str, Any]:
    if stage not in {"calibration", "final"}:
        raise KeyError(stage)
    if candidate != "direct" and not (
        candidate.startswith("mcv_") and int(candidate.split("_")[1]) in TA_PATHS
    ):
        raise ValueError(candidate)
    if candidate == "direct" and replicate != 0:
        raise ValueError("direct Actor evaluation is deterministic conditional on training")
    if candidate != "direct" and not 0 <= replicate < TA_INFERENCE_REPLICATES:
        raise ValueError(replicate)
    gpu = _require_h100(device)
    bank = load_bank(stage, require_truth=True)
    implementation = _configure()
    (
        run_id,
        training_path,
        selected,
        benchmark,
        actor,
        _,
        checkpoint_load_seconds,
    ) = _load_selected(seed, device)
    if candidate == "direct":
        prediction, metadata = _direct_predict(
            implementation, actor, benchmark, bank["t"], bank["x"], device
        )
    else:
        paths = int(candidate.split("_")[1])
        output, metadata = source_ta._evaluate_value_mcv_policy(
            benchmark,
            actor,
            bank["t"],
            bank["x"],
            paths=paths,
            time_steps=TA_STEPS,
            seed=stable_seed("ta", stage, seed, candidate, replicate),
            device=device,
            parent_batch=8,
            path_chunk=min(32, paths),
        )
        prediction = np.asarray(output["mcv"], dtype=np.float64)
    error = prediction - bank["reference"]
    result_id = f"ta__{stage}__s{seed}__{candidate}__r{replicate}"
    prediction_path = PREDICTION_DIR / stage / "tangent_actor" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        t=bank["t"],
        x=bank["x"],
        prediction=prediction,
        reference=bank["reference"],
    )
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": "complete" if np.all(np.isfinite(prediction)) else "numerical_failure",
        "stage": stage,
        "method": "tangent_actor",
        "training_seed": int(seed),
        "candidate": candidate,
        "inference_replicate": int(replicate),
        "selected_actor": selected,
        "source_training_run_id": run_id,
        "source_training_artifact": str(training_path),
        "source_training_artifact_sha256": file_sha256(training_path),
        "points": len(prediction),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
        "inference_seconds": float(metadata["wall_seconds"]),
        "checkpoint_load_seconds": checkpoint_load_seconds,
        "selection_time_seconds": float(metadata["wall_seconds"]),
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "peak_gpu_bytes": metadata.get("peak_gpu_bytes"),
        "metadata": metadata,
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"], bank["x"]),
        "gpu": gpu,
        "truth_used_for_training_or_actor_selection": False,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    output_path = RAW_DIR / stage / "tangent_actor" / f"{result_id}.json"
    write_json_atomic(output_path, record)
    return record


def timing(
    seed: int,
    candidate: str,
    device: torch.device,
    *,
    max_active_seconds: float = 3600.0,
) -> dict[str, Any]:
    """Measure actual incremental timing-bank prefixes; never read truth."""
    gpu = _require_h100(device)
    bank = load_bank("timing", require_truth=False)
    if "reference" in bank:
        raise RuntimeError("truth leaked into timing bank")
    implementation = _configure()
    (
        run_id,
        training_path,
        selected,
        benchmark,
        actor,
        _,
        checkpoint_load_seconds,
    ) = _load_selected(seed, device)
    training_record = json.loads(training_path.read_text(encoding="utf-8"))
    training_seconds = float(training_record["invocation_seconds"])
    predictions: list[np.ndarray] = []
    prefix_rows: list[dict[str, Any]] = []
    inference_elapsed = 0.0
    previous = 0
    peak = 0
    for count in TIMING_COUNTS:
        fixed_seconds = training_seconds + checkpoint_load_seconds
        # Avoid knowingly starting a large doubling interval that cannot fit
        # in the frozen per-replicate active-compute cap.  This forecast uses
        # only timing already observed in this truth-free run; it never
        # fabricates or reports the unexecuted prefix.
        if not _next_prefix_fits(
            fixed_seconds=fixed_seconds,
            inference_seconds=inference_elapsed,
            previous_queries=previous,
            next_queries=count,
            max_active_seconds=max_active_seconds,
        ):
            break
        local_t = bank["t"][previous:count]
        local_x = bank["x"][previous:count]
        if candidate == "direct":
            local, metadata = _direct_predict(
                implementation, actor, benchmark, local_t, local_x, device
            )
        else:
            paths = int(candidate.split("_")[1])
            output, metadata = source_ta._evaluate_value_mcv_policy(
                benchmark,
                actor,
                local_t,
                local_x,
                paths=paths,
                time_steps=TA_STEPS,
                seed=stable_seed("ta", "timing", seed, candidate, previous, count),
                device=device,
                parent_batch=8,
                path_chunk=min(32, paths),
            )
            local = np.asarray(output["mcv"], dtype=np.float64)
        predictions.append(np.asarray(local, dtype=np.float64))
        inference_elapsed += float(metadata["wall_seconds"])
        peak = max(peak, int(metadata.get("peak_gpu_bytes") or 0))
        prefix_rows.append({
            "queries": count,
            "increment_queries": count - previous,
            "inference_seconds": inference_elapsed,
            "checkpoint_load_seconds": checkpoint_load_seconds,
            "training_seconds": training_seconds,
            "total_seconds": training_seconds + checkpoint_load_seconds + inference_elapsed,
        })
        previous = count
    prediction = np.concatenate(predictions) if predictions else np.empty(0)
    status = "complete" if previous == TIMING_COUNTS[-1] else "TIME_CAP"
    result_id = f"ta__timing__s{seed}__{candidate}"
    prediction_path = PREDICTION_DIR / "timing" / "tangent_actor" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        prediction=prediction,
        prefix_counts=np.asarray([row["queries"] for row in prefix_rows], dtype=np.int64),
    )
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": status,
        "method": "tangent_actor",
        "training_seed": int(seed),
        "candidate": candidate,
        "selected_actor": selected,
        "source_training_run_id": run_id,
        "source_training_artifact": str(training_path),
        "prefix_rows": prefix_rows,
        "completed_queries": previous,
        "training_seconds": training_seconds,
        "checkpoint_load_seconds": checkpoint_load_seconds,
        "inference_seconds": inference_elapsed,
        "active_seconds": training_seconds + checkpoint_load_seconds + inference_elapsed,
        "peak_gpu_bytes": peak,
        "nonfinite_count": int((~np.isfinite(prediction)).sum()),
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"][:previous], bank["x"][:previous]),
        "gpu": gpu,
        "truth_loaded": False,
        "actual_prefix_execution": True,
    }
    write_json_atomic(RAW_DIR / "timing" / "tangent_actor" / f"{result_id}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    train_parser = subparsers.add_parser("train")
    train_parser.add_argument("--seed", type=int, required=True)
    eval_parser = subparsers.add_parser("evaluate")
    eval_parser.add_argument("--stage", choices=("calibration", "final"), required=True)
    eval_parser.add_argument("--seed", type=int, required=True)
    eval_parser.add_argument("--candidate", required=True)
    eval_parser.add_argument("--replicate", type=int, required=True)
    timing_parser = subparsers.add_parser("timing")
    timing_parser.add_argument("--seed", type=int, required=True)
    timing_parser.add_argument("--candidate", required=True)
    timing_parser.add_argument("--max-active-seconds", type=float, default=3600.0)
    for item in (train_parser, eval_parser, timing_parser):
        item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if args.action == "train":
        result = train(args.seed, device)
    elif args.action == "evaluate":
        result = evaluate(
            args.stage, args.seed, args.candidate, args.replicate, device
        )
    else:
        result = timing(
            args.seed,
            args.candidate,
            device,
            max_active_seconds=args.max_active_seconds,
        )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "rmse", "active_seconds", "completed_queries",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

