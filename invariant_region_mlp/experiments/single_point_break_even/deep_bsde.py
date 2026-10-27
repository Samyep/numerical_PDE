"""Native scalar-Y0 Deep BSDE calibration, final accuracy, and timing."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import random
import time
from typing import Any, Iterable

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    seed_everything,
    synchronize,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port import (
    DeepBSDEConfig,
    _sample_dw,
    _validation,
    clipped_terminal_loss,
)
from invariant_region_mlp.experiments.single_point_32_ice.benchmarks import (
    make_benchmark,
)
from invariant_region_mlp.experiments.single_point_32_ice.deep_bsde import (
    FixedInitialDeepBSDE,
)

from .protocol import (
    CHECKPOINT_DIR,
    DB_CANDIDATE_UPDATES,
    DB_TIMING_COUNTS,
    DIMENSION,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    array_sha256,
    ensure_directories,
    file_sha256,
    load_bank,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


CONFIG = DeepBSDEConfig(
    benchmark="allen_cahn_d100",
    total_time=0.3,
    intervals=20,
    hidden_widths=(110, 110),
    learning_rate=5e-4,
    iterations=8000,
    batch_size=64,
    validation_size=256,
    logging_frequency=2000,
    y_init_range=(0.3, 0.5),
    dtype="float64",
)


def _require_h100(device: torch.device) -> str:
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("formal break-even computation requires CUDA")
    name = torch.cuda.get_device_name(device)
    if "H100" not in name.upper():
        raise RuntimeError(f"formal timing requires NVIDIA H100, found {name!r}")
    return name


def _rng_state() -> dict[str, Any]:
    payload: dict[str, Any] = {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.random.get_rng_state(),
    }
    if torch.cuda.is_available():
        payload["torch_cuda"] = torch.cuda.get_rng_state_all()
    return payload


def _restore_rng(payload: dict[str, Any]) -> None:
    random.setstate(payload["python"])
    np.random.set_state(payload["numpy"])
    torch.random.set_rng_state(payload["torch_cpu"])
    if torch.cuda.is_available() and "torch_cuda" in payload:
        torch.cuda.set_rng_state_all(payload["torch_cuda"])


def _atomic_torch(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def _train_one(
    x0: np.ndarray,
    *,
    full_seed: int,
    capture_updates: Iterable[int],
    max_updates: int,
    device: torch.device,
    resume_path: Path | None,
    checkpoint_path: Path | None,
) -> dict[str, Any]:
    capture = tuple(sorted({int(value) for value in capture_updates}))
    if not capture or capture[-1] > max_updates:
        raise ValueError((capture, max_updates))
    benchmark = make_benchmark("allen_cahn", DIMENSION)
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    seed_everything(full_seed)
    model = FixedInitialDeepBSDE(
        benchmark, CONFIG, torch.as_tensor(x0, dtype=torch.float64)
    ).to(device=device, dtype=torch.float64)
    optimizer = torch.optim.Adam(model.parameters(), lr=CONFIG.learning_rate, eps=1e-8)
    train_generator = torch.Generator(device=device).manual_seed(full_seed + 1)
    valid_generator = torch.Generator(device=device).manual_seed(full_seed + 2)
    validation_dw = _sample_dw(
        CONFIG,
        DIMENSION,
        CONFIG.validation_size,
        generator=valid_generator,
        device=device,
    )
    completed = 0
    candidates: list[dict[str, Any]] = []
    nonfinite = 0
    elapsed_previous = 0.0
    if resume_path is not None and resume_path.is_file():
        state = torch.load(resume_path, map_location=device, weights_only=False)
        if state["revision"] != REVISION or state["config"] != asdict(CONFIG):
            raise RuntimeError("Deep BSDE resume metadata mismatch")
        if state["x0_sha256"] != array_sha256(x0):
            raise RuntimeError("Deep BSDE resume point mismatch")
        model.load_state_dict(state["model"])
        optimizer.load_state_dict(state["optimizer"])
        train_generator.set_state(state["train_generator"])
        validation_dw = state["validation_dw"].to(device=device, dtype=torch.float64)
        _restore_rng(state["rng"])
        completed = int(state["completed_updates"])
        candidates = list(state["candidates"])
        nonfinite = int(state["nonfinite_count"])
        elapsed_previous = float(state["elapsed_seconds"])
        synchronize(device)
        started = time.perf_counter()
    validation_hash = array_sha256(validation_dw.detach().cpu().numpy())
    capture_set = set(capture)
    while completed < max_updates:
        model.train()
        dw = _sample_dw(
            CONFIG,
            DIMENSION,
            CONFIG.batch_size,
            generator=train_generator,
            device=device,
        )
        terminal_y, terminal_x = model.rollout(dw)
        loss = clipped_terminal_loss(terminal_y - benchmark.terminal(terminal_x))
        if not torch.isfinite(loss):
            nonfinite += 1
            break
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        completed += 1
        if completed in capture_set:
            validation_loss, bad = _validation(model, validation_dw)
            nonfinite += bad
            synchronize(device)
            elapsed = elapsed_previous + time.perf_counter() - started
            candidates.append({
                "updates": completed,
                "y0": float(model.y_init.detach()),
                "validation_terminal_clipped_mse": validation_loss,
                "active_seconds": elapsed,
                "nonfinite_count_cumulative": nonfinite,
            })
            if resume_path is not None:
                _atomic_torch(resume_path, {
                    "revision": REVISION,
                    "config": asdict(CONFIG),
                    "x0_sha256": array_sha256(x0),
                    "model": model.state_dict(),
                    "optimizer": optimizer.state_dict(),
                    "train_generator": train_generator.get_state(),
                    "validation_dw": validation_dw.detach().cpu(),
                    "validation_dw_sha256": validation_hash,
                    "rng": _rng_state(),
                    "completed_updates": completed,
                    "candidates": candidates,
                    "nonfinite_count": nonfinite,
                    "elapsed_seconds": elapsed,
                })
    synchronize(device)
    active_seconds = elapsed_previous + time.perf_counter() - started
    status = "complete" if completed == max_updates and nonfinite == 0 else "numerical_failure"
    if checkpoint_path is not None:
        _atomic_torch(checkpoint_path, {
            "revision": REVISION,
            "config": asdict(CONFIG),
            "x0": np.asarray(x0, dtype=np.float64),
            "completed_updates": completed,
            "candidates": candidates,
            "state_dict": {
                key: value.detach().cpu() for key, value in model.state_dict().items()
            },
        })
    return {
        "status": status,
        "completed_updates": completed,
        "candidates": candidates,
        "active_seconds": active_seconds,
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "validation_dw_sha256": validation_hash,
    }


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
    ensure_directories()
    gpu = _require_h100(device)
    bank = load_bank(stage, require_truth=True)
    if not 0 <= point < len(bank["x"]):
        raise ValueError(point)
    if stage == "calibration" and (
        max_updates != 8000 or tuple(capture_updates) != DB_CANDIDATE_UPDATES
    ):
        raise ValueError("calibration must run the frozen 2k/4k/8k trajectory")
    capture_tag = "-".join(map(str, capture_updates))
    result_id = f"db__{stage}__q{point:03d}__s{seed}__u{max_updates}__c{capture_tag}"
    output_path = RAW_DIR / stage / "deep_bsde" / f"{result_id}.json"
    if output_path.is_file():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        if existing.get("status") in {"complete", "numerical_failure"}:
            return existing
    resume_path = CHECKPOINT_DIR / stage / "deep_bsde" / f"{result_id}_resume.pt"
    checkpoint_path = CHECKPOINT_DIR / stage / "deep_bsde" / f"{result_id}.pt"
    result = _train_one(
        bank["x"][point],
        full_seed=stable_seed("deep_bsde", stage, point, seed),
        capture_updates=capture_updates,
        max_updates=max_updates,
        device=device,
        resume_path=resume_path,
        checkpoint_path=checkpoint_path,
    )
    candidates = []
    reference = float(bank["reference"][point])
    for row in result["candidates"]:
        local = dict(row)
        local["reference"] = reference
        local["error"] = local["y0"] - reference
        local["absolute_error"] = abs(local["error"])
        candidates.append(local)
    prediction_path = PREDICTION_DIR / stage / "deep_bsde" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        updates=np.asarray([row["updates"] for row in candidates], dtype=np.int64),
        prediction=np.asarray([row["y0"] for row in candidates], dtype=np.float64),
        reference=np.asarray(reference, dtype=np.float64),
    )
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": result["status"],
        "stage": stage,
        "method": "deep_bsde_native_scalar_y0",
        "point": point,
        "seed": seed,
        "max_updates": max_updates,
        "capture_updates": list(capture_updates),
        "configuration": asdict(CONFIG),
        "candidates": candidates,
        "active_seconds": result["active_seconds"],
        "completed_updates": result["completed_updates"],
        "nonfinite_count": result["nonfinite_count"],
        "peak_gpu_bytes": result["peak_gpu_bytes"],
        "parameter_count": result["parameter_count"],
        "validation_dw_sha256": result["validation_dw_sha256"],
        "checkpoint": str(checkpoint_path),
        "checkpoint_sha256": file_sha256(checkpoint_path),
        "resume_checkpoint": str(resume_path),
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(
            np.asarray([row["y0"] for row in candidates], dtype=np.float64)
        ),
        "query_input_sha256": array_sha256(
            bank["t"][point:point + 1], bank["x"][point:point + 1]
        ),
        "gpu": gpu,
        "native_scalar_y0": True,
        "amortized_decoder": False,
        "fresh_model_per_query": True,
        "truth_loaded_only_after_frozen_training_schedule": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(output_path, record)
    return record


def timing(
    updates: int,
    replicate: int,
    device: torch.device,
    *,
    max_active_seconds: float = 3600.0,
) -> dict[str, Any]:
    gpu = _require_h100(device)
    bank = load_bank("timing", require_truth=False)
    if "reference" in bank:
        raise RuntimeError("truth leaked into timing bank")
    predictions: list[float] = []
    point_seconds: list[float] = []
    prefix_rows: list[dict[str, Any]] = []
    active = 0.0
    next_prefix = 0
    nonfinite = 0
    peak = 0
    for point in range(DB_TIMING_COUNTS[-1]):
        if active >= max_active_seconds:
            break
        local = _train_one(
            bank["x"][point],
            full_seed=stable_seed("deep_bsde", "timing", updates, replicate, point),
            capture_updates=(updates,),
            max_updates=updates,
            device=device,
            resume_path=None,
            checkpoint_path=None,
        )
        active += float(local["active_seconds"])
        point_seconds.append(float(local["active_seconds"]))
        peak = max(peak, int(local["peak_gpu_bytes"]))
        nonfinite += int(local["nonfinite_count"])
        if local["candidates"]:
            predictions.append(float(local["candidates"][-1]["y0"]))
        else:
            predictions.append(float("nan"))
        completed = point + 1
        if completed == DB_TIMING_COUNTS[next_prefix]:
            prefix_rows.append({
                "queries": completed,
                "training_seconds": active,
                "total_seconds": active,
            })
            next_prefix += 1
            if next_prefix == len(DB_TIMING_COUNTS):
                break
    completed = len(predictions)
    prediction = np.asarray(predictions, dtype=np.float64)
    status = "complete" if completed == DB_TIMING_COUNTS[-1] else "TIME_CAP"
    if nonfinite or not np.all(np.isfinite(prediction)):
        status = "numerical_failure"
    result_id = f"db__timing__u{updates}__r{replicate}"
    prediction_path = PREDICTION_DIR / "timing" / "deep_bsde" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        prediction=prediction,
        point_seconds=np.asarray(point_seconds, dtype=np.float64),
    )
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": status,
        "method": "deep_bsde_native_scalar_y0",
        "updates": updates,
        "replicate": replicate,
        "prefix_rows": prefix_rows,
        "point_seconds": point_seconds,
        "completed_queries": completed,
        "training_seconds": active,
        "active_seconds": active,
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": peak,
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"][:completed], bank["x"][:completed]),
        "gpu": gpu,
        "truth_loaded": False,
        "fresh_model_per_query": True,
        "actual_prefix_execution": True,
    }
    write_json_atomic(RAW_DIR / "timing" / "deep_bsde" / f"{result_id}.json", record)
    return record


def _parse_updates(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item)


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--stage", choices=("calibration", "final"), required=True)
    evaluate_parser.add_argument("--point", type=int, required=True)
    evaluate_parser.add_argument("--seed", type=int, required=True)
    evaluate_parser.add_argument("--max-updates", type=int, required=True)
    evaluate_parser.add_argument("--capture-updates", type=_parse_updates, required=True)
    timing_parser = subparsers.add_parser("timing")
    timing_parser.add_argument("--updates", type=int, required=True)
    timing_parser.add_argument("--replicate", type=int, required=True)
    timing_parser.add_argument("--max-active-seconds", type=float, default=3600.0)
    for item in (evaluate_parser, timing_parser):
        item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if args.action == "evaluate":
        result = evaluate(
            args.stage,
            args.point,
            args.seed,
            device,
            max_updates=args.max_updates,
            capture_updates=args.capture_updates,
        )
    else:
        result = timing(
            args.updates,
            args.replicate,
            device,
            max_active_seconds=args.max_active_seconds,
        )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "active_seconds", "completed_queries",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

