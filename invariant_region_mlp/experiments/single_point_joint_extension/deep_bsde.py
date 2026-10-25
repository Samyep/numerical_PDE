"""Isolated adapter for the completed suite's checkpointable Deep BSDE."""

from __future__ import annotations

import argparse
import json
import math
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.single_point_32_ice import deep_bsde as _frozen

from .benchmark import make_benchmark
from .protocol import (
    CHECKPOINT_DIR,
    PDE,
    PREDICTION_DIR,
    QUERY_COUNT,
    RAW_DIR,
    REVISION,
    TRAINING_SEEDS,
    array_sha256,
    deep_bsde_configuration,
    ensure_directories,
    file_sha256,
    query_path,
    require_preflight,
    require_smoke,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


FixedInitialDeepBSDE = _frozen.FixedInitialDeepBSDE


def _json_safe(value: Any) -> Any:
    """Replace non-finite scalar diagnostics without hiding numerical failure."""
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _safe_deep_record(payload: Any) -> Any:
    """Make a Deep-BSDE record JSON-safe and classify terminal nonfinites."""
    safe = _json_safe(payload)
    if isinstance(payload, dict) and "reported_y0" in payload:
        estimate = payload.get("reported_y0")
        if not isinstance(estimate, (int, float)) or not math.isfinite(float(estimate)):
            safe["status"] = "numerical_failure"
            safe["failure_reason"] = "nonfinite_prediction"
    return safe


def _write_deep_record(path, payload: Any) -> None:
    """Persist a failed Deep-BSDE run as data instead of raising on NaN."""
    safe = _safe_deep_record(payload)
    write_json_atomic(path, safe)


def _configure() -> None:
    _frozen.make_benchmark = lambda pde, d: make_benchmark(d)
    _frozen.CHECKPOINT_DIR = CHECKPOINT_DIR
    _frozen.PREDICTION_DIR = PREDICTION_DIR
    _frozen.QUERY_COUNT = QUERY_COUNT
    _frozen.RAW_DIR = RAW_DIR
    _frozen.REVISION = REVISION
    _frozen.TRAINING_SEEDS = TRAINING_SEEDS
    _frozen.array_sha256 = array_sha256
    _frozen.deep_bsde_configuration = (
        lambda pde, d, pilot=False: deep_bsde_configuration(d, pilot=pilot)
    )
    _frozen.ensure_directories = ensure_directories
    _frozen.file_sha256 = file_sha256
    _frozen.query_path = lambda pde, d: query_path(d)
    _frozen.require_stage_a = require_preflight
    _frozen.require_stage_b = require_smoke
    _frozen.stable_seed = stable_seed
    _frozen.write_json_atomic = _write_deep_record
    _frozen.write_npz_atomic = write_npz_atomic


def _recover_completed_nonfinite(
    dimension: int,
    point: int,
    seed: int,
    *,
    pilot: bool,
) -> dict[str, Any] | None:
    """Recover a terminal NaN record from already completed frozen artifacts."""
    tag = "pilot" if pilot else "formal"
    run_id = (
        f"deep_bsde__{PDE}__d{dimension}__q{point:02d}__s{seed}__{tag}"
    )
    output = RAW_DIR / tag / "deep_bsde" / f"{run_id}.json"
    resume = CHECKPOINT_DIR / tag / "deep_bsde" / f"{run_id}_resume.pt"
    final = CHECKPOINT_DIR / tag / "deep_bsde" / f"{run_id}.pt"
    prediction = PREDICTION_DIR / tag / "deep_bsde" / f"{run_id}.npz"
    if not (resume.is_file() and final.is_file() and prediction.is_file()):
        return None

    state = torch.load(resume, map_location="cpu", weights_only=False)
    config = dict(state["config"])
    if int(state.get("completed_updates", -1)) < int(config["iterations"]):
        return None
    model_state = state.get("model", {})
    y_init = model_state.get("y_init")
    if not isinstance(y_init, torch.Tensor):
        return None
    estimate = float(y_init.detach().cpu().reshape(-1)[0])
    if math.isfinite(estimate):
        return None
    with np.load(prediction, allow_pickle=False) as stored:
        stored_prediction = float(np.asarray(stored["prediction"]).reshape(-1)[0])
    if math.isfinite(stored_prediction):
        raise RuntimeError(f"checkpoint/prediction nonfinite mismatch: {run_id}")

    x0 = np.load(query_path(dimension) / "x.npy", allow_pickle=False)[point]
    reference = float(
        np.load(query_path(dimension) / "reference.npy", allow_pickle=False)[point]
    )
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "method": "deep_bsde_official_pytorch_compatibility_port",
        "status": "numerical_failure",
        "failure_reason": "nonfinite_prediction",
        "classification": "faithful_transfer",
        "pde": PDE,
        "dimension": dimension,
        "point": point,
        "seed": seed,
        "x0_sha256": array_sha256(np.asarray(x0, dtype=np.float64)),
        "configuration": config,
        "reported_y0": None,
        "reference": reference,
        "error": None,
        "absolute_error": None,
        "history": _json_safe(state.get("history", [])),
        "completed_updates": int(state["completed_updates"]),
        "nonfinite_count": int(state.get("nonfinite_count", 0)),
        "training_seconds": float(state.get("elapsed_seconds", 0.0)),
        "parameter_count": None,
        "peak_gpu_bytes": None,
        "checkpoint": str(final),
        "checkpoint_sha256": file_sha256(final),
        "resume_checkpoint": str(resume),
        "resume_checkpoint_sha256": file_sha256(resume),
        "validation_dw_sha256": state.get("validation_dw_sha256"),
        "prediction_file": str(prediction),
        "prediction_sha256": file_sha256(prediction),
        "selection": "fixed final scheduled iterate; no truth used",
        "native_scalar_y0": True,
        "amortized_y0": False,
        "terminal_aware_features": False,
        "mcv": False,
        "full_training_cost_charged_per_initial_state": True,
        "recovered_from_completed_artifacts": True,
        "recovery_note": (
            "Training completed before the strict JSON writer rejected NaN; "
            "this record was reconstructed without rerunning optimization."
        ),
    }
    write_json_atomic(output, record)
    return record


def execute(
    dimension: int,
    point: int,
    seed: int,
    device: torch.device,
    *,
    pilot: bool = False,
    stop_after_updates: int | None = None,
):
    _configure()
    tag = "pilot" if pilot else "formal"
    run_id = (
        f"deep_bsde__{PDE}__d{int(dimension)}__q{int(point):02d}__"
        f"s{int(seed)}__{tag}"
    )
    output = RAW_DIR / tag / "deep_bsde" / f"{run_id}.json"
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") == "numerical_failure":
            return existing
    recovered = _recover_completed_nonfinite(
        int(dimension), int(point), int(seed), pilot=pilot
    )
    if recovered is not None:
        return recovered
    # The frozen runner loads every checkpoint tensor with ``map_location``.
    # PyTorch requires a Generator state to be a CPU ByteTensor even when the
    # generator itself is CUDA-backed.  Preserve every RNG byte while fixing
    # only that device mismatch in this isolated adapter.
    original_torch_load = torch.load

    def load_with_cpu_generator_state(*args, **kwargs):
        state = original_torch_load(*args, **kwargs)
        if isinstance(state, dict) and isinstance(state.get("train_generator"), torch.Tensor):
            state["train_generator"] = state["train_generator"].cpu()
        return state

    torch.load = load_with_cpu_generator_state
    try:
        result = _frozen.execute(
            PDE,
            int(dimension),
            int(point),
            int(seed),
            device,
            pilot=pilot,
            stop_after_updates=stop_after_updates,
        )
        return _safe_deep_record(result)
    finally:
        torch.load = original_torch_load


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--point", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--stop-after-updates", type=int)
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    result = execute(
        args.dimension,
        args.point,
        args.seed,
        device,
        pilot=args.pilot,
        stop_after_updates=args.stop_after_updates,
    )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "completed_updates", "reported_y0", "absolute_error"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["FixedInitialDeepBSDE", "execute"]
