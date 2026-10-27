"""Accuracy and actual-prefix timing for the audited Full-History MLP."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    synchronize,
)
from invariant_region_mlp.experiments.single_point_32_ice.benchmarks import (
    make_benchmark,
)
from invariant_region_mlp.experiments.single_point_32_ice.fh_mlp import (
    AuditedFullHistoryMLP,
    _UnprojectedDriver,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    MLPStats,
)

from .protocol import (
    DIMENSION,
    FH_CONFIGS,
    FH_REPLICATES,
    FH_TIMING_REPLICATES,
    PREDICTION_DIR,
    RAW_DIR,
    REVISION,
    TIMING_COUNTS,
    array_sha256,
    ensure_directories,
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


def _solve_point(
    benchmark,
    t: float,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    seed: int,
    device: torch.device,
) -> tuple[float, MLPStats]:
    solver = AuditedFullHistoryMLP(
        benchmark,
        M,
        device,
        seed=seed,
        dtype=torch.float64,
        time_beta_alpha=0.5,
    )
    with torch.no_grad():
        output = solver.solve(
            n,
            torch.as_tensor([t], dtype=torch.float64, device=device),
            torch.as_tensor(x[None, :], dtype=torch.float64, device=device),
        )
    return float(output[0, 0].detach().cpu()), solver.stats


def _merge(total: MLPStats, local: MLPStats) -> None:
    total.merge(local)


def evaluate(
    stage: str,
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
) -> dict[str, Any]:
    if stage not in {"calibration", "final"}:
        raise KeyError(stage)
    if (n, M) not in FH_CONFIGS:
        raise ValueError((n, M))
    if not 0 <= replicate < FH_REPLICATES:
        raise ValueError(replicate)
    ensure_directories()
    gpu = _require_h100(device)
    bank = load_bank(stage, require_truth=True)
    benchmark = _UnprojectedDriver(make_benchmark("allen_cahn", DIMENSION))
    prediction = np.full(len(bank["t"]), np.nan, dtype=np.float64)
    point_seconds = np.zeros(len(prediction), dtype=np.float64)
    total = MLPStats()
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for point in range(len(prediction)):
        point_started = time.perf_counter()
        prediction[point], stats = _solve_point(
            benchmark,
            float(bank["t"][point]),
            bank["x"][point],
            n=n,
            M=M,
            seed=stable_seed("fh", stage, n, M, replicate, point),
            device=device,
        )
        point_seconds[point] = time.perf_counter() - point_started
        _merge(total, stats)
    synchronize(device)
    elapsed = time.perf_counter() - started
    error = prediction - bank["reference"]
    result_id = f"fh__{stage}__n{n}_M{M}__r{replicate}"
    prediction_path = PREDICTION_DIR / stage / "fh_mlp" / f"{result_id}.npz"
    write_npz_atomic(
        prediction_path,
        prediction=prediction,
        reference=bank["reference"],
        point_seconds=point_seconds,
    )
    nonfinite = int((~np.isfinite(prediction)).sum()) + int(total.nonfinite_count)
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": "complete" if nonfinite == 0 else "numerical_failure",
        "stage": stage,
        "method": "fh_mlp_standard",
        "n": n,
        "M": M,
        "replicate": replicate,
        "points": len(prediction),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
        "inference_seconds": elapsed,
        "point_seconds": point_seconds.tolist(),
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "stats": asdict(total),
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"], bank["x"]),
        "gpu": gpu,
        "dtype": "float64",
        "time_beta_alpha": 0.5,
        "terminal_gradient_estimator": "pathwise sigma*grad(g)",
        "driver_gradient_estimator": "untruncated Bismut weight",
        "high_low_subtrees": "independent",
        "truth_used_inside_solver": False,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RAW_DIR / stage / "fh_mlp" / f"{result_id}.json", record)
    return record


def timing(
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
    *,
    max_active_seconds: float = 3600.0,
) -> dict[str, Any]:
    if (n, M) not in FH_CONFIGS:
        raise ValueError((n, M))
    if not 0 <= replicate < FH_TIMING_REPLICATES:
        raise ValueError(replicate)
    gpu = _require_h100(device)
    bank = load_bank("timing", require_truth=False)
    if "reference" in bank:
        raise RuntimeError("truth leaked into timing bank")
    benchmark = _UnprojectedDriver(make_benchmark("allen_cahn", DIMENSION))
    prediction = np.full(len(bank["t"]), np.nan, dtype=np.float64)
    total = MLPStats()
    rows: list[dict[str, Any]] = []
    next_prefix = 0
    completed = 0
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for point in range(len(prediction)):
        prediction[point], stats = _solve_point(
            benchmark,
            float(bank["t"][point]),
            bank["x"][point],
            n=n,
            M=M,
            seed=stable_seed("fh", "timing", n, M, replicate, point),
            device=device,
        )
        _merge(total, stats)
        completed = point + 1
        if completed == TIMING_COUNTS[next_prefix]:
            synchronize(device)
            active = time.perf_counter() - started
            rows.append({
                "queries": completed,
                "inference_seconds": active,
                "total_seconds": active,
            })
            next_prefix += 1
            if next_prefix == len(TIMING_COUNTS) or active >= max_active_seconds:
                break
    prediction = prediction[:completed]
    elapsed = rows[-1]["inference_seconds"] if rows else 0.0
    status = "complete" if completed == TIMING_COUNTS[-1] else "TIME_CAP"
    result_id = f"fh__timing__n{n}_M{M}__r{replicate}"
    prediction_path = PREDICTION_DIR / "timing" / "fh_mlp" / f"{result_id}.npz"
    write_npz_atomic(prediction_path, prediction=prediction)
    nonfinite = int((~np.isfinite(prediction)).sum()) + int(total.nonfinite_count)
    record = {
        "revision": REVISION,
        "run_id": result_id,
        "status": "numerical_failure" if nonfinite else status,
        "method": "fh_mlp_standard",
        "n": n,
        "M": M,
        "replicate": replicate,
        "prefix_rows": rows,
        "completed_queries": completed,
        "inference_seconds": elapsed,
        "active_seconds": elapsed,
        "nonfinite_count": nonfinite,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "stats": asdict(total),
        "prediction_file": str(prediction_path),
        "prediction_sha256": array_sha256(prediction),
        "query_input_sha256": array_sha256(bank["t"][:completed], bank["x"][:completed]),
        "gpu": gpu,
        "truth_loaded": False,
        "actual_prefix_execution": True,
    }
    write_json_atomic(RAW_DIR / "timing" / "fh_mlp" / f"{result_id}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    evaluate_parser = subparsers.add_parser("evaluate")
    evaluate_parser.add_argument("--stage", choices=("calibration", "final"), required=True)
    evaluate_parser.add_argument("--n", type=int, required=True)
    evaluate_parser.add_argument("--M", type=int, required=True)
    evaluate_parser.add_argument("--replicate", type=int, required=True)
    timing_parser = subparsers.add_parser("timing")
    timing_parser.add_argument("--n", type=int, required=True)
    timing_parser.add_argument("--M", type=int, required=True)
    timing_parser.add_argument("--replicate", type=int, required=True)
    timing_parser.add_argument("--max-active-seconds", type=float, default=3600.0)
    for item in (evaluate_parser, timing_parser):
        item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if args.action == "evaluate":
        result = evaluate(args.stage, args.n, args.M, args.replicate, device)
    else:
        result = timing(
            args.n,
            args.M,
            args.replicate,
            device,
            max_active_seconds=args.max_active_seconds,
        )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "rmse", "active_seconds", "completed_queries",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

