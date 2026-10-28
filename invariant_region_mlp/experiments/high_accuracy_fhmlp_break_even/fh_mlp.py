"""Isolated adapter for the audited unprojected Full-History MLP."""

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
from invariant_region_mlp.experiments.single_point_break_even import fh_mlp as _source
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    MLPStats,
)

from .benchmark import make_benchmark
from .protocol import (
    CALIBRATION_TASK_CAP_SECONDS,
    DIMENSION,
    FH_ACCURACY_REPLICATES,
    FH_CONFIGS,
    FH_TIMING_REPLICATES,
    PREDICTION_DIR,
    PROFILE_QUERIES,
    RAW_DIR,
    REVISION,
    array_sha256,
    ensure_directories,
    load_bank,
    require_cost_gate,
    require_frozen_timing_grid,
    require_preflight,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


AuditedFullHistoryMLP = _source.AuditedFullHistoryMLP
UnprojectedDriver = _source._UnprojectedDriver


def _configure(*, timing: bool = False) -> None:
    _source.DIMENSION = DIMENSION
    _source.FH_CONFIGS = FH_CONFIGS
    _source.FH_REPLICATES = FH_ACCURACY_REPLICATES
    _source.FH_TIMING_REPLICATES = FH_TIMING_REPLICATES
    _source.PREDICTION_DIR = PREDICTION_DIR
    _source.RAW_DIR = RAW_DIR
    _source.REVISION = REVISION
    if timing:
        frozen = require_frozen_timing_grid()
        _source.TIMING_COUNTS = tuple(int(value) for value in frozen["timing_counts"])
    _source.array_sha256 = array_sha256
    _source.ensure_directories = ensure_directories
    _source.load_bank = load_bank
    _source.make_benchmark = lambda _pde, dimension: make_benchmark(dimension)
    _source.stable_seed = stable_seed
    _source.write_json_atomic = write_json_atomic
    _source.write_npz_atomic = write_npz_atomic


def evaluate(
    stage: str,
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
) -> dict[str, Any]:
    require_preflight()
    gate = require_cost_gate()
    key = f"n{int(n)}_M{int(M)}"
    if key not in gate["eligible_fh_configs"]:
        raise RuntimeError(f"configuration is not cost-gate eligible: {key}")
    _configure()
    return _source.evaluate(stage, int(n), int(M), int(replicate), device)


def timing(
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
    *,
    max_active_seconds: float,
) -> dict[str, Any]:
    require_preflight()
    _configure(timing=True)
    return _source.timing(
        int(n), int(M), int(replicate), device,
        max_active_seconds=float(max_active_seconds),
    )


def profile(n: int, M: int, replicate: int, device: torch.device) -> dict[str, Any]:
    """Truth-free actual execution used only for resource forecasting."""
    require_preflight()
    _configure()
    bank = load_bank("timing", require_truth=False)
    benchmark = UnprojectedDriver(make_benchmark(DIMENSION))
    predictions = np.full(PROFILE_QUERIES, np.nan, dtype=np.float64)
    total = MLPStats()
    gpu = _source._require_h100(device)
    synchronize(device)
    torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    point_seconds: list[float] = []
    for point in range(PROFILE_QUERIES):
        local_started = time.perf_counter()
        predictions[point], stats = _source._solve_point(
            benchmark,
            float(bank["t"][point]),
            bank["x"][point],
            n=int(n),
            M=int(M),
            seed=stable_seed("fh", "profile", n, M, replicate, point),
            device=device,
        )
        synchronize(device)
        point_seconds.append(time.perf_counter() - local_started)
        total.merge(stats)
    elapsed = time.perf_counter() - started
    nonfinite = int((~np.isfinite(predictions)).sum()) + int(total.nonfinite_count)
    run_id = f"fh__profile__n{n}_M{M}__r{replicate}"
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "status": "complete" if nonfinite == 0 else "numerical_failure",
        "method": "fh_mlp_standard",
        "n": int(n),
        "M": int(M),
        "replicate": int(replicate),
        "queries": PROFILE_QUERIES,
        "active_seconds": elapsed,
        "point_seconds": point_seconds,
        "seconds_per_query": elapsed / PROFILE_QUERIES,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)),
        "nonfinite_count": nonfinite,
        "stats": asdict(total),
        "prediction_sha256": array_sha256(predictions),
        "query_input_sha256": array_sha256(
            bank["t"][:PROFILE_QUERIES], bank["x"][:PROFILE_QUERIES]
        ),
        "gpu": gpu,
        "truth_loaded": False,
        "exact_audited_recursion": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(RAW_DIR / "profile" / "fh_mlp" / f"{run_id}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    item = sub.add_parser("profile")
    item.add_argument("--n", type=int, required=True)
    item.add_argument("--M", type=int, required=True)
    item.add_argument("--replicate", type=int, required=True)
    item = sub.add_parser("evaluate")
    item.add_argument("--stage", choices=("calibration", "final"), required=True)
    item.add_argument("--n", type=int, required=True)
    item.add_argument("--M", type=int, required=True)
    item.add_argument("--replicate", type=int, required=True)
    item = sub.add_parser("timing")
    item.add_argument("--n", type=int, required=True)
    item.add_argument("--M", type=int, required=True)
    item.add_argument("--replicate", type=int, required=True)
    item.add_argument("--max-active-seconds", type=float, required=True)
    for item in sub.choices.values():
        item.add_argument("--device", default="cuda")
    args = parser.parse_args()
    device = torch.device(args.device)
    if args.action == "profile":
        result = profile(args.n, args.M, args.replicate, device)
    elif args.action == "evaluate":
        result = evaluate(args.stage, args.n, args.M, args.replicate, device)
    else:
        result = timing(
            args.n, args.M, args.replicate, device,
            max_active_seconds=args.max_active_seconds,
        )
    print(json.dumps({
        key: result.get(key)
        for key in ("run_id", "status", "rmse", "active_seconds", "completed_queries")
    }, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["AuditedFullHistoryMLP", "UnprojectedDriver", "evaluate", "profile", "timing"]
