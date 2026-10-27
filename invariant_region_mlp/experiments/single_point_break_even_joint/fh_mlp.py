"""Isolated adapter for the audited Full-History MLP recursion."""

from __future__ import annotations

import argparse
import json

import torch

from invariant_region_mlp.experiments.single_point_break_even import fh_mlp as _source

from .benchmark import make_benchmark
from .protocol import (
    ACTIVE_CAP_SECONDS,
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
    require_preflight,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


AuditedFullHistoryMLP = _source.AuditedFullHistoryMLP
UnprojectedDriver = _source._UnprojectedDriver


def _configure() -> None:
    _source.DIMENSION = DIMENSION
    _source.FH_CONFIGS = FH_CONFIGS
    _source.FH_REPLICATES = FH_REPLICATES
    _source.FH_TIMING_REPLICATES = FH_TIMING_REPLICATES
    _source.PREDICTION_DIR = PREDICTION_DIR
    _source.RAW_DIR = RAW_DIR
    _source.REVISION = REVISION
    _source.TIMING_COUNTS = TIMING_COUNTS
    _source.array_sha256 = array_sha256
    _source.ensure_directories = ensure_directories
    _source.load_bank = load_bank
    _source.make_benchmark = lambda _pde, dimension: make_benchmark(dimension)
    _source.stable_seed = stable_seed
    _source.write_json_atomic = write_json_atomic
    _source.write_npz_atomic = write_npz_atomic


def evaluate(stage: str, n: int, M: int, replicate: int, device: torch.device):
    require_preflight()
    _configure()
    return _source.evaluate(stage, int(n), int(M), int(replicate), device)


def timing(n: int, M: int, replicate: int, device: torch.device, *, max_active_seconds: float = ACTIVE_CAP_SECONDS):
    require_preflight()
    _configure()
    return _source.timing(int(n), int(M), int(replicate), device, max_active_seconds=max_active_seconds)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("evaluate"); p.add_argument("--stage", choices=("calibration", "final"), required=True); p.add_argument("--n", type=int, required=True); p.add_argument("--M", type=int, required=True); p.add_argument("--replicate", type=int, required=True)
    p = sub.add_parser("timing"); p.add_argument("--n", type=int, required=True); p.add_argument("--M", type=int, required=True); p.add_argument("--replicate", type=int, required=True); p.add_argument("--max-active-seconds", type=float, default=ACTIVE_CAP_SECONDS)
    for item in sub.choices.values(): item.add_argument("--device", default="cuda")
    args = parser.parse_args(); device = torch.device(args.device)
    if args.action == "evaluate": result = evaluate(args.stage, args.n, args.M, args.replicate, device)
    else: result = timing(args.n, args.M, args.replicate, device, max_active_seconds=args.max_active_seconds)
    print(json.dumps({key: result.get(key) for key in ("run_id", "status", "rmse", "active_seconds", "completed_queries")}, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["AuditedFullHistoryMLP", "UnprojectedDriver", "evaluate", "timing"]
