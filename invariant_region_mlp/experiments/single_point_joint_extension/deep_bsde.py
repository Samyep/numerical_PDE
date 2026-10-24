"""Isolated adapter for the completed suite's checkpointable Deep BSDE."""

from __future__ import annotations

import argparse
import json

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
    _frozen.write_json_atomic = write_json_atomic
    _frozen.write_npz_atomic = write_npz_atomic


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
    return _frozen.execute(
        PDE,
        int(dimension),
        int(point),
        int(seed),
        device,
        pilot=pilot,
        stop_after_updates=stop_after_updates,
    )


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
