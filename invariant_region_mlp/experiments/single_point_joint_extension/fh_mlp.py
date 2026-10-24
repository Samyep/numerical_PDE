"""Isolated adapter for the completed suite's audited FH-MLP solver."""

from __future__ import annotations

import argparse
import json

import torch

from invariant_region_mlp.experiments.single_point_32_ice import fh_mlp as _frozen

from .benchmark import make_benchmark
from .protocol import (
    FH_CONFIGURATIONS,
    FH_REPLICATES,
    PDE,
    PREDICTION_DIR,
    QUERY_COUNT,
    RAW_DIR,
    REVISION,
    array_sha256,
    ensure_directories,
    query_path,
    require_preflight,
    require_smoke,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


AuditedFullHistoryMLP = _frozen.AuditedFullHistoryMLP
UnprojectedDriver = _frozen._UnprojectedDriver


def _configure() -> None:
    # The completed implementation reads these names from its own module at
    # call time.  Binding them here preserves the audited recursion while
    # isolating all artifacts and gates in the extension directory.
    _frozen.make_benchmark = lambda pde, d: make_benchmark(d)
    _frozen.FH_CONFIGURATIONS = FH_CONFIGURATIONS
    _frozen.FH_REPLICATES = FH_REPLICATES
    _frozen.PREDICTION_DIR = PREDICTION_DIR
    _frozen.QUERY_COUNT = QUERY_COUNT
    _frozen.RAW_DIR = RAW_DIR
    _frozen.REVISION = REVISION
    _frozen.array_sha256 = array_sha256
    _frozen.ensure_directories = ensure_directories
    _frozen.query_path = lambda pde, d: query_path(d)
    _frozen.require_stage_a = require_preflight
    _frozen.require_stage_b = require_smoke
    _frozen.stable_seed = stable_seed
    _frozen.write_json_atomic = write_json_atomic
    _frozen.write_npz_atomic = write_npz_atomic


def execute(
    dimension: int,
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
    *,
    pilot: bool = False,
):
    _configure()
    return _frozen.execute(
        PDE, int(dimension), int(n), int(M), int(replicate), device, pilot=pilot
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--M", type=int, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    result = execute(
        args.dimension, args.n, args.M, args.replicate, device, pilot=args.pilot
    )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "finite_count", "nonfinite_count", "inference_seconds"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["AuditedFullHistoryMLP", "UnprojectedDriver", "execute"]
