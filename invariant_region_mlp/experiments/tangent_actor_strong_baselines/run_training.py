"""CLI for calibration and compute-budget training runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from .common import BENCHMARKS, ensure_dirs
from .deep_picard import train_deep_picard
from .neural_baselines import train_deep_bsde_timed, train_sdgd_pinn_timed
from .tangent_baseline import train_tangent_method


METHODS = (
    "tangent_actor",
    "standard_tangent_actor",
    "deep_bsde",
    "sdgd_pinn",
    "deep_picard",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--benchmark", choices=tuple(BENCHMARKS), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--calibration-seconds", type=float)
    parser.add_argument("--calibration-rounds", type=int)
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("a CUDA allocation is required")
    if args.method == "deep_bsde":
        result = train_deep_bsde_timed(
            args.benchmark,
            args.seed,
            device,
            calibration_only_seconds=args.calibration_seconds,
        )
    elif args.method == "sdgd_pinn":
        result = train_sdgd_pinn_timed(
            args.benchmark,
            args.seed,
            device,
            calibration_only_seconds=args.calibration_seconds,
        )
    elif args.method == "deep_picard":
        result = train_deep_picard(
            args.benchmark,
            args.seed,
            device,
            calibration_rounds=args.calibration_rounds,
        )
    elif args.method in {"tangent_actor", "standard_tangent_actor"}:
        iterations = args.calibration_rounds or 3
        result = train_tangent_method(
            args.method, args.benchmark, args.seed, device, iterations=iterations
        )
    else:  # pragma: no cover
        raise KeyError(args.method)
    print(
        json.dumps(
            {
                "method": args.method,
                "benchmark": args.benchmark,
                "seed": args.seed,
                "status": "complete",
                "records": len(result.get("checkpoints", result.get("phases", []))),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
