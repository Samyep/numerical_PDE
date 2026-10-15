"""CLI entry point for staged high-budget training."""

from __future__ import annotations

import argparse
import json

import torch

from .actor_scaling import BRANCHES, train_branch, train_pilot_suite
from .common import BENCHMARKS, ensure_dirs
from .original_baselines import train_deep_bsde_original, train_sdgd_pinn_original


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--method",
        choices=("actor_pilot", "actor_branch", "deep_bsde_original", "sdgd_pinn_original"),
        required=True,
    )
    parser.add_argument("--benchmark", choices=tuple(BENCHMARKS), required=True)
    parser.add_argument("--seed", type=int, choices=(0, 1), required=True)
    parser.add_argument("--branch", choices=tuple(BRANCHES))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    ensure_dirs()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("a CUDA allocation is required")
    if args.method == "actor_pilot":
        result = train_pilot_suite(args.benchmark, args.seed, device)
    elif args.method == "actor_branch":
        if args.branch is None:
            parser.error("--branch is required for actor_branch")
        result = train_branch(args.benchmark, args.seed, args.branch, device)
    elif args.method == "deep_bsde_original":
        result = train_deep_bsde_original(args.benchmark, args.seed, device)
    else:
        result = train_sdgd_pinn_original(args.benchmark, args.seed, device)
    print(
        json.dumps(
            {
                "status": "complete", "method": args.method,
                "benchmark": args.benchmark, "seed": args.seed,
                "records": len(result) if isinstance(result, list)
                else len(result.get("phases", result.get("checkpoints", []))),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

