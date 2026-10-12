"""Run or aggregate the prescribed Stage B grid."""

from __future__ import annotations

import argparse

from .aggregate import aggregate
from .common import load_config
from .runner import METHODS, run_method


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--d", type=int)
    parser.add_argument("--k", type=int)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.aggregate_only:
        aggregate("stage_b")
        return
    config = load_config()
    if args.all:
        # Preserve the specified method-wise execution order so every MLP
        # cell is checkpointed before actor training, followed by DPI, SDGD,
        # and the amortized Deep BSDE baseline.
        for method in METHODS:
            for d in config["pde"]["dimensions"]:
                for k in config["pde"]["effective_dimensions"]:
                    run_method("stage_b", method, d, k, force=args.force)
    elif args.method is not None and args.d is not None and args.k is not None:
        run_method("stage_b", args.method, args.d, args.k, force=args.force)
    else:
        parser.error("use --all or provide --method, --d, and --k")
    aggregate("stage_b")


if __name__ == "__main__":
    main()
