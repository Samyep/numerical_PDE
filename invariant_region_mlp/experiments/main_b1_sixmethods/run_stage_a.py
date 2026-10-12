"""Run or aggregate the prescribed Stage A grid."""

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
        aggregate("stage_a")
        return
    if args.all:
        for d, k in load_config()["stage_a"]["cells"]:
            for method in METHODS:
                run_method("stage_a", method, d, k, force=args.force)
    elif args.method is not None and args.d is not None and args.k is not None:
        run_method("stage_a", args.method, args.d, args.k, force=args.force)
    else:
        parser.error("use --all or provide --method, --d, and --k")
    aggregate("stage_a")


if __name__ == "__main__":
    main()
