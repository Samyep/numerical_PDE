"""Run the data-only enrichment ablation with the frozen HCFL training code."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys

import torch


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "euler_2d_normal6_fixed_transverse"
SPEC = importlib.util.spec_from_file_location(
    "hcfl_normal6_fixed_training", BASE / "run_experiment.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(BASE / "run_experiment.py")
E = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = E
SPEC.loader.exec_module(E)


RESULTS = HERE / "results"
E.RESULTS = RESULTS


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--baseline-data-dir", type=Path, required=True)
    commands = parser.add_subparsers(dest="command", required=True)
    train = commands.add_parser("train")
    train.add_argument("--seed", type=int, required=True)
    train.add_argument("--max-updates", type=int, default=50000)
    train.add_argument("--validation-interval", type=int, default=500)
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    args = parser.parse_args()

    data = args.data_dir.resolve()
    baseline_data = args.baseline_data_dir.resolve()
    for name in ("train.npz", "validation.npz", "test.npz"):
        if not (data / name).is_file():
            raise FileNotFoundError(data / name)
    if not (baseline_data / "train.npz").is_file():
        raise FileNotFoundError(baseline_data / "train.npz")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(
        {"command": args.command, "device": str(device), "data": str(data)},
        flush=True,
    )
    if args.command == "train":
        E.train_model(
            args.seed,
            args.max_updates,
            args.validation_interval,
            device,
            data,
            statistics_archive=baseline_data / "train.npz",
        )
    else:
        E.evaluate(
            args.seeds,
            device,
            data,
            metric_scale_archive=baseline_data / "train.npz",
        )


if __name__ == "__main__":
    main()
