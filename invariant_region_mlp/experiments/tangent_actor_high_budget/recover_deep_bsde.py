"""Recover metadata for completed Deep-BSDE checkpoints after JSON failure.

The original-form solver can become non-finite during its rollout while its
raw amortized initial-value network remains finite.  The timed checkpoint files
are written before the final strict-JSON serialization.  This utility records
those immutable checkpoints without retraining and represents unavailable
non-finite validation objectives as JSON null.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import torch

from .common import BUDGETS, CHECKPOINT_DIR, RAW_DIR, atomic_json


def _partial_diagnostic(benchmark: str, seed: int) -> dict[str, Any]:
    pattern = re.compile(
        rf'"benchmark": "{re.escape(benchmark)}".*?'
        rf'"nonfinite_count": (?P<nonfinite>\d+).*?'
        rf'"peak_gpu_bytes": (?P<peak>\d+).*?'
        rf'"seed": {seed}.*?'
        rf'"training_seconds": (?P<seconds>[0-9.eE+-]+)',
        re.DOTALL,
    )
    for path in sorted(RAW_DIR.glob("tmp*")):
        text = path.read_text(encoding="utf-8", errors="replace")
        match = pattern.search(text)
        if match:
            return {
                "partial_json_path": str(path),
                "first_checkpoint_nonfinite_count": int(match.group("nonfinite")),
                "peak_gpu_bytes": int(match.group("peak")),
                "first_checkpoint_training_seconds": float(match.group("seconds")),
            }
    raise FileNotFoundError(
        f"no failed strict-JSON fragment found for {benchmark} seed {seed}"
    )


def recover(benchmark: str, seed: int) -> dict[str, Any]:
    diagnostic = _partial_diagnostic(benchmark, seed)
    paths = [
        CHECKPOINT_DIR / f"deep_bsde_original_{benchmark}_s{seed}_b{int(target)}.pt"
        for target in BUDGETS
    ]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing completed checkpoint(s): {missing}")
    first_seconds = diagnostic["first_checkpoint_training_seconds"]
    inferred_start = paths[0].stat().st_mtime - first_seconds
    checkpoints = []
    for target, path in zip(BUDGETS, paths):
        payload = torch.load(path, map_location="cpu", weights_only=False)
        state = payload["state_dict"]
        parameter_count = sum(tensor.numel() for tensor in state.values())
        elapsed = max(float(target) + 1e-6, path.stat().st_mtime - inferred_start)
        learning_rate = 5e-4 * (
            1.0 if elapsed < 480.0 else (0.2 if elapsed < 800.0 else 0.04)
        )
        checkpoints.append(
            {
                "method": "deep_bsde_original", "benchmark": benchmark,
                "seed": seed, "target_budget_seconds": target,
                "training_seconds": elapsed, "step": int(payload["step"]),
                "validation_objective": None, "checkpoint": str(path),
                "learning_rate": learning_rate, "time_steps": 16,
                "batch_size": 128, "parameter_count": parameter_count,
                "peak_gpu_bytes": diagnostic["peak_gpu_bytes"],
                "nonfinite_count": None,
                "selection_criterion": "predeclared elapsed-time checkpoint",
                "final_test_used_for_selection": False,
                "metadata_recovered_without_retraining": True,
                "training_seconds_resolution_seconds": 1.0 if target != BUDGETS[0] else 0.0,
            }
        )
    result = {
        "method": "deep_bsde_original", "benchmark": benchmark, "seed": seed,
        "status": "completed_checkpoints_recovered_after_strict_json_failure",
        "configuration": {
            "adaptation": "raw amortized Y0(x0); official-style nonshared BN/ReLU Z_n",
            "terminal_aware_features": False, "time_steps": 16,
            "batch_size": 128, "learning_rate": 5e-4,
            "z_width": int(torch.load(paths[0], map_location="cpu", weights_only=False)["architecture"]["width"]),
            "validation_size": 512,
            "official_source_commit": "e76fed80b1995daf4a0dd6d3e2cf64a0931dabda",
        },
        "failure_diagnostic": {
            "failure_stage": "final strict JSON serialization",
            "cause": "non-finite Deep-BSDE rollout diagnostics",
            "training_repeated": False,
            "validation_objective_available": False,
            "first_checkpoint_nonfinite_count_lower_bound": diagnostic["first_checkpoint_nonfinite_count"],
            "checkpoint_training_seconds_recovered_from_file_times": True,
            "later_checkpoint_time_resolution_seconds": 1.0,
            "partial_json_path": diagnostic["partial_json_path"],
        },
        "history": [], "checkpoints": checkpoints,
    }
    output = RAW_DIR / f"deep_bsde_original_{benchmark}_s{seed}.json"
    atomic_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--seed", type=int, choices=(0, 1), required=True)
    args = parser.parse_args()
    result = recover(args.benchmark, args.seed)
    print(json.dumps({
        "status": result["status"], "benchmark": args.benchmark,
        "seed": args.seed, "checkpoints": len(result["checkpoints"]),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
