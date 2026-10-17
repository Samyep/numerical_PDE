"""Instrumentation-only launcher for the pinned official DPI repository."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import random
import time

import numpy as np
import torch

from .common import OFFICIAL_WORK_DIR, RAW_DIR, atomic_json, ensure_dirs


DPI_SOURCE = Path(
    "/work/nvme/bgvf/ywang151/original_pde_source_audit/DeepPicardIteration"
)
CONFIGS = {
    "dpi": DPI_SOURCE / "scripts/hjb/base_100d_T1.0_w0.1_0.yaml",
    "pinn_hte": DPI_SOURCE / "scripts/hjb/pinn_100d_T1.0_v16.yaml",
    "d_dbsde": DPI_SOURCE / "scripts/hjb/diffusion_100d_T1.0.yaml",
}


def run_official(method: str, seed: int, calibration: bool) -> dict[str, object]:
    """Run upstream classes without changing any loss or network operation."""
    ensure_dirs()
    # Imports intentionally occur here so PYTHONPATH can point at the pinned
    # checkout and its isolated compatibility dependencies in the Slurm job.
    import lightning as pl
    from picard import load_cfg, PicardRunner
    import picard.solution as solutions

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    pl.seed_everything(seed, workers=True)
    suffix = "calibration" if calibration else "official"
    exp_dir = OFFICIAL_WORK_DIR / "dpi_runs" / f"{method}_s{seed}_{suffix}"
    capture = exp_dir.parent / f"{exp_dir.name}_final_state.pt"
    exp_dir.parent.mkdir(parents=True, exist_ok=True)

    # PINN and Diffusion training are upstream hand-written loops.  Upstream
    # saves periodic PINN checkpoints but no final Diffusion checkpoint.  The
    # wrappers below call the untouched loop and then export the final dnn;
    # they do not alter optimization or data generation.
    original_pinn_train = solutions.PINNSolution.train
    original_diffusion_train = solutions.DiffusionSolution.train

    def pinn_train_and_capture(self):
        result = original_pinn_train(self)
        torch.save(self.dnn.state_dict(), capture)
        return result

    def diffusion_train_and_capture(self):
        result = original_diffusion_train(self)
        torch.save(self.dnn.state_dict(), capture)
        return result

    solutions.PINNSolution.train = pinn_train_and_capture
    solutions.DiffusionSolution.train = diffusion_train_and_capture

    overrides = ["NAME", str(exp_dir), "FORCE", "True"]
    if calibration:
        if method == "dpi":
            overrides += [
                "PICARD.N", "1",
                "DATA.DATA_SIZE", "256",
                "DATA.kwargs.n_estimate_terminal", "128",
                "DATA.kwargs.n_estimate_integral", "128",
                "TRAIN.N_EPOCHS", "2",
                "TRAIN.BATCH_SIZE", "128",
                "EVAL.FREQ", "None",
            ]
        else:
            overrides += ["TRAIN.N_EPOCHS", "100", "EVAL.FREQ", "50"]
    cfg = load_cfg(str(CONFIGS[method]), overrides)
    cfg.freeze()
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
    started = time.perf_counter()
    runner = PicardRunner(cfg)
    runner.run()
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - started
    if method == "dpi":
        final_model = exp_dir / f"model_{cfg.PICARD.N}.ckpt"
    else:
        final_model = capture
    record: dict[str, object] = {
        "method": method,
        "seed": seed,
        "calibration": calibration,
        "official_source_commit": "e8616c2f531c9f95299b246a13c1180b8f02822a",
        "config_source": str(CONFIGS[method]),
        "expanded_config": str(exp_dir / "config.yaml"),
        "experiment_directory": str(exp_dir),
        "final_model": str(final_model),
        "training_seconds": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated())
        if torch.cuda.is_available()
        else None,
        "instrumentation": (
            "fixed RNG seeds, wall/memory accounting, and final state_dict export only; "
            "upstream architecture, terminal construction, loss, optimizer, labels, and iterations unchanged"
        ),
        "overrides": overrides,
        "terminal_aware_feature_added_by_us": False,
        "mcv_added_by_us": False,
        "status": "complete",
    }
    atomic_json(RAW_DIR / f"official_dpi_{method}_s{seed}_{suffix}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=tuple(CONFIGS), required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--calibration", action="store_true")
    args = parser.parse_args()
    record = run_official(args.method, args.seed, args.calibration)
    print(json.dumps({"status": "complete", "seconds": record["training_seconds"]}))


if __name__ == "__main__":
    main()

