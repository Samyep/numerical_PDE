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


def run_official(
    method: str,
    seed: int,
    calibration: bool,
    *,
    full_data_one_iteration: bool = False,
) -> dict[str, object]:
    """Run upstream classes without changing any loss or network operation."""
    if calibration and full_data_one_iteration:
        raise ValueError("calibration modes are mutually exclusive")
    if full_data_one_iteration and method != "dpi":
        raise ValueError("full-data one-iteration timing is only defined for DPI")
    ensure_dirs()
    # Imports intentionally occur here so PYTHONPATH can point at the pinned
    # checkout and its isolated compatibility dependencies in the Slurm job.
    import lightning as pl
    from lightning.pytorch.callbacks import RichProgressBar
    from picard import load_cfg, PicardRunner
    import picard.solution as solutions

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    pl.seed_everything(seed, workers=True)
    suffix = (
        "full_data_one_iteration"
        if full_data_one_iteration
        else "calibration"
        if calibration
        else "official"
    )
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

    # The upstream custom Rich progress callback crashes before the first
    # training batch under a noninteractive Slurm stream with the installed
    # Rich/Lightning combination (``Console.clear_live(): pop from empty
    # list``).  Remove only that display callback.  Evaluation callbacks,
    # model summaries, loggers, data, losses, and optimizer steps are kept.
    original_reporting_callbacks = PicardRunner.get_reporting_callbacks

    def reporting_callbacks_without_rich(self):
        callbacks, trainer_kwargs = original_reporting_callbacks(self)
        callbacks = [
            callback
            for callback in callbacks
            if not isinstance(callback, RichProgressBar)
        ]
        trainer_kwargs["enable_progress_bar"] = False
        return callbacks, trainer_kwargs

    PicardRunner.get_reporting_callbacks = reporting_callbacks_without_rich

    overrides = ["NAME", str(exp_dir), "FORCE", "True"]
    launcher_correction = None
    if method == "d_dbsde":
        # The derived Diffusion YAML inherits PICARD.N=40 from the DPI base,
        # although DiffusionSolution is a standalone 15,000-epoch method and
        # run_one() does not update Picard state.  Leaving N=40 trains 40 fresh
        # independent D-DBSDE networks and discards the first 39.  The README
        # invokes this YAML as one D-DBSDE run, so execute that native solver
        # once.  This changes no loss, architecture, epoch count, or output.
        overrides += ["PICARD.N", "1"]
        launcher_correction = (
            "set inherited PICARD.N from 40 to 1 for standalone DiffusionSolution; "
            "otherwise upstream PicardRunner redundantly trains and discards 39 full networks"
        )
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
    elif full_data_one_iteration:
        # Runtime calibration with every official data, MC-label, network,
        # optimizer, and evaluation setting intact.  Only the outer Picard
        # count is shortened from 40 to 1 so that a faithful per-iteration
        # cost is measured before reserving the full run.
        overrides += ["PICARD.N", "1"]
    # Upstream derived YAML files store ``BASE`` as a path relative to the
    # official scripts/hjb directory.  Resolve the configuration under that
    # native working directory, then return to the prepared GMM-tensor
    # directory before constructing the equation.  No config value changes.
    input_directory = Path.cwd()
    os.chdir(CONFIGS[method].parent)
    try:
        cfg = load_cfg(str(CONFIGS[method]), overrides)
    finally:
        os.chdir(input_directory)
    cfg.freeze()
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
    started = time.perf_counter()
    runner = PicardRunner(cfg)
    # The repository was released before PyTorch 2.6, when torch.load used
    # ``weights_only=False`` by default.  Picard checkpoints intentionally
    # contain the equation object, so iteration 2 cannot read iteration 1
    # under the new default.  Restore the legacy default only for calls where
    # upstream omitted the keyword.  These checkpoints are generated locally
    # in this same pinned-source process.
    original_torch_load = torch.load

    def torch_load_with_legacy_default(*args, **kwargs):
        kwargs.setdefault("weights_only", False)
        return original_torch_load(*args, **kwargs)

    torch.load = torch_load_with_legacy_default
    try:
        runner.run()
    finally:
        torch.load = original_torch_load
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
        "full_data_one_iteration": full_data_one_iteration,
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
            "fixed RNG seeds, wall/memory accounting, final state_dict export, and "
            "non-numerical Rich progress display disabled; PyTorch >=2.6 "
            "checkpoint loading restored to its legacy weights_only=False default; "
            "upstream architecture, terminal construction, loss, optimizer, labels, and iterations unchanged"
        ),
        "launcher_compatibility_correction": launcher_correction,
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
    parser.add_argument("--full-data-one-iteration", action="store_true")
    args = parser.parse_args()
    record = run_official(
        args.method,
        args.seed,
        args.calibration,
        full_data_one_iteration=args.full_data_one_iteration,
    )
    print(json.dumps({"status": "complete", "seconds": record["training_seconds"]}))


if __name__ == "__main__":
    main()

