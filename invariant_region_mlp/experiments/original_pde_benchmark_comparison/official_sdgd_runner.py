"""Run the pinned SDGD HJB-Log script verbatim and export its final model."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import runpy
import sys
import time

import torch

from .common import OFFICIAL_WORK_DIR, RAW_DIR, atomic_json, ensure_dirs


SOURCE = Path(
    "/work/nvme/bgvf/ywang151/original_pde_source_audit/SDGD_PINN/Section_5.3/HJB_Log.py"
)


def run(seed: int, calibration: bool) -> dict[str, object]:
    ensure_dirs()
    work = OFFICIAL_WORK_DIR / "sdgd"
    work.mkdir(parents=True, exist_ok=True)
    (work / "saved_loss_l2").mkdir(exist_ok=True)
    epochs = 100 if calibration else 10_000
    old_argv = sys.argv
    sys.argv = [
        str(SOURCE),
        "--SEED", str(seed),
        "--epochs", str(epochs),
        "--device", "cuda",
    ]
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
    started = time.perf_counter()
    previous = Path.cwd()
    os.chdir(work)
    try:
        namespace = runpy.run_path(str(SOURCE), run_name="__main__")
    finally:
        os.chdir(previous)
        sys.argv = old_argv
    if torch.cuda.is_available():
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - started
    model = namespace["model"].u_net
    suffix = "calibration" if calibration else "official"
    checkpoint = work / f"sdgd_hjb_log_s{seed}_{suffix}.pt"
    torch.save(model.state_dict(), checkpoint)
    spreadsheets = sorted(
        (work / "saved_loss_l2").glob(f"*SEED={seed}.xlsx"),
        key=lambda path: path.stat().st_mtime,
    )
    record: dict[str, object] = {
        "method": "sdgd_pinn_official",
        "benchmark": "hjb_log_d100",
        "seed": seed,
        "calibration": calibration,
        "epochs": epochs,
        "training_seconds_including_official_per_epoch_test_reporting": elapsed,
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated())
        if torch.cuda.is_available()
        else None,
        "checkpoint": str(checkpoint),
        "trajectory": str(spreadsheets[-1]) if spreadsheets else None,
        "official_source": str(SOURCE),
        "official_source_commit": "c597cc4fe9caac0b3ca118ea84453405f25ff833",
        "source_executed_with_runpy_unmodified": True,
        "postrun_instrumentation": "export final u_net state_dict and wall/memory metadata",
        "published_defaults": {
            "width": 1024,
            "hidden_layers": 3,
            "activation": "tanh",
            "collocation_points": 100,
            "pgd_steps_per_epoch": 20,
            "pgd_step_size": 0.05,
            "dimension_samples": 10,
            "optimizer": "Adam",
            "learning_rate": 0.001,
            "linear_lr_schedule": True,
            "hard_terminal_ansatz": True,
            "method": 2,
            "adversarial_method": 2,
        },
        "terminal_aware_feature_added_by_us": False,
        "mcv_added_by_us": False,
        "status": "complete",
    }
    atomic_json(RAW_DIR / f"official_sdgd_s{seed}_{suffix}.json", record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--calibration", action="store_true")
    args = parser.parse_args()
    result = run(args.seed, args.calibration)
    print(json.dumps({"status": "complete", "seconds": result[
        "training_seconds_including_official_per_epoch_test_reporting"
    ]}))


if __name__ == "__main__":
    main()

