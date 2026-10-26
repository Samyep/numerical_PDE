"""Post-freeze RNG replay of the historical Dense-Fourier DPI callback.

This is deliberately not a frontier trajectory.  The historical reporting
callback sampled 10,000 states from the global Torch RNG before evaluating the
exact solution.  Here we replay only those random draws, at the same callback
frequency, without evaluating exact truth.  This isolates whether the changed
RNG trajectory explains the failed 80-round historical-consistency check.  Its
prediction is never merged into ``checkpoint_metrics.csv`` or any paper
frontier.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import platform
import time
from types import MethodType

import numpy as np
import torch

from .baseline_frontier import (
    _bank,
    _configuration,
    _fresh_runner,
    _model_for_checkpoint,
    _seed_all,
)
from .common import file_sha256, save_npz_atomic, synchronize, value_metrics, write_json_atomic
from .prepare import historical_results_root
from .protocol import RESULTS_DIR, REVISION, TRAINING_SEED, WORK_DIR, ensure_directories


DIAGNOSTIC_DIR = RESULTS_DIR / "diagnostics"
RESULT_PATH = DIAGNOSTIC_DIR / "dense_dpi_callback_replay.json"
PREDICTION_PATH = DIAGNOSTIC_DIR / "dense_dpi_callback_replay_predictions.npz"


def _predict(model, t: np.ndarray, x: np.ndarray, device: torch.device) -> np.ndarray:
    pieces: list[np.ndarray] = []
    synchronize(device)
    with torch.no_grad():
        for begin in range(0, len(t), 128):
            end = min(begin + 128, len(t))
            tt = torch.as_tensor(t[begin:end, None], dtype=torch.float64, device=device)
            xx = torch.as_tensor(x[begin:end], dtype=torch.float64, device=device)
            pieces.append(
                model(torch.cat((tt, xx), dim=-1))[:, 0].detach().cpu().numpy()
            )
    synchronize(device)
    return np.concatenate(pieces).astype(np.float64)


def run() -> dict:
    ensure_directories()
    DIAGNOSTIC_DIR.mkdir(parents=True, exist_ok=True)
    if RESULT_PATH.is_file():
        existing = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
        if existing.get("status") == "COMPLETE":
            return existing

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("the callback replay requires the formal CUDA environment")

    historical_root = historical_results_root()
    historical_raw_path = (
        historical_root / "raw"
        / "dpi_dense_fourier_d100_semilinear_paper_official_s0.json"
    )
    historical = json.loads(historical_raw_path.read_text(encoding="utf-8"))
    historical_prediction_path = Path(historical["prediction_file"])
    historical_config_path = Path(historical["expanded_config"])
    historical_checkpoint_path = Path(historical["checkpoint"])
    for required in (
        historical_prediction_path, historical_config_path, historical_checkpoint_path,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)

    _seed_all(TRAINING_SEED)
    cfg, _source, _overrides, _profile = _configuration("dpi", "dense_fourier", "C0")
    cfg = cfg.clone()
    cfg.defrost()
    audit_work = WORK_DIR / "diagnostic_dense_dpi_rng_replay_epoch_hook_s0"
    if audit_work.exists():
        raise FileExistsError(
            f"refusing to overwrite an incomplete diagnostic work directory: {audit_work}"
        )
    cfg.NAME = str(audit_work)
    cfg.FORCE = True
    # Keep validation disabled.  Re-enabling the upstream reporting loop adds
    # 1,280 Lightning validation loops and confounds the RNG test with a large,
    # software-version-dependent control-flow cost.  The equivalent RNG draw
    # is installed below at the same boundary between training epochs.
    cfg.EVAL.FREQ = None
    cfg.PICARD.N = 80
    cfg.freeze()

    capture = DIAGNOSTIC_DIR / "unused_pinn_capture.pt"
    runner = _fresh_runner(cfg, capture)

    # The historical upstream EvalCallback ran once after each group of eight
    # training batches.  There are eight batches per epoch, so this is one draw
    # between consecutive epochs.  Replaying just equation.sample_x(t) in an
    # epoch-end hook preserves the global Torch draw tape without opening
    # nonterminal truth or constructing Lightning validation loops.
    import lightning.pytorch as pl

    callback_counter = {"calls": 0}

    class RngReplayAtEpochEnd(pl.Callback):
        def __init__(self, equation):
            super().__init__()
            self.equation = equation

        def on_train_epoch_end(self, trainer, pl_module):
            del trainer
            self.equation = self.equation.to(device=pl_module.device)
            t = torch.linspace(
                0.0, self.equation.T, 10_000, device=pl_module.device
            ).reshape(10_000, 1)
            with torch.no_grad():
                self.equation.sample_x(t)
            callback_counter["calls"] += 1

    base_callbacks = runner.get_reporting_callbacks

    def callbacks_with_rng_replay(self):
        del self
        callbacks, kwargs = base_callbacks()
        callbacks.append(RngReplayAtEpochEnd(runner.equation))
        return callbacks, kwargs

    runner.get_reporting_callbacks = MethodType(callbacks_with_rng_replay, runner)

    synchronize(device)
    if torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()
    started = time.perf_counter()
    old_cwd = Path.cwd()
    os.chdir(WORK_DIR)
    try:
        runner.run()
    finally:
        os.chdir(old_cwd)
    synchronize(device)
    replay_training_seconds = time.perf_counter() - started

    replay_checkpoint = audit_work / "model_80.ckpt"
    if not replay_checkpoint.is_file():
        raise FileNotFoundError(replay_checkpoint)
    model, _equation = _model_for_checkpoint("dpi", cfg, replay_checkpoint, device)
    t, x, truth, bank_path = _bank("dense_fourier")
    replay_prediction = _predict(model, t, x, device)
    with np.load(historical_prediction_path, allow_pickle=False) as payload:
        historical_prediction = np.asarray(payload["prediction"], dtype=np.float64)
    if replay_prediction.shape != historical_prediction.shape:
        raise RuntimeError(
            f"prediction shape mismatch: {replay_prediction.shape} vs {historical_prediction.shape}"
        )

    metrics = value_metrics(replay_prediction, truth)
    prediction_delta = replay_prediction - historical_prediction
    save_npz_atomic(
        PREDICTION_PATH,
        t=t,
        x=x,
        truth=truth,
        replay_prediction=replay_prediction,
        historical_prediction=historical_prediction,
    )
    historical_rmse = float(historical["value_rmse"])
    replay_rmse = float(metrics["rmse"])
    result = {
        "revision": REVISION,
        "status": "COMPLETE",
        "classification": "POST_FREEZE_DIAGNOSTIC_EXCLUDED_FROM_FRONTIER",
        "pde": "dense_fourier",
        "dimension": 100,
        "method": "dpi",
        "seed": TRAINING_SEED,
        "picard_iterations": 80,
        "epochs_per_picard": 16,
        "historical_eval_frequency": 8,
        "frontier_eval_frequency": None,
        "diagnostic_eval_frequency": None,
        "exact_nonterminal_truth_callback_enabled": False,
        "rng_only_callback_replay_enabled": True,
        "rng_replay_schedule": "one 10000-state sample_x draw after each training epoch",
        "rng_only_callback_calls": callback_counter["calls"],
        "expected_rng_only_callback_calls": 80 * 16,
        "ground_truth_used_for_model_selection_or_stopping": False,
        "nonterminal_ground_truth_used_during_training": False,
        "historical_ground_truth_callback_advances_global_torch_rng": True,
        "historical_rmse": historical_rmse,
        "replay_rmse": replay_rmse,
        "relative_rmse_difference": abs(replay_rmse - historical_rmse) / abs(historical_rmse),
        "historical_training_seconds": float(historical["training_seconds"]),
        "replay_training_seconds": replay_training_seconds,
        "prediction_rmse_vs_historical": float(np.sqrt(np.mean(prediction_delta**2))),
        "prediction_max_abs_vs_historical": float(np.max(np.abs(prediction_delta))),
        "prediction_file": str(PREDICTION_PATH),
        "prediction_sha256": file_sha256(PREDICTION_PATH),
        "historical_prediction_file": str(historical_prediction_path),
        "historical_prediction_sha256": file_sha256(historical_prediction_path),
        "test_bank": str(bank_path),
        "replay_checkpoint": str(replay_checkpoint),
        "replay_checkpoint_sha256": file_sha256(replay_checkpoint),
        "historical_checkpoint": str(historical_checkpoint_path),
        "historical_checkpoint_sha256": file_sha256(historical_checkpoint_path),
        "historical_config": str(historical_config_path),
        "historical_config_sha256": file_sha256(historical_config_path),
        "replay_config": str(audit_work / "config.yaml"),
        "replay_config_sha256": file_sha256(audit_work / "config.yaml"),
        "semantic_config_difference_from_frontier_C0": {
            "EVAL.FREQ": {"frontier_and_replay": None, "historical": 8},
            "NAME": "output-directory-only",
            "diagnostic_instrumentation": "epoch-end RNG-only draw replay",
        },
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated()),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "hostname": platform.node(),
        "torch_version": torch.__version__,
        "cuda_device": torch.cuda.get_device_name(0),
        "main_frontier_artifacts_modified": False,
        **{f"replay_{key}": value for key, value in metrics.items()},
    }
    write_json_atomic(RESULT_PATH, result)
    return result


def main() -> None:
    print(json.dumps(run(), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
