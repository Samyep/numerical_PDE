"""Common-bank evaluation of the frozen official SDGD-PINN network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import torch
from torch import nn

from .common import (
    OFFICIAL_WORK_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
)


class OfficialSDGDMLP(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        layers: list[nn.Module] = []
        dimensions = [101, 1024, 1024, 1024, 1]
        for index in range(len(dimensions) - 1):
            layers.append(nn.Linear(dimensions[index], dimensions[index + 1]))
            if index != len(dimensions) - 2:
                layers.append(nn.Tanh())
        self.nn = nn.Sequential(*layers)

    def forward(self, tx: torch.Tensor) -> torch.Tensor:
        x, t = tx[:, :-1], tx[:, -1:]
        terminal = torch.log1p(x.square().sum(dim=1, keepdim=True)) - np.log(2.0)
        return self.nn(tx) * (1.0 - t) + terminal


def evaluate(seed: int, device: torch.device) -> dict[str, object]:
    ensure_dirs()
    record = json.loads((RAW_DIR / f"official_sdgd_s{seed}_official.json").read_text())
    model = OfficialSDGDMLP().to(device)
    model.load_state_dict(
        torch.load(record["checkpoint"], map_location=device, weights_only=False)
    )
    model.eval()
    bank = torch.load(
        OFFICIAL_WORK_DIR / "sdgd/HJB_Log_100_grad.pkl",
        map_location="cpu",
        weights_only=False,
    )
    tx = bank["X"].numpy()
    reference = bank["Y"].numpy().reshape(-1)
    predictions: list[np.ndarray] = []
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
        torch.cuda.synchronize(device)
    started = time.perf_counter()
    with torch.no_grad():
        for begin in range(0, len(tx), 1024):
            local = torch.as_tensor(tx[begin : begin + 1024], device=device)
            predictions.append(model(local)[:, 0].cpu().numpy())
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    inference_seconds = time.perf_counter() - started
    prediction = np.concatenate(predictions).astype(np.float64)
    error = prediction - reference.astype(np.float64)
    with np.load(RESULTS_DIR / "evaluation_bank_hjb_log_d100.npz") as common:
        common_indices = common["test_indices"]
        common_reference = common["value"][common_indices]
    common_error = prediction[common_indices] - common_reference
    output = PREDICTION_DIR / f"official_sdgd_s{seed}.npz"
    np.savez_compressed(output, prediction=prediction)
    result: dict[str, object] = {
        "method": "sdgd_pinn_official",
        "benchmark": "hjb_log_d100",
        "seed": seed,
        "states": len(reference),
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(
            np.linalg.norm(error) / np.linalg.norm(reference)
        ),
        "value_bias": float(np.mean(error)),
        "common_960_value_rmse": float(
            np.sqrt(np.mean(np.square(common_error)))
        ),
        "common_960_value_relative_l2": float(
            np.linalg.norm(common_error) / np.linalg.norm(common_reference)
        ),
        "common_960_value_bias": float(np.mean(common_error)),
        "training_seconds": record[
            "training_seconds_including_official_per_epoch_test_reporting"
        ],
        "inference_seconds": inference_seconds,
        "total_seconds": float(
            record["training_seconds_including_official_per_epoch_test_reporting"]
        )
        + inference_seconds,
        "peak_gpu_bytes_inference": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "prediction_file": str(output),
        "identical_bank_for_ours_required": True,
        "test_labels_used_by_official_code_for_reporting_each_epoch": True,
        "test_labels_used_in_training_loss": False,
    }
    atomic_json(RAW_DIR / f"official_sdgd_s{seed}_evaluation.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = evaluate(args.seed, torch.device(args.device))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

