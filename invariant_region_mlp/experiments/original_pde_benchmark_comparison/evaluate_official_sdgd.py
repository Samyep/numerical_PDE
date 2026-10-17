"""Common-bank evaluation of the frozen official SDGD-PINN network."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
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
from .references import hjb_log_reference


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


def evaluate(
    seed: int, device: torch.device, *, calibration: bool = False
) -> dict[str, object]:
    ensure_dirs()
    suffix = "calibration" if calibration else "official"
    record = json.loads((RAW_DIR / f"official_sdgd_s{seed}_{suffix}.json").read_text())
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
    with np.load(RESULTS_DIR / "evaluation_bank_hjb_log_d100.npz") as common:
        common_indices = common["test_indices"]
        common_reference = common["value"][common_indices]
        common_reference_z = common["z"][common_indices]

    # Measure the shared 960-state task directly.  The full official bank has
    # 20,000 states, so charging its forward time to the common task would
    # overstate amortized inference cost by roughly 20.8x.
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
        torch.cuda.synchronize(device)
    common_started = time.perf_counter()
    with torch.no_grad():
        common_prediction = (
            model(torch.as_tensor(tx[common_indices], device=device))[:, 0]
            .cpu()
            .numpy()
            .astype(np.float64)
        )
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    common_inference_seconds = time.perf_counter() - common_started
    common_error = common_prediction - common_reference

    predictions: list[np.ndarray] = []
    if device.type == "cuda":
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
    gradient_chunks: list[np.ndarray] = []
    gradient_started = time.perf_counter()
    for begin in range(0, len(common_indices), 256):
        local_indices = common_indices[begin : begin + 256]
        local = torch.as_tensor(tx[local_indices], device=device).requires_grad_(True)
        value = model(local)[:, 0]
        gradient = torch.autograd.grad(value.sum(), local)[0][:, :-1]
        gradient_chunks.append(
            (np.sqrt(2.0) * gradient.detach().cpu().numpy()).astype(np.float64)
        )
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    gradient_evaluation_seconds = time.perf_counter() - gradient_started
    common_z = np.concatenate(gradient_chunks)
    common_z_error = common_z - common_reference_z

    native_input = torch.zeros((1, 101), dtype=torch.float32, device=device)
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    native_started = time.perf_counter()
    with torch.no_grad():
        native_prediction = float(model(native_input)[0, 0].cpu().item())
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    native_inference_seconds = time.perf_counter() - native_started

    native_tx = torch.zeros((1, 101), dtype=torch.float32, device=device).requires_grad_(True)
    native_value_tensor = model(native_tx)[:, 0]
    native_gradient = torch.autograd.grad(native_value_tensor.sum(), native_tx)[0][0, :-1]
    native_reference = float(
        hjb_log_reference(np.asarray([0.0]), np.zeros((1, 100)))[0]
    )
    output = PREDICTION_DIR / f"official_sdgd_s{seed}_{suffix}.npz"
    np.savez_compressed(
        output,
        prediction=prediction,
        common_prediction=common_prediction,
        common_indices=common_indices,
    )
    training_seconds = float(
        record["training_seconds_including_official_per_epoch_test_reporting"]
    )
    nonfinite_count = int((~np.isfinite(prediction)).sum())
    nonfinite_count += int((~np.isfinite(common_prediction)).sum())
    nonfinite_count += int((~np.isfinite(common_z)).sum())
    nonfinite_count += int(not np.isfinite(native_prediction))
    nonfinite_count += int((~torch.isfinite(native_gradient)).sum().item())
    evaluation_repository_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(__file__).resolve().parents[3],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    result: dict[str, object] = {
        "method": "sdgd_pinn_official",
        "benchmark": "hjb_log_d100",
        "seed": seed,
        "calibration": calibration,
        "evaluation_repository_commit": evaluation_repository_commit,
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
        "common_960_gradient_component_rmse": float(
            np.sqrt(np.mean(np.square(common_z_error)))
        ),
        "common_960_gradient_vector_rmse": float(
            np.sqrt(np.mean(np.square(common_z_error).sum(axis=1)))
        ),
        "common_960_gradient_relative_l2": float(
            np.linalg.norm(common_z_error) / np.linalg.norm(common_reference_z)
        ),
        "gradient_evaluation_seconds_reported_separately": gradient_evaluation_seconds,
        "native": {
            "point": "t=0,x=0",
            "prediction": native_prediction,
            "reference": native_reference,
            "absolute_error": abs(native_prediction - native_reference),
            "gradient_vector_error": float(
                np.sqrt(2.0) * torch.linalg.vector_norm(native_gradient).item()
            ),
            "exact_gradient": "zero by radial symmetry",
        },
        "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "training_seconds": training_seconds,
        "common_960_inference_seconds": common_inference_seconds,
        "common_960_total_seconds": training_seconds + common_inference_seconds,
        "native_inference_seconds": native_inference_seconds,
        "native_total_seconds": training_seconds + native_inference_seconds,
        "inference_seconds": inference_seconds,
        "total_seconds": training_seconds + inference_seconds,
        "nonfinite_count": nonfinite_count,
        "peak_gpu_bytes_inference": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "prediction_file": str(output),
        "identical_bank_for_ours_required": True,
        "test_labels_used_by_official_code_for_reporting_each_epoch": True,
        "test_labels_used_in_training_loss": False,
    }
    atomic_json(RAW_DIR / f"official_sdgd_s{seed}_{suffix}_evaluation.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--calibration", action="store_true")
    args = parser.parse_args()
    result = evaluate(
        args.seed, torch.device(args.device), calibration=args.calibration
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

