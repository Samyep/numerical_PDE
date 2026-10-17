"""Evaluate frozen official DPI/PINN-HTE/D-DBSDE models on one fixed bank."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
import torch

from .common import PREDICTION_DIR, RAW_DIR, RESULTS_DIR, atomic_json, ensure_dirs
from .official_dpi_runner import CONFIGS


def evaluate(
    method: str, seed: int, device: torch.device, *, calibration: bool = False
) -> dict[str, object]:
    from picard import load_cfg, PicardRunner
    from picard.solution_jac import PicardSolutionGradientWrapper

    ensure_dirs()
    suffix = "calibration" if calibration else "official"
    record_path = RAW_DIR / f"official_dpi_{method}_s{seed}_{suffix}.json"
    record = json.loads(record_path.read_text())
    cfg = load_cfg(record["expanded_config"])
    equation = PicardRunner.get_equation(cfg).to(device=device)

    class ReconstructionRunner(SimpleNamespace):
        def get_solution_plain(self):
            solution_cls = PicardRunner.get_solution_cls(self.cfg)
            return solution_cls.construct_solution(self)

    dummy = ReconstructionRunner(equation=equation, cfg=cfg)
    # All three official methods first call PicardRunner.get_solution().  The
    # selected equation contains a gradient term, so the trained network is
    # the authors' gradient-supervised wrapper around the PISGradNet.
    model = PicardSolutionGradientWrapper.construct_solution(dummy).to(device=device)
    payload = torch.load(record["final_model"], map_location=device, weights_only=False)
    state_dict = payload["state_dict"] if "state_dict" in payload else payload
    model.load_state_dict(state_dict)
    model.eval()
    with np.load(RESULTS_DIR / "evaluation_bank_dpi_ou_gmm_d100.npz") as bank:
        indices = bank["test_indices"]
        t = bank["t"][indices]
        x = bank["x"][indices]
        reference = bank["value"][indices]
        reference_z = bank["z"][indices]
    predictions: list[np.ndarray] = []
    gradients: list[np.ndarray] = []
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
        torch.cuda.synchronize(device)
    started = time.perf_counter()
    for begin in range(0, len(t), 128):
        end = min(begin + 128, len(t))
        local_t = torch.as_tensor(t[begin:end], dtype=torch.float64, device=device)
        local_x = torch.as_tensor(
            x[begin:end], dtype=torch.float64, device=device
        ).requires_grad_(True)
        value = model(torch.cat((local_t[:, None], local_x), dim=1))[:, 0]
        gradient = torch.autograd.grad(value.sum(), local_x)[0]
        predictions.append(value.detach().cpu().numpy())
        gradients.append(gradient.detach().cpu().numpy())
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    inference_seconds = time.perf_counter() - started
    prediction = np.concatenate(predictions)
    gradient = np.concatenate(gradients)
    error = prediction - reference
    gradient_error = gradient - reference_z  # alpha=1, so z=grad(u)
    output = PREDICTION_DIR / f"official_dpi_{method}_s{seed}_{suffix}.npz"
    np.savez_compressed(output, prediction=prediction, gradient=gradient)
    result: dict[str, object] = {
        "method": method,
        "benchmark": "dpi_ou_gmm_d100",
        "seed": seed,
        "calibration": calibration,
        "states": len(reference),
        "value_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "value_relative_l2": float(
            np.linalg.norm(error) / np.linalg.norm(reference)
        ),
        "value_bias": float(np.mean(error)),
        "gradient_component_rmse": float(
            np.sqrt(np.mean(np.square(gradient_error)))
        ),
        "gradient_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient_error).sum(axis=1)))
        ),
        "inference_seconds": inference_seconds,
        "training_seconds": record["training_seconds"],
        "total_seconds": float(record["training_seconds"]) + inference_seconds,
        "peak_gpu_bytes_inference": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "prediction_file": str(output),
        "test_bank_used_for_training_or_selection": False,
    }
    atomic_json(
        RAW_DIR / f"official_dpi_{method}_s{seed}_{suffix}_evaluation.json",
        result,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=tuple(CONFIGS), required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--calibration", action="store_true")
    args = parser.parse_args()
    result = evaluate(
        args.method,
        args.seed,
        torch.device(args.device),
        calibration=args.calibration,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

