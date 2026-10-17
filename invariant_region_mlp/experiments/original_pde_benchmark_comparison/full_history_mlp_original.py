"""Fixed-budget Full-History MLP control on supported original PDEs."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    GenericFullHistoryMLP,
    MLPStats,
    recommended_chunk,
)

from .benchmarks import AllenCahnBenchmark, HJBLogBenchmark, OriginalBenchmark
from .common import PREDICTION_DIR, RAW_DIR, RESULTS_DIR, atomic_json, stable_seed, synchronize
from .references import hjb_log_reference


FIXED_CONFIGURATIONS = ((2, 16), (2, 32), (3, 8), (3, 12), (3, 16), (3, 24))


class OriginalPDEMLPDriver:
    """Minimal interface for the existing generic recursive estimator."""

    def __init__(self, benchmark: OriginalBenchmark) -> None:
        self.benchmark = benchmark
        self.d = benchmark.d
        self.T = benchmark.T
        self.sigma = benchmark.sigma
        self.has_gradient = benchmark.actor_kind == "gradient"
        # The original-PDE control is raw: no invariant-region correction is
        # inserted into a competitor.  Infinite bounds make the inherited
        # clamp exactly the identity on every finite recursive state.
        self.r_low = -float("inf")
        self.r_high = float("inf")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.benchmark.terminal(x)

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        return self.benchmark.terminal_gradient(x)

    def driver(
        self, value: torch.Tensor, z: torch.Tensor | None
    ) -> torch.Tensor:
        if isinstance(self.benchmark, HJBLogBenchmark):
            if z is None:
                raise RuntimeError("HJB-Log Full-History MLP requires z")
            return -0.5 * z.square().sum(dim=-1)
        if isinstance(self.benchmark, AllenCahnBenchmark):
            return value - value.pow(3)
        raise TypeError("the zero-drift Full-History implementation does not support this PDE")

    @staticmethod
    def working_to_original(value: np.ndarray) -> np.ndarray:
        return np.asarray(value)


def _metrics(prediction: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    error = np.asarray(prediction) - np.asarray(reference)
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "relative_l2": float(np.linalg.norm(error) / np.linalg.norm(reference)),
        "bias": float(np.mean(error)),
        "mae": float(np.mean(np.abs(error))),
    }


@torch.no_grad()
def _predict(
    benchmark: OriginalBenchmark,
    t: np.ndarray,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, dict[str, Any]]:
    driver = OriginalPDEMLPDriver(benchmark)
    times = np.asarray(t, dtype=np.float32)
    states = np.asarray(x, dtype=np.float32)
    prediction = np.empty(len(states), dtype=np.float64)
    chunk = recommended_chunk(benchmark.d, n, M)
    total = MLPStats()
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(states), chunk)):
        end = min(begin + chunk, len(states))
        solver = GenericFullHistoryMLP(
            driver,
            M,
            device,
            seed=stable_seed(
                "original_pde_full_history_mlp", benchmark.name, seed, n, M, chunk_index
            ),
        )
        state = solver.solve(
            n,
            torch.as_tensor(times[begin:end], dtype=torch.float32, device=device),
            torch.as_tensor(states[begin:end], dtype=torch.float32, device=device),
        )
        prediction[begin:end] = state[:, 0].double().cpu().numpy()
        total.merge(solver.stats)
    synchronize(device)
    elapsed = time.perf_counter() - started
    if total.nonfinite_count or not np.all(np.isfinite(prediction)):
        raise FloatingPointError(
            f"Full-History MLP produced {total.nonfinite_count} nonfinite generator values"
        )
    return prediction, {
        "wall_seconds": elapsed,
        "chunk_size": chunk,
        "stats": asdict(total),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
    }


def run(
    benchmark_name: str,
    n: int,
    M: int,
    seed: int,
    device: torch.device,
) -> dict[str, Any]:
    if (n, M) not in FIXED_CONFIGURATIONS:
        raise ValueError("configuration was not in the frozen Full-History MLP grid")
    benchmark: OriginalBenchmark = (
        HJBLogBenchmark()
        if benchmark_name == "hjb_log_d100"
        else AllenCahnBenchmark()
    )
    with np.load(RESULTS_DIR / f"evaluation_bank_{benchmark_name}.npz") as bank:
        indices = bank["test_indices"]
        t = bank["t"][indices]
        x = bank["x"][indices]
        reference = bank["value"][indices]
    function_prediction, function_meta = _predict(
        benchmark, t, x, n=n, M=M, seed=seed, device=device
    )
    native_t = np.asarray([0.0])
    native_x = np.zeros((1, benchmark.d), dtype=np.float64)
    native_prediction, native_meta = _predict(
        benchmark,
        native_t,
        native_x,
        n=n,
        M=M,
        seed=stable_seed("native", seed),
        device=device,
    )
    if benchmark_name == "hjb_log_d100":
        native_reference = float(
            hjb_log_reference(np.asarray([0.0]), native_x)[0]
        )
    else:
        audit = json.loads((RESULTS_DIR / "reference_accuracy_audit.json").read_text())
        native_reference = float(audit["allen_cahn"]["native_u_0_0_extrapolated"])
    prediction_path = PREDICTION_DIR / (
        f"full_history_mlp_{benchmark_name}_n{n}_M{M}_s{seed}.npz"
    )
    np.savez_compressed(
        prediction_path,
        prediction=function_prediction,
        test_indices=indices,
        native=native_prediction,
    )
    total_seconds = float(function_meta["wall_seconds"] + native_meta["wall_seconds"])
    result = {
        "method": "full_history_mlp_raw",
        "benchmark": benchmark_name,
        "seed": seed,
        "n": n,
        "M": M,
        "configuration_grid_frozen_before_results": True,
        "function": {
            "states": len(reference),
            "metrics": _metrics(function_prediction, reference),
            **function_meta,
        },
        "native": {
            "point": "t=0,x=0",
            "prediction": float(native_prediction[0]),
            "reference": native_reference,
            "absolute_error": abs(float(native_prediction[0]) - native_reference),
            **native_meta,
        },
        "training_seconds": 0.0,
        "parameter_count": 0,
        "inference_seconds": total_seconds,
        "total_seconds": total_seconds,
        "pathwise_terminal_gradient": benchmark.actor_kind == "gradient",
        "recursive_value_projection": False,
        "terminal_aware_neural_feature": False,
        "mcv": False,
        "dtype": "float32",
        "prediction_file": str(prediction_path),
        "ground_truth_used_for_configuration_selection": False,
    }
    atomic_json(
        RAW_DIR / f"full_history_mlp_{benchmark_name}_n{n}_M{M}_s{seed}.json",
        result,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--benchmark", choices=("hjb_log_d100", "allen_cahn_d100"), required=True
    )
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--M", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = run(args.benchmark, args.n, args.M, args.seed, torch.device(args.device))
    print(
        json.dumps(
            {
                "status": "complete",
                "rmse": result["function"]["metrics"]["rmse"],
                "seconds": result["total_seconds"],
            }
        )
    )


if __name__ == "__main__":
    main()
