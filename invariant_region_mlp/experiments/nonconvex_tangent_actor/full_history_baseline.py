"""Raw Full-History MLP adaptation for the new semilinear drivers."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    GenericFullHistoryMLP,
    MLPStats,
    recommended_chunk,
)

from .common import (
    PREDICTION_DIR,
    RAW_DIR,
    RESULTS_DIR,
    atomic_json,
    ensure_dirs,
    environment_payload,
    metrics,
    stable_seed,
    synchronize,
)
from .configuration import ExperimentSpec


FIXED_GRID = (
    (1, 64),
    (1, 256),
    (2, 8),
    (2, 16),
    (2, 32),
    (2, 64),
    (3, 4),
    (3, 8),
    (3, 12),
    (3, 16),
    (3, 24),
    (3, 32),
)


class RawDriver:
    """Disable invariant-region projection in the competitor."""

    def __init__(self, driver) -> None:
        self.original = driver
        self.d, self.T = driver.d, driver.T
        self.sigma, self.has_gradient = driver.sigma, True
        self.r_low, self.r_high = -float("inf"), float("inf")

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.original.terminal(x)

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        return self.original.terminal_gradient_torch(x)

    def driver(self, y: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        return self.original.driver(y, z)

    @staticmethod
    def working_to_original(value: np.ndarray) -> np.ndarray:
        return np.asarray(value)


@torch.no_grad()
def _predict(
    raw_driver: RawDriver,
    t: np.ndarray,
    x: np.ndarray,
    *,
    n: int,
    M: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, dict[str, Any]]:
    states = np.asarray(x, dtype=np.float32)
    times = np.asarray(t, dtype=np.float32)
    prediction = np.empty(len(states), dtype=np.float64)
    chunk = recommended_chunk(raw_driver.d, n, M)
    total = MLPStats()
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    for chunk_index, begin in enumerate(range(0, len(states), chunk)):
        end = min(begin + chunk, len(states))
        solver = GenericFullHistoryMLP(
            raw_driver,
            M,
            device,
            seed=stable_seed(
                "nonconvex_full_history", seed, n, M, chunk_index
            ),
        )
        output = solver.solve(
            n,
            torch.as_tensor(times[begin:end], dtype=torch.float32, device=device),
            torch.as_tensor(states[begin:end], dtype=torch.float32, device=device),
        )
        prediction[begin:end] = output[:, 0].double().cpu().numpy()
        total.merge(solver.stats)
    synchronize(device)
    elapsed = time.perf_counter() - started
    return prediction, {
        "wall_seconds": elapsed,
        "chunk_size": chunk,
        "stats": asdict(total),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "finite": bool(
            total.nonfinite_count == 0 and np.all(np.isfinite(prediction))
        ),
    }


def run(
    spec: ExperimentSpec,
    *,
    n: int,
    M: int,
    baseline_seed: int,
    actor_seed: int,
    device: torch.device,
) -> dict[str, Any]:
    if (n, M) not in FIXED_GRID:
        raise ValueError("configuration is outside the frozen Full-History grid")
    ensure_dirs()
    actor_spec = ExperimentSpec(
        family=spec.family,
        d=spec.d,
        lam=spec.lam,
        eta=spec.eta,
        seed=actor_seed,
        budget="main",
        estimator="mcv",
    )
    bank_path = PREDICTION_DIR / f"post_freeze_{actor_spec.run_id}.npz"
    if not bank_path.exists():
        raise FileNotFoundError(
            "the frozen common test bank must exist before baseline execution"
        )
    with np.load(bank_path) as bank:
        t, x, truth = bank["test_t"], bank["test_x"], bank["truth"]
    raw_driver = RawDriver(spec.driver())
    prediction, metadata = _predict(
        raw_driver,
        t,
        x,
        n=n,
        M=M,
        seed=baseline_seed,
        device=device,
    )
    accuracy = (
        metrics(prediction, truth)
        if metadata["finite"]
        else {
            "rmse": None,
            "relative_l2": None,
            "bias": None,
            "mae": None,
            "max_abs": None,
        }
    )
    prediction_path = PREDICTION_DIR / (
        f"full_history_{actor_spec.run_id}_n{n}_M{M}_s{baseline_seed}.npz"
    )
    np.savez_compressed(prediction_path, prediction=prediction)
    result = {
        "method": "full_history_mlp_raw_adaptation",
        "spec": asdict(spec),
        "actor_test_bank_run_id": actor_spec.run_id,
        "baseline_seed": baseline_seed,
        "n": n,
        "M": M,
        "states": len(t),
        "metrics": accuracy,
        "runtime": metadata,
        "training_seconds": 0.0,
        "total_seconds": metadata["wall_seconds"],
        "value_projection": False,
        "terminal_aware_feature": False,
        "mcv": False,
        "oracle_projection": False,
        "configuration_grid_frozen": True,
        "ground_truth_used_for_algorithm": False,
        "prediction_path": str(prediction_path),
    }
    output = RAW_DIR / (
        f"baseline_full_history_{actor_spec.run_id}_n{n}_M{M}_s{baseline_seed}.json"
    )
    atomic_json(output, result)
    print(
        json.dumps(
            {
                "method": result["method"],
                "run_id": actor_spec.run_id,
                "n": n,
                "M": M,
                "seed": baseline_seed,
                "rmse": accuracy["rmse"],
                "seconds": metadata["wall_seconds"],
                "finite": metadata["finite"],
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--family", choices=("primary", "coupled"), default="primary")
    parser.add_argument("--d", type=int, default=100)
    parser.add_argument("--lam", type=float, default=0.2)
    parser.add_argument("--eta", type=float, default=1.0)
    parser.add_argument("--actor-seed", type=int, default=0)
    parser.add_argument("--baseline-seed", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--M", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    spec = ExperimentSpec(
        family=args.family,
        d=args.d,
        lam=args.lam,
        eta=args.eta if args.family == "coupled" else 0.0,
        seed=args.actor_seed,
        budget="main",
        estimator="mcv",
    )
    ensure_dirs()
    atomic_json(
        RAW_DIR
        / f"environment_full_history_{spec.run_id}_n{args.n}_M{args.M}_s{args.baseline_seed}.json",
        environment_payload("full_history_baseline"),
    )
    run(
        spec,
        n=args.n,
        M=args.M,
        baseline_seed=args.baseline_seed,
        actor_seed=args.actor_seed,
        device=torch.device(args.device),
    )


if __name__ == "__main__":
    main()
