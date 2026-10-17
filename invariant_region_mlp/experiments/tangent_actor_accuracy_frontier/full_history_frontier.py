"""Unmodified raw Full-History recursion on the frozen comparison bank."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import time

import numpy as np
import torch

from invariant_region_mlp.experiments.nonconvex_tangent_actor.full_history_baseline import RawDriver
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    GenericFullHistoryMLP,
    MLPStats,
    recommended_chunk,
)

from .common import (
    OLD_PREDICTION_DIR,
    OLD_REFERENCE_DIR,
    PREDICTION_DIR,
    RAW_DIR,
    REFERENCE_DIR,
    array_sha256,
    atomic_json,
    ensure_dirs,
    environment_payload,
    metrics,
    stable_seed,
    synchronize,
)
from .configuration import ExperimentSpec
from .reference2d import load_reference


class MemoryTiledFullHistoryMLP(GenericFullHistoryMLP):
    """The same estimator with terminal Monte Carlo accumulated in blocks."""

    max_terminal_coordinates = 64_000_000

    def _terminal_estimate(self, n, t, x):
        b, d = len(t), self.driver.d
        samples = self.M**n
        tau = torch.clamp(self.driver.T - t, min=0.0)
        parent_terminal = self.driver.terminal(x)
        value_sum = torch.zeros(b, dtype=self.dtype, device=self.device)
        gradient_sum = torch.zeros(
            (b, d), dtype=self.dtype, device=self.device
        ) if self.driver.has_gradient else None
        sample_chunk = max(
            1, min(samples, self.max_terminal_coordinates // max(b * d, 1))
        )
        for begin in range(0, samples, sample_chunk):
            count = min(sample_chunk, samples - begin)
            normal = self._normal((b, count, d))
            terminal_x = x[:, None, :] + self.driver.sigma * torch.sqrt(tau)[:, None, None] * normal
            terminal = self.driver.terminal(terminal_x.reshape(-1, d)).reshape(b, count)
            value_sum += (terminal - parent_terminal[:, None]).sum(dim=1)
            if gradient_sum is not None:
                terminal_z = self.driver.sigma * self.driver.terminal_gradient_torch(
                    terminal_x.reshape(-1, d)
                ).reshape(b, count, d)
                gradient_sum += terminal_z.sum(dim=1)
            del normal, terminal_x, terminal
            if gradient_sum is not None:
                del terminal_z
        output = torch.zeros((b, self.state_dim), dtype=self.dtype, device=self.device)
        output[:, 0] = parent_terminal + value_sum / samples
        if gradient_sum is not None:
            output[:, 1:] = gradient_sum / samples
        self.stats.terminal_samples += b * samples
        return output


@torch.no_grad()
def _predict_tiled(raw_driver, t, x, *, n, M, seed, device):
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
        solver = MemoryTiledFullHistoryMLP(
            raw_driver, M, device,
            seed=stable_seed("nonconvex_full_history", seed, n, M, chunk_index),
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
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "finite": bool(total.nonfinite_count == 0 and np.all(np.isfinite(prediction))),
        "terminal_sample_memory_tiling": True,
        "terminal_tiling_max_coordinates": MemoryTiledFullHistoryMLP.max_terminal_coordinates,
        "estimator_formula_changed": False,
    }


def _bank(spec: ExperimentSpec) -> tuple[np.ndarray, np.ndarray]:
    if spec.family == "primary":
        path = OLD_PREDICTION_DIR / (
            "post_freeze_primary_d100_k2_lam0p2_eta0_main_mcv_s0.npz"
        )
    else:
        path = PREDICTION_DIR / f"post_freeze_{spec.run_id}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    with np.load(path) as archive:
        return archive["test_t"].astype(np.float64), archive["test_x"].astype(np.float64)


def _truth(spec: ExperimentSpec, t: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, str]:
    driver = spec.driver()
    refined = REFERENCE_DIR / f"{driver.name}_refined.npz"
    reference_path = refined if refined.exists() else OLD_REFERENCE_DIR / f"{driver.name}.npz"
    reference = load_reference(reference_path)
    return reference.value(t, x @ driver.B), str(reference_path)


def run(
    spec: ExperimentSpec,
    *,
    n: int,
    M: int,
    baseline_seed: int,
    role: str,
    device: torch.device,
) -> dict:
    if n not in (3, 4) or M < 2 or M > 512:
        raise ValueError((n, M))
    ensure_dirs()
    t, x = _bank(spec)
    truth, reference_path = _truth(spec, t, x)
    prediction, runtime = _predict_tiled(
        RawDriver(spec.driver()),
        t,
        x,
        n=n,
        M=M,
        seed=baseline_seed,
        device=device,
    )
    accuracy = metrics(prediction, truth) if runtime["finite"] else {
        "rmse": None, "relative_l2": None, "bias": None, "mae": None, "max_abs": None
    }
    prediction_path = PREDICTION_DIR / (
        f"full_history_{spec.family}_d{spec.d}_n{n}_M{M}_s{baseline_seed}.npz"
    )
    np.savez_compressed(prediction_path, prediction=prediction)
    result = {
        "method": "full_history_mlp_raw_adaptation",
        "spec": asdict(spec),
        "baseline_seed": baseline_seed,
        "n": n,
        "M": M,
        "configuration_role": role,
        "states": len(t),
        "test_bank_sha256": array_sha256(t, x),
        "reference_path": reference_path,
        "metrics": accuracy,
        "runtime": runtime,
        "training_seconds": 0.0,
        "total_seconds": runtime["wall_seconds"],
        "value_projection": False,
        "terminal_aware_feature": False,
        "mcv": False,
        "oracle_projection": False,
        "ground_truth_used_for_algorithm_or_configuration_selection": False,
        "prediction_path": str(prediction_path),
    }
    output = RAW_DIR / (
        f"full_history_{spec.family}_d{spec.d}_n{n}_M{M}_s{baseline_seed}.json"
    )
    atomic_json(output, result)
    print(json.dumps({
        "family": spec.family, "n": n, "M": M, "seed": baseline_seed,
        "rmse": accuracy["rmse"], "seconds": runtime["wall_seconds"],
        "finite": runtime["finite"], "role": role,
    }, sort_keys=True), flush=True)
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
    parser.add_argument("--role", choices=("profile", "frontier", "endpoint"), required=True)
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
        RAW_DIR / f"environment_fh_{args.family}_n{args.n}_M{args.M}_s{args.baseline_seed}.json",
        environment_payload("full_history_frontier"),
    )
    run(
        spec, n=args.n, M=args.M, baseline_seed=args.baseline_seed,
        role=args.role, device=torch.device(args.device),
    )


if __name__ == "__main__":
    main()
