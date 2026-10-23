"""Standard full-history MLP for value-only, gradient, and joint drivers."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.original_pde_benchmark_comparison.common import (
    synchronize,
)
from invariant_region_mlp.experiments.tangent_actor_strong_baselines.full_history_mlp import (
    GenericFullHistoryMLP,
    MLPStats,
)

from .benchmarks import make_benchmark
from .protocol import (
    FH_CONFIGURATIONS,
    FH_REPLICATES,
    PREDICTION_DIR,
    QUERY_COUNT,
    RAW_DIR,
    REVISION,
    array_sha256,
    ensure_directories,
    query_path,
    require_stage_a,
    require_stage_b,
    stable_seed,
    write_json_atomic,
    write_npz_atomic,
)


class _UnprojectedDriver:
    """Delegate a PDE while disabling metadata-only value projection guards."""

    r_low = -math.inf
    r_high = math.inf

    def __init__(self, implementation) -> None:
        self.impl = implementation

    def __getattr__(self, name: str):
        return getattr(self.impl, name)


class AuditedFullHistoryMLP(GenericFullHistoryMLP):
    """Generic recursion with the exact Bismut denominator in float64.

    The inherited implementation clamps ``delta`` at 1e-7 in the gradient
    weight.  That is useful as a robustness heuristic but is not the standard
    estimator.  This formal runner uses only the dtype's positive floor, so the
    preregistered standard estimator is not silently truncated.
    """

    def solve(self, n: int, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        t = t.to(dtype=self.dtype, device=self.device).reshape(-1)
        x = x.to(dtype=self.dtype, device=self.device).reshape(-1, self.driver.d)
        batch, d = len(t), self.driver.d
        self.stats.recursive_calls += 1
        self.stats.recursive_states += batch
        if n <= 0:
            return torch.zeros(
                (batch, self.state_dim), dtype=self.dtype, device=self.device
            )
        output = self._terminal_estimate(n, t, x)
        remaining = torch.clamp(self.driver.T - t, min=0.0)
        alpha = self.time_beta_alpha
        tiny = torch.finfo(self.dtype).tiny
        for level in range(1, n):
            siblings = self.M ** (n - level)
            uniform = torch.clamp(
                torch.rand(
                    (batch, siblings), dtype=self.dtype, device=self.device,
                    generator=self.generator,
                ),
                min=tiny,
            )
            fraction = uniform ** (1.0 / alpha)
            delta = remaining[:, None] * fraction
            normal = self._normal((batch, siblings, d))
            child_x = x[:, None, :] + self.driver.sigma * torch.sqrt(delta)[:, :, None] * normal
            child_t = t[:, None] + delta
            flat_t = child_t.reshape(-1)
            flat_x = child_x.reshape(-1, d)
            high = self.solve(level, flat_t, flat_x).reshape(
                batch, siblings, self.state_dim
            )
            difference = self._driver(high)
            self.stats.driver_evaluations += batch * siblings
            if level > 1:
                low = self.solve(level - 1, flat_t, flat_x).reshape(
                    batch, siblings, self.state_dim
                )
                difference = difference - self._driver(low)
                self.stats.driver_evaluations += batch * siblings
            importance = remaining[:, None] * fraction ** (1.0 - alpha) / alpha
            output[:, 0] += torch.mean(importance * difference, dim=1)
            if self.driver.has_gradient:
                bismut = normal / torch.sqrt(torch.clamp(delta, min=tiny))[:, :, None]
                output[:, 1:] += torch.mean(
                    importance[:, :, None] * difference[:, :, None] * bismut,
                    dim=1,
                )
            self.stats.transition_samples += batch * siblings
        return output


def _merge_stats(total: MLPStats, local: MLPStats) -> None:
    total.merge(local)


def execute(
    pde: str,
    d: int,
    n: int,
    M: int,
    replicate: int,
    device: torch.device,
    *,
    pilot: bool = False,
) -> dict[str, Any]:
    if (n, M) not in FH_CONFIGURATIONS and not pilot:
        raise ValueError((n, M))
    if not 0 <= replicate < FH_REPLICATES and not pilot:
        raise ValueError(replicate)
    require_stage_a()
    if not pilot:
        require_stage_b()
    ensure_directories()
    tag = "pilot" if pilot else "formal"
    run_id = f"fh_mlp_standard__{pde}__d{d}__n{n}_M{M}__r{replicate:02d}__{tag}"
    output = RAW_DIR / tag / "fh_mlp" / f"{run_id}.json"
    if output.is_file():
        existing = json.loads(output.read_text(encoding="utf-8"))
        if existing.get("status") in {"complete", "pilot_complete", "numerical_failure"}:
            return existing
    x = np.load(query_path(pde, d) / "x.npy", allow_pickle=False)
    t = np.load(query_path(pde, d) / "t.npy", allow_pickle=False)
    if pilot:
        x, t = x[:4], t[:4]
    benchmark = _UnprojectedDriver(make_benchmark(pde, d))
    prediction = np.full(len(x), np.nan, dtype=np.float64)
    states = np.full((len(x), 1 + (d if benchmark.has_gradient else 0)), np.nan, dtype=np.float64)
    point_seconds = np.zeros(len(x), dtype=np.float64)
    total = MLPStats()
    synchronize(device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    with torch.no_grad():
        for point in range(len(x)):
            point_started = time.perf_counter()
            solver = AuditedFullHistoryMLP(
                benchmark,
                M,
                device,
                seed=stable_seed("fh_mlp", pde, d, n, M, replicate, point, tag),
                dtype=torch.float64,
                time_beta_alpha=0.5,
            )
            result = solver.solve(
                n,
                torch.as_tensor(t[point:point + 1], dtype=torch.float64, device=device),
                torch.as_tensor(x[point:point + 1], dtype=torch.float64, device=device),
            )
            states[point] = result[0].detach().cpu().numpy()
            prediction[point] = states[point, 0]
            point_seconds[point] = time.perf_counter() - point_started
            _merge_stats(total, solver.stats)
    synchronize(device)
    elapsed = time.perf_counter() - started
    finite = np.isfinite(prediction)
    reference = np.load(query_path(pde, d) / "reference.npy", allow_pickle=False)[:len(x)]
    prediction_path = PREDICTION_DIR / tag / "fh_mlp" / f"{run_id}.npz"
    write_npz_atomic(
        prediction_path,
        t=t,
        x=x,
        prediction=prediction,
        state=states,
        reference=reference,
        point_seconds=point_seconds,
    )
    status = "pilot_complete" if pilot else "complete"
    if not np.all(finite) or total.nonfinite_count:
        status = "numerical_failure"
    record = {
        "revision": REVISION,
        "run_id": run_id,
        "method": "fh_mlp_standard",
        "status": status,
        "pde": pde,
        "dimension": d,
        "n": n,
        "M": M,
        "replicate": replicate,
        "points": len(x),
        "prediction_file": str(prediction_path),
        "query_input_sha256": array_sha256(t, x),
        "prediction_sha256": array_sha256(prediction),
        "training_seconds": 0.0,
        "inference_seconds": elapsed,
        "point_seconds": point_seconds.tolist(),
        "stats": asdict(total),
        "finite_count": int(finite.sum()),
        "nonfinite_count": int((~finite).sum()),
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else None,
        "dtype": "float64",
        "time_beta_alpha": 0.5,
        "terminal_gradient_estimator": "pathwise sigma*grad(g)",
        "driver_gradient_estimator": "untruncated Bismut weight",
        "multilevel_coupling": "common evaluation time/state; independent recursive high/low subtrees",
        "ground_truth_used_for_configuration_selection": False,
        "main_configuration": bool((n, M) == (3, 16)),
        "nonlinear_plugin_unbiased_claimed": False,
    }
    write_json_atomic(output, record)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pde", required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--M", type=int, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA device requested but unavailable")
    result = execute(
        args.pde, args.dimension, args.n, args.M, args.replicate,
        device, pilot=args.pilot,
    )
    print(json.dumps({key: result.get(key) for key in (
        "run_id", "status", "finite_count", "nonfinite_count", "inference_seconds"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()

