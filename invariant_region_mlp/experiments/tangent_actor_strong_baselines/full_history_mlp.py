"""GPU full-history recursive multilevel Picard estimator for generic f(u,z).

This is an inference-only Monte Carlo solver.  It therefore has zero training
cost and is compared on total compute, not on the training-only frontier.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import time
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    ConcaveSemilinearDriver,
)

from .common import RAW_DIR, atomic_json, driver_for, stable_seed, synchronize


@dataclass
class MLPStats:
    terminal_samples: int = 0
    transition_samples: int = 0
    driver_evaluations: int = 0
    recursive_states: int = 0
    recursive_calls: int = 0
    value_projection_count: int = 0
    nonfinite_count: int = 0

    def merge(self, other: "MLPStats") -> None:
        for name in self.__dataclass_fields__:
            setattr(self, name, getattr(self, name) + getattr(other, name))


class GenericFullHistoryMLP:
    def __init__(
        self,
        driver: ConcaveSemilinearDriver,
        M: int,
        device: torch.device,
        *,
        seed: int,
        dtype: torch.dtype = torch.float32,
        time_beta_alpha: float = 0.5,
    ) -> None:
        self.driver, self.M, self.device = driver, int(M), device
        self.dtype, self.time_beta_alpha = dtype, float(time_beta_alpha)
        self.generator = torch.Generator(device=device).manual_seed(seed)
        self.stats = MLPStats()

    def _normal(self, shape: tuple[int, ...]) -> torch.Tensor:
        return torch.randn(
            shape, dtype=self.dtype, device=self.device, generator=self.generator
        )

    def _terminal_estimate(
        self, n: int, t: torch.Tensor, x: torch.Tensor
    ) -> torch.Tensor:
        b, d = len(t), self.driver.d
        samples = self.M**n
        tau = torch.clamp(self.driver.T - t, min=0.0)
        normal = self._normal((b, samples, d))
        terminal_x = x[:, None, :] + self.driver.sigma * torch.sqrt(tau)[:, None, None] * normal
        terminal = self.driver.terminal(terminal_x.reshape(-1, d)).reshape(b, samples)
        parent_terminal = self.driver.terminal(x)
        output = torch.zeros((b, d + 1), dtype=self.dtype, device=self.device)
        # Centered terminal value estimator preserves its expectation.
        output[:, 0] = parent_terminal + torch.mean(
            terminal - parent_terminal[:, None], dim=1
        )
        terminal_z = self.driver.sigma * self.driver.terminal_gradient_torch(
            terminal_x.reshape(-1, d)
        ).reshape(b, samples, d)
        output[:, 1:] = terminal_z.mean(dim=1)
        self.stats.terminal_samples += b * samples
        return output

    def _driver(self, state: torch.Tensor) -> torch.Tensor:
        raw_value = state[..., 0]
        value = torch.clamp(raw_value, self.driver.r_low, self.driver.r_high)
        self.stats.value_projection_count += int((value != raw_value).sum().item())
        z = state[..., 1:] if self.driver.has_gradient else None
        result = self.driver.driver(value, z)
        self.stats.nonfinite_count += int((~torch.isfinite(result)).sum().item())
        return result

    def solve(self, n: int, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        t = t.to(dtype=self.dtype, device=self.device).reshape(-1)
        x = x.to(dtype=self.dtype, device=self.device).reshape(-1, self.driver.d)
        b, d = len(t), self.driver.d
        self.stats.recursive_calls += 1
        self.stats.recursive_states += b
        if n <= 0:
            return torch.zeros((b, d + 1), dtype=self.dtype, device=self.device)
        output = self._terminal_estimate(n, t, x)
        remaining = torch.clamp(self.driver.T - t, min=0.0)
        alpha = self.time_beta_alpha
        tiny = torch.finfo(self.dtype).tiny
        for level in range(1, n):
            siblings = self.M ** (n - level)
            uniform = torch.clamp(
                torch.rand(
                    (b, siblings),
                    dtype=self.dtype,
                    device=self.device,
                    generator=self.generator,
                ),
                min=tiny,
            )
            fraction = uniform ** (1.0 / alpha)
            delta = remaining[:, None] * fraction
            normal = self._normal((b, siblings, d))
            child_x = x[:, None, :] + self.driver.sigma * torch.sqrt(delta)[:, :, None] * normal
            child_t = t[:, None] + delta
            flat_t = child_t.reshape(-1)
            flat_x = child_x.reshape(-1, d)
            high = self.solve(level, flat_t, flat_x).reshape(b, siblings, d + 1)
            difference = self._driver(high)
            self.stats.driver_evaluations += b * siblings
            if level > 1:
                low = self.solve(level - 1, flat_t, flat_x).reshape(b, siblings, d + 1)
                difference = difference - self._driver(low)
                self.stats.driver_evaluations += b * siblings
            importance = remaining[:, None] * fraction ** (1.0 - alpha) / alpha
            output[:, 0] += torch.mean(importance * difference, dim=1)
            bismut = normal / torch.sqrt(torch.clamp(delta, min=1e-7))[:, :, None]
            output[:, 1:] += torch.mean(
                importance[:, :, None] * difference[:, :, None] * bismut, dim=1
            )
            self.stats.transition_samples += b * siblings
        return output


def recommended_chunk(d: int, n: int, M: int) -> int:
    root_coordinates = max(1, d * (M**n))
    return max(1, min(32, 8_000_000 // root_coordinates))


@torch.no_grad()
def predict_full_history_mlp(
    benchmark: str,
    x: np.ndarray,
    device: torch.device,
    *,
    n: int,
    M: int,
    seed: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    driver = driver_for(benchmark)
    states = np.asarray(x, dtype=np.float32)
    prediction = np.empty(len(states), dtype=np.float64)
    chunk = recommended_chunk(driver.d, n, M)
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
            seed=stable_seed("full_history_mlp", benchmark, seed, n, M, chunk_index),
        )
        local_x = torch.as_tensor(states[begin:end], dtype=torch.float32, device=device)
        local_t = torch.zeros(end - begin, dtype=torch.float32, device=device)
        state = solver.solve(n, local_t, local_x)
        prediction[begin:end] = (
            driver.working_to_original(state[:, 0].cpu().numpy()).astype(np.float64)
        )
        total.merge(solver.stats)
    synchronize(device)
    elapsed = time.perf_counter() - started
    if total.nonfinite_count or not np.all(np.isfinite(prediction)):
        raise FloatingPointError(
            f"full-history MLP produced {total.nonfinite_count} nonfinite driver values"
        )
    metadata = {
        "method": "full_history_mlp",
        "benchmark": benchmark,
        "seed": seed,
        "n": n,
        "M": M,
        "chunk_size": chunk,
        "training_seconds": 0.0,
        "inference_seconds": elapsed,
        "total_seconds": elapsed,
        "stats": asdict(total),
        "pathwise_terminal_gradient": True,
        "value_domain_truncation": [driver.r_low, driver.r_high],
        "dtype": "float32",
        "peak_gpu_bytes": int(torch.cuda.max_memory_allocated(device))
        if device.type == "cuda"
        else None,
        "ground_truth_used_for_configuration_selection": False,
    }
    return prediction, metadata


DEFAULT_CONFIGURATIONS = (
    (1, 16),
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
    (3, 40),
    (3, 48),
    (3, 64),
)


__all__ = [
    "DEFAULT_CONFIGURATIONS", "GenericFullHistoryMLP", "MLPStats",
    "predict_full_history_mlp", "recommended_chunk",
]
