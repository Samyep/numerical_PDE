"""Adapters for the four preregistered equations."""

from __future__ import annotations

from functools import lru_cache
import math
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.main_highdim_expansion.benchmarks import (
    exact_u as highdim_exact_u,
    make_b6,
    make_dense_fourier,
    metadata as highdim_metadata,
)
from invariant_region_mlp.experiments.original_pde_benchmark_comparison.references import (
    hjb_log_reference,
)
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.benchmarks import (
    AllenCahnBenchmark,
    HJBLogBenchmark,
)


class StandardAdapter:
    """Supply the common FH-MLP/Deep-BSDE interface without changing a PDE."""

    reconstruction_sign = 1.0
    r_low = -math.inf
    r_high = math.inf

    def __init__(self, implementation: Any, pde: str) -> None:
        self.impl = implementation
        self.pde = pde
        self.name = implementation.name
        self.d = implementation.d
        self.T = implementation.T
        self.sigma = implementation.sigma
        self.actor_kind = implementation.actor_kind
        self.has_gradient = implementation.actor_kind == "gradient"
        if pde == "allen_cahn":
            # Metadata only: the formal standard FH recursion is not projected.
            self.r_low, self.r_high = -math.inf, math.inf

    def terminal(self, x: torch.Tensor) -> torch.Tensor:
        return self.impl.terminal(x)

    def terminal_gradient_torch(self, x: torch.Tensor) -> torch.Tensor:
        if hasattr(self.impl, "terminal_gradient_torch"):
            return self.impl.terminal_gradient_torch(x)
        return self.impl.terminal_gradient(x)

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_gradient_torch(x)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        return self.impl.base_drift(t, x)

    def driver(self, value: torch.Tensor, z: torch.Tensor | None) -> torch.Tensor:
        return self.impl.driver(value, z)

    def tangent_coefficients(self, value, z):
        return self.impl.tangent_coefficients(value, z)

    def sample_states(self, *args, **kwargs):
        return self.impl.sample_states(*args, **kwargs)

    @staticmethod
    def working_to_original(value: np.ndarray) -> np.ndarray:
        return np.asarray(value)


@lru_cache(maxsize=None)
def make_benchmark(pde: str, dimension: int):
    d = int(dimension)
    if pde == "hjb_log":
        return StandardAdapter(HJBLogBenchmark(d), pde)
    if pde == "allen_cahn":
        return StandardAdapter(AllenCahnBenchmark(d), pde)
    if pde == "b6_fullrank":
        return StandardAdapter(make_b6(d), pde)
    if pde == "dense_fourier":
        # DenseFourierDriver already implements the complete common interface.
        return make_dense_fourier(d)
    raise KeyError(pde)


def exact_values(pde: str, dimension: int, t: np.ndarray, x: np.ndarray) -> np.ndarray:
    if pde == "hjb_log":
        return np.asarray(hjb_log_reference(t, x, order=256), dtype=np.float64)
    if pde == "allen_cahn":
        raise RuntimeError("Allen-Cahn references require the convergence-audited radial solver")
    benchmark = make_benchmark(pde, dimension)
    implementation = benchmark.impl if isinstance(benchmark, StandardAdapter) else benchmark
    return np.asarray(highdim_exact_u(implementation, t, x), dtype=np.float64)


def benchmark_metadata(pde: str, dimension: int) -> dict[str, Any]:
    benchmark = make_benchmark(pde, dimension)
    if pde in {"b6_fullrank", "dense_fourier"}:
        implementation = benchmark.impl if isinstance(benchmark, StandardAdapter) else benchmark
        return highdim_metadata(implementation)
    return {
        "pde": pde,
        "name": benchmark.name,
        "dimension": benchmark.d,
        "T": benchmark.T,
        "sigma": benchmark.sigma,
        "actor_kind": benchmark.actor_kind,
        "oracle_projection_used": False,
    }


__all__ = [
    "StandardAdapter", "benchmark_metadata", "exact_values", "make_benchmark",
]

