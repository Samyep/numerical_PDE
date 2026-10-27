"""Genuinely nonseparable joint value-gradient benchmark adapter."""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
import math
from typing import Any

import numpy as np
import torch

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.drivers import (
    NonseparableDriver,
)

from .protocol import BASE_SEED, PDE


def frozen_basis(dimension: int) -> np.ndarray:
    d = int(dimension)
    if d < 2:
        raise ValueError(d)
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, d, 17]))
    raw = rng.standard_normal((d, 2))
    basis, triangular = np.linalg.qr(raw, mode="reduced")
    signs = np.where(np.diag(triangular) < 0.0, -1.0, 1.0)
    # Own writable storage: ``torch.as_tensor`` warns (and has undefined
    # mutation semantics) for NumPy views marked read-only.
    return np.array(basis * signs[None, :], dtype=np.float64, copy=True)


class JointNonseparableBenchmark(NonseparableDriver):
    """f(y,z)=y-y^3-0.5|z+y v|^2 with a frozen dense v."""

    actor_kind = "joint"
    input_scale = 1.0
    value_bounds = (0.0, 1.0)

    def __init__(self, dimension: int) -> None:
        d = int(dimension)
        super().__init__(d=d, T=0.5, eta=1.0, basis_seed=BASE_SEED)
        basis = frozen_basis(d)
        self.B = basis
        self.beta = basis[:, 0]
        self.v = basis[:, 0]
        self.pde = PDE
        self.family = PDE
        self.name = f"{PDE}_d{d}"

    def terminal_gradient(self, x: torch.Tensor) -> torch.Tensor:
        return self.terminal_gradient_torch(x)

    def base_drift(self, t: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        del t
        return torch.zeros_like(x)

    def sample_states(
        self, count: int, *, seed: int, purpose: str = "function"
    ) -> tuple[np.ndarray, np.ndarray]:
        if purpose not in {"function", "official"}:
            raise KeyError(purpose)
        rng = np.random.default_rng(int(seed))
        t = rng.uniform(0.0, np.nextafter(self.T, 0.0), int(count))
        x = rng.standard_normal((int(count), self.d))
        return t.astype(np.float64), x.astype(np.float64)

    def metadata(self) -> dict[str, Any]:
        gram = self.B.T @ self.B
        return {
            "pde": PDE,
            "name": self.name,
            "family": self.family,
            "dimension": self.d,
            "effective_dimension": 2,
            "T": self.T,
            "sigma": self.sigma,
            "eta": self.eta,
            "actor_kind": self.actor_kind,
            "driver": "y-y^3-0.5*|z+y*v|^2",
            "mixed_derivative": "nabla_z f_y=-v",
            "basis_seed_sequence": [BASE_SEED, self.d, 17],
            "basis_sha256": sha256(np.ascontiguousarray(self.B).tobytes()).hexdigest(),
            "v_sha256": sha256(np.ascontiguousarray(self.v).tobytes()).hexdigest(),
            "basis_orthogonality_max_abs": float(
                np.max(np.abs(gram - np.eye(2)))
            ),
            "terminal_range": [0.2, 0.8],
            "tangent_domain": [self.r_low, self.r_high],
            "actor_input": "full ambient x",
            "oracle_projection_used": False,
        }


@lru_cache(maxsize=None)
def make_benchmark(dimension: int) -> JointNonseparableBenchmark:
    return JointNonseparableBenchmark(int(dimension))


def sample_states(
    dimension: int, count: int, *, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    return make_benchmark(int(dimension)).sample_states(
        int(count), seed=int(seed), purpose="function"
    )


def reduced_to_ambient_z(
    dimension: int, reduced_z: np.ndarray
) -> np.ndarray:
    return np.asarray(reduced_z, dtype=np.float64) @ frozen_basis(dimension).T


__all__ = [
    "JointNonseparableBenchmark",
    "frozen_basis",
    "make_benchmark",
    "reduced_to_ambient_z",
    "sample_states",
]
