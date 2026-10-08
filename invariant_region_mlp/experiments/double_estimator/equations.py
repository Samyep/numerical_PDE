"""Equation registry and fresh fixed points for the double-estimator study.

Only this new experiment directory is extended.  The established equation
implementations are imported without modification; the multi-ridge equation
is reproduced here so importing it does not monkey-patch the shared
mechanism-suite projection functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import math
from pathlib import Path
from typing import Any

import numpy as np

from invariant_region_mlp.experiments.candidate_screen.candidates import LQGameHJB
from invariant_region_mlp.experiments.mechanism_suite.equations import (
    ExactEquation,
    RidgeLSEHJB,
    _as_points,
    stable_seed,
)
from invariant_region_mlp.experiments.mechanism_suite.norm_hjb import NormDriverHJB


BASE_SEED = 20261201
POINT_COUNT = 1200
SQRT2 = math.sqrt(2.0)
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
P4_REFERENCE = (
    PACKAGE_ROOT
    / "results"
    / "mechanism_suite"
    / "reference_cache"
    / "norm_hjb_beta2_lambda1_T0p5.npz"
)


def _logsumexp(values: np.ndarray) -> np.ndarray:
    maximum = np.max(values, axis=-1, keepdims=True)
    return maximum[..., 0] + np.log(np.sum(np.exp(values - maximum), axis=-1))


@dataclass(frozen=True)
class MultiRidgeHJB(ExactEquation):
    """The fixed k=10, scale=4 effective-dimension benchmark."""

    d: int
    k: int = 10
    scale: float = 4.0
    T: float = 0.1
    seed: int = 0
    n_grad_samples: int = 4000
    sigma: float = SQRT2
    name: str = "MR_multiridge_lse_hjb"
    family: str = "multiridge"

    def __post_init__(self) -> None:
        if self.d <= self.k:
            raise ValueError("multi-ridge requires d > k")
        # Keep the PDE itself identical to effdim_mlp.py.  BASE_SEED is used
        # for all newly sampled points and MLP trees, not to redefine the PDE.
        rng = np.random.default_rng(
            np.random.SeedSequence([20261008, self.k, self.d, self.seed])
        )
        vertices = rng.standard_normal((2 * self.k, self.k), dtype=np.float64)
        vertices /= np.linalg.norm(vertices, axis=1, keepdims=True)
        vertices *= self.scale * rng.uniform(0.5, 1.5, (2 * self.k, 1))
        basis, _ = np.linalg.qr(rng.standard_normal((self.d, self.k), dtype=np.float64))
        directions = vertices @ basis.T
        object.__setattr__(self, "B", vertices)
        object.__setattr__(self, "Q", basis)
        object.__setattr__(self, "A", directions)
        object.__setattr__(self, "mu", 0.0)
        object.__setattr__(self, "logc", -math.log(len(vertices)) * np.ones(len(vertices)))
        object.__setattr__(self, "n2", np.sum(directions * directions, axis=1))

        # Reproduce the data-only subspace estimate and its certified hull.
        samples = np.random.default_rng(99).uniform(
            -1.0, 1.0, (self.n_grad_samples, self.d)
        )
        logits = self.logc + samples @ directions.T
        weights = np.exp(logits - _logsumexp(logits)[:, None])
        gradients = -(weights @ directions)
        left, singular_values, _ = np.linalg.svd(gradients.T, full_matrices=False)
        estimated_basis = left[:, : self.k]
        cosines = np.linalg.svd(basis.T @ estimated_basis, compute_uv=False)
        coordinates = -SQRT2 * vertices @ (basis.T @ estimated_basis)
        ambient_vertices = -SQRT2 * directions
        object.__setattr__(self, "Qh", estimated_basis)
        object.__setattr__(
            self,
            "max_principal_angle",
            float(np.arccos(np.clip(np.min(cosines), -1.0, 1.0))),
        )
        gap = (
            float(singular_values[self.k - 1] / max(singular_values[self.k], 1e-300))
            if len(singular_values) > self.k
            else float("inf")
        )
        object.__setattr__(self, "sv_gap", gap)
        object.__setattr__(self, "sub_lo", np.min(coordinates, axis=0))
        object.__setattr__(self, "sub_hi", np.max(coordinates, axis=0))
        object.__setattr__(self, "_box_lo", np.min(ambient_vertices, axis=0))
        object.__setattr__(self, "_box_hi", np.max(ambient_vertices, axis=0))

    def _logits(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        return self.logc + state @ self.A.T + (self.T - time)[..., None] * self.n2

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return -_logsumexp(self._logits(t, x))

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        weights = np.exp(logits - _logsumexp(logits)[..., None])
        return -SQRT2 * (weights @ self.A)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return -_logsumexp(self.logc + np.asarray(x, dtype=np.float64) @ self.A.T)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    def fzero_reference(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        # Not needed for verdicts; f_zero is evaluated by the same terminal
        # Monte Carlo estimator as all other methods.
        raise NotImplementedError

    @property
    def box_low(self) -> np.ndarray:
        return self._box_lo

    @property
    def box_high(self) -> np.ndarray:
        return self._box_hi

    @property
    def z_ball_radius(self) -> float:
        return SQRT2 * float(np.sqrt(np.max(self.n2)))

    @property
    def z_center(self) -> np.ndarray:
        # This is the established data-free ambient hull centre used by the
        # effective-dimension screen.
        return 0.5 * (self._box_lo + self._box_hi)

    @property
    def sub_box_center(self) -> np.ndarray:
        return (0.5 * (self.sub_lo + self.sub_hi)) @ self.Qh.T

    @property
    def dose_scale(self) -> float:
        return float(np.mean(0.5 * (self._box_hi - self._box_lo)))

    @property
    def picard_proxy(self) -> float:
        return self.z_ball_radius * self.T

    def certificate_labels(self) -> dict[str, str]:
        return {
            "sub_box": "PDE-derived Hopf-Cole hull in a data-estimated subspace",
            "box": "ambient coordinatewise hull of -sqrt(2) a_j",
        }

    def certificate_parameters(self) -> dict[str, Any]:
        return {
            "k": self.k,
            "scale": self.scale,
            "seed": self.seed,
            "max_principal_angle": self.max_principal_angle,
            "sv_gap": self.sv_gap,
            "A_sha256": hashlib.sha256(self.A.tobytes()).hexdigest(),
        }


C2_CONFIGS: dict[str, dict[str, float]] = {
    "C2-convex": {"a": 2.0, "b": 0.0, "frac_A": 0.5},
    "C2-cancel": {"a": 3.0, "b": 1.0, "frac_A": 0.25},
    "C2-flip": {"a": 3.0, "b": 1.0, "frac_A": 0.125},
}


@lru_cache(maxsize=None)
def make_equation(pde_id: str, d: int) -> ExactEquation:
    """Return a cached immutable equation instance."""

    if pde_id == "P1":
        return RidgeLSEHJB(d=d)
    if pde_id == "MR":
        return MultiRidgeHJB(d=d)
    if pde_id in C2_CONFIGS:
        config = C2_CONFIGS[pde_id]
        return LQGameHJB(
            d=d,
            a=config["a"],
            b=config["b"],
            frac_A=config["frac_A"],
            direction_seed=20261207,
            name=pde_id,
        )
    if pde_id == "P4":
        if not P4_REFERENCE.is_file():
            raise FileNotFoundError(f"missing frozen P4 reference: {P4_REFERENCE}")
        return NormDriverHJB(
            d=d,
            reference_path=str(P4_REFERENCE),
            T=0.5,
        )
    raise ValueError(f"unknown PDE {pde_id!r}")


@lru_cache(maxsize=None)
def fixed_points(pde_id: str, d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create the new 1,200-point sample and fixed 20% validation split."""

    equation = make_equation(pde_id, d)
    rng = np.random.default_rng(
        np.random.SeedSequence(
            [BASE_SEED, d, POINT_COUNT, stable_seed(pde_id)]
        )
    )
    t = rng.uniform(0.0, np.nextafter(equation.T, 0.0), POINT_COUNT).astype(np.float64)
    x = rng.uniform(-1.0, 1.0, (POINT_COUNT, d)).astype(np.float64)
    is_validation = np.zeros(POINT_COUNT, dtype=bool)
    selected = rng.choice(POINT_COUNT, size=POINT_COUNT // 5, replace=False)
    is_validation[selected] = True
    # Cached arrays are treated as immutable by callers.
    t.setflags(write=False)
    x.setflags(write=False)
    is_validation.setflags(write=False)
    return t, x, is_validation


def equation_kind(equation: Any) -> str:
    if hasattr(equation, "dA") and hasattr(equation, "a") and hasattr(equation, "b"):
        return "game"
    if equation.family == "norm_hjb":
        return "norm"
    if equation.family in {"ridge_lse", "multiridge"}:
        return "quadratic"
    raise ValueError(f"double generator is unavailable for family {equation.family!r}")

