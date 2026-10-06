"""Published high-dimensional gradient-dependent nonlinear benchmark.

The published SCaSML benchmark uses ``sigma=sqrt(2)``.  The public SCaSML
repository currently uses ``sigma=0.25``; that value is deliberately *not*
the default here and must only be requested explicitly for a diagnostic.

The MLP state convention throughout this experiment is ``(u, z)`` with
``z = sigma * grad(u)``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np


PUBLISHED_SIGMA = math.sqrt(2.0)
PUBLIC_REPOSITORY_SIGMA = 0.25
SCASML_REPOSITORY_URL = "https://github.com/Francis-Fan-create/SCaSML"
SCASML_AUDITED_COMMIT = "c2412b942b1720084accf9ac9651c6085e073919"
SCASML_REPORT_URL = "https://2prime.github.io/scasml_techreport.pdf"


def stable_sigmoid(value: np.ndarray | float) -> np.ndarray:
    """Numerically stable logistic sigmoid in float64."""

    value = np.asarray(value, dtype=np.float64)
    return np.exp(-np.logaddexp(0.0, -value))


@dataclass(frozen=True)
class ViscousBurgersEquation:
    """Published semilinear PDE and its exact solution."""

    d: int
    T: float = 0.5
    sigma: float = PUBLISHED_SIGMA

    def __post_init__(self) -> None:
        if self.d <= 0:
            raise ValueError("d must be positive")
        if self.T <= 0.0:
            raise ValueError("T must be positive")
        if self.sigma <= 0.0:
            raise ValueError("sigma must be positive")

    @property
    def mu(self) -> float:
        return -(1.0 / self.d + 0.5 * self.sigma**2)

    @property
    def z_upper(self) -> float:
        return self.sigma / 4.0

    @property
    def z_ball_radius(self) -> float:
        return self.sigma * math.sqrt(self.d) / 4.0

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return stable_sigmoid(np.asarray(t, dtype=np.float64) + np.sum(x, axis=-1))

    def exact_q(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        return u * (1.0 - u)

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        q = self.exact_q(t, x)
        return np.broadcast_to((self.sigma * q)[..., None], np.asarray(x).shape).copy()

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        z = self.exact_z(t, x)
        return np.concatenate((u[..., None], z), axis=-1)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return self.exact_u(self.T, x)

    def generator(self, u: np.ndarray, z: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        return self.sigma * u * np.sum(z, axis=-1)

    def exact_generator(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        u = self.exact_u(t, x)
        # Every exact coordinate equals sigma*u*(1-u).
        return self.sigma**2 * self.d * u**2 * (1.0 - u)

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        """Evaluate the PDE residual using analytic derivatives of the exact u."""

        u = self.exact_u(t, x)
        q = u * (1.0 - u)
        u_t = q
        sum_grad = self.d * q
        laplacian = self.d * q * (1.0 - 2.0 * u)
        z = np.broadcast_to((self.sigma * q)[..., None], np.asarray(x).shape)
        return (
            u_t
            + self.mu * sum_grad
            + 0.5 * self.sigma**2 * laplacian
            + self.generator(u, z)
        )

    def fzero_reference(
        self,
        t: np.ndarray,
        x: np.ndarray,
        quadrature_order: int = 80,
    ) -> np.ndarray:
        """High-accuracy solution after deleting f, via Gauss-Hermite quadrature.

        The terminal condition depends only on the coordinate sum, so the
        d-dimensional Gaussian expectation reduces exactly to one dimension.
        ``hermgauss`` integrates against exp(-s^2); the sqrt(2) rescaling below
        converts it to a standard-normal expectation.
        """

        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        tau = np.maximum(self.T - t, 0.0)
        mean_argument = self.T + np.sum(x, axis=1) + self.d * self.mu * tau
        std_argument = self.sigma * np.sqrt(self.d * tau)
        arguments = mean_argument[:, None] + math.sqrt(2.0) * std_argument[:, None] * nodes[None, :]
        return np.sum(stable_sigmoid(arguments) * weights[None, :], axis=1) / math.sqrt(math.pi)

    def provenance(self) -> dict[str, Any]:
        return {
            "benchmark": "published SCaSML gradient-dependent nonlinear / viscous-Burgers-like PDE",
            "dimension": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "mu": self.mu,
            "state_convention": "z = sigma * grad(u)",
            "published_report": SCASML_REPORT_URL,
            "public_repository": SCASML_REPOSITORY_URL,
            "public_repository_audited_commit": SCASML_AUDITED_COMMIT,
            "public_repository_sigma": PUBLIC_REPOSITORY_SIGMA,
        }


def make_test_points(
    d: int,
    *,
    n_interior: int = 1000,
    n_boundary: int = 200,
    seed: int = 20261006,
    validation_fraction: float = 0.2,
) -> dict[str, np.ndarray]:
    """Create the fixed published-domain test geometry and stratified split."""

    if not 0.0 <= validation_fraction < 1.0:
        raise ValueError("validation_fraction must lie in [0,1)")
    rng = np.random.default_rng(np.random.SeedSequence([seed, d, n_interior, n_boundary]))
    xi = rng.uniform(-0.5, 0.5, size=(n_interior, d))
    xb = rng.uniform(-0.5, 0.5, size=(n_boundary, d))
    if n_boundary:
        coordinate = rng.integers(0, d, size=n_boundary)
        sign = rng.choice(np.array([-0.5, 0.5]), size=n_boundary)
        xb[np.arange(n_boundary), coordinate] = sign
    x = np.concatenate((xi, xb), axis=0)
    t = rng.uniform(0.0, 0.5, size=len(x))
    is_boundary = np.concatenate(
        (np.zeros(n_interior, dtype=bool), np.ones(n_boundary, dtype=bool))
    )
    is_validation = np.zeros(len(x), dtype=bool)
    for indices in (np.arange(n_interior), np.arange(n_interior, len(x))):
        count = int(round(validation_fraction * len(indices)))
        if count:
            chosen = rng.choice(indices, size=count, replace=False)
            is_validation[chosen] = True
    return {
        "t": t,
        "x": x,
        "is_boundary": is_boundary,
        "is_validation": is_validation,
    }


def distribution_summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "q10": float(np.quantile(values, 0.1)),
        "q90": float(np.quantile(values, 0.9)),
        "max": float(np.max(values)),
    }


def nonlinear_signal_summary(
    d: int,
    *,
    n_samples: int = 100_000,
    seed: int = 20261006,
    sigma: float = PUBLISHED_SIGMA,
) -> dict[str, Any]:
    """Evaluate the exact nonlinear signal and wrong-PDE bias on a fixed sample."""

    equation = ViscousBurgersEquation(d=d, sigma=sigma)
    rng = np.random.default_rng(np.random.SeedSequence([seed, d, n_samples]))
    x = rng.uniform(-0.5, 0.5, size=(n_samples, d))
    t = rng.uniform(0.0, equation.T, size=n_samples)
    u = equation.exact_u(t, x)
    z = equation.exact_z(t, x)
    z_norm = np.linalg.norm(z, axis=1)
    f = equation.generator(u, z)
    fzero = equation.fzero_reference(t, x)
    error = fzero - u
    return {
        "d": d,
        "n_samples": n_samples,
        "seed": seed,
        "sigma": sigma,
        "u": distribution_summary(u),
        "z_l2": distribution_summary(z_norm),
        "f": distribution_summary(f),
        "abs_f": distribution_summary(np.abs(f)),
        "f_zero_error": distribution_summary(np.abs(f)),
        "fzero_relative_l2": float(np.linalg.norm(error) / np.linalg.norm(u)),
        "fzero_mae": float(np.mean(np.abs(error))),
        "fzero_max_abs": float(np.max(np.abs(error))),
        "pde_residual_max_abs": float(np.max(np.abs(equation.analytic_residual(t, x)))),
    }
