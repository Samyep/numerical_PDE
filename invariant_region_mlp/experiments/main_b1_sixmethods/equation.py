"""Coupled multi-ridge quadratic HJB and exact Hopf--Cole reference.

The exact reference is used only for correctness, reference-only calibration,
and held-out evaluation.  Training code in this package receives only the
terminal condition and its explicit derivative.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Any

import numpy as np
from scipy.special import logsumexp


T_FINAL = 0.25
SIGMA = math.sqrt(2.0)
BASIS_SEED = 20261006  # keeps k=1 identical to the established P1 direction
COEFFICIENT_SEED = 20261012
DEFAULT_COEFFICIENT_SCALE = 2.0


def _canonical_qr(matrix: np.ndarray) -> np.ndarray:
    """Reduced QR with deterministic column signs."""

    q, r = np.linalg.qr(matrix, mode="reduced")
    signs = np.sign(np.diag(r))
    signs[signs == 0.0] = 1.0
    return q * signs[None, :]


def _as_points(
    t: np.ndarray | float, x: np.ndarray, d: int
) -> tuple[np.ndarray, np.ndarray]:
    state = np.asarray(x, dtype=np.float64)
    if state.shape[-1] != d:
        raise ValueError(f"last coordinate must have length {d}, got {state.shape}")
    time = np.asarray(t, dtype=np.float64)
    return np.broadcast_to(time, state.shape[:-1]), state


def effective_rank(eigenvalues: np.ndarray, *, tolerance: float = 1e-14) -> float:
    """Entropy effective rank of nonnegative eigenvalues."""

    values = np.maximum(np.asarray(eigenvalues, dtype=np.float64), 0.0)
    values = values[values > tolerance * max(float(values.max(initial=0.0)), 1.0)]
    if not len(values) or float(values.sum()) == 0.0:
        return 0.0
    probabilities = values / values.sum()
    return float(np.exp(-np.sum(probabilities * np.log(probabilities))))


@dataclass(frozen=True)
class B1Equation:
    """B1 with dense active basis and dense affine-spanning coefficients."""

    d: int
    k: int
    coefficient_scale: float = DEFAULT_COEFFICIENT_SCALE
    basis_seed: int = BASIS_SEED
    coefficient_seed: int = COEFFICIENT_SEED
    T: float = T_FINAL
    sigma: float = SIGMA
    name: str = "B1_coupled_multiridge_quadratic_hjb"
    family: str = "multiridge"

    def __post_init__(self) -> None:
        if self.d <= 0 or self.k <= 0 or self.k > self.d:
            raise ValueError("require 1 <= k <= d")
        if self.coefficient_scale <= 0.0:
            raise ValueError("coefficient_scale must be positive")

        if self.k == 1:
            # Byte-for-byte generation rule used by RidgeLSEHJB/P1.
            rng = np.random.default_rng(
                np.random.SeedSequence([self.basis_seed, self.d, 1])
            )
            basis = rng.standard_normal((self.d, 1), dtype=np.float64)
            basis /= np.linalg.norm(basis, axis=0, keepdims=True)
            coefficients = np.array([[0.5], [1.5], [3.0]], dtype=np.float64)
        else:
            basis_rng = np.random.default_rng(
                np.random.SeedSequence([self.basis_seed, self.d, self.k])
            )
            basis = _canonical_qr(
                basis_rng.standard_normal((self.d, self.k), dtype=np.float64)
            )

            # Coefficients depend on k, not d, so ambient-dimension comparisons
            # retain exactly the same active-coordinate terminal geometry.
            coeff_rng = np.random.default_rng(
                np.random.SeedSequence([self.coefficient_seed, self.k, 2 * self.k + 1])
            )
            coefficients = coeff_rng.standard_normal(
                (2 * self.k + 1, self.k), dtype=np.float64
            )
            coefficients -= coefficients.mean(axis=0, keepdims=True)
            rms_norm = float(np.sqrt(np.mean(np.sum(coefficients**2, axis=1))))
            coefficients *= self.coefficient_scale / rms_norm

        affine_rank = int(
            np.linalg.matrix_rank(coefficients[1:] - coefficients[:1], tol=1e-12)
        )
        if affine_rank != self.k:
            raise RuntimeError(
                f"coefficient construction lost affine rank: {affine_rank} != {self.k}"
            )

        directions = coefficients @ basis.T
        object.__setattr__(self, "B", basis)
        object.__setattr__(self, "b", coefficients)
        object.__setattr__(self, "K", len(coefficients))
        # Compatibility names consumed by the established GPU MLP adapter.
        object.__setattr__(self, "A", directions)
        object.__setattr__(self, "logc", np.full(len(coefficients), -math.log(len(coefficients))))
        object.__setattr__(self, "n2", np.sum(directions * directions, axis=1))
        object.__setattr__(self, "mu", 0.0)
        object.__setattr__(self, "affine_rank", affine_rank)

    def _logits(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        time, state = _as_points(t, x, self.d)
        return (
            self.logc
            + np.einsum("...d,md->...m", state, self.A)
            + (self.T - time)[..., None] * self.n2
        )

    def terminal_logits(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        return self.logc + np.einsum("...d,md->...m", state, self.A)

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return -logsumexp(self._logits(t, x), axis=-1)

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return -logsumexp(self.terminal_logits(x), axis=-1)

    def weights(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        logits = self._logits(t, x)
        return np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))

    def terminal_weights(self, x: np.ndarray) -> np.ndarray:
        logits = self.terminal_logits(x)
        return np.exp(logits - logsumexp(logits, axis=-1, keepdims=True))

    def exact_grad(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return -self.weights(t, x) @ self.A

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.exact_grad(t, x)

    def terminal_grad(self, x: np.ndarray) -> np.ndarray:
        return -self.terminal_weights(x) @ self.A

    def terminal_z(self, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.terminal_grad(x)

    def exact_state(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        value = self.exact_u(t, x)
        return np.concatenate((value[..., None], self.exact_z(t, x)), axis=-1)

    @staticmethod
    def generator(_u: np.ndarray, z: np.ndarray) -> np.ndarray:
        vector = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(vector * vector, axis=-1)

    def exact_generator(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.generator(self.exact_u(t, x), self.exact_z(t, x))

    @property
    def box_low(self) -> np.ndarray:
        return np.min(-self.sigma * self.A, axis=0)

    @property
    def box_high(self) -> np.ndarray:
        return np.max(-self.sigma * self.A, axis=0)

    @property
    def z_ball_radius(self) -> float:
        # Every exact z is a convex combination of these vertices.
        return self.sigma * float(np.sqrt(np.max(self.n2)))

    @property
    def z_center(self) -> np.ndarray:
        return 0.5 * (self.box_low + self.box_high)

    @property
    def dose_scale(self) -> float:
        return float(np.mean(0.5 * (self.box_high - self.box_low)))

    def exact_active_hessian(
        self, t: np.ndarray | float, x: np.ndarray
    ) -> np.ndarray:
        probabilities = self.weights(t, x)
        mean = probabilities @ self.b
        second = np.einsum("...m,mi,mj->...ij", probabilities, self.b, self.b)
        covariance = second - np.einsum("...i,...j->...ij", mean, mean)
        return -covariance

    def analytic_residual(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        probabilities = self.weights(t, x)
        mean_a = probabilities @ self.A
        mean_norm2 = probabilities @ self.n2
        u_t = mean_norm2
        laplacian = -(mean_norm2 - np.sum(mean_a * mean_a, axis=-1))
        z = -self.sigma * mean_a
        return u_t + laplacian + self.generator(np.zeros_like(u_t), z)

    def heat_only_mc(
        self,
        t: np.ndarray | float,
        x: np.ndarray,
        *,
        n_samples: int,
        seed: int,
        chunk_size: int = 256,
    ) -> np.ndarray:
        """Antithetic MC for E[g(X_T)] in active coordinates."""

        time, state = _as_points(t, x, self.d)
        flat_x = state.reshape(-1, self.d)
        flat_t = time.reshape(-1)
        if n_samples < 2:
            raise ValueError("n_samples must be at least two")
        half = (n_samples + 1) // 2
        rng = np.random.default_rng(seed)
        output = np.empty(len(flat_x), dtype=np.float64)
        for start in range(0, len(flat_x), chunk_size):
            stop = min(start + chunk_size, len(flat_x))
            y = flat_x[start:stop] @ self.B
            normal = rng.standard_normal((stop - start, half, self.k), dtype=np.float64)
            normal = np.concatenate((normal, -normal), axis=1)[:, :n_samples]
            tau = np.maximum(self.T - flat_t[start:stop], 0.0)
            terminal_y = y[:, None, :] + np.sqrt(2.0 * tau)[:, None, None] * normal
            logits = self.logc + np.einsum("bnk,mk->bnm", terminal_y, self.b)
            values = -logsumexp(logits, axis=-1)
            output[start:stop] = values.mean(axis=1)
        return output.reshape(state.shape[:-1])

    def provenance(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "family": self.family,
            "d": self.d,
            "k": self.k,
            "K": self.K,
            "T": self.T,
            "sigma": self.sigma,
            "mu": self.mu,
            "coefficient_scale": self.coefficient_scale if self.k > 1 else None,
            "basis_seed": self.basis_seed,
            "coefficient_seed": self.coefficient_seed if self.k > 1 else None,
            "affine_rank": self.affine_rank,
            "basis_sha256": hashlib.sha256(self.B.tobytes()).hexdigest(),
            "coefficient_sha256": hashlib.sha256(self.b.tobytes()).hexdigest(),
            "direction_sha256": hashlib.sha256(self.A.tobytes()).hexdigest(),
            "coefficient_rule": (
                "P1 fixed (0.5,1.5,3)" if self.k == 1 else
                "centered Gaussian Kxk matrix, globally rescaled to requested RMS row norm"
            ),
        }


def make_equation(d: int, k: int, coefficient_scale: float = DEFAULT_COEFFICIENT_SCALE) -> B1Equation:
    return B1Equation(d=d, k=k, coefficient_scale=coefficient_scale)
