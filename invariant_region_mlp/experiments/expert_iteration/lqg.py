"""LQG benchmark, deterministic audit reference, and preregistered MC gate.

The experiment itself uses the antithetic Monte Carlo reference required by
PREREGISTRATION.md.  A scaled Gauss--Laguerre representation of the Hopf--Cole
formula is also provided for states created inside the recursive solver, where
running an adaptive Monte Carlo reference at every child would be prohibitive.
The deterministic evaluator is independently checked against the MC reference.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy.special import roots_laguerre


BASE_SEED = 20261101
STUDY_A = 1
PROBLEM_LQG = 1
SQRT2 = math.sqrt(2.0)


def seed_sequence(d: int, seed: int, component: int, *extra: int) -> np.random.SeedSequence:
    """Return a component-isolated seed following the frozen scheme."""

    return np.random.SeedSequence(
        [BASE_SEED, STUDY_A, PROBLEM_LQG, int(d), int(seed), int(component), *map(int, extra)]
    )


def uniform_unit_ball(rng: np.random.Generator, count: int, d: int) -> np.ndarray:
    direction = rng.standard_normal((count, d), dtype=np.float64)
    direction /= np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-300)
    radius = rng.random(count, dtype=np.float64) ** (1.0 / d)
    return direction * radius[:, None]


@dataclass
class LQGEquation:
    """The preregistered SCaSML LQG interpretation."""

    d: int
    c1: np.ndarray
    c2: np.ndarray
    T: float = 0.5
    lam: float = 1.0
    quadrature_order: int = 64

    family: str = "lqg"
    name: str = "SCaSML_LQG_interpretation"
    mu: float = 0.0
    sigma: float = SQRT2

    def __post_init__(self) -> None:
        self.c1 = np.asarray(self.c1, dtype=np.float64)
        self.c2 = np.asarray(self.c2, dtype=np.float64)
        if self.c1.shape != (self.d - 1,) or self.c2.shape != (self.d - 1,):
            raise ValueError("c1 and c2 must have shape (d-1,)")
        # Symmetric tridiagonal matrix A with q(x)=1+x^T A x.
        diagonal = np.zeros(self.d, dtype=np.float64)
        diagonal[:-1] += self.c1
        diagonal[1:] += self.c1 + self.c2
        off = -self.c1
        A = np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1)
        eigenvalues, eigenvectors = np.linalg.eigh(A)
        if np.min(eigenvalues) <= 0.0:
            raise ValueError("terminal quadratic form must be positive definite")
        self.A = A
        self.trace_A = float(np.trace(A))
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors
        nodes, weights = roots_laguerre(self.quadrature_order)
        self.quad_nodes = np.asarray(nodes, dtype=np.float64)
        self.quad_weights = np.asarray(weights, dtype=np.float64)

    @classmethod
    def create(cls, d: int, *, quadrature_order: int = 64) -> "LQGEquation":
        rng = np.random.default_rng(seed_sequence(d, 0, 1))
        c1 = rng.uniform(0.5, 1.5, d - 1)
        c2 = rng.uniform(0.5, 1.5, d - 1)
        return cls(d=d, c1=c1, c2=c2, quadrature_order=quadrature_order)

    def quadratic(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        return np.sum(
            self.c1 * (x[..., :-1] - x[..., 1:]) ** 2
            + self.c2 * x[..., 1:] ** 2,
            axis=-1,
        )

    def terminal(self, x: np.ndarray) -> np.ndarray:
        return np.log((1.0 + self.quadratic(x)) / 2.0)

    def apply_A(self, x: np.ndarray) -> np.ndarray:
        """Apply the terminal tridiagonal quadratic form in O(d) work."""

        x = np.asarray(x, dtype=np.float64)
        result = np.zeros_like(x)
        difference = x[..., :-1] - x[..., 1:]
        result[..., :-1] += self.c1 * difference
        result[..., 1:] -= self.c1 * difference
        result[..., 1:] += self.c2 * x[..., 1:]
        return result

    def grad_terminal(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        q = 1.0 + self.quadratic(x)
        return 2.0 * self.apply_A(x) / q[..., None]

    def sigma_grad_terminal(self, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.grad_terminal(x)

    @staticmethod
    def generator(u: np.ndarray, z: np.ndarray) -> np.ndarray:
        del u
        return -0.5 * np.sum(np.asarray(z, dtype=np.float64) ** 2, axis=-1)

    def hopf_cole(self, t: np.ndarray, x: np.ndarray, *, chunk_size: int = 256) -> tuple[np.ndarray, np.ndarray]:
        """Evaluate (u,z) by a scaled one-dimensional Laplace integral.

        Scaling by 1+E[X^T A X] keeps Gauss--Laguerre nodes in the mass of the
        integral even at d=160.  At h=0 the formula reduces exactly to g and
        its gradient (up to roundoff).
        """

        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.d)
        out_u = np.empty(len(t), dtype=np.float64)
        out_z = np.empty((len(t), self.d), dtype=np.float64)
        eig = self.eigenvalues
        Q = self.eigenvectors
        nodes = self.quad_nodes
        log_weights = np.log(self.quad_weights)
        for begin in range(0, len(t), chunk_size):
            end = min(begin + chunk_size, len(t))
            tb = t[begin:end]
            xb = x[begin:end]
            h = np.maximum(self.T - tb, 0.0)
            terminal = h <= 8.0 * np.finfo(np.float64).eps
            if np.any(terminal):
                out_u[begin:end][terminal] = self.terminal(xb[terminal])
                out_z[begin:end][terminal] = self.sigma_grad_terminal(xb[terminal])
            active = ~terminal
            if not np.any(active):
                continue
            xa = xb[active]
            ha = h[active]
            y = xa @ Q
            mean_q = np.sum(eig[None, :] * y * y, axis=1) + 2.0 * ha * self.trace_A
            rate = 1.0 + mean_q
            s = nodes[None, :] / rate[:, None]
            denom = 1.0 + 4.0 * ha[:, None, None] * s[:, :, None] * eig[None, None, :]
            log_phi = -0.5 * np.sum(np.log(denom), axis=2)
            log_phi -= s * np.sum(
                eig[None, None, :] * y[:, None, :] ** 2 / denom,
                axis=2,
            )
            # Change of variable r=rate*s inside the Laguerre integral.
            log_terms = (
                log_weights[None, :]
                - np.log(rate)[:, None]
                + nodes[None, :] * (1.0 - 1.0 / rate[:, None])
                + log_phi
            )
            maximum = np.max(log_terms, axis=1, keepdims=True)
            terms = np.exp(log_terms - maximum)
            integral = np.exp(maximum[:, 0]) * np.sum(terms, axis=1)
            weights_normalized = terms / np.sum(terms, axis=1, keepdims=True)
            coefficient = 2.0 * s[:, :, None] * eig[None, None, :] / denom
            grad_eig = np.sum(
                weights_normalized[:, :, None] * coefficient * y[:, None, :],
                axis=1,
            )
            ua = -np.log(2.0 * integral) / self.lam
            za = self.sigma * (grad_eig @ Q.T) / self.lam
            block_u = out_u[begin:end]
            block_z = out_z[begin:end]
            block_u[active] = ua
            block_z[active] = za
        return out_u, out_z

    def exact_state(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        u, z = self.hopf_cole(t, x)
        return np.concatenate([u[:, None], z], axis=1)

    def exact_generator(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        _, z = self.hopf_cole(np.asarray(t).reshape(-1), np.asarray(x).reshape(-1, self.d))
        return self.generator(np.zeros(len(z)), z).reshape(np.asarray(t).shape)

    def provenance(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "dimension": self.d,
            "T": self.T,
            "lambda": self.lam,
            "mu": self.mu,
            "sigma": self.sigma,
            "quadrature_order": self.quadrature_order,
            "coefficient_seed": seed_sequence(self.d, 0, 1).entropy,
            "interpretation": "lambda=1, no drift, T=0.5, unit-ball evaluation",
        }


def make_evaluation_points(eq: LQGEquation) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 200 validation then 1000 test points, fixed across all methods."""

    rng = np.random.default_rng(seed_sequence(eq.d, 0, 2))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), 1200)
    x = uniform_unit_ball(rng, 1200, eq.d)
    is_validation = np.zeros(1200, dtype=bool)
    is_validation[:200] = True
    return t, x, is_validation


@dataclass
class ReferenceAccumulator:
    count_pairs: int
    sum_w: np.ndarray
    sum_w2: np.ndarray
    sum_g: np.ndarray
    sum_grad_w: np.ndarray


def _mc_pair_batch(eq: LQGEquation, x: np.ndarray, h: np.ndarray, normal: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Common random numbers across evaluation points are legitimate marginal
    # MC samples and substantially reduce random-number and A-application work.
    # Expand Q(x +/- s*xi)=Q(x)+s^2 Q(xi)+/-2s<x,Axi>.
    scale = eq.sigma * np.sqrt(h)
    Ax = eq.apply_A(x)
    An = eq.apply_A(normal)
    xAx = eq.quadratic(x)[:, None]
    nAn = np.sum(normal * An, axis=-1)
    xAn = Ax @ normal.T
    aAa = scale[:, None] ** 2 * nAn[None, :]
    xAa = scale[:, None] * xAn
    q_plus = 1.0 + xAx + aAa + 2.0 * xAa
    q_minus = 1.0 + xAx + aAa - 2.0 * xAa
    w_plus = 2.0 / q_plus
    w_minus = 2.0 / q_minus
    w_pair = 0.5 * (w_plus + w_minus)
    g_pair = 0.5 * (np.log(q_plus / 2.0) + np.log(q_minus / 2.0))
    inverse_plus = 1.0 / q_plus**2
    inverse_minus = 1.0 / q_minus**2
    grad_pair_sum = 2.0 * (
        Ax * np.sum(inverse_plus + inverse_minus, axis=1)[:, None]
        + scale[:, None] * ((inverse_plus - inverse_minus) @ An)
    )
    return w_pair, g_pair, grad_pair_sum


def antithetic_reference(
    eq: LQGEquation,
    t: np.ndarray,
    x: np.ndarray,
    *,
    tolerance: float = 1e-4,
    initial_pairs: int = 16_384,
    batch_pairs: int = 16_384,
    max_pairs: int = 16_777_216,
    point_chunk: int = 8,
    progress: Any | None = None,
) -> dict[str, np.ndarray]:
    """Adaptive antithetic MC reference satisfying the per-point u SE gate."""

    t = np.asarray(t, dtype=np.float64).reshape(-1)
    x = np.asarray(x, dtype=np.float64).reshape(len(t), eq.d)
    u = np.empty(len(t), dtype=np.float64)
    z = np.empty((len(t), eq.d), dtype=np.float64)
    f_zero = np.empty(len(t), dtype=np.float64)
    se = np.empty(len(t), dtype=np.float64)
    pairs_used = np.empty(len(t), dtype=np.int64)
    h_all = np.maximum(eq.T - t, 0.0)
    for chunk_index, begin in enumerate(range(0, len(t), point_chunk)):
        end = min(begin + point_chunk, len(t))
        xb = x[begin:end]
        hb = h_all[begin:end]
        b = end - begin
        rng = np.random.default_rng(seed_sequence(eq.d, 0, 20, chunk_index))
        acc = ReferenceAccumulator(
            count_pairs=0,
            sum_w=np.zeros(b),
            sum_w2=np.zeros(b),
            sum_g=np.zeros(b),
            sum_grad_w=np.zeros((b, eq.d)),
        )
        current_se = np.full(b, np.inf)
        while acc.count_pairs < max_pairs and (
            acc.count_pairs < initial_pairs or np.any(current_se > tolerance)
        ):
            take = min(batch_pairs, max_pairs - acc.count_pairs)
            normal = rng.standard_normal((take, eq.d), dtype=np.float64)
            wp, gp, grad_pair_sum = _mc_pair_batch(eq, xb, hb, normal)
            acc.sum_w += np.sum(wp, axis=1)
            acc.sum_w2 += np.sum(wp * wp, axis=1)
            acc.sum_g += np.sum(gp, axis=1)
            acc.sum_grad_w += grad_pair_sum
            acc.count_pairs += take
            mean_w = acc.sum_w / acc.count_pairs
            variance_w = np.maximum(
                (acc.sum_w2 - acc.sum_w * acc.sum_w / acc.count_pairs)
                / max(acc.count_pairs - 1, 1),
                0.0,
            )
            current_se = np.sqrt(variance_w / acc.count_pairs) / np.maximum(mean_w, 1e-300)
        mean_w = acc.sum_w / acc.count_pairs
        u[begin:end] = -np.log(mean_w) / eq.lam
        z[begin:end] = eq.sigma * (acc.sum_grad_w / acc.count_pairs) / mean_w[:, None]
        f_zero[begin:end] = acc.sum_g / acc.count_pairs
        se[begin:end] = current_se
        pairs_used[begin:end] = acc.count_pairs
        if progress is not None:
            progress(end, len(t), float(np.max(current_se)), acc.count_pairs)
    return {
        "u": u,
        "z": z,
        "f_zero_u": f_zero,
        "u_se": se,
        "pairs": pairs_used,
        "terminal_samples": 2 * pairs_used,
    }


def relative_l2(prediction: np.ndarray, truth: np.ndarray) -> float:
    prediction = np.asarray(prediction, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    return float(np.linalg.norm(prediction - truth) / np.linalg.norm(truth))
