"""Independent numerical references and convergence audits."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any

import numpy as np
from scipy.integrate import quad
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import factorized
from scipy.special import logsumexp, roots_laguerre

from .benchmarks import DPIOUBenchmark


def hjb_log_reference(
    t: np.ndarray,
    x: np.ndarray,
    *,
    order: int = 128,
    return_gradient: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Hopf--Cole reference using scaled Gauss--Laguerre quadrature.

    ``exp(-u)=E[2/(1+||x+sqrt(2 tau) W||^2)]`` reduces exactly to a smooth
    one-dimensional integral.  The state enters only through ``||x||^2``;
    this fact is used by the reference calculator and is never exposed to an
    Actor.
    """
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(len(times), -1)
    tau = np.maximum(1.0 - times, 0.0)
    radius2 = np.square(states).sum(axis=1)
    nodes, weights = roots_laguerre(order)
    values = np.empty(len(times), dtype=np.float64)
    gradients = np.empty_like(states) if return_gradient else None
    terminal = tau <= 2e-14
    values[terminal] = np.log((1.0 + radius2[terminal]) / 2.0)
    if gradients is not None:
        gradients[terminal] = 2.0 * states[terminal] / (
            1.0 + radius2[terminal, None]
        )
    active = np.flatnonzero(~terminal)
    # Small batches avoid a large (#states, order) temporary for 20k-point PINN banks.
    for begin in range(0, len(active), 4096):
        index = active[begin : begin + 4096]
        local_tau = tau[index, None]
        local_r2 = radius2[index, None]
        scale = 1.0 + 2.0 * states.shape[1] * local_tau + local_r2
        s = nodes[None, :] / scale
        denominator = 1.0 + 4.0 * local_tau * s
        log_relative = (
            nodes[None, :]
            - s
            - 0.5 * states.shape[1] * np.log(denominator)
            - s * local_r2 / denominator
        )
        with np.errstate(divide="ignore"):
            log_terms = np.log(weights)[None, :] + log_relative
        log_sum = logsumexp(log_terms, axis=1)
        log_j = log_sum - np.log(scale[:, 0])
        values[index] = -math.log(2.0) - log_j
        if gradients is not None:
            normalized = np.exp(log_terms - log_sum[:, None])
            k_over_j = np.sum(normalized * s / denominator, axis=1)
            gradients[index] = 2.0 * states[index] * k_over_j[:, None]
    if gradients is None:
        return values
    return values, math.sqrt(2.0) * gradients


def hjb_log_quad_scalar(t: float, radius2: float, d: int = 100) -> float:
    tau = max(1.0 - t, 0.0)
    if tau == 0.0:
        return math.log((1.0 + radius2) / 2.0)

    def integrand(s: float) -> float:
        denominator = 1.0 + 4.0 * tau * s
        return math.exp(
            -s - 0.5 * d * math.log(denominator) - s * radius2 / denominator
        )

    integral, _ = quad(integrand, 0.0, np.inf, epsabs=2e-14, epsrel=2e-13, limit=500)
    return -math.log(2.0 * integral)


def dpi_exact_reference(
    benchmark: DPIOUBenchmark,
    t: np.ndarray,
    x: np.ndarray,
    *,
    return_gradient: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    times = np.asarray(t, dtype=np.float64).reshape(-1)
    states = np.asarray(x, dtype=np.float64).reshape(len(times), benchmark.d)
    tau = benchmark.T - times
    decay = np.exp(-benchmark.theta * tau)
    variance = benchmark.initial_variance * decay**2 + benchmark.alpha / (
        2.0 * benchmark.theta
    ) * (1.0 - decay**2)
    means = benchmark.mean.numpy()[None, :, :] * decay[:, None, None]
    difference = states[:, None, :] - means
    logits = (
        np.log(benchmark.pi.numpy())[None, :]
        - 0.5 * benchmark.d * np.log(2.0 * np.pi * variance[:, None])
        - 0.5 * np.square(difference).sum(axis=-1) / variance[:, None]
    )
    values = -logsumexp(logits, axis=1)
    if not return_gradient:
        return values
    responsibilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
    gradient = (
        responsibilities[:, :, None]
        * difference
        / variance[:, None, None]
    ).sum(axis=1)
    return values, benchmark.sigma * gradient


@dataclass(frozen=True)
class AllenRadialConfig:
    s_max: float = 300.0
    ds: float = 0.1
    dt: float = 0.0001
    dimension: int = 100
    horizon: float = 0.3


class AllenRadialReference:
    """Strang reaction / Crank--Nicolson squared-radius solver."""

    def __init__(self, config: AllenRadialConfig) -> None:
        self.config = config
        ratio = config.s_max / config.ds
        if abs(ratio - round(ratio)) > 1e-10:
            raise ValueError("s_max must be divisible by ds")
        self.s = np.linspace(0.0, config.s_max, int(round(ratio)) + 1)
        n = len(self.s)
        lower = np.zeros(n)
        diagonal = np.zeros(n)
        upper = np.zeros(n)
        # Generator of S=||X||^2 for dX=sqrt(2)dW:
        # L=4s*d_ss+2d*d_s.  A forward spatial difference for the drift
        # preserves the Markov-generator sign pattern at s=0 and near it.
        diffusion = 4.0 * self.s[:-1] / config.ds**2
        convection = 2.0 * config.dimension / config.ds
        diagonal[:-1] = -2.0 * diffusion - convection
        upper[:-1] = diffusion + convection
        lower[1:-1] = diffusion[1:]
        # Far boundary is fixed during the diffusion substep.
        diagonal[-1] = 0.0
        operator = diags(
            (lower[1:], diagonal, upper[:-1]), offsets=(-1, 0, 1), format="csc"
        )
        identity = diags(np.ones(n), format="csc")
        self._left = csc_matrix(identity - 0.5 * config.dt * operator)
        self._right = csc_matrix(identity + 0.5 * config.dt * operator)
        self._solve = factorized(self._left)

    @staticmethod
    def _reaction(value: np.ndarray, duration: float) -> np.ndarray:
        # Exact flow of y'=y-y^3 for nonnegative y.
        square = np.square(value)
        denominator = np.sqrt(square + (1.0 - square) * np.exp(-2.0 * duration))
        return value / denominator

    def _step(self, value: np.ndarray) -> np.ndarray:
        half = self._reaction(value, 0.5 * self.config.dt)
        rhs = self._right @ half
        rhs[-1] = half[-1]
        diffused = self._solve(rhs)
        return self._reaction(diffused, 0.5 * self.config.dt)

    def evaluate(self, t: np.ndarray, radius2: np.ndarray) -> np.ndarray:
        times = np.asarray(t, dtype=np.float64).reshape(-1)
        radii = np.asarray(radius2, dtype=np.float64).reshape(-1)
        if len(times) != len(radii):
            raise ValueError("time/state length mismatch")
        if np.min(times) < -1e-14 or np.max(times) > self.config.horizon + 1e-14:
            raise ValueError("time outside Allen-Cahn horizon")
        if np.min(radii) < 0.0 or np.max(radii) > self.config.s_max:
            raise ValueError("radius outside radial reference grid")
        tau = self.config.horizon - times
        order = np.argsort(tau)
        result = np.empty(len(times), dtype=np.float64)
        value = 0.5 / (1.0 + 0.2 * self.s)
        current_tau = 0.0
        query_offset = 0
        while query_offset < len(order) and tau[order[query_offset]] <= 1e-15:
            index = order[query_offset]
            result[index] = np.interp(radii[index], self.s, value)
            query_offset += 1
        while query_offset < len(order):
            previous = value
            value = self._step(value)
            next_tau = current_tau + self.config.dt
            while query_offset < len(order) and tau[order[query_offset]] <= next_tau + 1e-15:
                index = order[query_offset]
                fraction = np.clip(
                    (tau[index] - current_tau) / self.config.dt, 0.0, 1.0
                )
                spatial = (1.0 - fraction) * previous + fraction * value
                result[index] = np.interp(radii[index], self.s, spatial)
                query_offset += 1
            current_tau = next_tau
            if current_tau > self.config.horizon + self.config.dt:
                raise RuntimeError("radial integration overran horizon")
        return result

    def metadata(self) -> dict[str, Any]:
        return {
            **asdict(self.config),
            "coordinate": "s=||x||^2",
            "pde": "U_tau=4s U_ss+2d U_s+U-U^3",
            "reaction_integrator": "exact half steps",
            "diffusion_integrator": "Crank-Nicolson",
            "drift_spatial_difference": "forward Markov-generator stencil",
        }


def allen_reference(
    t: np.ndarray,
    x: np.ndarray,
    *,
    config: AllenRadialConfig | None = None,
) -> np.ndarray:
    states = np.asarray(x, dtype=np.float64)
    solver = AllenRadialReference(config or AllenRadialConfig())
    return solver.evaluate(t, np.square(states).sum(axis=1))


def allen_reference_extrapolated(
    t: np.ndarray,
    x: np.ndarray,
    *,
    s_max: float = 200.0,
    fine_ds: float = 0.0125,
) -> tuple[np.ndarray, dict[str, Any]]:
    """First-order Richardson limit with an observable discretization error.

    The monotone forward stencil is first order in ``ds``.  Time errors are
    made negligible by halving ``dt`` with ``ds``; the coarse/fine difference
    and the extrapolation correction are retained in the audit.
    """
    states = np.asarray(x, dtype=np.float64)
    radius2 = np.square(states).sum(axis=1)
    fine_config = AllenRadialConfig(
        s_max=s_max, ds=fine_ds, dt=fine_ds * 0.0005
    )
    coarse_config = AllenRadialConfig(
        s_max=s_max, ds=2.0 * fine_ds, dt=fine_ds * 0.001
    )
    coarse = AllenRadialReference(coarse_config).evaluate(t, radius2)
    fine = AllenRadialReference(fine_config).evaluate(t, radius2)
    extrapolated = 2.0 * fine - coarse
    return extrapolated, {
        "coarse": asdict(coarse_config),
        "fine": asdict(fine_config),
        "richardson_order": 1,
        "maximum_coarse_fine_difference": float(np.max(np.abs(fine - coarse))),
        "maximum_extrapolation_correction": float(
            np.max(np.abs(extrapolated - fine))
        ),
        "first_query_coarse": float(coarse[0]) if len(coarse) else None,
        "first_query_fine": float(fine[0]) if len(fine) else None,
        "first_query_extrapolated": float(extrapolated[0]) if len(extrapolated) else None,
    }


__all__ = [
    "AllenRadialConfig",
    "AllenRadialReference",
    "allen_reference",
    "allen_reference_extrapolated",
    "dpi_exact_reference",
    "hjb_log_quad_scalar",
    "hjb_log_reference",
]
