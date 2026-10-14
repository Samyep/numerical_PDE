"""The frozen B1 and HJB-Log equations and independent exact references."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import torch
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from main_b1_sixmethods.equation import B1Equation


SQRT2 = math.sqrt(2.0)


@dataclass(frozen=True)
class HJBLogEquation:
    """HJB-Log with a deterministic one-dimensional Hopf--Cole reference."""

    d: int = 400
    T: float = 1.0
    sigma: float = SQRT2
    name: str = "hjb_log"
    family: str = "quadratic_hjb"
    quadrature_nodes: int = 256
    quadrature_upper: float = 80.0
    k: int = 0

    def __post_init__(self) -> None:
        if self.d <= 0 or self.T <= 0:
            raise ValueError("d and T must be positive")
        nodes, weights = leggauss(self.quadrature_nodes)
        y = 0.5 * self.quadrature_upper * (nodes + 1.0)
        w = 0.5 * self.quadrature_upper * weights
        object.__setattr__(self, "_quad_y", y.astype(np.float64))
        object.__setattr__(self, "_quad_logw", np.log(w.astype(np.float64)))

    def terminal(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        return np.log1p(np.sum(state * state, axis=-1)) - math.log(2.0)

    def terminal_grad(self, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        return 2.0 * state / (1.0 + np.sum(state * state, axis=-1))[..., None]

    def terminal_z(self, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.terminal_grad(x)

    def _quadrature(self, t: np.ndarray | float, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        state = np.asarray(x, dtype=np.float64)
        if state.shape[-1] != self.d:
            raise ValueError(f"expected final dimension {self.d}, got {state.shape}")
        time = np.broadcast_to(np.asarray(t, dtype=np.float64), state.shape[:-1])
        flat_t = time.reshape(-1)
        flat_x = state.reshape(-1, self.d)
        r2 = np.sum(flat_x * flat_x, axis=1)
        values = np.empty(len(flat_t), dtype=np.float64)
        ratios = np.empty(len(flat_t), dtype=np.float64)
        terminal = np.isclose(flat_t, self.T, rtol=0.0, atol=2e-15)
        if np.any(terminal):
            values[terminal] = np.log1p(r2[terminal]) - math.log(2.0)
            ratios[terminal] = 1.0 / (1.0 + r2[terminal])
        indices = np.flatnonzero(~terminal)
        y = self._quad_y[None, :]
        logw = self._quad_logw[None, :]
        for start in range(0, len(indices), 1024):
            chosen = indices[start : start + 1024]
            tau = np.maximum(self.T - flat_t[chosen], 0.0)[:, None]
            denominator = 1.0 + 4.0 * tau * y / self.d
            exponent = (
                -y / self.d
                - 0.5 * self.d * np.log1p(4.0 * tau * y / self.d)
                - (y / self.d) * r2[chosen, None] / denominator
            )
            log_terms = logw + exponent
            maximum = np.max(log_terms, axis=1, keepdims=True)
            scaled = np.exp(log_terms - maximum)
            normalizer = np.sum(scaled, axis=1)
            log_integral = (
                math.log(2.0 / self.d)
                + maximum[:, 0]
                + np.log(normalizer)
            )
            values[chosen] = -log_integral
            ratios[chosen] = np.sum(
                scaled * (y / self.d) / denominator, axis=1
            ) / normalizer
        return values.reshape(state.shape[:-1]), ratios.reshape(state.shape[:-1])

    def exact_u(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self._quadrature(t, x)[0]

    def exact_grad(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        state = np.asarray(x, dtype=np.float64)
        ratio = self._quadrature(t, state)[1]
        return 2.0 * state * ratio[..., None]

    def exact_z(self, t: np.ndarray | float, x: np.ndarray) -> np.ndarray:
        return self.sigma * self.exact_grad(t, x)

    def adaptive_reference(self, t: float, x: np.ndarray) -> tuple[float, np.ndarray, float]:
        """Slow independent scipy.quad reference for audit points."""

        state = np.asarray(x, dtype=np.float64).reshape(self.d)
        tau = max(self.T - float(t), 0.0)
        r2 = float(np.dot(state, state))
        if tau == 0.0:
            value = math.log1p(r2) - math.log(2.0)
            grad = 2.0 * state / (1.0 + r2)
            return value, self.sigma * grad, 0.0

        # Integrate in y=d*s.  At high dimension the original s-integrand is
        # concentrated in an O(1/d) boundary layer that generic quadrature can
        # otherwise miss.
        def integrand(y: float) -> float:
            s = y / self.d
            denominator = 1.0 + 4.0 * tau * s
            return math.exp(
                -s - 0.5 * self.d * math.log(denominator) - s * r2 / denominator
            )

        integral_y, error = quad(
            integrand, 0.0, np.inf, epsabs=2e-13, epsrel=2e-13, limit=300
        )

        def moment_integrand(y: float) -> float:
            return integrand(y) * y / (1.0 + 4.0 * tau * y / self.d)

        moment, moment_error = quad(
            moment_integrand, 0.0, np.inf, epsabs=2e-13, epsrel=2e-13, limit=300
        )
        value = -math.log(2.0 * integral_y / self.d)
        grad = 2.0 * state * moment / (self.d * integral_y)
        return value, self.sigma * grad, max(error, moment_error)

    def analytic_residual(self, t: np.ndarray, x: np.ndarray, h: float = 2e-4) -> np.ndarray:
        """Finite-difference audit of the exact reference at a few points."""

        times = np.asarray(t, dtype=np.float64).reshape(-1)
        states = np.asarray(x, dtype=np.float64).reshape(-1, self.d)
        output = np.empty(len(times), dtype=np.float64)
        for row, (time, state) in enumerate(zip(times, states, strict=True)):
            u_plus = self.exact_u(time + h, state[None, :])[0]
            u_minus = self.exact_u(time - h, state[None, :])[0]
            u_t = (u_plus - u_minus) / (2.0 * h)
            center = self.exact_u(time, state[None, :])[0]
            laplacian = 0.0
            for coordinate in range(self.d):
                delta = np.zeros(self.d, dtype=np.float64)
                delta[coordinate] = h
                laplacian += (
                    self.exact_u(time, (state + delta)[None, :])[0]
                    - 2.0 * center
                    + self.exact_u(time, (state - delta)[None, :])[0]
                ) / h**2
            gradient = self.exact_grad(time, state[None, :])[0]
            output[row] = u_t + laplacian - float(np.dot(gradient, gradient))
        return output

    def provenance(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "d": self.d,
            "T": self.T,
            "sigma": self.sigma,
            "terminal": "log(1+||x||^2)-log(2)",
            "quadrature": {
                "change_of_variable": "y=d*s",
                "rule": "Gauss-Legendre",
                "nodes": self.quadrature_nodes,
                "upper_y": self.quadrature_upper,
                "dtype": "float64",
            },
        }


Equation = B1Equation | HJBLogEquation


def make_equation(name: str) -> Equation:
    if name == "b1":
        return B1Equation(d=400, k=40)
    if name in {"hjb_log", "log"}:
        return HJBLogEquation(d=400)
    raise ValueError(f"unknown equation {name}")


def equation_key(equation: Equation) -> str:
    return "b1_d400_k40" if isinstance(equation, B1Equation) else "hjb_log_d400"


def torch_terminal(equation: Equation, x: torch.Tensor) -> torch.Tensor:
    if isinstance(equation, B1Equation):
        directions = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
        log_weights = torch.as_tensor(equation.logc, dtype=x.dtype, device=x.device)
        return -torch.logsumexp(x @ directions.T + log_weights, dim=-1)
    return torch.log1p(torch.sum(x.square(), dim=-1)) - math.log(2.0)


def torch_terminal_z(equation: Equation, x: torch.Tensor) -> torch.Tensor:
    if isinstance(equation, B1Equation):
        directions = torch.as_tensor(equation.A, dtype=x.dtype, device=x.device)
        log_weights = torch.as_tensor(equation.logc, dtype=x.dtype, device=x.device)
        probabilities = torch.softmax(x @ directions.T + log_weights, dim=-1)
        return -equation.sigma * (probabilities @ directions)
    denominator = 1.0 + torch.sum(x.square(), dim=-1, keepdim=True)
    return equation.sigma * 2.0 * x / denominator
