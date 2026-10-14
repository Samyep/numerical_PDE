"""Independent Riccati and one-dimensional finite-difference references."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import solve_banded

from .problems import BoundedControlProblem, LQGProblem, QuarticControlProblem


class RiccatiReference:
    def __init__(self, problem: LQGProblem) -> None:
        self.problem = problem

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            del t
            P = y[0]
            return np.array([
                2.0 * problem.gamma * P + P**2 / problem.R - problem.Q,
                -problem.d * P,
            ])

        self.solution = solve_ivp(
            rhs,
            (problem.T, 0.0),
            np.array([problem.G, 0.0]),
            rtol=1e-11,
            atol=1e-13,
            dense_output=True,
        )
        if not self.solution.success:
            raise RuntimeError(self.solution.message)

    def coefficients(self, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = self.solution.sol(np.asarray(t, dtype=np.float64))
        return values[0], values[1]

    def value(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=np.float64).reshape(-1)
        x = np.asarray(x, dtype=np.float64).reshape(len(t), self.problem.d)
        P, c = self.coefficients(t)
        return 0.5 * P * np.square(x).sum(axis=-1) + c

    def action(self, t: np.ndarray, x: np.ndarray) -> np.ndarray:
        P, _ = self.coefficients(np.asarray(t, dtype=np.float64).reshape(-1))
        return -(P / self.problem.R)[:, None] * np.asarray(x, dtype=np.float64)


@dataclass
class ScalarReference:
    problem_name: str
    times: np.ndarray
    space: np.ndarray
    values: np.ndarray
    gradients: np.ndarray
    actions: np.ndarray
    dx: float
    dt: float

    def _interpolate(self, array: np.ndarray, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        points = np.column_stack((np.asarray(t).reshape(-1), np.asarray(s).reshape(-1)))
        interpolator = RegularGridInterpolator(
            (self.times, self.space), array, bounds_error=False, fill_value=None
        )
        return np.asarray(interpolator(points), dtype=np.float64)

    def value(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.values, t, s)

    def gradient(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.gradients, t, s)

    def action(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.actions, t, s)


def solve_scalar_reference(
    problem: QuarticControlProblem | BoundedControlProblem,
    *,
    half_width: float = 10.0,
    dx: float = 0.025,
    dt: float = 0.0005,
) -> ScalarReference:
    """Semi-implicit diffusion and explicit Hamiltonian integration in tau=T-t."""
    nx = int(round(2.0 * half_width / dx)) + 1
    nt = int(round(problem.T / dt))
    dx = 2.0 * half_width / (nx - 1)
    dt = problem.T / nt
    space = np.linspace(-half_width, half_width, nx)
    tau = np.linspace(0.0, problem.T, nt + 1)
    values_tau = np.empty((nt + 1, nx), dtype=np.float64)
    values_tau[0] = problem.terminal_1d(space)
    ratio = dt / dx**2
    n_inner = nx - 2
    band = np.zeros((3, n_inner), dtype=np.float64)
    band[0, 1:] = -ratio
    band[1] = 1.0 + 2.0 * ratio
    band[2, :-1] = -ratio

    for n in range(nt):
        current = values_tau[n]
        gradient = np.empty_like(current)
        gradient[1:-1] = (current[2:] - current[:-2]) / (2.0 * dx)
        gradient[0] = (current[1] - current[0]) / dx
        gradient[-1] = (current[-1] - current[-2]) / dx
        action = problem.greedy_numpy(gradient)
        hamiltonian = problem.hamiltonian_numpy(action, gradient)
        left = current[0] + dt * hamiltonian[0]
        right = current[-1] + dt * hamiltonian[-1]
        rhs = current[1:-1] + dt * hamiltonian[1:-1]
        rhs[0] += ratio * left
        rhs[-1] += ratio * right
        following = np.empty_like(current)
        following[0], following[-1] = left, right
        following[1:-1] = solve_banded((1, 1), band, rhs, check_finite=False)
        if not np.all(np.isfinite(following)):
            raise FloatingPointError(f"nonfinite {problem.name} reference at step {n}")
        values_tau[n + 1] = following

    gradients_tau = np.gradient(values_tau, dx, axis=1, edge_order=2)
    actions_tau = problem.greedy_numpy(gradients_tau)
    times = problem.T - tau[::-1]
    return ScalarReference(
        problem_name=problem.name,
        times=times,
        space=space,
        values=values_tau[::-1].copy(),
        gradients=gradients_tau[::-1].copy(),
        actions=actions_tau[::-1].copy(),
        dx=dx,
        dt=dt,
    )


def validate_scalar_reference(
    problem: QuarticControlProblem | BoundedControlProblem,
) -> tuple[ScalarReference, dict[str, Any]]:
    coarse = solve_scalar_reference(problem, dx=0.05, dt=0.001)
    fine = solve_scalar_reference(problem, dx=0.025, dt=0.0005)
    t = np.repeat(np.linspace(0.0, problem.T, 5), 321)
    s = np.tile(np.linspace(-4.0, 4.0, 321), 5)
    value_delta = fine.value(t, s) - coarse.value(t, s)
    action_delta = fine.action(t, s) - coarse.action(t, s)
    terminal_delta = fine.value(
        np.full(321, problem.T), np.linspace(-4.0, 4.0, 321)
    ) - problem.terminal_1d(np.linspace(-4.0, 4.0, 321))
    metrics = {
        "problem": problem.name,
        "coarse_dx": coarse.dx,
        "coarse_dt": coarse.dt,
        "fine_dx": fine.dx,
        "fine_dt": fine.dt,
        "audit_domain": [-4.0, 4.0],
        "value_refinement_rmse": float(np.sqrt(np.mean(value_delta**2))),
        "value_refinement_max_abs": float(np.max(np.abs(value_delta))),
        "action_refinement_rmse": float(np.sqrt(np.mean(action_delta**2))),
        "action_refinement_max_abs": float(np.max(np.abs(action_delta))),
        "terminal_max_abs": float(np.max(np.abs(terminal_delta))),
        "boundary_half_width": 10.0,
    }
    return fine, metrics

