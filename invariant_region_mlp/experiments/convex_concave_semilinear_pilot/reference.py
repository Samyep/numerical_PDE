"""Independent periodic spectral references and exact tangent policies."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import torch
from torch import nn

from .drivers import ConcaveSemilinearDriver


@dataclass
class PeriodicReference:
    name: str
    times: np.ndarray
    space: np.ndarray
    values: np.ndarray
    derivatives: np.ndarray
    dt: float
    n: int

    def _interpolate(self, array: np.ndarray, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=np.float64).reshape(-1)
        s = np.mod(np.asarray(s, dtype=np.float64).reshape(-1), 2.0 * np.pi)
        ft = np.clip((t - self.times[0]) / (self.times[1] - self.times[0]), 0, len(self.times) - 1)
        fs = s / (2.0 * np.pi) * self.n
        it = np.floor(ft).astype(int)
        it = np.minimum(it, len(self.times) - 2)
        ix = np.floor(fs).astype(int) % self.n
        jt = ft - it
        jx = fs - np.floor(fs)
        ix1 = (ix + 1) % self.n
        lower = array[it, ix] * (1.0 - jx) + array[it, ix1] * jx
        upper = array[it + 1, ix] * (1.0 - jx) + array[it + 1, ix1] * jx
        return lower * (1.0 - jt) + upper * jt

    def value(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.values, t, s)

    def derivative(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.derivatives, t, s)


def _spectral_derivative(value: np.ndarray, order: int = 1) -> np.ndarray:
    n = value.shape[-1]
    wave = np.fft.fftfreq(n, d=1.0 / n)
    multiplier = (1j * wave) ** order
    return np.fft.ifft(np.fft.fft(value, axis=-1) * multiplier, axis=-1).real


def _dealias(value: np.ndarray) -> np.ndarray:
    n = value.shape[-1]
    wave = np.fft.fftfreq(n, d=1.0 / n)
    spectrum = np.fft.fft(value)
    spectrum[np.abs(wave) > n / 3] = 0.0
    return np.fft.ifft(spectrum).real


def solve_periodic_reference(
    driver: ConcaveSemilinearDriver,
    *,
    n: int = 512,
    dt: float = 0.000125,
    store_dt: float = 0.0005,
) -> PeriodicReference:
    """Second-order diffusion splitting in forward time tau=T-t."""
    steps = int(round(driver.T / dt))
    dt = driver.T / steps
    store_every = int(round(store_dt / dt))
    if steps % store_every:
        raise ValueError("store interval must divide the time grid")
    space = np.arange(n, dtype=np.float64) * (2.0 * np.pi / n)
    beta_state = space[:, None] * driver.beta[None, :]
    value = driver.terminal_numpy(beta_state)
    wave = np.fft.fftfreq(n, d=1.0 / n)
    diffusion_half = np.exp(-wave**2 * dt / 2.0)
    stored = [value.copy()]

    def nonlinear(current: np.ndarray) -> np.ndarray:
        derivative = _spectral_derivative(current)
        z_scalar = math.sqrt(2.0) * derivative
        return _dealias(driver.driver_numpy(current, z_scalar))

    for step in range(1, steps + 1):
        value = np.fft.ifft(np.fft.fft(value) * diffusion_half).real
        first = nonlinear(value)
        midpoint = value + 0.5 * dt * first
        value = value + dt * nonlinear(midpoint)
        value = np.fft.ifft(np.fft.fft(value) * diffusion_half).real
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"nonfinite {driver.name} reference at step {step}")
        if step % store_every == 0:
            stored.append(value.copy())
    tau_values = np.asarray(stored)
    tau_times = np.linspace(0.0, driver.T, len(tau_values))
    derivatives_tau = _spectral_derivative(tau_values)
    return PeriodicReference(
        name=driver.name,
        times=(driver.T - tau_times[::-1]).copy(),
        space=space,
        values=tau_values[::-1].copy(),
        derivatives=derivatives_tau[::-1].copy(),
        dt=dt,
        n=n,
    )


def validate_periodic_reference(
    driver: ConcaveSemilinearDriver,
) -> tuple[PeriodicReference, dict[str, Any]]:
    coarse = solve_periodic_reference(driver, n=256, dt=0.00025, store_dt=0.001)
    fine = solve_periodic_reference(driver, n=512, dt=0.000125, store_dt=0.0005)
    times = np.repeat(np.linspace(0.0, driver.T, 11), 513)
    space = np.tile(np.linspace(0.0, 2.0 * np.pi, 513, endpoint=False), 11)
    value_delta = fine.value(times, space) - coarse.value(times, space)
    derivative_delta = fine.derivative(times, space) - coarse.derivative(times, space)
    terminal = fine.value(np.full(513, driver.T), space[:513])
    terminal_state = space[:513, None] * driver.beta[None, :]
    terminal_truth = driver.terminal_numpy(terminal_state)

    # Residual in forward tau, using central differences on the stored fine grid.
    tau_values = fine.values[::-1]
    tau_dt = fine.times[1] - fine.times[0]
    time_derivative = (tau_values[2:] - tau_values[:-2]) / (2.0 * tau_dt)
    center = tau_values[1:-1]
    spatial_first = _spectral_derivative(center)
    spatial_second = _spectral_derivative(center, order=2)
    rhs = spatial_second + driver.driver_numpy(center, math.sqrt(2.0) * spatial_first)
    residual = time_derivative - rhs
    domain_min = float(fine.values.min())
    domain_max = float(fine.values.max())
    return fine, {
        "problem": driver.name,
        "coarse_n": coarse.n,
        "coarse_dt": coarse.dt,
        "fine_n": fine.n,
        "fine_dt": fine.dt,
        "value_refinement_rmse": float(np.sqrt(np.mean(value_delta**2))),
        "value_refinement_max_abs": float(np.max(np.abs(value_delta))),
        "gradient_refinement_rmse": float(np.sqrt(np.mean(derivative_delta**2))),
        "gradient_refinement_max_abs": float(np.max(np.abs(derivative_delta))),
        "terminal_max_abs": float(np.max(np.abs(terminal - terminal_truth))),
        "pde_residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "pde_residual_max_abs": float(np.max(np.abs(residual))),
        "working_value_min": domain_min,
        "working_value_max": domain_max,
        "inside_tangent_domain": bool(
            domain_min >= driver.r_low - 2e-6 and domain_max <= driver.r_high + 2e-6
        ),
    }


class ReferenceTangentActor(nn.Module):
    """Post-training exact tangent lookup; never used for supervision."""

    def __init__(self, driver: ConcaveSemilinearDriver, reference: PeriodicReference) -> None:
        super().__init__()
        self.driver = driver
        self.dt = float(reference.times[1] - reference.times[0])
        self.register_buffer("values", torch.as_tensor(reference.values, dtype=torch.float32))
        self.register_buffer("derivatives", torch.as_tensor(reference.derivatives, dtype=torch.float32))
        self.register_buffer("beta", torch.as_tensor(driver.beta, dtype=torch.float32))

    def _lookup(self, table: torch.Tensor, t: torch.Tensor, s: torch.Tensor) -> torch.Tensor:
        ft = torch.clamp(t / self.dt, 0.0, table.shape[0] - 1.0)
        fs = torch.remainder(s, 2.0 * math.pi) / (2.0 * math.pi) * table.shape[1]
        it = torch.floor(ft).long().clamp(max=table.shape[0] - 2)
        ix = torch.floor(fs).long() % table.shape[1]
        ix1 = (ix + 1) % table.shape[1]
        wt, wx = ft - it, fs - torch.floor(fs)
        lower = table[it, ix] * (1.0 - wx) + table[it, ix1] * wx
        upper = table[it + 1, ix] * (1.0 - wx) + table[it + 1, ix1] * wx
        return lower * (1.0 - wt) + upper * wt

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        del T
        beta = self.beta.to(dtype=x.dtype)
        s = x @ beta
        r = self._lookup(self.values, t, s)
        if not self.driver.has_gradient:
            return r, None
        derivative = self._lookup(self.derivatives, t, s)
        q = math.sqrt(2.0) * derivative[:, None] * beta[None, :]
        return r, q
