"""Independent two-dimensional periodic spectral ground truth."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
import torch
from torch import nn

from .drivers import ConcaveSemilinearDriver


def _waves(n: int) -> np.ndarray:
    return np.fft.fftfreq(n, d=1.0 / n)


def _spectral_fields(value: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = value.shape[-1]
    wave = _waves(n)
    spectrum = np.fft.fftn(value, axes=(-2, -1))
    d1 = np.fft.ifftn(spectrum * (1j * wave)[:, None], axes=(-2, -1)).real
    d2 = np.fft.ifftn(spectrum * (1j * wave)[None, :], axes=(-2, -1)).real
    laplace = np.fft.ifftn(
        spectrum * (-(wave[:, None] ** 2 + wave[None, :] ** 2)), axes=(-2, -1)
    ).real
    return d1, d2, laplace


def _dealias(value: np.ndarray) -> np.ndarray:
    n = value.shape[-1]
    wave = _waves(n)
    keep = (np.abs(wave)[:, None] <= n / 3) & (np.abs(wave)[None, :] <= n / 3)
    spectrum = np.fft.fftn(value, axes=(-2, -1))
    spectrum *= keep
    return np.fft.ifftn(spectrum, axes=(-2, -1)).real


@dataclass
class PeriodicReference2D:
    name: str
    times: np.ndarray
    space: np.ndarray
    values: np.ndarray
    derivative_1: np.ndarray
    derivative_2: np.ndarray
    dt: float
    store_dt: float
    n: int

    def _interpolate(self, table: np.ndarray, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        times = np.asarray(t, dtype=np.float64).reshape(-1)
        coordinates = np.mod(np.asarray(s, dtype=np.float64).reshape(-1, 2), 2.0 * np.pi)
        ft = np.clip(times / self.store_dt, 0.0, len(self.times) - 1.0)
        it = np.floor(ft).astype(int)
        it = np.minimum(it, len(self.times) - 2)
        wt = ft - it
        fs = coordinates / (2.0 * np.pi) * self.n
        ix = np.floor(fs).astype(int) % self.n
        wx = fs - np.floor(fs)
        ix1 = (ix + 1) % self.n

        def spatial(time_index: np.ndarray) -> np.ndarray:
            a00 = table[time_index, ix[:, 0], ix[:, 1]]
            a10 = table[time_index, ix1[:, 0], ix[:, 1]]
            a01 = table[time_index, ix[:, 0], ix1[:, 1]]
            a11 = table[time_index, ix1[:, 0], ix1[:, 1]]
            low = a00 * (1.0 - wx[:, 0]) + a10 * wx[:, 0]
            high = a01 * (1.0 - wx[:, 0]) + a11 * wx[:, 0]
            return low * (1.0 - wx[:, 1]) + high * wx[:, 1]

        lower = spatial(it)
        upper = spatial(it + 1)
        return lower * (1.0 - wt) + upper * wt

    def value(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return self._interpolate(self.values, t, s)

    def gradient(self, t: np.ndarray, s: np.ndarray) -> np.ndarray:
        return np.stack(
            (
                self._interpolate(self.derivative_1, t, s),
                self._interpolate(self.derivative_2, t, s),
            ),
            axis=-1,
        )


def solve_periodic_reference_2d(
    driver: ConcaveSemilinearDriver,
    *,
    n: int = 128,
    dt: float = 0.00025,
    store_dt: float = 0.002,
) -> PeriodicReference2D:
    """Strang diffusion splitting with midpoint nonlinear integration in tau=T-t."""
    if driver.k != 2:
        raise ValueError("the numerical reference is deliberately restricted to k=2")
    steps = int(round(driver.T / dt))
    dt = driver.T / steps
    store_every = int(round(store_dt / dt))
    if store_every < 1 or steps % store_every:
        raise ValueError("store interval must divide time grid")
    store_dt = store_every * dt
    space = np.arange(n, dtype=np.float64) * (2.0 * np.pi / n)
    s1, s2 = np.meshgrid(space, space, indexing="ij")
    coordinates = np.stack((s1, s2), axis=-1)
    value = np.asarray(driver.terminal_from_coordinates_numpy(coordinates), dtype=np.float64)
    wave = _waves(n)
    diffusion_half = np.exp(
        -(wave[:, None] ** 2 + wave[None, :] ** 2) * dt / 2.0
    )
    stored = [value.copy()]

    def diffuse(current: np.ndarray) -> np.ndarray:
        return np.fft.ifftn(
            np.fft.fftn(current, axes=(-2, -1)) * diffusion_half, axes=(-2, -1)
        ).real

    def nonlinear(current: np.ndarray) -> np.ndarray:
        d1, d2, _ = _spectral_fields(current)
        z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
        return _dealias(driver.driver_reduced_numpy(current, z))

    for step in range(1, steps + 1):
        value = diffuse(value)
        first = nonlinear(value)
        value = value + dt * nonlinear(value + 0.5 * dt * first)
        value = diffuse(value)
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"nonfinite {driver.name} reference at step {step}")
        if step % store_every == 0:
            stored.append(value.copy())
    tau_values = np.asarray(stored)
    derivative_1, derivative_2, _ = _spectral_fields(tau_values)
    tau_times = np.arange(len(tau_values), dtype=np.float64) * store_dt
    return PeriodicReference2D(
        name=driver.name,
        times=(driver.T - tau_times[::-1]).copy(),
        space=space,
        values=tau_values[::-1].copy(),
        derivative_1=derivative_1[::-1].copy(),
        derivative_2=derivative_2[::-1].copy(),
        dt=dt,
        store_dt=store_dt,
        n=n,
    )


def validate_periodic_reference_2d(
    driver: ConcaveSemilinearDriver,
) -> tuple[PeriodicReference2D, dict[str, Any]]:
    coarse = solve_periodic_reference_2d(driver, n=64, dt=0.0005, store_dt=0.004)
    fine = solve_periodic_reference_2d(driver, n=128, dt=0.00025, store_dt=0.002)
    rng = np.random.default_rng(77831)
    query_t = rng.uniform(0.0, driver.T, 32768)
    query_s = rng.uniform(0.0, 2.0 * np.pi, (32768, 2))
    value_delta = fine.value(query_t, query_s) - coarse.value(query_t, query_s)
    gradient_delta = fine.gradient(query_t, query_s) - coarse.gradient(query_t, query_s)
    terminal = fine.value(np.full(len(query_t), driver.T), query_s)
    terminal_truth = driver.terminal_from_coordinates_numpy(query_s)

    # A centered residual on a subset of stored fine-time slices avoids using
    # the splitting update itself as the validation equation.
    tau_values = fine.values[::-1]
    indices = np.arange(1, len(tau_values) - 1, 8)
    center = tau_values[indices]
    time_derivative = (tau_values[indices + 1] - tau_values[indices - 1]) / (
        2.0 * fine.store_dt
    )
    d1, d2, laplace = _spectral_fields(center)
    z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
    residual = time_derivative - laplace - driver.driver_reduced_numpy(center, z)
    domain_min = float(fine.values.min())
    domain_max = float(fine.values.max())
    return fine, {
        "problem": driver.family,
        "reference_dimension": 2,
        "solver": "periodic Fourier Strang splitting with explicit midpoint nonlinear step",
        "coarse_n": coarse.n,
        "coarse_dt": coarse.dt,
        "coarse_store_dt": coarse.store_dt,
        "fine_n": fine.n,
        "fine_dt": fine.dt,
        "fine_store_dt": fine.store_dt,
        "query_points": len(query_t),
        "value_refinement_rmse": float(np.sqrt(np.mean(value_delta**2))),
        "value_refinement_max_abs": float(np.max(np.abs(value_delta))),
        "gradient_refinement_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient_delta).sum(axis=1)))
        ),
        "gradient_refinement_max_norm": float(
            np.max(np.linalg.norm(gradient_delta, axis=1))
        ),
        "terminal_interpolation_rmse": float(np.sqrt(np.mean((terminal - terminal_truth) ** 2))),
        "terminal_interpolation_max_abs": float(np.max(np.abs(terminal - terminal_truth))),
        "pde_residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "pde_residual_max_abs": float(np.max(np.abs(residual))),
        "working_value_min": domain_min,
        "working_value_max": domain_max,
        "inside_tangent_domain": bool(
            domain_min >= driver.r_low - 2e-6 and domain_max <= driver.r_high + 2e-6
        ),
    }


class ReferenceTangentActor2D(nn.Module):
    """Post-training lookup of the independent 2D reference tangent."""

    def __init__(self, driver: ConcaveSemilinearDriver, reference: PeriodicReference2D) -> None:
        super().__init__()
        self.driver = driver
        self.store_dt = reference.store_dt
        self.register_buffer("values", torch.as_tensor(reference.values, dtype=torch.float32))
        self.register_buffer("d1", torch.as_tensor(reference.derivative_1, dtype=torch.float32))
        self.register_buffer("d2", torch.as_tensor(reference.derivative_2, dtype=torch.float32))
        self.register_buffer("basis", torch.as_tensor(driver.B[:, :2], dtype=torch.float32))

    def _lookup(self, table: torch.Tensor, t: torch.Tensor, s: torch.Tensor) -> torch.Tensor:
        ft = torch.clamp(t / self.store_dt, 0.0, table.shape[0] - 1.0)
        it = torch.floor(ft).long().clamp(max=table.shape[0] - 2)
        wt = ft - it
        fs = torch.remainder(s, 2.0 * math.pi) / (2.0 * math.pi) * table.shape[1]
        ix = torch.floor(fs).long() % table.shape[1]
        ix1 = (ix + 1) % table.shape[1]
        wx = fs - torch.floor(fs)

        def spatial(time_index: torch.Tensor) -> torch.Tensor:
            a00 = table[time_index, ix[:, 0], ix[:, 1]]
            a10 = table[time_index, ix1[:, 0], ix[:, 1]]
            a01 = table[time_index, ix[:, 0], ix1[:, 1]]
            a11 = table[time_index, ix1[:, 0], ix1[:, 1]]
            low = a00 * (1.0 - wx[:, 0]) + a10 * wx[:, 0]
            high = a01 * (1.0 - wx[:, 0]) + a11 * wx[:, 0]
            return low * (1.0 - wx[:, 1]) + high * wx[:, 1]

        return spatial(it) * (1.0 - wt) + spatial(it + 1) * wt

    def predict(
        self, t: torch.Tensor, x: torch.Tensor, T: float
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        del T
        basis = self.basis.to(dtype=x.dtype)
        coordinates = x @ basis
        r = self._lookup(self.values, t, coordinates)
        if not self.driver.has_gradient:
            return r, None
        reduced_gradient = torch.stack(
            (self._lookup(self.d1, t, coordinates), self._lookup(self.d2, t, coordinates)),
            dim=-1,
        )
        q = math.sqrt(2.0) * (reduced_gradient @ basis.T)
        return r, q
