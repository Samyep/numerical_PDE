"""Independent periodic 2D references and nonconvexity audits."""

from __future__ import annotations

import math
from pathlib import Path
import time
from typing import Any, Literal

import numpy as np

from invariant_region_mlp.experiments.convex_concave_semilinear_extended.reference2d import (
    PeriodicReference2D,
)

from .drivers import SineQuadraticDriver


def _waves(n: int) -> np.ndarray:
    return np.fft.fftfreq(n, d=1.0 / n)


def _spectral_fields(value: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = value.shape[-1]
    wave = _waves(n)
    spectrum = np.fft.fftn(value, axes=(-2, -1))
    d1 = np.fft.ifftn(spectrum * (1j * wave)[:, None], axes=(-2, -1)).real
    d2 = np.fft.ifftn(spectrum * (1j * wave)[None, :], axes=(-2, -1)).real
    laplace = np.fft.ifftn(
        spectrum * (-(wave[:, None] ** 2 + wave[None, :] ** 2)),
        axes=(-2, -1),
    ).real
    return d1, d2, laplace


def _dealias(value: np.ndarray) -> np.ndarray:
    n = value.shape[-1]
    wave = _waves(n)
    keep = (np.abs(wave)[:, None] <= n / 3) & (
        np.abs(wave)[None, :] <= n / 3
    )
    spectrum = np.fft.fftn(value, axes=(-2, -1))
    spectrum *= keep
    return np.fft.ifftn(spectrum, axes=(-2, -1)).real


def _reference_from_tau(
    driver: SineQuadraticDriver,
    stored: list[np.ndarray],
    *,
    n: int,
    dt: float,
    store_dt: float,
) -> PeriodicReference2D:
    tau_values = np.asarray(stored, dtype=np.float64)
    d1, d2, _ = _spectral_fields(tau_values)
    tau_times = np.arange(len(tau_values), dtype=np.float64) * store_dt
    return PeriodicReference2D(
        name=driver.name,
        times=(driver.T - tau_times[::-1]).copy(),
        space=np.arange(n, dtype=np.float64) * (2.0 * np.pi / n),
        values=tau_values[::-1].copy(),
        derivative_1=d1[::-1].copy(),
        derivative_2=d2[::-1].copy(),
        dt=dt,
        store_dt=store_dt,
        n=n,
    )


def solve_strang_reference(
    driver: SineQuadraticDriver,
    *,
    n: int,
    dt: float,
    store_dt: float = 0.002,
    nonlinear_scheme: Literal["midpoint", "heun"] = "midpoint",
    dealias: bool = True,
) -> PeriodicReference2D:
    """Fourier Strang splitting with an explicit second-order nonlinear step."""
    steps = int(round(driver.T / dt))
    dt = driver.T / steps
    store_every = int(round(store_dt / dt))
    if store_every < 1 or steps % store_every:
        raise ValueError("store_dt must divide the reference time grid")
    store_dt = store_every * dt
    space = np.arange(n, dtype=np.float64) * (2.0 * np.pi / n)
    s1, s2 = np.meshgrid(space, space, indexing="ij")
    coordinates = np.stack((s1, s2), axis=-1)
    value = np.asarray(
        driver.terminal_from_coordinates_numpy(coordinates), dtype=np.float64
    )
    wave = _waves(n)
    diffusion_half = np.exp(
        -(wave[:, None] ** 2 + wave[None, :] ** 2) * dt / 2.0
    )
    stored = [value.copy()]

    def diffuse(current: np.ndarray) -> np.ndarray:
        return np.fft.ifftn(
            np.fft.fftn(current, axes=(-2, -1)) * diffusion_half,
            axes=(-2, -1),
        ).real

    def nonlinear(current: np.ndarray) -> np.ndarray:
        d1, d2, _ = _spectral_fields(current)
        z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
        result = driver.driver_reduced_numpy(current, z)
        return _dealias(result) if dealias else result

    for step in range(1, steps + 1):
        value = diffuse(value)
        first = nonlinear(value)
        if nonlinear_scheme == "midpoint":
            value = value + dt * nonlinear(value + 0.5 * dt * first)
        elif nonlinear_scheme == "heun":
            value = value + 0.5 * dt * (first + nonlinear(value + dt * first))
        else:
            raise KeyError(nonlinear_scheme)
        value = diffuse(value)
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"nonfinite reference at time step {step}")
        if step % store_every == 0:
            stored.append(value.copy())
    return _reference_from_tau(
        driver, stored, n=n, dt=dt, store_dt=store_dt
    )


def solve_mol_rk4_reference(
    driver: SineQuadraticDriver,
    *,
    n: int = 48,
    dt: float = 0.00005,
    store_dt: float = 0.002,
) -> PeriodicReference2D:
    """Independent method-of-lines Fourier discretization with classical RK4."""
    steps = int(round(driver.T / dt))
    dt = driver.T / steps
    store_every = int(round(store_dt / dt))
    if store_every < 1 or steps % store_every:
        raise ValueError("store_dt must divide the RK4 time grid")
    store_dt = store_every * dt
    space = np.arange(n, dtype=np.float64) * (2.0 * np.pi / n)
    s1, s2 = np.meshgrid(space, space, indexing="ij")
    coordinates = np.stack((s1, s2), axis=-1)
    value = np.asarray(
        driver.terminal_from_coordinates_numpy(coordinates), dtype=np.float64
    )
    stored = [value.copy()]

    def rhs(current: np.ndarray) -> np.ndarray:
        d1, d2, laplace = _spectral_fields(current)
        z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
        return laplace + _dealias(driver.driver_reduced_numpy(current, z))

    for step in range(1, steps + 1):
        k1 = rhs(value)
        k2 = rhs(value + 0.5 * dt * k1)
        k3 = rhs(value + 0.5 * dt * k2)
        k4 = rhs(value + dt * k3)
        value = value + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"nonfinite RK4 reference at time step {step}")
        if step % store_every == 0:
            stored.append(value.copy())
    return _reference_from_tau(
        driver, stored, n=n, dt=dt, store_dt=store_dt
    )


def _comparison(
    candidate: PeriodicReference2D,
    reference: PeriodicReference2D,
    t: np.ndarray,
    coordinates: np.ndarray,
) -> dict[str, float]:
    value_delta = candidate.value(t, coordinates) - reference.value(t, coordinates)
    gradient_delta = candidate.gradient(t, coordinates) - reference.gradient(
        t, coordinates
    )
    return {
        "value_rmse": float(np.sqrt(np.mean(np.square(value_delta)))),
        "value_max_abs": float(np.max(np.abs(value_delta))),
        "gradient_vector_rmse": float(
            np.sqrt(np.mean(np.square(gradient_delta).sum(axis=1)))
        ),
        "gradient_max_norm": float(
            np.max(np.linalg.norm(gradient_delta, axis=1))
        ),
    }


def curvature_coverage(
    driver: SineQuadraticDriver,
    reference: PeriodicReference2D,
    *,
    samples: int = 32768,
    seed: int = 20261017,
) -> dict[str, Any]:
    """Measure curvature signs on the interior evaluation distribution."""
    rng = np.random.default_rng(seed + driver.d)
    t = rng.uniform(0.0, 0.8 * driver.T, samples)
    # B^T x is N(0,I_2) for x~N(0,I_d); constructing only these coordinates is
    # permitted for post-hoc reference auditing, never Actor training.
    coordinates = rng.standard_normal((samples, 2))
    value = reference.value(t, coordinates)
    f_uu = (
        -4.0 * math.pi**2 * driver.lam * np.sin(2.0 * math.pi * value)
        - driver.eta**2
    )
    eigenvalues = driver.reduced_joint_hessian_eigenvalues(value)
    positive = f_uu > 0.0
    negative = f_uu < 0.0
    strong_threshold = 0.25
    return {
        "distribution": "t~Uniform(0,0.8T), x~N(0,I_d), reference only",
        "samples": samples,
        "time_max": float(0.8 * driver.T),
        "value_min": float(np.min(value)),
        "value_max": float(np.max(value)),
        "value_mean": float(np.mean(value)),
        "value_std": float(np.std(value)),
        "value_quantiles": np.quantile(value, [0.01, 0.1, 0.5, 0.9, 0.99]).tolist(),
        "f_uu_positive_fraction": float(np.mean(positive)),
        "f_uu_negative_fraction": float(np.mean(negative)),
        "strong_positive_fraction": float(np.mean(f_uu > strong_threshold)),
        "strong_negative_fraction": float(np.mean(f_uu < -strong_threshold)),
        "f_uu_quantiles": np.quantile(
            f_uu, [0.01, 0.1, 0.5, 0.9, 0.99]
        ).tolist(),
        "joint_hessian_indefinite_fraction": float(
            np.mean((eigenvalues[:, 0] < 0.0) & (eigenvalues[:, -1] > 0.0))
        ),
        "joint_hessian_negative_definite_fraction": float(
            np.mean(eigenvalues[:, -1] < 0.0)
        ),
        "meaningful_sign_change_gate": bool(
            np.mean(positive) >= 0.10
            and np.mean(negative) >= 0.10
            and np.mean(f_uu > strong_threshold) >= 0.05
            and np.mean(f_uu < -strong_threshold) >= 0.05
        ),
    }


def validate_reference(
    driver: SineQuadraticDriver,
    *,
    query_points: int = 16384,
) -> tuple[PeriodicReference2D, dict[str, Any]]:
    """Four-level convergence plus scheme, residual, and aliasing checks."""
    started = time.perf_counter()
    coarse = solve_strang_reference(driver, n=48, dt=0.001)
    medium = solve_strang_reference(driver, n=64, dt=0.0005)
    fine = solve_strang_reference(driver, n=96, dt=0.00025)
    production = solve_strang_reference(driver, n=128, dt=0.000125)
    no_dealias = solve_strang_reference(
        driver, n=96, dt=0.00025, dealias=False
    )
    heun = solve_strang_reference(
        driver, n=64, dt=0.00025, nonlinear_scheme="heun"
    )
    rk4 = solve_mol_rk4_reference(driver, n=48, dt=0.00005)
    rng = np.random.default_rng(48719 + driver.d + int(1000 * driver.lam))
    t = rng.uniform(0.0, driver.T, query_points)
    coordinates = rng.uniform(0.0, 2.0 * np.pi, (query_points, 2))

    tau_values = production.values[::-1]
    indices = np.arange(1, len(tau_values) - 1, 4)
    center = tau_values[indices]
    time_derivative = (tau_values[indices + 1] - tau_values[indices - 1]) / (
        2.0 * production.store_dt
    )
    d1, d2, laplace = _spectral_fields(center)
    z = math.sqrt(2.0) * np.stack((d1, d2), axis=-1)
    residual = time_derivative - laplace - driver.driver_reduced_numpy(center, z)
    terminal_value = production.value(np.full(query_points, driver.T), coordinates)
    terminal_gradient = production.gradient(
        np.full(query_points, driver.T), coordinates
    )
    exact_terminal = driver.terminal_from_coordinates_numpy(coordinates)
    exact_terminal_gradient = driver.terminal_coordinate_gradient_numpy(coordinates)

    audit: dict[str, Any] = {
        "driver": driver.metadata(),
        "solver": "periodic Fourier Strang splitting + explicit midpoint",
        "production": {
            "n": production.n,
            "dt": production.dt,
            "store_dt": production.store_dt,
        },
        "levels": [
            {"name": "coarse", "n": coarse.n, "dt": coarse.dt},
            {"name": "medium", "n": medium.n, "dt": medium.dt},
            {"name": "fine", "n": fine.n, "dt": fine.dt},
            {"name": "production", "n": production.n, "dt": production.dt},
        ],
        "coarse_vs_production": _comparison(coarse, production, t, coordinates),
        "medium_vs_production": _comparison(medium, production, t, coordinates),
        "fine_vs_production": _comparison(fine, production, t, coordinates),
        "aliasing_n96_dealiased_vs_raw": _comparison(
            no_dealias, fine, t, coordinates
        ),
        "heun_n64_vs_midpoint_production": _comparison(
            heun, production, t, coordinates
        ),
        "independent_mol_rk4_n48_vs_split_n48": _comparison(
            rk4, coarse, t, coordinates
        ),
        "terminal_value_rmse": float(
            np.sqrt(np.mean(np.square(terminal_value - exact_terminal)))
        ),
        "terminal_gradient_vector_rmse": float(
            np.sqrt(
                np.mean(
                    np.square(terminal_gradient - exact_terminal_gradient).sum(axis=1)
                )
            )
        ),
        "pde_residual_rmse": float(np.sqrt(np.mean(np.square(residual)))),
        "pde_residual_max_abs": float(np.max(np.abs(residual))),
        "solution_grid_min": float(np.min(production.values)),
        "solution_grid_max": float(np.max(production.values)),
        "inside_invariant_interval": bool(
            np.min(production.values) >= -2e-8
            and np.max(production.values) <= 1.0 + 2e-8
        ),
        "reference_error_proxy_value_rmse": _comparison(
            fine, production, t, coordinates
        )["value_rmse"],
        "query_points": query_points,
        "wall_seconds": time.perf_counter() - started,
    }
    audit["curvature_coverage"] = curvature_coverage(driver, production)
    return production, audit


def save_reference(path: Path, reference: PeriodicReference2D) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.npz")
    np.savez_compressed(
        temporary,
        name=np.asarray(reference.name),
        times=reference.times,
        space=reference.space,
        values=reference.values,
        derivative_1=reference.derivative_1,
        derivative_2=reference.derivative_2,
        dt=np.asarray(reference.dt),
        store_dt=np.asarray(reference.store_dt),
        n=np.asarray(reference.n),
    )
    temporary.replace(path)


def load_reference(path: Path) -> PeriodicReference2D:
    with np.load(path) as archive:
        return PeriodicReference2D(
            name=str(archive["name"].item()),
            times=archive["times"].copy(),
            space=archive["space"].copy(),
            values=archive["values"].copy(),
            derivative_1=archive["derivative_1"].copy(),
            derivative_2=archive["derivative_2"].copy(),
            dt=float(archive["dt"]),
            store_dt=float(archive["store_dt"]),
            n=int(archive["n"]),
        )


__all__ = [
    "curvature_coverage",
    "load_reference",
    "save_reference",
    "solve_mol_rk4_reference",
    "solve_strang_reference",
    "validate_reference",
]
