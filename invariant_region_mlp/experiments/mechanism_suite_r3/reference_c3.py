"""Production one-dimensional reference for the round-3 cubic HJ equation.

The spatial Hamiltonian is discretised with the monotone Godunov flux used by
the candidate screen.  Diffusion is advanced implicitly, while the Godunov
term is explicit with ``dt / ds = 1/32``.  The four registered spatial grids
therefore also halve the time step at every refinement.  All Richardson
comparisons are made at shared nodes; no interpolated values enter the
refinement or domain gates.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.linalg import solve_banded


DS_LEVELS = (0.02, 0.01, 0.005, 0.0025)
DOMAINS = (6.0, 8.0)
OUTPUT_STEPS = 400
REFERENCE_SCHEMA = 1


def terminal_profile(s: np.ndarray, *, A: float, beta: float) -> np.ndarray:
    s = np.asarray(s, dtype=np.float64)
    return A * (
        np.logaddexp(beta * s, -beta * s) - math.log(2.0)
    ) / beta


def solve_monotone_godunov(
    *,
    kappa: float,
    A: float,
    beta: float,
    T: float,
    L: float,
    ds: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Solve ``psi_tau = psi_ss - kappa |psi_s|^3`` on ``[-L,L]``.

    The asymptotically exact affine boundary value
    ``g(+-L) - kappa A^3 tau`` is imposed at both ends.  The implicit
    diffusion matrix is an M-matrix and the explicit Godunov step obeys the
    Hamiltonian CFL bound for ``|psi_s| <= A``.
    """

    intervals = int(round(2.0 * L / ds))
    if not math.isclose(intervals * ds, 2.0 * L, abs_tol=1e-12):
        raise ValueError("L and ds must define an aligned uniform grid")
    s = np.linspace(-L, L, intervals + 1, dtype=np.float64)
    tau = np.linspace(0.0, T, OUTPUT_STEPS + 1, dtype=np.float64)
    # dt/ds=1/32 gives 3*kappa*A^2*dt/ds=0.1875 for the registered
    # (kappa,A)=(0.5,2).  The conservative choice also leaves a sufficiently
    # fine common time grid for the independent residual check.
    requested_dt = ds / 32.0
    substeps = int(round((T / OUTPUT_STEPS) / requested_dt))
    if substeps < 1 or not math.isclose(
        substeps * requested_dt, T / OUTPUT_STEPS, rel_tol=0.0, abs_tol=1e-14
    ):
        raise ValueError("registered output spacing must be divisible by ds/32")
    dt = (T / OUTPUT_STEPS) / substeps
    cfl = 3.0 * kappa * A * A * dt / ds
    if cfl > 1.0 + 1e-14:
        raise ValueError(f"Hamiltonian CFL violated: {cfl}")

    value = terminal_profile(s, A=A, beta=beta)
    snapshots = np.empty((len(tau), len(s)), dtype=np.float64)
    snapshots[0] = value
    interior = len(s) - 2
    ratio = dt / (ds * ds)
    band = np.zeros((3, interior), dtype=np.float64)
    band[0, 1:] = -ratio
    band[1, :] = 1.0 + 2.0 * ratio
    band[2, :-1] = -ratio

    step_index = 0
    for output_index in range(1, len(tau)):
        for _ in range(substeps):
            step_index += 1
            backward = (value[1:-1] - value[:-2]) / ds
            forward = (value[2:] - value[1:-1]) / ds
            godunov_slope = np.maximum(
                np.maximum(backward, 0.0), np.maximum(-forward, 0.0)
            )
            hamiltonian = kappa * godunov_slope**3
            rhs = value[1:-1] - dt * hamiltonian
            next_tau = step_index * dt
            left = terminal_profile(np.asarray([-L]), A=A, beta=beta)[0]
            right = terminal_profile(np.asarray([L]), A=A, beta=beta)[0]
            boundary_shift = kappa * A**3 * next_tau
            left -= boundary_shift
            right -= boundary_shift
            rhs[0] += ratio * left
            rhs[-1] += ratio * right
            next_value = np.empty_like(value)
            next_value[0] = left
            next_value[-1] = right
            next_value[1:-1] = solve_banded(
                (1, 1), band, rhs, overwrite_ab=False, overwrite_b=True,
                check_finite=False,
            )
            value = next_value
        snapshots[output_index] = value

    metadata = {
        "L": L,
        "ds": ds,
        "dt": dt,
        "space_points": len(s),
        "time_steps": step_index,
        "output_time_points": len(tau),
        "hamiltonian_cfl": cfl,
        "scheme": "implicit diffusion plus explicit monotone Godunov Hamiltonian",
        "boundary": "asymptotic affine Dirichlet g(+-L)-kappa*A^3*tau",
    }
    return tau, s, snapshots, metadata


def _inner_max(value: np.ndarray, s: np.ndarray, radius: float = 3.0) -> float:
    return float(np.max(np.abs(value[:, np.abs(s) <= radius])))


def _hierarchy(
    solutions: list[tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]]
) -> dict[str, Any]:
    tau = solutions[0][0]
    grids = [item[1] for item in solutions]
    values = [item[2] for item in solutions]
    if not all(np.array_equal(tau, item[0]) for item in solutions):
        raise RuntimeError("reference time grids are not shared")

    raw_differences = []
    first_richardson = []
    for index in range(3):
        fine_on_coarse = values[index + 1][:, ::2]
        raw_differences.append(
            _inner_max(fine_on_coarse - values[index], grids[index])
        )
        first_richardson.append(2.0 * fine_on_coarse - values[index])

    first_differences = [
        _inner_max(
            first_richardson[index + 1][:, ::2] - first_richardson[index],
            grids[index],
        )
        for index in range(2)
    ]
    second_coarse = (
        4.0 * first_richardson[1][:, ::2] - first_richardson[0]
    ) / 3.0
    second_fine = (
        4.0 * first_richardson[2][:, ::2] - first_richardson[1]
    ) / 3.0
    finest_difference = _inner_max(
        second_fine[:, ::2] - second_coarse, grids[0]
    )
    decreasing = bool(
        raw_differences[1] < raw_differences[0]
        and raw_differences[2] < raw_differences[1]
        and first_differences[1] < first_differences[0]
    )
    return {
        "tau": tau,
        "grids": grids,
        "values": values,
        "first_richardson": first_richardson,
        "second_coarse": second_coarse,
        "second_fine": second_fine,
        "raw_successive_differences": raw_differences,
        "first_richardson_successive_differences": first_differences,
        "finest_richardson_difference": finest_difference,
        "successive_differences_decrease": decreasing,
    }


def _independent_fd_residual(
    spline: RectBivariateSpline,
    *,
    kappa: float,
    T: float,
    seed: int,
    n_points: int = 4096,
    time_step: float = 5.0e-4,
    spatial_step: float = 1.0e-2,
) -> dict[str, Any]:
    """Fourth-order check at the stored reference's resolved mesh scales."""

    rng = np.random.default_rng(np.random.SeedSequence([seed, 303, 1]))
    tau = rng.uniform(3.0 * time_step, T - 3.0 * time_step, n_points)
    s = rng.uniform(-3.0 + 3.0 * spatial_step, 3.0 - 3.0 * spatial_step, n_points)

    def ev(tt: np.ndarray, ss: np.ndarray) -> np.ndarray:
        return spline.ev(tt, ss)

    psi_tau = (
        -ev(tau + 2.0 * time_step, s)
        + 8.0 * ev(tau + time_step, s)
        - 8.0 * ev(tau - time_step, s)
        + ev(tau - 2.0 * time_step, s)
    ) / (12.0 * time_step)
    psi_s = (
        -ev(tau, s + 2.0 * spatial_step)
        + 8.0 * ev(tau, s + spatial_step)
        - 8.0 * ev(tau, s - spatial_step)
        + ev(tau, s - 2.0 * spatial_step)
    ) / (12.0 * spatial_step)
    psi_ss = (
        -ev(tau, s + 2.0 * spatial_step)
        + 16.0 * ev(tau, s + spatial_step)
        - 30.0 * ev(tau, s)
        + 16.0 * ev(tau, s - spatial_step)
        - ev(tau, s - 2.0 * spatial_step)
    ) / (12.0 * spatial_step**2)
    residual = -psi_tau + psi_ss - kappa * np.abs(psi_s) ** 3
    return {
        "kind": "independent fourth-order finite-difference residual of the interpolated reference",
        "n_points": n_points,
        "time_step": time_step,
        "spatial_step": spatial_step,
        "step_rationale": "fixed to the resolved storage scales (dtau=0.000625, ds=0.01), not tuned per point",
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rmse_residual": float(np.sqrt(np.mean(residual**2))),
    }


def build_and_audit_reference(
    path: Path,
    *,
    kappa: float = 0.5,
    A: float = 2.0,
    beta: float = 2.0,
    T: float = 0.25,
    seed: int = 20261207,
) -> dict[str, Any]:
    """Build the registered eight solves, audit them, and save atomically."""

    hierarchies: dict[str, dict[str, Any]] = {}
    solve_metadata: list[dict[str, Any]] = []
    for L in DOMAINS:
        solves = []
        for ds in DS_LEVELS:
            solved = solve_monotone_godunov(
                kappa=kappa, A=A, beta=beta, T=T, L=L, ds=ds
            )
            solves.append(solved)
            solve_metadata.append(solved[3])
        hierarchies[str(int(L))] = _hierarchy(solves)

    h6 = hierarchies["6"]
    h8 = hierarchies["8"]
    s6 = h6["grids"][1]
    s8 = h8["grids"][1]
    offset = int(round((8.0 - 6.0) / DS_LEVELS[1]))
    l8_on_l6 = h8["second_fine"][:, offset : offset + len(s6)]
    domain_difference = _inner_max(l8_on_l6 - h6["second_fine"], s6)

    tau = h8["tau"]
    s = s8
    value = h8["second_fine"]
    spline = RectBivariateSpline(tau, s, value, kx=3, ky=3, s=0.0)
    fd = _independent_fd_residual(
        spline, kappa=kappa, T=T, seed=seed
    )
    finest = float(h8["finest_richardson_difference"])
    decreasing = bool(
        h6["successive_differences_decrease"]
        and h8["successive_differences_decrease"]
    )
    passed = bool(
        decreasing
        and finest <= 1.0e-6
        and domain_difference <= 1.0e-7
        and fd["max_abs_residual"] <= 1.0e-5
    )
    metadata = {
        "schema_version": REFERENCE_SCHEMA,
        "parameters": {"kappa": kappa, "A": A, "beta": beta, "T": T},
        "registered_domains": list(DOMAINS),
        "registered_ds": list(DS_LEVELS),
        "stored_reference": "L=8 second Richardson estimate on ds=0.01 shared nodes",
        "solves": solve_metadata,
        "hierarchies": {
            key: {
                "raw_successive_differences": item["raw_successive_differences"],
                "first_richardson_successive_differences": item[
                    "first_richardson_successive_differences"
                ],
                "finest_richardson_difference": item[
                    "finest_richardson_difference"
                ],
                "successive_differences_decrease": item[
                    "successive_differences_decrease"
                ],
            }
            for key, item in hierarchies.items()
        },
        "finest_richardson_difference_L8": finest,
        "finest_richardson_threshold": 1.0e-6,
        "domain_L6_L8_max_difference": domain_difference,
        "domain_threshold": 1.0e-7,
        "independent_fd_residual": fd,
        "fd_residual_threshold": 1.0e-5,
        "successive_differences_decrease": decreasing,
        "passed": passed,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        tau=tau,
        s=s,
        value=value,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    temporary.replace(path)
    return metadata


def load_reference(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    with np.load(path, allow_pickle=False) as data:
        tau = np.asarray(data["tau"], dtype=np.float64)
        s = np.asarray(data["s"], dtype=np.float64)
        value = np.asarray(data["value"], dtype=np.float64)
        metadata = json.loads(str(data["metadata_json"]))
    return tau, s, value, metadata
