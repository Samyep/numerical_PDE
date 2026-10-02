"""Controlled training-distribution ablation for direct-vector 1D Euler HCFL.

The only scientific factor changed is the final 100 trajectories in the
580-trajectory broad training set. The baseline uses the existing random-
extreme generator; the ablation replaces those 100 trajectories with balanced,
parameter-randomized Euler wave regimes. Architecture, initialization,
normalization, optimizer, update budget, safety layers, and evaluation are
identical between arms.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import torch


HERE = Path(__file__).resolve().parent
HCFL_ROOT = HERE.parents[1]
LOCAL_RUN = HCFL_ROOT / "local_run"
if str(LOCAL_RUN) not in sys.path:
    sys.path.insert(0, str(LOCAL_RUN))

import euler_ablation_runner as base  # noqa: E402


ARMS = ("broad_random", "wave_coverage")
CANONICAL_CASES = (
    "sod",
    "lax",
    "collision",
    "strong_pressure",
    "near_vacuum_expansion",
)
RHO_FLOOR = 1.0e-5
PRESSURE_FLOOR = 1.0e-5
CANONICAL_NSNAP = 64
ENTROPY_RESIDUAL_TOLERANCE = 1.0e-5
FD_ENTROPY_TOLERANCE = 1.0e-8


def strict_rollout_reference(initial: np.ndarray, nsnap: int = base.NSNAP) -> np.ndarray:
    """Generate the existing fine-grid reference without silent floor repair."""
    U = np.asarray(initial, dtype=np.float64).copy()
    ntraj = U.shape[0]
    dx = 1.0 / base.NREF

    def restrict(V: np.ndarray) -> np.ndarray:
        return V.reshape(ntraj, base.NCOARSE, base.FACTOR, 3).mean(axis=2).astype(np.float32)

    snapshots = [restrict(U)]
    for _ in range(nsnap - 1):
        remaining = base.DT_SNAPSHOT
        while remaining > 1.0e-14:
            rho = U[..., 0]
            pressure = base.np_pressure(U)
            velocity = U[..., 1] / rho
            sound_speed = np.sqrt(base.GAMMA * pressure / rho)
            max_speed = float(np.max(np.abs(velocity) + sound_speed))
            dt = min(remaining, 0.25 * dx / max(max_speed, 1.0e-12))
            U = base.np_ssprk2(U, dt, dx)
            if not np.isfinite(U).all() or (U[..., 0] <= 0).any() or (base.np_pressure(U) <= 0).any():
                raise RuntimeError("Strict reference solver lost admissibility; no state repair was applied.")
            remaining -= dt
        snapshots.append(restrict(U))
    return np.stack(snapshots, axis=1)


def make_baseline_training_data(seed: int) -> torch.Tensor:
    """Reproduce the current 220/260/100 broad baseline with strict references."""
    ordinary = strict_rollout_reference(base.generate_ic(220, base.NREF, 6000 + seed, ood=False))
    broad = strict_rollout_reference(base.generate_ic(260, base.NREF, 7000 + seed, ood=True))
    extreme = strict_rollout_reference(base.generate_extreme_ic(100, base.NREF, 8000 + seed))
    return torch.from_numpy(np.concatenate([ordinary, broad, extreme], axis=0))


def _two_state_periodic(
    N: int,
    cut: int,
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    shift: int,
) -> np.ndarray:
    rho = np.r_[np.full(cut, left[0]), np.full(N - cut, right[0])]
    velocity = np.r_[np.full(cut, left[1]), np.full(N - cut, right[1])]
    pressure = np.r_[np.full(cut, left[2]), np.full(N - cut, right[2])]
    return base.prim_to_cons(
        np.roll(rho, shift),
        np.roll(velocity, shift),
        np.roll(pressure, shift),
    )


def generate_wave_coverage_ic(ntraj: int, N: int, seed: int) -> np.ndarray:
    """Balanced contacts/compressions/expansions/pressure jumps/collisions."""
    rng = np.random.default_rng(seed)
    regimes = np.resize(np.array([
        "contact",
        "compression",
        "expansion",
        "pressure_jump",
        "collision",
    ], dtype=object), ntraj)
    rng.shuffle(regimes)
    out = np.empty((ntraj, N, 3), dtype=np.float64)

    for j, regime in enumerate(regimes):
        cut = int(rng.integers(N // 3, 2 * N // 3))
        shift = int(rng.integers(0, N))

        if regime == "contact":
            rho_low = rng.uniform(0.15, 0.75)
            rho_high = min(rho_low * rng.uniform(2.0, 5.0), 3.2)
            rho_left, rho_right = (rho_low, rho_high) if rng.random() < 0.5 else (rho_high, rho_low)
            common_velocity = rng.uniform(-1.1, 1.1)
            common_pressure = rng.uniform(0.18, 4.2)
            left = (rho_left, common_velocity, common_pressure)
            right = (rho_right, common_velocity, common_pressure)

        elif regime == "compression":
            center_velocity = rng.uniform(-0.35, 0.35)
            velocity_gap = rng.uniform(0.6, 2.8)
            left = (
                rng.uniform(0.18, 3.0),
                center_velocity + 0.5 * velocity_gap,
                rng.uniform(0.15, 4.2),
            )
            right = (
                rng.uniform(0.18, 3.0),
                center_velocity - 0.5 * velocity_gap,
                rng.uniform(0.15, 4.2),
            )

        elif regime == "expansion":
            center_velocity = rng.uniform(-0.35, 0.35)
            velocity_gap = rng.uniform(0.6, 3.0)
            left = (
                rng.uniform(0.18, 3.0),
                center_velocity - 0.5 * velocity_gap,
                rng.uniform(0.15, 4.2),
            )
            right = (
                rng.uniform(0.18, 3.0),
                center_velocity + 0.5 * velocity_gap,
                rng.uniform(0.15, 4.2),
            )

        elif regime == "pressure_jump":
            pressure_low = rng.uniform(0.08, 0.30)
            pressure_high = min(pressure_low * rng.uniform(8.0, 20.0), 5.0)
            pressure_left, pressure_right = (
                (pressure_low, pressure_high)
                if rng.random() < 0.5
                else (pressure_high, pressure_low)
            )
            left = (rng.uniform(0.20, 3.0), rng.uniform(-0.45, 0.45), pressure_left)
            right = (rng.uniform(0.20, 3.0), rng.uniform(-0.45, 0.45), pressure_right)

        elif regime == "collision":
            center_velocity = rng.uniform(-0.20, 0.20)
            left = (
                rng.uniform(0.30, 2.5),
                center_velocity + rng.uniform(0.8, 2.2),
                rng.uniform(0.10, 1.2),
            )
            right = (
                rng.uniform(0.30, 2.5),
                center_velocity - rng.uniform(0.8, 2.2),
                rng.uniform(0.10, 1.2),
            )
        else:  # pragma: no cover - guarded by the fixed regime list
            raise ValueError(regime)

        out[j] = _two_state_periodic(N, cut, left, right, shift)

    return out


def make_wave_coverage_training_data(baseline: torch.Tensor, seed: int) -> torch.Tensor:
    """Keep the first 480 baseline trajectories and replace only the final 100."""
    coverage_ic = generate_wave_coverage_ic(100, base.NREF, 9000 + seed)
    coverage = torch.from_numpy(strict_rollout_reference(coverage_ic))
    return torch.cat([baseline[:480].clone(), coverage], dim=0)


def generate_moderate_ood_ic(ntraj: int, N: int, seed: int) -> np.ndarray:
    """Moderate smooth frequency shift: modes 4--6 versus training modes 1--3."""
    rng = np.random.default_rng(seed)
    x = (np.arange(N) + 0.5) / N
    out = np.empty((ntraj, N, 3), dtype=np.float64)
    for j in range(ntraj):
        rho0 = rng.uniform(0.90, 1.55)
        pressure0 = rng.uniform(0.90, 1.55)
        rho = rho0 + rng.uniform(0.06, 0.24) * np.sin(
            2 * np.pi * rng.integers(4, 7) * x + rng.uniform(0, 2 * np.pi)
        )
        pressure = pressure0 + rng.uniform(0.06, 0.24) * np.sin(
            2 * np.pi * rng.integers(4, 7) * x + rng.uniform(0, 2 * np.pi)
        )
        velocity = rng.uniform(0.15, 0.65) * np.sin(
            2 * np.pi * rng.integers(4, 7) * x + rng.uniform(0, 2 * np.pi)
        )
        out[j] = base.prim_to_cons(rho, velocity, pressure)
    return out


def canonical_initial_conditions() -> dict[str, np.ndarray]:
    """Fixed canonical tests, none of which is inserted into training."""
    states = {
        "sod": ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
        "lax": ((0.445, 0.698, 3.528), (0.5, 0.0, 0.571)),
        "collision": ((1.0, 2.0, 1.0), (1.0, -2.0, 1.0)),
        "strong_pressure": ((1.0, 0.0, 5.0), (1.0, 0.0, 0.05)),
        "near_vacuum_expansion": ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4)),
    }
    return {
        name: _two_state_periodic(base.NREF, base.NREF // 2, left, right, 0)[None, ...]
        for name, (left, right) in states.items()
    }


def make_evaluation_suite(seed: int) -> dict[str, torch.Tensor]:
    suite = {
        "ordinary_id": torch.from_numpy(strict_rollout_reference(
            base.generate_ic(90, base.NREF, 2000 + seed, ood=False)
        )),
        "broad_random_in_support": torch.from_numpy(strict_rollout_reference(
            base.generate_ic(90, base.NREF, 3000 + seed, ood=True)
        )),
        "moderate_ood_high_frequency": torch.from_numpy(strict_rollout_reference(
            generate_moderate_ood_ic(90, base.NREF, 4000 + seed)
        )),
    }
    for name, initial in canonical_initial_conditions().items():
        suite[name] = torch.from_numpy(strict_rollout_reference(initial, nsnap=CANONICAL_NSNAP))
    return suite


def admissible(
    U: torch.Tensor,
    rho_floor: float = RHO_FLOOR,
    pressure_floor: float = PRESSURE_FLOOR,
) -> torch.Tensor:
    return (U[..., 0] >= rho_floor) & (base.t_pressure(U) >= pressure_floor)


def total_entropy(U: torch.Tensor) -> torch.Tensor:
    return base.entropy(U.double()).sum(dim=-1)


def entropy_residual64(Fh: torch.Tensor, U: torch.Tensor) -> torch.Tensor:
    """Evaluate the Tadmor residual in float64, including the state transform."""
    return base.entropy_residual(Fh.double(), U.double())


def strict_entropy_projection(Fh: torch.Tensor, U: torch.Tensor) -> torch.Tensor:
    """Project with a roundoff backoff and verify after conversion to float32.

    The base model intentionally remains unchanged. This function closes a
    finite-precision leak that appears only for very large near-vacuum entropy
    variables: an exact double projection can land a few float32 ulps outside
    the half-space after casting.
    """
    output = Fh
    eps = torch.finfo(Fh.dtype).eps
    U64 = U.double()
    right = torch.roll(U64, -1, dims=-2)
    normal = base.entropy_variables(right) - base.entropy_variables(U64)
    bound = base.entropy_potential(right) - base.entropy_potential(U64)
    norm_squared = (normal * normal).sum(dim=-1)

    for _ in range(4):
        flux64 = output.double()
        residual = (normal * flux64).sum(dim=-1) - bound
        scale = bound.abs() + (normal * flux64).abs().sum(dim=-1)
        target = -8.0 * eps * scale
        mask = (residual > target) & (norm_squared > 1.0e-14)
        if not bool(mask.any()):
            break
        correction = torch.zeros_like(residual)
        correction[mask] = (residual[mask] - target[mask]) / norm_squared[mask]
        output = (flux64 - correction[..., None] * normal).to(Fh.dtype)

    return output


def local_admissibility_limiter(
    U: torch.Tensor,
    F_hi: torch.Tensor,
    F_lo: torch.Tensor,
    lam: float,
    max_outer: int = 12,
    n_bisect: int = 30,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Existing shared-interface local convex limiter, made self-contained."""
    delta_flux = F_hi - F_lo
    U_lo = U - lam * (F_lo - torch.roll(F_lo, 1, dims=-2))
    if not bool(admissible(U_lo).all()):
        raise RuntimeError("Low-order state is not admissible; reduce the time step.")

    batch, cells, _ = U.shape
    alpha = torch.ones(batch, cells, dtype=U.dtype, device=U.device)

    for _ in range(max_outer):
        flux = F_lo + alpha[..., None] * delta_flux
        U_current = U - lam * (flux - torch.roll(flux, 1, dims=-2))
        good = admissible(U_current)
        if bool(good.all()):
            return flux, alpha

        bad = ~good
        delta_state = U_current - U_lo
        count = int(bad.sum().item())
        lo = torch.zeros(count, dtype=torch.float64, device=U.device)
        hi = torch.ones(count, dtype=torch.float64, device=U.device)
        U0 = U_lo[bad].double()
        direction = delta_state[bad].double()

        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            okay = admissible(U0 + mid[:, None] * direction)
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)

        theta = torch.ones(batch, cells, dtype=U.dtype, device=U.device)
        theta[bad] = torch.clamp(lo.to(U.dtype) - 1.0e-5, min=0.0)
        interface_factor = torch.minimum(theta, torch.roll(theta, -1, dims=-1))
        alpha = alpha * interface_factor

    flux = F_lo + alpha[..., None] * delta_flux
    U_current = U - lam * (flux - torch.roll(flux, 1, dims=-2))
    if not bool(admissible(U_current).all()):
        bad_trajectory = (~admissible(U_current)).any(dim=-1)
        alpha[bad_trajectory] = 0.0
        flux = F_lo + alpha[..., None] * delta_flux
    if not bool(admissible(U - lam * (flux - torch.roll(flux, 1, dims=-2))).all()):
        raise RuntimeError("Local admissibility limiter failed after conservative fallback.")
    return flux, alpha


def global_entropy_limiter(
    U: torch.Tensor,
    F_hi: torch.Tensor,
    F_lo: torch.Tensor,
    lam: float,
    n_bisect: int = 40,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Existing per-trajectory fully-discrete total-entropy line search."""
    U_lo = U - lam * (F_lo - torch.roll(F_lo, 1, dims=-2))
    U_hi = U - lam * (F_hi - torch.roll(F_hi, 1, dims=-2))
    entropy_before = total_entropy(U)
    entropy_hi = total_entropy(U_hi)
    beta = torch.ones(U.shape[0], dtype=U.dtype, device=U.device)
    need = entropy_hi > entropy_before

    if bool(need.any()):
        count = int(need.sum().item())
        lo = torch.zeros(count, dtype=torch.float64, device=U.device)
        hi = torch.ones(count, dtype=torch.float64, device=U.device)
        U0 = U_lo[need].double()
        direction = (U_hi - U_lo)[need].double()
        entropy_target = entropy_before[need]

        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            okay = total_entropy(U0 + mid[:, None, None] * direction) <= entropy_target
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)

        beta[need] = torch.clamp(lo.to(U.dtype) - 1.0e-6, min=0.0)

    # Recheck the actual reconstructed float32 flux/update. The double-state
    # interpolation above can otherwise choose a boundary coefficient that
    # rounds just outside the entropy set when converted back to float32.
    delta_flux = F_hi - F_lo
    flux = F_lo + beta[:, None, None] * delta_flux
    U_check = U - lam * (flux - torch.roll(flux, 1, dims=-2))
    rounded_bad = total_entropy(U_check) > entropy_before

    if bool(rounded_bad.any()):
        count = int(rounded_bad.sum().item())
        lo = torch.zeros(count, dtype=torch.float64, device=U.device)
        hi = beta[rounded_bad].double()
        U_bad = U[rounded_bad]
        F_lo_bad = F_lo[rounded_bad]
        delta_bad = delta_flux[rounded_bad]
        target_bad = entropy_before[rounded_bad]

        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            flux_mid = F_lo_bad + mid.to(U.dtype)[:, None, None] * delta_bad
            state_mid = U_bad - lam * (flux_mid - torch.roll(flux_mid, 1, dims=-2))
            okay = total_entropy(state_mid) <= target_bad
            lo = torch.where(okay, mid, lo)
            hi = torch.where(okay, hi, mid)

        beta[rounded_bad] = torch.clamp(lo.to(U.dtype) - 1.0e-6, min=0.0)
        flux = F_lo + beta[:, None, None] * delta_flux

    # Float32 evaluation can be locally non-monotone at a few ulps. A final
    # checked fallback to the already-safe low-order endpoint closes that leak.
    U_final = U - lam * (flux - torch.roll(flux, 1, dims=-2))
    still_bad = total_entropy(U_final) > entropy_before + FD_ENTROPY_TOLERANCE
    if bool(still_bad.any()):
        beta[still_bad] = 0.0
        flux = F_lo + beta[:, None, None] * delta_flux

    return flux, beta


@torch.no_grad()
def advance_safe_snapshot(
    model: base.Solver,
    U: torch.Tensor,
    cfl: float = 0.42,
) -> tuple[torch.Tensor, dict[str, float]]:
    """Advance one saved interval using the fixed full HCFL safety stack."""
    dx = 1.0 / base.NCOARSE
    remaining = base.DT_SNAPSHOT
    stats = {
        "local_active": 0.0,
        "local_total": 0.0,
        "fd_active": 0.0,
        "fd_total": 0.0,
        "fd_beta_sum": 0.0,
        "fd_beta_min": 1.0,
        "entropy_violations": 0.0,
        "entropy_total": 0.0,
        "max_entropy_residual": -float("inf"),
        "max_total_entropy_change": -float("inf"),
        "substeps": 0.0,
    }

    while remaining > 1.0e-14:
        rho = U[..., 0].clamp_min(1.0e-10)
        pressure = base.t_pressure(U).clamp_min(1.0e-10)
        velocity = U[..., 1] / rho
        sound_speed = torch.sqrt(base.GAMMA * pressure / rho)
        max_speed = float((torch.abs(velocity) + sound_speed).max())
        dt = min(remaining, cfl * dx / max(max_speed, 1.0e-12))

        # Enforce the checked low-order premise for both positivity and the
        # fully-discrete total-entropy line search.
        for _ in range(24):
            lam = dt / dx
            F_lo = strict_entropy_projection(base.t_rusanov(U), U)
            U_lo = U - lam * (F_lo - torch.roll(F_lo, 1, dims=-2))
            entropy_ok = bool((total_entropy(U_lo) <= total_entropy(U) + 1.0e-10).all())
            if bool(admissible(U_lo).all()) and entropy_ok:
                break
            dt *= 0.5
        else:
            raise RuntimeError("Could not establish the low-order safety premise.")

        F_hi = strict_entropy_projection(model.flux(U), U)
        F_local, alpha = local_admissibility_limiter(U, F_hi, F_lo, lam)
        F_final, beta = global_entropy_limiter(U, F_local, F_lo, lam)

        residual = entropy_residual64(F_final, U)
        U_next = U - lam * (F_final - torch.roll(F_final, 1, dims=-2))
        entropy_change = total_entropy(U_next) - total_entropy(U)

        if not torch.isfinite(U_next).all() or not bool(admissible(U_next).all()):
            raise RuntimeError("Safe HCFL update produced a non-admissible state.")

        stats["local_active"] += float((alpha < 1.0 - 1.0e-7).sum())
        stats["local_total"] += float(alpha.numel())
        stats["fd_active"] += float((beta < 1.0 - 1.0e-7).sum())
        stats["fd_total"] += float(beta.numel())
        stats["fd_beta_sum"] += float(beta.sum())
        stats["fd_beta_min"] = min(stats["fd_beta_min"], float(beta.min()))
        stats["entropy_violations"] += float((residual > ENTROPY_RESIDUAL_TOLERANCE).sum())
        stats["entropy_total"] += float(residual.numel())
        stats["max_entropy_residual"] = max(
            stats["max_entropy_residual"], float(residual.max())
        )
        stats["max_total_entropy_change"] = max(
            stats["max_total_entropy_change"], float(entropy_change.max())
        )
        stats["substeps"] += 1.0

        U = U_next
        remaining -= dt

    return U, stats


def merge_step_stats(total: dict[str, float], step: dict[str, float]) -> None:
    for key in (
        "local_active",
        "local_total",
        "fd_active",
        "fd_total",
        "fd_beta_sum",
        "entropy_violations",
        "entropy_total",
        "substeps",
    ):
        total[key] += step[key]
    total["fd_beta_min"] = min(total["fd_beta_min"], step["fd_beta_min"])
    total["max_entropy_residual"] = max(
        total["max_entropy_residual"], step["max_entropy_residual"]
    )
    total["max_total_entropy_change"] = max(
        total["max_total_entropy_change"], step["max_total_entropy_change"]
    )


@torch.no_grad()
def evaluate_dataset(
    model: base.Solver,
    data: torch.Tensor,
    state_std: torch.Tensor,
    split: str,
) -> dict[str, Any]:
    started = time.perf_counter()
    U = data[:, 0].clone()
    errors: list[float] = []
    min_rho = float(U[..., 0].min())
    min_pressure = float(base.t_pressure(U).min())
    totals = defaultdict(float)
    totals["fd_beta_min"] = 1.0
    totals["max_entropy_residual"] = -float("inf")
    totals["max_total_entropy_change"] = -float("inf")

    for snapshot in range(1, data.shape[1]):
        U, step_stats = advance_safe_snapshot(model, U)
        merge_step_stats(totals, step_stats)
        min_rho = min(min_rho, float(U[..., 0].min()))
        min_pressure = min(min_pressure, float(base.t_pressure(U).min()))
        errors.append(float((((U - data[:, snapshot]) / state_std) ** 2).mean()))

    return {
        "split": split,
        "ntraj": int(data.shape[0]),
        "rollout_nrmse": float(np.sqrt(np.mean(errors))),
        "min_rho": min_rho,
        "min_pressure": min_pressure,
        "entropy_violation_rate": totals["entropy_violations"] / max(totals["entropy_total"], 1.0),
        "max_entropy_residual": totals["max_entropy_residual"],
        "local_limiter_intervention_rate": totals["local_active"] / max(totals["local_total"], 1.0),
        "fd_entropy_intervention_rate": totals["fd_active"] / max(totals["fd_total"], 1.0),
        "mean_fd_beta": totals["fd_beta_sum"] / max(totals["fd_total"], 1.0),
        "min_fd_beta": totals["fd_beta_min"],
        "max_total_entropy_change": totals["max_total_entropy_change"],
        "mean_substeps_per_snapshot": totals["substeps"] / (data.shape[1] - 1),
        "eval_seconds": time.perf_counter() - started,
    }


def train_model(
    train_data: torch.Tensor,
    mean: np.ndarray,
    std: np.ndarray,
    state_std: torch.Tensor,
    seed: int,
    iterations: int,
    width: int,
    batch_size: int,
    learning_rate: float,
    arm: str,
) -> tuple[base.Solver, list[dict[str, Any]], float]:
    """Train one arm using an identical initialization and minibatch schedule."""
    torch.manual_seed(12000 + seed)
    model = base.Solver("direct", mean, std, width=width)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    generator = torch.Generator().manual_seed(13000 + seed)
    curve: list[dict[str, Any]] = []
    started = time.perf_counter()

    model.train()
    for iteration in range(iterations):
        indices = torch.randint(0, train_data.shape[0], (batch_size,), generator=generator)
        times = torch.randint(0, train_data.shape[1] - 1, (batch_size,), generator=generator)
        U = train_data[indices, times]
        target = train_data[indices, times + 1]
        prediction = model.one_step(U)
        loss = (((prediction - target) / state_std) ** 2).mean()

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        if iteration % 100 == 0 or iteration == iterations - 1:
            value = float(loss.detach())
            curve.append({"seed": seed, "arm": arm, "iteration": iteration, "loss": value})
            print(json.dumps({"seed": seed, "arm": arm, "iteration": iteration, "loss": value}))

    elapsed = time.perf_counter() - started
    model.eval()
    return model, curve, elapsed


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_comparison(metrics: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {(row["arm"], row["split"]): row for row in metrics}
    rows: list[dict[str, Any]] = []
    split_order = (
        "ordinary_id",
        "broad_random_in_support",
        "moderate_ood_high_frequency",
        *CANONICAL_CASES,
    )
    for split in split_order:
        baseline = indexed[("broad_random", split)]
        coverage = indexed[("wave_coverage", split)]
        old = float(baseline["rollout_nrmse"])
        new = float(coverage["rollout_nrmse"])
        rows.append({
            "seed": baseline["seed"],
            "split": split,
            "baseline_nrmse": old,
            "wave_coverage_nrmse": new,
            "relative_nrmse_change_percent": 100.0 * (new / old - 1.0),
            "baseline_local_limiter_rate": baseline["local_limiter_intervention_rate"],
            "wave_coverage_local_limiter_rate": coverage["local_limiter_intervention_rate"],
            "baseline_fd_entropy_rate": baseline["fd_entropy_intervention_rate"],
            "wave_coverage_fd_entropy_rate": coverage["fd_entropy_intervention_rate"],
        })
    return rows


def screening_decision(metrics: list[dict[str, Any]]) -> dict[str, Any]:
    indexed = {(row["arm"], row["split"]): row for row in metrics}

    def relative_change(split: str) -> float:
        old = float(indexed[("broad_random", split)]["rollout_nrmse"])
        new = float(indexed[("wave_coverage", split)]["rollout_nrmse"])
        return new / old - 1.0

    ordinary_change = relative_change("ordinary_id")
    moderate_change = relative_change("moderate_ood_high_frequency")
    canonical_changes = {case: relative_change(case) for case in CANONICAL_CASES}
    old_mean = float(np.mean([
        indexed[("broad_random", case)]["rollout_nrmse"] for case in CANONICAL_CASES
    ]))
    new_mean = float(np.mean([
        indexed[("wave_coverage", case)]["rollout_nrmse"] for case in CANONICAL_CASES
    ]))
    improved_count = sum(change < 0.0 for change in canonical_changes.values())
    safe = all(
        float(row["min_rho"]) >= RHO_FLOOR
        and float(row["min_pressure"]) >= PRESSURE_FLOOR
        and float(row["entropy_violation_rate"]) == 0.0
        and float(row["max_total_entropy_change"]) <= FD_ENTROPY_TOLERANCE
        for row in metrics
    )
    passed = (
        ordinary_change <= 0.05
        and moderate_change <= 0.05
        and improved_count >= 3
        and (new_mean / old_mean - 1.0) <= -0.05
        and max(canonical_changes.values()) <= 0.10
        and safe
    )
    return {
        "passed": passed,
        "recommendation": "confirm_with_more_seeds" if passed else "stop_or_modify",
        "ordinary_id_relative_change": ordinary_change,
        "moderate_ood_relative_change": moderate_change,
        "canonical_relative_changes": canonical_changes,
        "canonical_improved_count": improved_count,
        "baseline_canonical_mean_nrmse": old_mean,
        "wave_coverage_canonical_mean_nrmse": new_mean,
        "canonical_mean_relative_change": new_mean / old_mean - 1.0,
        "all_rollouts_safe": safe,
    }


def self_test() -> None:
    initial = generate_wave_coverage_ic(5, base.NREF, 123)
    assert initial.shape == (5, base.NREF, 3)
    assert np.all(initial[..., 0] > 0)
    assert np.all(base.np_pressure(initial) > 0)
    trajectory = strict_rollout_reference(initial)
    train = torch.from_numpy(trajectory)
    primitive = base.primitive(train)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    model = base.Solver("direct", mean, std, width=16).eval()
    updated, stats = advance_safe_snapshot(model, train[:, 0].clone())
    assert bool(admissible(updated).all())
    assert stats["entropy_violations"] == 0.0
    print("self-test passed")


def run(args: argparse.Namespace) -> None:
    output = Path(args.outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    total_started = time.perf_counter()

    print("Generating strict baseline training data...")
    data_started = time.perf_counter()
    baseline_data = make_baseline_training_data(args.seed)
    coverage_data = make_wave_coverage_training_data(baseline_data, args.seed)
    evaluation_suite = make_evaluation_suite(args.seed)
    data_seconds = time.perf_counter() - data_started

    # Deliberately derive normalization only from the baseline and reuse it in
    # both arms so preprocessing/loss scaling cannot confound the data ablation.
    primitive = base.primitive(baseline_data)
    mean = primitive.mean(dim=(0, 1, 2)).numpy()
    std = primitive.std(dim=(0, 1, 2)).numpy()
    state_std = torch.tensor([
        float(baseline_data[..., channel].std()) for channel in range(3)
    ])

    datasets = {
        "broad_random": baseline_data,
        "wave_coverage": coverage_data,
    }
    metrics: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    run_rows: list[dict[str, Any]] = []

    for arm in ARMS:
        model, curve, train_seconds = train_model(
            datasets[arm],
            mean,
            std,
            state_std,
            args.seed,
            args.iters,
            args.width,
            args.batch_size,
            args.lr,
            arm,
        )
        curves.extend(curve)
        parameter_count = sum(parameter.numel() for parameter in model.parameters())

        eval_started = time.perf_counter()
        for split, data in evaluation_suite.items():
            row = evaluate_dataset(model, data, state_std, split)
            row.update({
                "seed": args.seed,
                "arm": arm,
                "parameter_count": parameter_count,
                "train_seconds": train_seconds,
            })
            metrics.append(row)
            print(json.dumps({
                "seed": args.seed,
                "arm": arm,
                "split": split,
                "nrmse": row["rollout_nrmse"],
            }))
        run_rows.append({
            "seed": args.seed,
            "arm": arm,
            "parameter_count": parameter_count,
            "train_seconds": train_seconds,
            "evaluation_seconds": time.perf_counter() - eval_started,
            "final_logged_loss": curve[-1]["loss"],
        })

    comparison = make_comparison(metrics)
    decision = screening_decision(metrics)
    write_csv(output / f"metrics_seed{args.seed}.csv", metrics)
    write_csv(output / f"comparison_seed{args.seed}.csv", comparison)
    write_csv(output / f"training_curve_seed{args.seed}.csv", curves)
    write_csv(output / f"runtime_seed{args.seed}.csv", run_rows)

    metadata = {
        "seed": args.seed,
        "single_scientific_factor": "training_data_distribution",
        "baseline_training_composition": {
            "ordinary": 220,
            "broad_random": 260,
            "random_extreme": 100,
        },
        "wave_coverage_training_composition": {
            "ordinary": 220,
            "broad_random": 260,
            "balanced_wave_regimes": 100,
        },
        "model": "direct-vector HLLC-HCFL",
        "iterations": args.iters,
        "width": args.width,
        "batch_size": args.batch_size,
        "learning_rate": args.lr,
        "reference": "512-cell Rusanov + SSP-RK2, strict/no repair",
        "ordinary_and_moderate_nsnap": base.NSNAP,
        "canonical_nsnap": CANONICAL_NSNAP,
        "entropy_residual_tolerance": ENTROPY_RESIDUAL_TOLERANCE,
        "fully_discrete_entropy_tolerance": FD_ENTROPY_TOLERANCE,
        "pyclaw_available": False,
        "data_generation_seconds": data_seconds,
        "total_wall_seconds": time.perf_counter() - total_started,
        "screening_decision": decision,
    }
    (output / f"metadata_seed{args.seed}.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(decision, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--iters", type=int, default=1100)
    parser.add_argument("--width", type=int, default=72)
    parser.add_argument("--batch-size", type=int, default=56)
    parser.add_argument("--lr", type=float, default=3.0e-4)
    parser.add_argument("--outdir", default=str(HERE / "results"))
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    parsed = parse_args()
    if parsed.self_test:
        self_test()
    else:
        run(parsed)
