"""Phase-2 safety layers for 1D shallow-water HCFL.

Adds three fully discrete safeguards around the learned entropy-projected flux:
1) adaptive CFL substepping so the Rusanov low-order step is positivity preserving;
2) a conservative interface correction limiter that guarantees h >= h_floor;
3) an optional per-trajectory global convex line search that guarantees
   non-increase of the periodic-domain total physical entropy for the step.

The entropy half-space guarantee is preserved because every final interface flux
is a convex combination of two entropy-feasible fluxes (Rusanov and the HCFL
candidate). Conservation is preserved because one shared interface flux is used
by adjacent cells.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Tuple

import numpy as np
import torch

import swe_1d_hcfl as base

DX = 1.0 / base.NCOARSE
DT_SNAPSHOT = base.DT_SNAPSHOT


def fv_step_lambda(U: torch.Tensor, flux: torch.Tensor, lam: float) -> torch.Tensor:
    return U - lam * (flux - torch.roll(flux, 1, dims=-2))


def total_entropy(U: torch.Tensor) -> torch.Tensor:
    """Periodic-domain physical entropy, accumulated in float64."""
    V = U.double()
    h = V[..., 0].clamp_min(1e-14)
    m = V[..., 1]
    eta = 0.5 * m * m / h + 0.5 * base.G * h * h
    return eta.sum(dim=-1)


def positivity_limited_flux(
    U: torch.Tensor,
    candidate_flux: torch.Tensor,
    lam: float,
    h_floor: float = 1e-4,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """Blend the entropy-feasible candidate with Rusanov to preserve h >= h_floor.

    Let F = F_L + alpha_{i+1/2}(F_H-F_L), where F_L is Rusanov and
    alpha is shared by the two cells adjacent to an interface.  A cellwise
    budget bounds the sum of negative height corrections.  Setting each
    interface alpha to the minimum of its two neighboring cell factors gives a
    simple sufficient positivity condition.

    This guarantee assumes the low-order Rusanov forward-Euler step is itself
    positive; ``advance_snapshot_safe`` enforces a conservative CFL condition
    before calling this routine.
    """
    Flo = base.t_rusanov(U)
    dF = candidate_flux - Flo
    Ulo = fv_step_lambda(U, Flo, lam)

    dFh = dF[..., 0]
    q_right = -lam * dFh
    q_left = lam * torch.roll(dFh, 1, dims=-1)

    dangerous = torch.relu(-q_right) + torch.relu(-q_left)
    budget = torch.clamp(Ulo[..., 0] - h_floor, min=0.0)

    theta = torch.ones_like(budget)
    active = dangerous > budget
    theta[active] = budget[active] / (dangerous[active] + 1e-30)
    theta = theta.clamp(0.0, 1.0)

    # Interface i+1/2 affects cells i and i+1.
    alpha = torch.minimum(theta, torch.roll(theta, -1, dims=-1))
    flux = Flo + alpha[..., None] * dF
    return flux, alpha, Flo, Ulo


def global_entropy_limiter(
    U: torch.Tensor,
    safe_flux: torch.Tensor,
    low_flux: torch.Tensor,
    lam: float,
    n_bisect: int = 44,
    safety: float = 1e-5,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Hard-enforce one-step global entropy non-increase on a periodic domain.

    We further blend
        F(beta) = F_L + beta (F_safe - F_L), beta in [0,1].
    Since the update is affine in beta and the physical entropy is convex for
    h>0, the sublevel set {beta : E(U^{n+1}(beta)) <= E(U^n)} is an interval
    containing beta=0 whenever the Rusanov step is entropy dissipative.  A
    per-trajectory bisection selects the largest feasible beta to numerical
    precision.  The tiny ``safety`` offset protects against float32 roundoff.
    """
    Ulo = fv_step_lambda(U, low_flux, lam)
    Usafe = fv_step_lambda(U, safe_flux, lam)
    E0 = total_entropy(U)
    Elo = total_entropy(Ulo)
    E1 = total_entropy(Usafe)

    # The global limiter needs a feasible low-order endpoint.  Under the CFL
    # used by ``advance_snapshot_safe`` the Rusanov step is expected to be
    # entropy dissipative; fail loudly rather than silently claiming a
    # fully-discrete guarantee if that premise is violated.
    if torch.any(Elo > E0 + 1e-10):
        raise RuntimeError("Low-order Rusanov step increased total entropy; reduce CFL.")

    beta = torch.ones(U.shape[0], dtype=U.dtype, device=U.device)
    need = E1 > E0
    if need.any():
        count = int(need.sum())
        lo = torch.zeros(count, dtype=torch.float64, device=U.device)
        hi = torch.ones(count, dtype=torch.float64, device=U.device)
        Ulo_n = Ulo[need].double()
        dU = (Usafe - Ulo)[need].double()
        E0_n = E0[need]

        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            Emid = total_entropy(Ulo_n + mid[:, None, None] * dU)
            ok = Emid <= E0_n
            lo = torch.where(ok, mid, lo)
            hi = torch.where(ok, hi, mid)

        beta[need] = torch.clamp(lo - safety, min=0.0).float()

    flux = low_flux + beta[:, None, None] * (safe_flux - low_flux)
    return flux, beta


@torch.no_grad()
def advance_snapshot_safe(
    model: base.Solver,
    U: torch.Tensor,
    cfl: float = 0.45,
    h_floor: float = 1e-4,
    enforce_global_entropy: bool = True,
):
    """Advance one saved-data interval using adaptive CFL substeps."""
    remaining = DT_SNAPSHOT
    positivity_activity = []
    entropy_activity = []
    nsub = 0

    while remaining > 1e-14:
        h = U[..., 0].clamp_min(1e-10)
        u = U[..., 1] / h
        max_speed = float((torch.abs(u) + torch.sqrt(base.G * h)).max())
        dt = min(remaining, cfl * DX / max(max_speed, 1e-12))
        lam = dt / DX

        # The learned model already applies the Tadmor hard projection.
        Fent = model.flux(U)
        Fpos, alpha, Flo, Ulo = positivity_limited_flux(U, Fent, lam, h_floor)

        # The low-order premise must hold for the positivity proof.
        if float(Ulo[..., 0].min()) < -1e-7:
            raise RuntimeError("Rusanov positivity premise failed; reduce CFL.")

        if enforce_global_entropy:
            Ffinal, beta = global_entropy_limiter(U, Fpos, Flo, lam)
        else:
            Ffinal = Fpos
            beta = torch.ones(U.shape[0], dtype=U.dtype, device=U.device)

        U = fv_step_lambda(U, Ffinal, lam)
        positivity_activity.append(float((alpha < 1.0 - 1e-8).float().mean()))
        entropy_activity.append(float((beta < 1.0 - 1e-7).float().mean()))
        remaining -= dt
        nsub += 1

    return U, float(np.mean(positivity_activity)), float(np.mean(entropy_activity)), nsub


def generate_near_dry_trajectory(ntraj: int, seed: int):
    """Fine-grid reference trajectories well outside the training depth range."""
    rng = np.random.default_rng(seed)
    N = base.NREF
    x = (np.arange(N) + 0.5) / N
    U = np.zeros((ntraj, N, 2), dtype=np.float64)

    for j in range(ntraj):
        typ = rng.choice(["dam", "piecewise", "smooth"], p=[0.55, 0.30, 0.15])
        if typ == "dam":
            if rng.random() < 0.5:
                hL, hR = rng.uniform(1.0, 2.5), rng.uniform(0.05, 0.20)
            else:
                hR, hL = rng.uniform(1.0, 2.5), rng.uniform(0.05, 0.20)
            uL, uR = rng.uniform(-0.8, 0.8, 2)
            cut = int(rng.integers(N // 3, 2 * N // 3))
            h = np.r_[np.full(cut, hL), np.full(N - cut, hR)]
            u = np.r_[np.full(cut, uL), np.full(N - cut, uR)]
        elif typ == "piecewise":
            K = int(rng.integers(3, 6))
            cuts = sorted(rng.choice(np.arange(1, N), K - 1, replace=False))
            h, u = np.empty(N), np.empty(N)
            start = 0
            for k in range(K):
                end = cuts[k] if k < K - 1 else N
                h[start:end] = rng.uniform(0.05, 2.5)
                u[start:end] = rng.uniform(-1.3, 1.3)
                start = end
        else:
            h0, amp = rng.uniform(0.25, 0.60), rng.uniform(0.10, 0.20)
            h = np.maximum(
                h0 + amp * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi)),
                0.05,
            )
            u = rng.uniform(0.5, 1.3) * np.sin(
                2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi)
            )

        U[j, :, 0] = h
        U[j, :, 1] = h * u

    dx = 1.0 / base.NREF

    def restrict(V):
        return V.reshape(ntraj, base.NCOARSE, base.FACTOR, 2).mean(2).astype(np.float32)

    snapshots = [restrict(U)]
    for _ in range(base.NSNAP - 1):
        remaining = DT_SNAPSHOT
        while remaining > 1e-14:
            h = U[..., 0]
            u = U[..., 1] / h
            max_speed = np.max(np.abs(u) + np.sqrt(base.G * h))
            dt = min(remaining, 0.28 * dx / max_speed)
            U = base.np_ssprk2(U, dt, dx)
            if U[..., 0].min() <= 0:
                raise RuntimeError("Fine-grid near-dry reference lost positivity.")
            remaining -= dt
        snapshots.append(restrict(U))

    return np.stack(snapshots, axis=1)
