"""Self-contained 1D Euler ML-ablation runner for HCFL.

This file intentionally has no dependency on the older experiment scripts.
It reproduces the current broad-data comparison between:
  direct          : HLLC + learned raw flux correction
  invariant       : HLLC + dimensionless/Galilean-aware raw correction
  characteristic  : HLLC + Roe-basis correction
  dissipation     : HLLC + learned characteristic-wise viscosity correction
  conv            : HLLC + raw correction with a Conv1d backbone

Every learned proposal is passed through the same hard Tadmor projection.

Examples
--------
python euler_ablation_runner.py --model direct --seed 0 --iters 1100
python euler_ablation_runner.py --model invariant --seed 0 --iters 1100
python euler_ablation_runner.py --model characteristic --seed 0 --iters 1100
python euler_ablation_runner.py --model dissipation --seed 0 --iters 1100
python euler_ablation_runner.py --model conv --seed 0 --iters 1100
"""

import argparse
import json
from pathlib import Path
import math
import numpy as np
import torch
from torch import nn
import torch.nn.functional as Fnn

torch.set_num_threads(2)

GAMMA = 1.4
NREF = 512
NCOARSE = 64
FACTOR = NREF // NCOARSE
DT_SNAPSHOT = 4e-4
LAMBDA = DT_SNAPSHOT / (1.0 / NCOARSE)
NSNAP = 16


# ---------------------------------------------------------------------
# Reference-data generator used in the exploratory experiments
# ---------------------------------------------------------------------
def np_pressure(U):
    rho = U[..., 0]
    m = U[..., 1]
    E = U[..., 2]
    return (GAMMA - 1.0) * (E - 0.5 * m * m / rho)


def np_flux(U):
    rho = U[..., 0]
    m = U[..., 1]
    E = U[..., 2]
    u = m / rho
    p = np_pressure(U)
    return np.stack([m, m * u + p, u * (E + p)], axis=-1)


def np_rusanov_pair(UL, UR):
    fL, fR = np_flux(UL), np_flux(UR)
    rhoL, rhoR = UL[..., 0], UR[..., 0]
    uL, uR = UL[..., 1] / rhoL, UR[..., 1] / rhoR
    pL, pR = np_pressure(UL), np_pressure(UR)
    cL = np.sqrt(GAMMA * pL / rhoL)
    cR = np.sqrt(GAMMA * pR / rhoR)
    a = np.maximum(np.abs(uL) + cL, np.abs(uR) + cR)
    return 0.5 * (fL + fR) - 0.5 * a[..., None] * (UR - UL)


def np_rhs(U, dx):
    Fh = np_rusanov_pair(U, np.roll(U, -1, axis=-2))
    return -(Fh - np.roll(Fh, 1, axis=-2)) / dx


def np_ssprk2(U, dt, dx):
    U1 = U + dt * np_rhs(U, dx)
    U2 = U1 + dt * np_rhs(U1, dx)
    return 0.5 * U + 0.5 * U2


def prim_to_cons(rho, u, p):
    return np.stack(
        [rho, rho * u, p / (GAMMA - 1.0) + 0.5 * rho * u * u],
        axis=-1,
    )


def generate_ic(ntraj, N, seed, ood=False):
    rng = np.random.default_rng(seed)
    x = (np.arange(N) + 0.5) / N
    U = np.zeros((ntraj, N, 3), dtype=np.float64)

    for j in range(ntraj):
        typ = rng.choice(["riemann", "piecewise", "smooth"], p=[0.45, 0.30, 0.25])

        if typ == "riemann":
            if ood:
                rhoL, rhoR = rng.uniform(0.25, 2.8, 2)
                uL, uR = rng.uniform(-1.3, 1.3, 2)
                pL, pR = rng.uniform(0.25, 3.0, 2)
            else:
                rhoL, rhoR = rng.uniform(0.65, 1.8, 2)
                uL, uR = rng.uniform(-0.65, 0.65, 2)
                pL, pR = rng.uniform(0.65, 1.8, 2)

            cut = int(rng.integers(N // 4, 3 * N // 4))
            rho = np.r_[np.full(cut, rhoL), np.full(N - cut, rhoR)]
            u = np.r_[np.full(cut, uL), np.full(N - cut, uR)]
            p = np.r_[np.full(cut, pL), np.full(N - cut, pR)]
            shift = int(rng.integers(0, N))
            rho, u, p = np.roll(rho, shift), np.roll(u, shift), np.roll(p, shift)

        elif typ == "piecewise":
            K = int(rng.integers(3, 6))
            cuts = sorted(rng.choice(np.arange(1, N), K - 1, replace=False))
            rho, u, p = np.empty(N), np.empty(N), np.empty(N)
            start = 0
            for k in range(K):
                end = cuts[k] if k < K - 1 else N
                if ood:
                    rhov = rng.uniform(0.25, 2.8)
                    uv = rng.uniform(-1.3, 1.3)
                    pv = rng.uniform(0.25, 3.0)
                else:
                    rhov = rng.uniform(0.65, 1.8)
                    uv = rng.uniform(-0.65, 0.65)
                    pv = rng.uniform(0.65, 1.8)
                rho[start:end], u[start:end], p[start:end] = rhov, uv, pv
                start = end

        else:
            if ood:
                rho0, p0 = rng.uniform(0.7, 1.8), rng.uniform(0.7, 1.8)
                ar, ap, au = rng.uniform(0.10, 0.45), rng.uniform(0.10, 0.45), rng.uniform(0.2, 1.0)
            else:
                rho0, p0 = rng.uniform(0.85, 1.4), rng.uniform(0.85, 1.4)
                ar, ap, au = rng.uniform(0.05, 0.25), rng.uniform(0.05, 0.25), rng.uniform(0.05, 0.5)

            rho = rho0 + ar * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi))
            p = p0 + ap * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi))
            u = au * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi))
            rho = np.maximum(rho, 0.2)
            p = np.maximum(p, 0.2)

        U[j] = prim_to_cons(rho, u, p)

    return U


def generate_extreme_ic(ntraj, N, seed):
    rng = np.random.default_rng(seed)
    x = (np.arange(N) + 0.5) / N
    U = np.zeros((ntraj, N, 3), dtype=np.float64)

    for j in range(ntraj):
        typ = rng.choice(["riemann", "piecewise", "smooth"], p=[0.55, 0.30, 0.15])

        if typ == "riemann":
            rhoL, rhoR = rng.uniform(0.10, 3.2, 2)
            pL, pR = rng.uniform(0.08, 5.0, 2)
            uL, uR = rng.uniform(-2.2, 2.2, 2)
            cut = int(rng.integers(N // 4, 3 * N // 4))
            rho = np.r_[np.full(cut, rhoL), np.full(N - cut, rhoR)]
            u = np.r_[np.full(cut, uL), np.full(N - cut, uR)]
            p = np.r_[np.full(cut, pL), np.full(N - cut, pR)]
            shift = int(rng.integers(0, N))
            rho, u, p = np.roll(rho, shift), np.roll(u, shift), np.roll(p, shift)

        elif typ == "piecewise":
            K = int(rng.integers(3, 6))
            cuts = sorted(rng.choice(np.arange(1, N), K - 1, replace=False))
            rho, u, p = np.empty(N), np.empty(N), np.empty(N)
            start = 0
            for k in range(K):
                end = cuts[k] if k < K - 1 else N
                rho[start:end] = rng.uniform(0.10, 3.2)
                u[start:end] = rng.uniform(-2.2, 2.2)
                p[start:end] = rng.uniform(0.08, 5.0)
                start = end

        else:
            rho0, p0 = rng.uniform(0.3, 2.0), rng.uniform(0.3, 2.0)
            rho = np.maximum(
                rho0 + rng.uniform(0.05, 0.25) * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi)),
                0.10,
            )
            p = np.maximum(
                p0 + rng.uniform(0.05, 0.25) * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi)),
                0.08,
            )
            u = rng.uniform(0.5, 2.0) * np.sin(2 * np.pi * rng.integers(1, 4) * x + rng.uniform(0, 2 * np.pi))

        U[j] = prim_to_cons(rho, u, p)

    return U


def _rollout_reference(U):
    ntraj = U.shape[0]
    dx = 1.0 / NREF

    def restrict(V):
        return V.reshape(ntraj, NCOARSE, FACTOR, 3).mean(axis=2).astype(np.float32)

    snaps = [restrict(U)]

    for _ in range(NSNAP - 1):
        rem = DT_SNAPSHOT
        while rem > 1e-14:
            rho = U[..., 0]
            p = np_pressure(U)
            u = U[..., 1] / rho
            c = np.sqrt(GAMMA * p / rho)
            max_speed = np.max(np.abs(u) + c)
            dt = min(rem, 0.25 * dx / max_speed)
            U = np_ssprk2(U, dt, dx)

            # Conservative floor repair for rare pathological random states.
            rho = U[..., 0]
            m = U[..., 1]
            kinetic = 0.5 * m * m / np.maximum(rho, 1e-12)
            bad = (rho <= 0) | (np_pressure(U) <= 0)
            if np.any(bad):
                U[..., 0] = np.maximum(U[..., 0], 0.05)
                U[..., 2] = np.maximum(U[..., 2], kinetic + 0.08 / (GAMMA - 1.0))

            rem -= dt

        snaps.append(restrict(U))

    return np.stack(snaps, axis=1)


def generate_trajectory(ntraj, seed, ood=False):
    return _rollout_reference(generate_ic(ntraj, NREF, seed, ood=ood))


def generate_extreme_trajectory(ntraj, seed):
    return _rollout_reference(generate_extreme_ic(ntraj, NREF, seed))


# ---------------------------------------------------------------------
# Torch Euler physics + hard entropy projection
# ---------------------------------------------------------------------
def t_pressure(U):
    rho = U[..., 0].clamp_min(1e-8)
    m = U[..., 1]
    E = U[..., 2]
    return (GAMMA - 1.0) * (E - 0.5 * m * m / rho)


def primitive(U):
    rho = U[..., 0].clamp_min(1e-8)
    return torch.stack([rho, U[..., 1] / rho, t_pressure(U)], dim=-1)


def t_flux(U):
    rho = U[..., 0].clamp_min(1e-8)
    m = U[..., 1]
    E = U[..., 2]
    u = m / rho
    p = t_pressure(U)
    return torch.stack([m, m * u + p, u * (E + p)], dim=-1)


def entropy(U):
    rho = U[..., 0].clamp_min(1e-12)
    p = t_pressure(U).clamp_min(1e-12)
    s = torch.log(p) - GAMMA * torch.log(rho)
    return -rho * s / (GAMMA - 1.0)


def entropy_variables(U):
    rho = U[..., 0].clamp_min(1e-12)
    m = U[..., 1]
    u = m / rho
    p = t_pressure(U).clamp_min(1e-12)
    s = torch.log(p) - GAMMA * torch.log(rho)
    v1 = (GAMMA - s) / (GAMMA - 1.0) - rho * u * u / (2.0 * p)
    v2 = rho * u / p
    v3 = -rho / p
    return torch.stack([v1, v2, v3], dim=-1)


def entropy_potential(U):
    return U[..., 1]


def t_rusanov(U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    fL, fR = t_flux(UL), t_flux(UR)
    rhoL, rhoR = UL[..., 0].clamp_min(1e-8), UR[..., 0].clamp_min(1e-8)
    uL, uR = UL[..., 1] / rhoL, UR[..., 1] / rhoR
    pL, pR = t_pressure(UL).clamp_min(1e-8), t_pressure(UR).clamp_min(1e-8)
    cL = torch.sqrt(GAMMA * pL / rhoL)
    cR = torch.sqrt(GAMMA * pR / rhoR)
    a = torch.maximum(torch.abs(uL) + cL, torch.abs(uR) + cR)
    return 0.5 * (fL + fR) - 0.5 * a[..., None] * (UR - UL)


def entropy_residual(Fh, U):
    UR = torch.roll(U, -1, dims=-2)
    a = entropy_variables(UR) - entropy_variables(U)
    b = entropy_potential(UR) - entropy_potential(U)
    return (a * Fh).sum(dim=-1) - b


def hard_entropy_projection(Fh, U):
    UR = torch.roll(U, -1, dims=-2)
    a = entropy_variables(UR) - entropy_variables(U)
    b = entropy_potential(UR) - entropy_potential(U)
    r = (a * Fh).sum(dim=-1) - b
    n2 = (a * a).sum(dim=-1)
    alpha = torch.zeros_like(r)
    mask = (r > 0) & (n2 > 1e-14)
    alpha[mask] = r[mask] / n2[mask]
    return Fh - alpha[..., None] * a


def fv_step(U, Fh):
    return U - LAMBDA * (Fh - torch.roll(Fh, 1, dims=-2))


# ---------------------------------------------------------------------
# HLLC base
# ---------------------------------------------------------------------
def hllc_pair(UL, UR):
    rhoL, rhoR = UL[..., 0].clamp_min(1e-10), UR[..., 0].clamp_min(1e-10)
    mL, mR = UL[..., 1], UR[..., 1]
    EL, ER = UL[..., 2], UR[..., 2]
    uL, uR = mL / rhoL, mR / rhoR
    pL, pR = t_pressure(UL).clamp_min(1e-10), t_pressure(UR).clamp_min(1e-10)
    cL, cR = torch.sqrt(GAMMA * pL / rhoL), torch.sqrt(GAMMA * pR / rhoR)

    SL = torch.minimum(uL - cL, uR - cR)
    SR = torch.maximum(uL + cL, uR + cR)
    den = rhoL * (SL - uL) - rhoR * (SR - uR)
    SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (den + 1e-14)

    rhoSL = rhoL * (SL - uL) / (SL - SM + 1e-14)
    rhoSR = rhoR * (SR - uR) / (SR - SM + 1e-14)
    ESL = rhoSL * (EL / rhoL + (SM - uL) * (SM + pL / (rhoL * (SL - uL) + 1e-14)))
    ESR = rhoSR * (ER / rhoR + (SM - uR) * (SM + pR / (rhoR * (SR - uR) + 1e-14)))

    USL = torch.stack([rhoSL, rhoSL * SM, ESL], dim=-1)
    USR = torch.stack([rhoSR, rhoSR * SM, ESR], dim=-1)

    fL, fR = t_flux(UL), t_flux(UR)
    FSL = fL + SL[..., None] * (USL - UL)
    FSR = fR + SR[..., None] * (USR - UR)

    return torch.where(
        (SL >= 0)[..., None],
        fL,
        torch.where((SM >= 0)[..., None], FSL, torch.where((SR > 0)[..., None], FSR, fR)),
    )


def t_hllc(U):
    return hllc_pair(U, torch.roll(U, -1, dims=-2))


# ---------------------------------------------------------------------
# Model variants
# ---------------------------------------------------------------------
class FullFlux(nn.Module):
    """Predict the complete interface flux from a five-cell primitive stencil.

    Unlike the correction architectures below, this model has no HLLC base,
    jump gate, or analytic output scale.  The shared ``Solver`` still applies
    the same hard Tadmor half-space projection to the network output.
    """

    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([
            (torch.roll(P, s, dims=-2) - self.mean) / self.std
            for s in [2, 1, 0, -1, -2]
        ], dim=-1)
        return self.net(feats)


class CentralConsistentFlux(nn.Module):
    """Central physical flux plus an exactly jump-gated learned correction.

    The correction is identically zero whenever the two interface states are
    equal, even if the outer stencil is not constant.  Consequently the raw
    proposal satisfies F_hat(U, U) = F(U) by construction.
    """

    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))
        self.register_buffer(
            "scale", torch.tensor([0.6, 1.2, 2.5], dtype=torch.float32)
        )

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([
            (torch.roll(P, s, dims=-2) - self.mean) / self.std
            for s in [2, 1, 0, -1, -2]
        ], dim=-1)
        correction = torch.tanh(self.net(feats))
        right_state = torch.roll(U, -1, dims=-2)
        right_primitive = torch.roll(P, -1, dims=-2)
        jump = torch.linalg.vector_norm(
            (right_primitive - P) / self.std, dim=-1
        )
        central_flux = 0.5 * (t_flux(U) + t_flux(right_state))
        return (
            central_flux
            + 0.18 * jump[..., None] * correction * self.scale
        )


class DirectFlux(nn.Module):
    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))
        self.register_buffer("scale", torch.tensor([0.6, 1.2, 2.5], dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([(torch.roll(P, s, dims=-2) - self.mean) / self.std for s in [2,1,0,-1,-2]], dim=-1)
        corr = torch.tanh(self.net(feats))
        PR = torch.roll(P, -1, dims=-2)
        jump = torch.sqrt((((PR - P) / self.std) ** 2).sum(dim=-1) + 1e-12)
        return t_hllc(U) + 0.18 * jump[..., None] * corr * self.scale


def invariant_features(U):
    P = primitive(U)
    rho = P[..., 0].clamp_min(1e-8)
    vel = P[..., 1]
    pres = P[..., 2].clamp_min(1e-8)
    PR = torch.roll(P, -1, dims=-2)
    rhoR, velR, presR = PR[..., 0].clamp_min(1e-8), PR[..., 1], PR[..., 2].clamp_min(1e-8)

    rho_ref = torch.sqrt(rho * rhoR)
    p_ref = torch.sqrt(pres * presR)
    u_ref = 0.5 * (vel + velR)
    c_ref = torch.sqrt(GAMMA * p_ref / rho_ref).clamp_min(1e-8)

    feats = []
    for s in [2,1,0,-1,-2]:
        Ps = torch.roll(P, s, dims=-2)
        feats += [
            torch.log(Ps[...,0].clamp_min(1e-8) / rho_ref),
            (Ps[...,1] - u_ref) / c_ref,
            torch.log(Ps[...,2].clamp_min(1e-8) / p_ref),
        ]
    z = torch.stack(feats, dim=-1)

    zL = torch.stack([torch.log(rho/rho_ref), (vel-u_ref)/c_ref, torch.log(pres/p_ref)], dim=-1)
    zR = torch.stack([torch.log(rhoR/rho_ref), (velR-u_ref)/c_ref, torch.log(presR/p_ref)], dim=-1)
    jump = torch.sqrt(((zR-zL)**2).sum(dim=-1) + 1e-12)
    flux_scale = torch.stack([rho_ref*c_ref, p_ref, p_ref*c_ref], dim=-1)
    return z, jump, flux_scale


class InvariantFlux(nn.Module):
    def __init__(self, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)

    def forward(self, U):
        z, jump, flux_scale = invariant_features(U)
        corr = torch.tanh(self.net(z))
        return t_hllc(U) + 0.18 * jump[..., None] * corr * flux_scale


def roe_basis(U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    rhoL, rhoR = UL[...,0].clamp_min(1e-8), UR[...,0].clamp_min(1e-8)
    uL, uR = UL[...,1]/rhoL, UR[...,1]/rhoR
    pL, pR = t_pressure(UL).clamp_min(1e-8), t_pressure(UR).clamp_min(1e-8)
    HL, HR = (UL[...,2]+pL)/rhoL, (UR[...,2]+pR)/rhoR

    srL, srR = torch.sqrt(rhoL), torch.sqrt(rhoR)
    den = (srL + srR).clamp_min(1e-8)
    u = (srL*uL + srR*uR) / den
    H = (srL*HL + srR*HR) / den
    c = torch.sqrt(((GAMMA-1.0)*(H-0.5*u*u)).clamp_min(1e-8))

    r1 = torch.stack([torch.ones_like(u), u-c, H-u*c], dim=-1)
    r2 = torch.stack([torch.ones_like(u), u, 0.5*u*u], dim=-1)
    r3 = torch.stack([torch.ones_like(u), u+c, H+u*c], dim=-1)
    R = torch.stack([r1,r2,r3], dim=-1)
    return R, u, c


def entropy_fixed_roe_waves(U):
    """Return Roe eigenvectors, jumps, and sonic-fixed wave speeds."""
    R, u, c = roe_basis(U)
    dU = torch.roll(U, -1, dims=-2) - U
    alpha = torch.linalg.solve(R, dU.unsqueeze(-1)).squeeze(-1)
    acoustic_delta = (0.1 * c).clamp_min(1e-6)

    def acoustic_abs(lam):
        magnitude = torch.abs(lam)
        return torch.where(
            magnitude >= acoustic_delta,
            magnitude,
            0.5 * (lam * lam / acoustic_delta + acoustic_delta),
        )

    speeds = torch.stack([
        acoustic_abs(u - c),
        torch.abs(u),
        acoustic_abs(u + c),
    ], dim=-1)
    return R, alpha, speeds


class RoeCompleteFlux(nn.Module):
    """Predict the complete flux in a local normalized Roe eigenbasis."""

    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([
            (torch.roll(P, s, dims=-2) - self.mean) / self.std
            for s in [2, 1, 0, -1, -2]
        ], dim=-1)
        coefficients = self.net(feats)
        R, _, _ = roe_basis(U)
        normalized_R = R / torch.linalg.vector_norm(
            R, dim=-2
        ).clamp_min(1e-8)[..., None, :]
        right_state = torch.roll(U, -1, dims=-2)
        left_flux = t_flux(U)
        right_flux = t_flux(right_state)
        flux_scale = torch.sqrt(
            0.5 * (left_flux.square() + right_flux.square()).sum(dim=-1)
        ).clamp_min(1e-4)
        return torch.einsum(
            "...ij,...j->...i",
            normalized_R,
            flux_scale[..., None] * coefficients,
        )


class CentralRoeSignedFlux(nn.Module):
    """Central flux with learned Roe multipliers that may reverse a wave."""

    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([
            (torch.roll(P, s, dims=-2) - self.mean) / self.std
            for s in [2, 1, 0, -1, -2]
        ], dim=-1)
        multipliers = 1.0 + 2.0 * torch.tanh(self.net(feats))
        R, alpha, speeds = entropy_fixed_roe_waves(U)
        dissipation = torch.einsum(
            "...ij,...j->...i", R, multipliers * speeds * alpha
        )
        right_state = torch.roll(U, -1, dims=-2)
        central_flux = 0.5 * (t_flux(U) + t_flux(right_state))
        return central_flux - 0.5 * dissipation


class CentralRoeUpwindFlux(nn.Module):
    """Central flux with bounded nonnegative automatic-upwind multipliers."""

    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([
            (torch.roll(P, s, dims=-2) - self.mean) / self.std
            for s in [2, 1, 0, -1, -2]
        ], dim=-1)
        multipliers = 1.0 + torch.tanh(self.net(feats))
        R, alpha, speeds = entropy_fixed_roe_waves(U)
        dissipation = torch.einsum(
            "...ij,...j->...i", R, multipliers * speeds * alpha
        )
        right_state = torch.roll(U, -1, dims=-2)
        central_flux = 0.5 * (t_flux(U) + t_flux(right_state))
        return central_flux - 0.5 * dissipation


class CharacteristicFlux(nn.Module):
    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([(torch.roll(P,s,dims=-2)-self.mean)/self.std for s in [2,1,0,-1,-2]], dim=-1)
        delta = torch.tanh(self.net(feats))

        R, _, c = roe_basis(U)
        rho_ref = torch.sqrt(U[...,0].clamp_min(1e-8) * torch.roll(U[...,0], -1, dims=-1).clamp_min(1e-8))
        D = torch.stack([rho_ref*c, rho_ref*c*c, rho_ref*c*c*c], dim=-1).clamp_min(1e-8)
        Rn = R / D[..., :, None]
        Rn = Rn / torch.linalg.vector_norm(Rn, dim=-2).clamp_min(1e-8)[..., None, :]
        B = D[..., :, None] * Rn

        PR = torch.roll(P, -1, dims=-2)
        jump = torch.sqrt((((PR-P)/self.std)**2).sum(dim=-1) + 1e-12)
        corr = torch.einsum("...ij,...j->...i", B, delta)
        return t_hllc(U) + 0.18 * jump[..., None] * corr


class DissipationFlux(nn.Module):
    def __init__(self, mean, std, width=72):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(15, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 3),
        )
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat([(torch.roll(P,s,dims=-2)-self.mean)/self.std for s in [2,1,0,-1,-2]], dim=-1)
        d = torch.tanh(self.net(feats))

        R, u, c = roe_basis(U)
        dU = torch.roll(U, -1, dims=-2) - U
        alpha = torch.linalg.solve(R, dU.unsqueeze(-1)).squeeze(-1)
        lam = torch.stack([torch.abs(u-c), torch.abs(u), torch.abs(u+c)], dim=-1)
        corr = -0.5 * torch.einsum("...ij,...j->...i", R, d * lam * alpha)
        return t_hllc(U) + corr


class ConvFlux(nn.Module):
    def __init__(self, mean, std, width=72):
        super().__init__()
        self.conv5 = nn.Conv1d(3, width, 5)
        self.conv1 = nn.Conv1d(width, width, 1)
        self.out = nn.Conv1d(width, 3, 1)
        nn.init.zeros_(self.out.weight)
        nn.init.zeros_(self.out.bias)
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32).view(1,3,1))
        self.register_buffer("std", torch.tensor(std, dtype=torch.float32).view(1,3,1))
        self.register_buffer("scale", torch.tensor([0.6,1.2,2.5], dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        x = P.permute(0,2,1)
        x = (x-self.mean)/self.std
        x = Fnn.pad(x, (2,2), mode="circular")
        x = torch.tanh(self.conv5(x))
        x = torch.tanh(self.conv1(x))
        corr = torch.tanh(self.out(x)).permute(0,2,1)
        PR = torch.roll(P, -1, dims=-2)
        stdv = self.std.view(1,1,3)
        jump = torch.sqrt((((PR-P)/stdv)**2).sum(dim=-1)+1e-12)
        return t_hllc(U) + 0.18 * jump[...,None] * corr * self.scale


class Solver(nn.Module):
    def __init__(self, model_name, mean, std, width=72):
        super().__init__()
        if model_name == "full":
            self.flux_net = FullFlux(mean, std, width)
        elif model_name == "roe_complete":
            self.flux_net = RoeCompleteFlux(mean, std, width)
        elif model_name == "central_consistent":
            self.flux_net = CentralConsistentFlux(mean, std, width)
        elif model_name == "central_roe_signed":
            self.flux_net = CentralRoeSignedFlux(mean, std, width)
        elif model_name == "central_roe_upwind":
            self.flux_net = CentralRoeUpwindFlux(mean, std, width)
        elif model_name == "direct":
            self.flux_net = DirectFlux(mean, std, width)
        elif model_name == "invariant":
            self.flux_net = InvariantFlux(width)
        elif model_name == "characteristic":
            self.flux_net = CharacteristicFlux(mean, std, width)
        elif model_name == "dissipation":
            self.flux_net = DissipationFlux(mean, std, width)
        elif model_name == "conv":
            self.flux_net = ConvFlux(mean, std, width)
        else:
            raise ValueError(model_name)

    def flux(self, U):
        return hard_entropy_projection(self.flux_net(U), U)

    def one_step(self, U):
        return fv_step(U, self.flux(U))


def make_broad_training_data(seed):
    idata = torch.tensor(generate_trajectory(220, 6000 + seed, ood=False))
    odata = torch.tensor(generate_trajectory(260, 7000 + seed, ood=True))
    xdata = torch.tensor(generate_extreme_trajectory(100, 8000 + seed))
    return torch.cat([idata, odata, xdata], dim=0)


def evaluate_moderate(model, seed, state_std):
    rows = {}
    with torch.no_grad():
        for split, arr in [
            ("ID", generate_trajectory(90, 2000 + seed, ood=False)),
            ("OOD", generate_trajectory(90, 3000 + seed, ood=True)),
        ]:
            data = torch.tensor(arr)
            U = data[:,0].clone()
            errs = []
            violation = 0
            total = 0
            for t in range(1, data.shape[1]):
                Fh = model.flux(U)
                r = entropy_residual(Fh, U)
                violation += int((r > 1e-5).sum())
                total += r.numel()
                U = fv_step(U, Fh)
                errs.append(float((((U-data[:,t])/state_std)**2).mean()))
            rows[split] = {
                "nrmse": float(np.sqrt(np.mean(errs))),
                "entropy_violation_rate": violation / max(total, 1),
            }
    return rows


def train(args):
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    train_data = make_broad_training_data(args.seed)
    P = primitive(train_data)
    mean = P.mean(dim=(0,1,2)).numpy()
    std = P.std(dim=(0,1,2)).numpy()
    state_std = torch.tensor([float(train_data[...,j].std()) for j in range(3)])

    torch.manual_seed(12000 + args.seed)
    model = Solver(args.model, mean, std, width=args.width)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    gen = torch.Generator().manual_seed(13000 + args.seed)

    for it in range(args.iters):
        bs = args.batch_size
        inds = torch.randint(0, train_data.shape[0], (bs,), generator=gen)
        ts = torch.randint(0, train_data.shape[1]-1, (bs,), generator=gen)
        U = torch.stack([train_data[i,t] for i,t in zip(inds,ts)], dim=0)
        target = torch.stack([train_data[i,t+1] for i,t in zip(inds,ts)], dim=0)

        pred = model.one_step(U)
        loss = (((pred-target)/state_std)**2).mean()

        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()

        if args.log_every and (it % args.log_every == 0 or it == args.iters-1):
            print(json.dumps({"iter": it, "loss": float(loss.detach())}))

    model.eval()
    stem = f"{args.model}_seed{args.seed}"
    torch.save(model.state_dict(), out / f"{stem}.pt")
    np.savez(out / f"{stem}_stats.npz", mean=mean, std=std, state_std=state_std.numpy())

    metrics = evaluate_moderate(model, args.seed, state_std)
    (out / f"{stem}_metrics.json").write_text(json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--model", choices=["full","roe_complete","central_consistent","central_roe_signed","central_roe_upwind","direct","invariant","characteristic","dissipation","conv"], required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--iters", type=int, default=1100)
    p.add_argument("--width", type=int, default=72)
    p.add_argument("--batch-size", type=int, default=56)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--outdir", default="outputs/euler_ablation")
    p.add_argument("--log-every", type=int, default=100)
    train(p.parse_args())
