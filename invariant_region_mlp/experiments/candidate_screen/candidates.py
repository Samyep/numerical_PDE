"""Candidate PDEs for the next benchmark round (screening only).

All candidates reuse the mechanism-suite recursion (`MechanismMLP`) unchanged and
the state convention z = sigma * grad(u), sigma = sqrt(2), mu = 0.

C1  l1-control HJB      u_t + Lap u - lam ||grad u||_1 = 0,   lam = 1/||w||_1
    Ridge solution identical to P4 (beta=2, lambda_f=1, T=0.5); only the
    generator off the ridge differs (noise aggregated in l1 instead of l2).
C2  zero-sum LQ game     u_t + Lap u - a |grad_A u|^2 + b |grad_B u|^2 = 0
    (Isaacs; controller on block A, adversary on block B; nonconvex in z).
    w carries mass rho = 1/2 on each block and a - b = 2, so along the ridge the
    PDE is exactly P1's (kappa = 1): same closed-form solution, same certificate.
    One-step theory predicts the orthogonal Jensen bias
        (c/2M) * [ -a (d_A - 1/2) + b (d_B - 1/2) ],
    i.e. negative (convex), ~0 (cancelling) or positive (flipped) by design.
C3  cubic viscous HJ     u_t + Lap u - kappa |grad u|^3 = 0, g = A log cosh(beta s)/beta
    One-dimensional reference by a monotone (Godunov) scheme with Richardson
    extrapolation; certificate |psi_s| <= A from Lipschitz preservation.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline

HERE = Path(__file__).resolve().parent
SUITE = HERE.parent / "mechanism_suite"
for p in (str(SUITE), str(SUITE.parent)):
    if p not in sys.path:
        sys.path.insert(0, p)

from mechanism_suite.equations import RidgeLSEHJB, ExactEquation, _as_points  # noqa: E402
from mechanism_suite.norm_hjb import NormDriverHJB  # noqa: E402

REF_P4 = (HERE.parent.parent / "results" / "mechanism_suite" / "reference_cache"
          / "norm_hjb_beta2_lambda1_T0p5.npz")
SQ2 = math.sqrt(2.0)


# --------------------------------------------------------------------------- C1
@dataclass(frozen=True)
class L1ControlHJB(NormDriverHJB):
    reference_path: str = str(REF_P4)
    T: float = 0.5
    name: str = "C1_l1_control_hjb"

    def __post_init__(self) -> None:
        super().__post_init__()
        rng = np.random.default_rng(np.random.SeedSequence([self.direction_seed, self.d, 101]))
        w = rng.standard_normal(self.d); w /= np.linalg.norm(w)
        object.__setattr__(self, "w", w)
        object.__setattr__(self, "lam_l1", 1.0 / float(np.sum(np.abs(w))))

    def generator(self, u, z):
        del u
        return -(self.lam_l1 / self.sigma) * np.sum(np.abs(np.asarray(z, dtype=np.float64)), axis=-1)


# --------------------------------------------------------------------------- C2
@dataclass(frozen=True)
class LQGameHJB(RidgeLSEHJB):
    a: float = 2.0
    b: float = 0.0
    frac_A: float = 0.5           # fraction of coordinates controlled by the minimiser
    name: str = "C2_lq_game"

    def __post_init__(self) -> None:
        super().__post_init__()
        if not math.isclose(0.5 * self.a - 0.5 * self.b, 1.0):
            raise ValueError("need a*rho - b*(1-rho) = 1 with rho = 1/2")
        dA = max(1, int(round(self.frac_A * self.d)))
        if dA >= self.d:
            raise ValueError("block B must be nonempty")
        rng = np.random.default_rng(np.random.SeedSequence([self.direction_seed, self.d, 202, dA]))
        wA = rng.standard_normal(dA); wB = rng.standard_normal(self.d - dA)
        w = np.concatenate([wA / np.linalg.norm(wA), wB / np.linalg.norm(wB)]) / SQ2
        object.__setattr__(self, "w", w)
        object.__setattr__(self, "dA", dA)

    def generator(self, u, z):
        del u
        z = np.asarray(z, dtype=np.float64)
        zA, zB = z[..., :self.dA], z[..., self.dA:]
        return -0.5 * self.a * np.sum(zA * zA, axis=-1) + 0.5 * self.b * np.sum(zB * zB, axis=-1)

    def predicted_orthogonal_bias_factor(self) -> float:
        """Coefficient K in  E f(z_hat) - f(E z_hat) = (c/M) * K  from the one-step theory."""
        return 0.5 * (-self.a * (self.dA - 0.5) + self.b * (self.d - self.dA - 0.5))


# --------------------------------------------------------------------------- C3
def _solve_cubic_1d(kappa, A, beta, T, L, ds):
    """psi_tau + kappa|psi_s|^3 = psi_ss on [-L,L], tau = T - t, Godunov flux for convex H >= 0."""
    s = np.arange(-L, L + ds / 2, ds)
    psi = A * (np.logaddexp(beta * s, -beta * s) - math.log(2.0)) / beta
    pmax = A
    dt = 0.45 / (2.0 / ds**2 + 3 * kappa * pmax**2 / ds)
    nt = int(math.ceil(T / dt)); dt = T / nt
    taus = [0.0]; snaps = [psi.copy()]
    save_every = max(1, nt // 400)
    for k in range(1, nt + 1):
        fwd = np.empty_like(psi); bwd = np.empty_like(psi)
        fwd[:-1] = (psi[1:] - psi[:-1]) / ds; fwd[-1] = A * math.tanh(beta * s[-1])
        bwd[1:] = (psi[1:] - psi[:-1]) / ds; bwd[0] = A * math.tanh(beta * s[0])
        H = kappa * np.maximum(np.maximum(bwd, 0.0), np.maximum(-fwd, 0.0)) ** 3
        lap = np.empty_like(psi)
        lap[1:-1] = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) / ds**2
        lap[0] = lap[1]; lap[-1] = lap[-2]
        psi = psi + dt * (lap - H)
        if k % save_every == 0 or k == nt:
            taus.append(k * dt); snaps.append(psi.copy())
    return s, np.array(taus), np.array(snaps)


@lru_cache(maxsize=4)
def cubic_reference(kappa, A, beta, T, L=6.0, ds=0.01):
    """Richardson-extrapolated (first-order scheme) value table on a common grid."""
    s1, t1, P1 = _solve_cubic_1d(kappa, A, beta, T, L, ds)
    s2, t2, P2 = _solve_cubic_1d(kappa, A, beta, T, L, ds / 2)
    # interpolate the fine solution onto the coarse (tau, s) grid
    fine = RectBivariateSpline(t2, s2, P2, kx=3, ky=3)
    P2c = fine(t1, s1)
    rich = 2.0 * P2c - P1
    spline = RectBivariateSpline(t1, s1, rich, kx=3, ky=3)
    inner = np.abs(s1) <= L - 1.0
    return spline, float(np.max(np.abs(P2c - P1)[:, inner]))


@dataclass(frozen=True)
class CubicHJB(ExactEquation):
    d: int
    kappa: float = 0.25
    A: float = 2.0
    beta: float = 2.0
    T: float = 0.25
    sigma: float = SQ2
    direction_seed: int = 20261006
    name: str = "C3_cubic_viscous_hj"
    family: str = "cubic_hj"

    def __post_init__(self) -> None:
        rng = np.random.default_rng(np.random.SeedSequence([self.direction_seed, self.d, 303]))
        w = rng.standard_normal(self.d); w /= np.linalg.norm(w)
        object.__setattr__(self, "w", w)
        object.__setattr__(self, "mu", 0.0)
        spline, gap = cubic_reference(self.kappa, self.A, self.beta, self.T)
        object.__setattr__(self, "_spline", spline)
        object.__setattr__(self, "reference_grid_gap", gap)

    def _tau_s(self, t, x):
        time, state = _as_points(t, x, self.d)
        s = np.einsum("...d,d->...", state, self.w)
        tau = np.broadcast_to(self.T - time, np.shape(s))
        return tau, s

    def exact_u(self, t, x):
        tau, s = self._tau_s(t, x)
        return self._spline.ev(tau.reshape(-1), s.reshape(-1)).reshape(np.shape(s))

    def exact_z(self, t, x):
        tau, s = self._tau_s(t, x)
        ps = self._spline.ev(tau.reshape(-1), s.reshape(-1), dx=0, dy=1).reshape(np.shape(s))
        return self.sigma * ps[..., None] * self.w

    def terminal(self, x):
        s = np.asarray(x) @ self.w
        return self.A * (np.logaddexp(self.beta * s, -self.beta * s) - math.log(2.0)) / self.beta

    def generator(self, u, z):
        del u
        n = np.linalg.norm(np.asarray(z, dtype=np.float64), axis=-1)
        return -self.kappa * (n / self.sigma) ** 3

    def fzero_reference(self, t, x, quadrature_order: int = 80):
        time, state = _as_points(t, x, self.d)
        nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
        s = state @ self.w; h = np.maximum(self.T - time, 0.0)
        sh = s[:, None] + 2.0 * np.sqrt(h)[:, None] * nodes[None, :]
        term = self.A * (np.logaddexp(self.beta * sh, -self.beta * sh) - math.log(2.0)) / self.beta
        return np.sum(term * weights, axis=1) / math.sqrt(math.pi)

    # certificate: z = sigma c w, |c| <= A   (Lipschitz preservation; z-only generator)
    @property
    def segment_coefficients(self):
        return (-self.A, self.A)

    @property
    def box_low(self):
        return -self.sigma * self.A * np.abs(self.w)

    @property
    def box_high(self):
        return self.sigma * self.A * np.abs(self.w)

    @property
    def z_ball_radius(self):
        return self.sigma * self.A

    @property
    def z_center(self):
        return np.zeros(self.d)

    @property
    def dose_scale(self):
        return self.sigma * self.A

    @property
    def picard_proxy(self):
        return 3.0 * self.kappa * self.A**2 / self.sigma * self.T

    def certificate_labels(self):
        return {"span": "PDE-derived: translation invariance orthogonal to w",
                "segment": "PDE-derived: Lipschitz preservation |psi_s| <= A",
                "box": "PDE-derived: coordinatewise hull of the segment",
                "ball": "PDE-derived: Euclidean hull of the segment"}

    def certificate_parameters(self) -> dict[str, Any]:
        return {"kappa": self.kappa, "A": self.A, "beta": self.beta,
                "reference_grid_gap": self.reference_grid_gap}
