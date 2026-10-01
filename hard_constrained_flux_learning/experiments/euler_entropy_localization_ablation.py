"""Ablations for fully-discrete entropy limiting in 1D Euler HCFL.

Requires:
  euler_1d_hcfl.py
  euler_hllc_hcfl.py
  euler_local_safe.py

This file records the two localization attempts studied after the HLLC/local-
admissibility phase:

1) strict cell-wise entropy budgets;
2) sequential interface-wise spending of a global total-entropy slack.

Both are mathematically safe constructions, but both were less accurate than
the current default: local admissibility limiting + a rare global scalar
fully-discrete entropy line search. The main paper method should therefore keep
the current global entropy scalar unless a better local construction is found.
"""
import numpy as np
import torch

import euler_1d_hcfl as E


def numerical_entropy_flux(Fh, U):
    """Tadmor numerical entropy flux Q = {v}^T F - {psi}."""
    UR = torch.roll(U, -1, dims=-2)
    vavg = 0.5 * (E.entropy_variables(U) + E.entropy_variables(UR))
    psiavg = 0.5 * (E.entropy_potential(U) + E.entropy_potential(UR))
    return (vavg * Fh).sum(dim=-1) - psiavg


def local_entropy_budget(U, Flo, lam):
    Q = numerical_entropy_flux(Flo, U)
    return E.entropy(U) - lam * (Q - torch.roll(Q, 1, dims=-1))


def total_entropy(U):
    return E.entropy(U.double()).sum(dim=-1)


def sequential_entropy_slack_limiter(
    U, Fhi, Flo, lam,
    rho_floor=1e-5, p_floor=1e-5,
    n_bisect=38, tol=1e-10, reverse=False,
):
    """Sequentially add interface corrections while spending global entropy slack.

    Start from a safe low-order state. At interface i, only cells i and i+1
    change. Choose the largest beta_i in [0,1] such that both cells remain
    admissible and the total entropy remains below its value at the beginning
    of the step.

    This is conservative and fully-discrete entropy safe, but order dependent.
    In the canonical Euler tests it was more restrictive than the default
    global scalar entropy line search.
    """
    U = U.double()
    Flo = E.hard_entropy_projection(Flo.double(), U)
    Fhi = E.hard_entropy_projection(Fhi.double(), U)

    Ucur = E.fv_step_lam(U, Flo, lam)
    if not E.admissible(Ucur, rho_floor, p_floor).all():
        return None

    E0 = total_entropy(U)
    Ecur = total_entropy(Ucur)
    if (Ecur > E0 + tol).any():
        return None

    slack = E0 - Ecur
    B, N, _ = U.shape
    dF = Fhi - Flo
    beta = torch.zeros(B, N, dtype=torch.float64)

    order = range(N - 1, -1, -1) if reverse else range(N)

    for i in order:
        j = (i + 1) % N
        L0 = Ucur[:, i, :].clone()
        R0 = Ucur[:, j, :].clone()
        d = dF[:, i, :]
        pair0 = E.entropy(L0) + E.entropy(R0)

        def feasible(c):
            VL = L0 - lam * c[:, None] * d
            VR = R0 + lam * c[:, None] * d
            adm = E.admissible(VL, rho_floor, p_floor) & E.admissible(VR, rho_floor, p_floor)
            de = E.entropy(VL) + E.entropy(VR) - pair0
            return adm & (de <= slack + tol), VL, VR, de

        one = torch.ones(B, dtype=torch.float64)
        good, _, _, _ = feasible(one)
        coeff = torch.ones(B, dtype=torch.float64)
        need = ~good

        if need.any():
            lo = torch.zeros(int(need.sum()), dtype=torch.float64)
            hi = torch.ones(int(need.sum()), dtype=torch.float64)
            L0n, R0n, dn = L0[need], R0[need], d[need]
            pair0n, slackn = pair0[need], slack[need]

            for _ in range(n_bisect):
                mid = 0.5 * (lo + hi)
                VL = L0n - lam * mid[:, None] * dn
                VR = R0n + lam * mid[:, None] * dn
                adm = E.admissible(VL, rho_floor, p_floor) & E.admissible(VR, rho_floor, p_floor)
                de = E.entropy(VL) + E.entropy(VR) - pair0n
                ok = adm & (de <= slackn + tol)
                lo = torch.where(ok, mid, lo)
                hi = torch.where(ok, hi, mid)

            coeff[need] = torch.clamp(lo - 1e-8, min=0.0)

        _, Lnew, Rnew, de = feasible(coeff)
        beta[:, i] = coeff
        Ucur[:, i, :] = Lnew
        Ucur[:, j, :] = Rnew
        slack = slack - de

    return Flo + beta[..., None] * dF, beta


if __name__ == "__main__":
    print("This module provides the fully-discrete entropy localization ablations.")
    print("See notes/euler_fully_discrete_entropy_ablation.md for the recorded results.")
