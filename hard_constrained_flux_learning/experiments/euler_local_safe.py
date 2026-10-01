
"""Local conservative admissibility limiting for 1D Euler HCFL.

This module assumes:
1) F_hi is already projected into the Tadmor entropy half-space.
2) F_lo is a robust low-order flux (we use entropy-projected Rusanov).
3) The low-order forward-Euler state is admissible under the chosen CFL.

The limiter keeps one shared alpha per interface, hence conservation is
preserved. Because both endpoint fluxes are entropy feasible and the Tadmor
constraint is affine in the interface flux, every convex interface blend also
remains entropy feasible.
"""
import torch
import euler_1d_hcfl as E

def local_admissibility_limiter(
    U, F_hi, F_lo, lam, rho_floor=1e-5, p_floor=1e-5,
    max_outer=12, n_bisect=30
):
    dF = F_hi - F_lo
    U_lo = E.fv_step_lam(U, F_lo, lam)
    if not E.admissible(U_lo, rho_floor, p_floor).all():
        raise RuntimeError("Low-order state is not admissible; reduce CFL.")

    B, N, _ = U.shape
    alpha = torch.ones(B, N, dtype=U.dtype, device=U.device)

    for _ in range(max_outer):
        F = F_lo + alpha[..., None] * dF
        U_cur = E.fv_step_lam(U, F, lam)
        good = E.admissible(U_cur, rho_floor, p_floor)
        if good.all():
            return F, alpha

        bad = ~good
        dU = U_cur - U_lo

        lo = torch.zeros(int(bad.sum()), dtype=torch.float64, device=U.device)
        hi = torch.ones(int(bad.sum()), dtype=torch.float64, device=U.device)
        U0 = U_lo[bad].double()
        dUb = dU[bad].double()

        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            ok = E.admissible(U0 + mid[:, None] * dUb, rho_floor, p_floor)
            lo = torch.where(ok, mid, lo)
            hi = torch.where(ok, hi, mid)

        theta = torch.ones(B, N, dtype=U.dtype, device=U.device)
        theta[bad] = torch.clamp(lo.float() - 1e-5, min=0.0)

        # Interface i is shared by cells i and i+1.
        factor = torch.minimum(theta, torch.roll(theta, -1, dims=-1))
        alpha = alpha * factor

    # Conservative fallback on any trajectory that still fails.
    F = F_lo + alpha[..., None] * dF
    U_cur = E.fv_step_lam(U, F, lam)
    if not E.admissible(U_cur, rho_floor, p_floor).all():
        bad_traj = (~E.admissible(U_cur, rho_floor, p_floor)).any(dim=-1)
        alpha[bad_traj] = 0.0
        F = F_lo + alpha[..., None] * dF

    return F, alpha
