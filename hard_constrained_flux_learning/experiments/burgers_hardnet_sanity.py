"""Minimal Burgers sanity check for finite hard entropy constraints.

The solver uses a conservative finite-volume update. A central-flux proposal is
projected onto either one quadratic-entropy constraint or a finite Kruzkov
family. Burgers is used as a diagnostic; the full Kruzkov family approaches the
classical Godunov construction for convex flux.
"""
import numpy as np


def phys_flux(u):
    return 0.5 * u**2


def central_flux(a, b):
    return 0.5 * (phys_flux(a) + phys_flux(b))


def godunov_flux(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    out = np.empty(np.broadcast(a, b).shape)
    rare = a <= b
    ar, br = a[rare], b[rare]
    rr = np.empty_like(ar)
    cross = (ar <= 0) & (br >= 0)
    rr[cross] = 0.0
    rr[ar > 0] = phys_flux(ar[ar > 0])
    rr[br < 0] = phys_flux(br[br < 0])
    out[rare] = rr
    sh = ~rare
    out[sh] = np.maximum(phys_flux(a[sh]), phys_flux(b[sh]))
    return out


def hard_quadratic_entropy(a, b):
    F = central_flux(a, b)
    threshold = (a*a + a*b + b*b) / 6.0
    F = np.where(b > a, np.minimum(F, threshold), F)
    F = np.where(b < a, np.maximum(F, threshold), F)
    return F


def make_hard_kruzkov(k_values):
    ks = np.asarray(k_values, dtype=float)
    fks = phys_flux(ks)

    def flux(a, b):
        a = np.asarray(a)
        b = np.asarray(b)
        F = central_flux(a, b).copy()
        A, B, K = a[:, None], b[:, None], ks[None, :]

        mask_r = (K > A) & (K < B)
        upper = np.min(np.where(mask_r, fks[None, :], np.inf), axis=1)
        idx = (a < b) & np.isfinite(upper)
        F[idx] = np.minimum(F[idx], upper[idx])

        mask_s = (K > B) & (K < A)
        lower = np.max(np.where(mask_s, fks[None, :], -np.inf), axis=1)
        idx = (a > b) & np.isfinite(lower)
        F[idx] = np.maximum(F[idx], lower[idx])
        return F

    return flux


def fv_step(u, dx, dt, flux):
    left = np.concatenate([[u[0]], u])
    right = np.concatenate([u, [u[-1]]])
    F = flux(left, right)
    return u - (dt / dx) * (F[1:] - F[:-1])


def run(u0, x, t_end, flux, cfl=0.15):
    u = u0.copy().astype(float)
    dx = x[1] - x[0]
    t = 0.0
    while t < t_end - 1e-14:
        dt = min(cfl * dx / max(np.max(np.abs(u)), 1e-10), t_end - t)
        u = fv_step(u, dx, dt, flux)
        t += dt
    return u


def main():
    n = 400
    x = np.linspace(-1, 1, n, endpoint=False) + 1/n
    u0 = np.where(x < 0, 2.0, -1.0)
    t_end = 0.18

    methods = {
        "central": central_flux,
        "Q": hard_quadratic_entropy,
        "K17": make_hard_kruzkov(np.linspace(-2.5, 2.5, 17)),
        "K257": make_hard_kruzkov(np.linspace(-2.5, 2.5, 257)),
        "godunov": godunov_flux,
    }

    for name, flux in methods.items():
        u = run(u0, x, t_end, flux)
        print(name, "min/max", float(u.min()), float(u.max()))


if __name__ == "__main__":
    main()
