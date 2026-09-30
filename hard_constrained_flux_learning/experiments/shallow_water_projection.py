"""Sanity check for the closed-form HardNet-style Tadmor projection.

This is not a full shallow-water solver. It verifies that a flux proposal can
be projected onto the interface entropy half-space to machine precision.
"""
import numpy as np


def shallow_water_flux(h, m, g=9.81):
    return np.stack([m, m * m / h + 0.5 * g * h * h], axis=-1)


def entropy_variables(h, m, g=9.81):
    u = m / h
    return np.stack([g * h - 0.5 * u * u, u], axis=-1)


def entropy_potential(h, m, g=9.81):
    return 0.5 * g * h * m


def hard_entropy_projection(proposal, hL, mL, hR, mR, g=9.81, eps=1e-14):
    vL = entropy_variables(hL, mL, g)
    vR = entropy_variables(hR, mR, g)
    a = vR - vL
    b = entropy_potential(hR, mR, g) - entropy_potential(hL, mL, g)

    residual = np.sum(a * proposal, axis=-1) - b
    norm2 = np.sum(a * a, axis=-1)
    alpha = np.where((residual > 0) & (norm2 > eps), residual / norm2, 0.0)
    projected = proposal - alpha[..., None] * a
    return projected


def main(seed=0, n=200_000):
    rng = np.random.default_rng(seed)
    g = 9.81

    hL = rng.uniform(0.2, 3.0, n)
    hR = rng.uniform(0.2, 3.0, n)
    uL = rng.uniform(-3.0, 3.0, n)
    uR = rng.uniform(-3.0, 3.0, n)
    mL, mR = hL * uL, hR * uR

    fL = shallow_water_flux(hL, mL, g)
    fR = shallow_water_flux(hR, mR, g)
    proposal = 0.5 * (fL + fR) + rng.normal(size=(n, 2)) * np.array([2.0, 8.0])

    vL = entropy_variables(hL, mL, g)
    vR = entropy_variables(hR, mR, g)
    a = vR - vL
    b = entropy_potential(hR, mR, g) - entropy_potential(hL, mL, g)

    r_before = np.sum(a * proposal, axis=-1) - b
    projected = hard_entropy_projection(proposal, hL, mL, hR, mR, g)
    r_after = np.sum(a * projected, axis=-1) - b

    print("fraction violating before:", np.mean(r_before > 0))
    print("max residual before:", np.max(r_before))
    print("fraction violating after > 1e-10:", np.mean(r_after > 1e-10))
    print("max residual after:", np.max(r_after))


if __name__ == "__main__":
    main()
