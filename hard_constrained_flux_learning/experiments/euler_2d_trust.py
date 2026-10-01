"""State-adaptive hard trust region for 2D Euler HCFL.

The learned candidate and the classical base are both first projected into the
Tadmor entropy half-space. Their convex blend is therefore also entropy
feasible.

Calibration uses only training trajectories. Let z(U) be the maximum over grid
cells of the standardized primitive-state norm. Let q be the 99th percentile
of z on training snapshots. At inference,

    tau(U) = clip((q / z(U))^2, tau_min, 1).

The realized face flux is

    F = F_base + tau(U) * (F_learned - F_base),

with F_base = entropy-projected HLLC.

This leaves in-distribution states unchanged and automatically falls back
toward HLLC under severe state-space shift.
"""
import torch

import euler_2d_hcfl as E2


def calibration_quantile(train_snapshots, mean, std, quantile=0.99):
    mean = torch.as_tensor(mean)
    std = torch.as_tensor(std)
    P = E2.primitive(train_snapshots)
    cell_score = torch.sqrt((((P - mean) / std) ** 2).mean(dim=-1))
    state_score = cell_score.amax(dim=(-2, -1))
    return float(torch.quantile(state_score.reshape(-1), quantile))


def state_score(U, mean, std):
    mean = torch.as_tensor(mean, dtype=U.dtype, device=U.device)
    std = torch.as_tensor(std, dtype=U.dtype, device=U.device)
    P = E2.primitive(U)
    return torch.sqrt((((P - mean) / std) ** 2).mean(dim=-1)).amax(dim=(-2, -1))


def trust_coefficient(U, mean, std, q99, tau_min=0.05, power=2.0):
    score = state_score(U, mean, std)
    tau = torch.clamp((q99 / (score + 1e-12)) ** power,
                      min=tau_min, max=1.0)
    return tau, score


def trusted_fluxes(model, U, mean, std, q99, tau_min=0.05, power=2.0):
    # Learned candidate is already hard-projected by model.fluxes.
    Fx_learn, Gy_learn = model.fluxes(U)

    # Project HLLC too, so both endpoints lie in the same convex Tadmor set.
    Fx_base = E2.hard_project(E2.hllc_flux(U, "x"), U, "x")
    Gy_base = E2.hard_project(E2.hllc_flux(U, "y"), U, "y")

    tau, score = trust_coefficient(U, mean, std, q99, tau_min, power)
    t = tau[:, None, None, None]

    Fx = Fx_base + t * (Fx_learn - Fx_base)
    Gy = Gy_base + t * (Gy_learn - Gy_base)
    return Fx, Gy, tau, score
