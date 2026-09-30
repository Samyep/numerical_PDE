"""Compact learned-flux Burgers prototype.

This file demonstrates the architecture, not the final training protocol:
local stencil -> DNN proposal flux -> optional hard entropy projection ->
conservative finite-volume update.

For the final research model, replace direct flux supervision with multi-step
trajectory-only supervision.
"""
import torch
from torch import nn


def burgers_flux(u):
    return 0.5 * u**2


class FluxNet(nn.Module):
    def __init__(self, width=32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(5, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 1),
        )

    def forward(self, stencil):
        return self.net(stencil).squeeze(-1)


def periodic_stencils(u):
    return torch.stack([
        torch.roll(u, 2, dims=-1),
        torch.roll(u, 1, dims=-1),
        u,
        torch.roll(u, -1, dims=-1),
        torch.roll(u, -2, dims=-1),
    ], dim=-1)


def hard_q(proposal, a, b):
    threshold = (a*a + a*b + b*b) / 6.0
    out = torch.where(b > a, torch.minimum(proposal, threshold), proposal)
    out = torch.where(b < a, torch.maximum(out, threshold), out)
    return out


def hard_kruzkov(proposal, a, b, ks):
    K = ks.view(*([1] * a.ndim), -1)
    A, B = a.unsqueeze(-1), b.unsqueeze(-1)
    FK = 0.5 * K**2

    upper = torch.where((K > A) & (K < B), FK, torch.inf).amin(dim=-1)
    lower = torch.where((K > B) & (K < A), FK, -torch.inf).amax(dim=-1)

    out = torch.where((a < b) & torch.isfinite(upper), torch.minimum(proposal, upper), proposal)
    out = torch.where((a > b) & torch.isfinite(lower), torch.maximum(out, lower), out)
    return out


def conservative_step(u, interface_flux, dt_dx):
    return u - dt_dx * (interface_flux - torch.roll(interface_flux, 1, dims=-1))


def one_step(model, u, dt_dx, mode="plain", ks=None):
    stencil = periodic_stencils(u)
    proposal = model(stencil)
    a, b = u, torch.roll(u, -1, dims=-1)

    if mode == "Q":
        flux = hard_q(proposal, a, b)
    elif mode == "K":
        flux = hard_kruzkov(proposal, a, b, ks)
    else:
        flux = proposal

    return conservative_step(u, flux, dt_dx)


if __name__ == "__main__":
    torch.manual_seed(0)
    model = FluxNet()
    u = torch.sin(torch.linspace(0, 2 * torch.pi, 128))
    ks = torch.linspace(-2.5, 2.5, 17)
    unext = one_step(model, u, dt_dx=0.02, mode="K", ks=ks)
    print(unext.shape)
