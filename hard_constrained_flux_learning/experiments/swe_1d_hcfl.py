import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
import torch.nn.functional as F
import matplotlib.pyplot as plt

torch.set_num_threads(2)

G = 9.81
NREF = 1024
NCOARSE = 64
FACTOR = NREF // NCOARSE
DT_SNAPSHOT = 5e-4
LAMBDA = DT_SNAPSHOT / (1.0 / NCOARSE)
NSNAP = 18


# ---------- reference solver ----------
def np_flux(U):
    h = U[..., 0]
    m = U[..., 1]
    u = m / h
    return np.stack([m, m * u + 0.5 * G * h * h], axis=-1)


def np_rusanov(UL, UR):
    fL, fR = np_flux(UL), np_flux(UR)
    hL, hR = UL[..., 0], UR[..., 0]
    uL, uR = UL[..., 1] / hL, UR[..., 1] / hR
    a = np.maximum(np.abs(uL) + np.sqrt(G * hL),
                   np.abs(uR) + np.sqrt(G * hR))
    return 0.5 * (fL + fR) - 0.5 * a[..., None] * (UR - UL)


def np_rhs(U, dx):
    Fh = np_rusanov(U, np.roll(U, -1, axis=-2))
    return -(Fh - np.roll(Fh, 1, axis=-2)) / dx


def np_ssprk2(U, dt, dx):
    U1 = U + dt * np_rhs(U, dx)
    return 0.5 * U + 0.5 * (U1 + dt * np_rhs(U1, dx))


def generate_ic(ntraj, N, seed, ood=False):
    rng = np.random.default_rng(seed)
    x = (np.arange(N) + 0.5) / N
    U = np.zeros((ntraj, N, 2), dtype=np.float64)

    for j in range(ntraj):
        typ = rng.choice(["riemann", "piecewise", "smooth"],
                         p=[0.40, 0.30, 0.30])

        if typ == "riemann":
            if ood:
                hL, hR = rng.uniform(0.45, 2.4, 2)
                uL, uR = rng.uniform(-1.2, 1.2, 2)
            else:
                hL, hR = rng.uniform(0.7, 1.7, 2)
                uL, uR = rng.uniform(-0.7, 0.7, 2)
            cut = int(rng.integers(N // 4, 3 * N // 4))
            h = np.r_[np.full(cut, hL), np.full(N - cut, hR)]
            u = np.r_[np.full(cut, uL), np.full(N - cut, uR)]
            shift = int(rng.integers(0, N))
            h, u = np.roll(h, shift), np.roll(u, shift)

        elif typ == "piecewise":
            K = int(rng.integers(3, 6))
            cuts = sorted(rng.choice(np.arange(1, N), K - 1, replace=False))
            h, u = np.empty(N), np.empty(N)
            start = 0
            for k in range(K):
                end = cuts[k] if k < K - 1 else N
                if ood:
                    hv, uv = rng.uniform(0.45, 2.4), rng.uniform(-1.2, 1.2)
                else:
                    hv, uv = rng.uniform(0.7, 1.7), rng.uniform(-0.7, 0.7)
                h[start:end], u[start:end] = hv, uv
                start = end

        else:
            h0 = rng.uniform(0.7, 1.8) if ood else rng.uniform(0.9, 1.4)
            ah = rng.uniform(0.1, 0.4) if ood else rng.uniform(0.05, 0.25)
            au = rng.uniform(0.2, 0.9) if ood else rng.uniform(0.05, 0.5)
            k1, k2 = int(rng.integers(1, 4)), int(rng.integers(1, 4))
            p1, p2 = rng.uniform(0, 2 * np.pi, 2)
            h = h0 + ah * np.sin(2 * np.pi * k1 * x + p1)
            if h.min() < 0.4:
                h += 0.4 - h.min()
            u = au * np.sin(2 * np.pi * k2 * x + p2) + rng.uniform(-0.2, 0.2)

        U[j, :, 0] = h
        U[j, :, 1] = h * u

    return U


def generate_trajectory(ntraj, seed, ood=False):
    U = generate_ic(ntraj, NREF, seed, ood)
    dx = 1.0 / NREF

    def restrict(V):
        return V.reshape(ntraj, NCOARSE, FACTOR, 2).mean(axis=2).astype(np.float32)

    snapshots = [restrict(U)]

    for _ in range(NSNAP - 1):
        remaining = DT_SNAPSHOT
        while remaining > 1e-14:
            h = U[..., 0]
            u = U[..., 1] / h
            max_speed = np.max(np.abs(u) + np.sqrt(G * h))
            dt = min(remaining, 0.35 * dx / max_speed)
            U = np_ssprk2(U, dt, dx)
            if U[..., 0].min() <= 0:
                raise RuntimeError("Reference solver produced non-positive water depth.")
            remaining -= dt
        snapshots.append(restrict(U))

    return np.stack(snapshots, axis=1)


# ---------- Torch fluxes / entropy ----------
def t_flux(U):
    h = U[..., 0].clamp_min(1e-6)
    m = U[..., 1]
    u = m / h
    return torch.stack([m, m * u + 0.5 * G * h * h], dim=-1)


def primitive(U):
    h = U[..., 0].clamp_min(1e-6)
    return torch.stack([h, U[..., 1] / h], dim=-1)


def entropy_variables(U):
    h = U[..., 0].clamp_min(1e-6)
    u = U[..., 1] / h
    return torch.stack([G * h - 0.5 * u * u, u], dim=-1)


def entropy_potential(U):
    h = U[..., 0].clamp_min(1e-6)
    m = U[..., 1]
    return 0.5 * G * h * m


def t_rusanov(U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    fL, fR = t_flux(UL), t_flux(UR)
    hL, hR = UL[..., 0].clamp_min(1e-6), UR[..., 0].clamp_min(1e-6)
    uL, uR = UL[..., 1] / hL, UR[..., 1] / hR
    a = torch.maximum(torch.abs(uL) + torch.sqrt(G * hL),
                      torch.abs(uR) + torch.sqrt(G * hR))
    return 0.5 * (fL + fR) - 0.5 * a[..., None] * (UR - UL)


def t_hll(U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    fL, fR = t_flux(UL), t_flux(UR)
    hL, hR = UL[..., 0].clamp_min(1e-6), UR[..., 0].clamp_min(1e-6)
    uL, uR = UL[..., 1] / hL, UR[..., 1] / hR
    cL, cR = torch.sqrt(G * hL), torch.sqrt(G * hR)
    sL = torch.minimum(uL - cL, uR - cR)
    sR = torch.maximum(uL + cL, uR + cR)
    mid = (sR[..., None] * fL - sL[..., None] * fR +
           (sL * sR)[..., None] * (UR - UL)) / (sR - sL).clamp_min(1e-8)[..., None]
    return torch.where((sL >= 0)[..., None], fL,
                       torch.where((sR <= 0)[..., None], fR, mid))


def entropy_residual(Fh, U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    a = entropy_variables(UR) - entropy_variables(UL)
    b = entropy_potential(UR) - entropy_potential(UL)
    return (a * Fh).sum(dim=-1) - b


def hard_entropy_projection(Fh, U):
    UL, UR = U, torch.roll(U, -1, dims=-2)
    a = entropy_variables(UR) - entropy_variables(UL)
    b = entropy_potential(UR) - entropy_potential(UL)
    residual = (a * Fh).sum(dim=-1) - b
    norm2 = (a * a).sum(dim=-1)

    # Important implementation detail: perform the division only on the active
    # mask. Eagerly forming residual/norm2 inside torch.where can create
    # unstable backward gradients on near-identical states even when that
    # branch is not selected.
    alpha = torch.zeros_like(residual)
    mask = (residual > 0) & (norm2 > 1e-14)
    alpha[mask] = residual[mask] / norm2[mask]
    return Fh - alpha[..., None] * a


def fv_step(U, Fh):
    return U - LAMBDA * (Fh - torch.roll(Fh, 1, dims=-2))


# ---------- learned flux ----------
class LearnedFlux(nn.Module):
    def __init__(self, h_mean, h_std, u_mean, u_std, width=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(10, width), nn.Tanh(),
            nn.Linear(width, width), nn.Tanh(),
            nn.Linear(width, 2),
        )
        # Safe starting point: coarse Rusanov.
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)

        self.register_buffer("mean", torch.tensor([h_mean, u_mean], dtype=torch.float32))
        self.register_buffer("std", torch.tensor([h_std, u_std], dtype=torch.float32))
        self.register_buffer("flux_scale", torch.tensor([0.8, 5.0], dtype=torch.float32))

    def forward(self, U):
        P = primitive(U)
        feats = torch.cat(
            [(torch.roll(P, s, dims=-2) - self.mean) / self.std
             for s in [2, 1, 0, -1, -2]],
            dim=-1,
        )
        correction = torch.tanh(self.net(feats))

        UL, UR = U, torch.roll(U, -1, dims=-2)
        PL, PR = primitive(UL), primitive(UR)
        jump = torch.sqrt((((PR - PL) / self.std) ** 2).sum(dim=-1) + 1e-12)

        # Five-point learned correction to a conservative Rusanov proposal.
        return t_rusanov(U) + 0.30 * jump[..., None] * correction * self.flux_scale


class Solver(nn.Module):
    def __init__(self, stats, mode="plain", soft_weight=1e-2):
        super().__init__()
        self.mode = mode
        self.soft_weight = soft_weight
        self.flux_net = LearnedFlux(*stats)

    def flux(self, U):
        raw = self.flux_net(U)
        if self.mode == "hard":
            return hard_entropy_projection(raw, U)
        return raw

    def one_step(self, U):
        raw = self.flux_net(U)
        penalty = torch.tensor(0.0, dtype=U.dtype, device=U.device)
        if self.mode == "hard":
            Fh = hard_entropy_projection(raw, U)
        else:
            Fh = raw
            if self.mode == "soft":
                penalty = torch.relu(entropy_residual(Fh, U)).pow(2).mean()
        return fv_step(U, Fh), penalty


def train_one_seed(seed, output_dir, iterations=1000):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_np = generate_trajectory(500, 1000 + seed, ood=False)
    val_np = generate_trajectory(100, 2000 + seed, ood=False)
    ood_np = generate_trajectory(100, 3000 + seed, ood=True)

    train = torch.tensor(train_np)
    val = torch.tensor(val_np)
    ood = torch.tensor(ood_np)

    hp = train[..., 0]
    up = train[..., 1] / hp
    stats = (
        float(hp.mean()), float(hp.std()),
        float(up.mean()), float(up.std()),
    )
    state_std = torch.tensor(
        [float(train[..., 0].std()), float(train[..., 1].std())],
        dtype=torch.float32,
    )

    methods = ["plain", "soft", "hard"]
    rows = []
    models = {}

    # Same neural initialization and minibatch schedule within each seed.
    torch.manual_seed(4000 + seed)
    base = Solver(stats, mode="plain")
    base_state = base.flux_net.state_dict()

    for mode in methods:
        torch.manual_seed(4000 + seed)
        model = Solver(stats, mode=mode)
        model.flux_net.load_state_dict(base_state)
        opt = torch.optim.Adam(model.parameters(), lr=5e-4)
        gen = torch.Generator().manual_seed(5000 + seed)

        for _ in range(iterations):
            B = 64
            inds = torch.randint(0, train.shape[0], (B,), generator=gen)
            ts = torch.randint(0, train.shape[1] - 1, (B,), generator=gen)
            U = torch.stack([train[i, t] for i, t in zip(inds, ts)], dim=0)
            target = torch.stack([train[i, t + 1] for i, t in zip(inds, ts)], dim=0)

            pred, penalty = model.one_step(U)
            traj_loss = (((pred - target) / state_std) ** 2).mean()
            loss = traj_loss + model.soft_weight * penalty

            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()

        models[mode] = model.eval()
        torch.save(model.state_dict(), output_dir / f"seed{seed}_{mode}.pt")

    @torch.no_grad()
    def eval_learned(model, data, split):
        U = data[:, 0].clone()
        mass0 = U[..., 0].sum(dim=-1)
        err2 = []
        min_h = float("inf")
        violation_count = 0
        violation_total = 0
        max_residual = 0.0
        max_mass_drift = 0.0
        failed = False

        for t in range(1, data.shape[1]):
            Fh = model.flux(U)
            res = entropy_residual(Fh, U)
            violation_count += int((res > 1e-5).sum())
            violation_total += res.numel()
            max_residual = max(max_residual, float(res.max()))

            U = fv_step(U, Fh)
            if not torch.isfinite(U).all():
                failed = True
                break

            min_h = min(min_h, float(U[..., 0].min()))
            max_mass_drift = max(
                max_mass_drift,
                float(torch.abs(U[..., 0].sum(dim=-1) - mass0).max()),
            )
            err2.append(float((((U - data[:, t]) / state_std) ** 2).mean()))

        rollout_nrmse = float(np.sqrt(np.mean(err2))) if err2 else np.inf
        return {
            "seed": seed,
            "split": split,
            "method": model.mode,
            "rollout_nrmse": rollout_nrmse,
            "min_h": min_h,
            "entropy_violation_rate": violation_count / max(1, violation_total),
            "max_entropy_residual": max_residual,
            "max_mass_drift": max_mass_drift,
            "failed": failed,
        }

    @torch.no_grad()
    def eval_baseline(name, data, split):
        U = data[:, 0].clone()
        mass0 = U[..., 0].sum(dim=-1)
        err2 = []
        min_h = float("inf")
        max_mass_drift = 0.0

        for t in range(1, data.shape[1]):
            Fh = t_rusanov(U) if name == "rusanov" else t_hll(U)
            U = fv_step(U, Fh)
            min_h = min(min_h, float(U[..., 0].min()))
            max_mass_drift = max(
                max_mass_drift,
                float(torch.abs(U[..., 0].sum(dim=-1) - mass0).max()),
            )
            err2.append(float((((U - data[:, t]) / state_std) ** 2).mean()))

        return {
            "seed": seed,
            "split": split,
            "method": name,
            "rollout_nrmse": float(np.sqrt(np.mean(err2))),
            "min_h": min_h,
            "entropy_violation_rate": 0.0,
            "max_entropy_residual": np.nan,
            "max_mass_drift": max_mass_drift,
            "failed": False,
        }

    for split, data in [("ID", val), ("OOD", ood)]:
        for mode in methods:
            rows.append(eval_learned(models[mode], data, split))
        rows.append(eval_baseline("rusanov", data, split))
        rows.append(eval_baseline("hll", data, split))

    result = pd.DataFrame(rows)
    result.to_csv(output_dir / f"seed{seed}_metrics.csv", index=False)

    # Representative OOD trajectory figure.
    with torch.no_grad():
        sample = ood[:1]
        x = (np.arange(NCOARSE) + 0.5) / NCOARSE
        final_ref = sample[0, -1, :, 0].numpy()
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(x, final_ref, linewidth=2.5, label="fine-grid reference")
        for mode in methods:
            U = sample[:, 0].clone()
            for _ in range(1, sample.shape[1]):
                U = fv_step(U, models[mode].flux(U))
            ax.plot(x, U[0, :, 0].numpy(), linewidth=1.2, label=mode)
        ax.set_xlabel("x")
        ax.set_ylabel("h")
        ax.set_title(f"1D SWE OOD rollout, seed {seed}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"seed{seed}_ood_height.png", dpi=170)
        plt.close(fig)

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=str, default="./swe_1d_runs")
    args = parser.parse_args()
    result = train_one_seed(args.seed, args.output)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
