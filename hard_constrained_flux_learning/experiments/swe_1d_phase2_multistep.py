"""Four-step trajectory fine-tuning for 1D SWE HCFL.

Start from a one-step-trained hard-entropy model produced by ``swe_1d_hcfl.py``.
Fine-tune through K conservative steps using the differentiable entropy projection
and the local positivity flux limiter.  Evaluation uses the full adaptive-CFL
safe rollout from ``swe_1d_phase2_safe.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import swe_1d_hcfl as base
import swe_1d_phase2_safe as safe


def train_multistep(
    seed: int,
    pretrained_dir: str,
    output_dir: str,
    iterations: int = 500,
    horizon: int = 4,
):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    train = torch.tensor(base.generate_trajectory(500, 1000 + seed, ood=False))
    val = torch.tensor(base.generate_trajectory(100, 2000 + seed, ood=False))
    ood = torch.tensor(base.generate_trajectory(100, 3000 + seed, ood=True))
    state_std = torch.tensor([float(train[..., 0].std()), float(train[..., 1].std())])

    checkpoint = Path(pretrained_dir) / f"seed{seed}_hard.pt"
    original = base.Solver((1, 1, 0, 1), mode="hard")
    original.load_state_dict(torch.load(checkpoint, map_location="cpu"))
    original.eval()

    model = base.Solver((1, 1, 0, 1), mode="hard")
    model.load_state_dict(original.state_dict())
    opt = torch.optim.Adam(model.parameters(), lr=2e-4)
    gen = torch.Generator().manual_seed(9000 + seed)

    # Ordinary training distribution satisfies the fixed saved-step CFL.  The
    # local limiter is differentiable almost everywhere and is enough to avoid
    # accidental positivity failures while backpropagating through the rollout.
    for _ in range(iterations):
        batch = 32
        inds = torch.randint(0, train.shape[0], (batch,), generator=gen)
        ts = torch.randint(0, train.shape[1] - horizon, (batch,), generator=gen)
        U = torch.stack([train[i, t] for i, t in zip(inds, ts)], 0)
        loss = 0.0

        for k in range(1, horizon + 1):
            target = torch.stack([train[i, t + k] for i, t in zip(inds, ts)], 0)
            Fent = model.flux(U)
            Fpos, _, _, _ = safe.positivity_limited_flux(
                U, Fent, base.LAMBDA, h_floor=1e-4
            )
            U = safe.fv_step_lambda(U, Fpos, base.LAMBDA)
            loss = loss + (((U - target) / state_std) ** 2).mean() / horizon

        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()

    model.eval()
    torch.save(model.state_dict(), out / f"seed{seed}_multistep.pt")

    @torch.no_grad()
    def evaluate(m, data, split, label):
        U = data[:, 0].clone()
        errors, min_h = [], float("inf")
        p_active, e_active, substeps = [], [], []
        max_dE = -float("inf")
        for t in range(1, data.shape[1]):
            E0 = safe.total_entropy(U)
            U, pa, ea, ns = safe.advance_snapshot_safe(m, U, cfl=0.45)
            errors.append(float((((U - data[:, t]) / state_std) ** 2).mean()))
            min_h = min(min_h, float(U[..., 0].min()))
            p_active.append(pa)
            e_active.append(ea)
            substeps.append(ns)
            max_dE = max(max_dE, float((safe.total_entropy(U) - E0).max()))
        return {
            "seed": seed,
            "split": split,
            "method": label,
            "rollout_nrmse": float(np.sqrt(np.mean(errors))),
            "min_h": min_h,
            "positivity_active": float(np.mean(p_active)),
            "fd_entropy_active": float(np.mean(e_active)),
            "mean_substeps": float(np.mean(substeps)),
            "max_total_entropy_change": max_dE,
        }

    rows = []
    for split, data in [("ID", val), ("OOD", ood)]:
        rows.append(evaluate(original, data, split, "HCFL-safe-one-step"))
        rows.append(evaluaue(model, data, split, "HCFL-safe-multistep"))

    df = pd.DataFrame(rows)
    df.to_csv(out / f"seed{seed}_metrics.csv", index=False)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--pretrained-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, default="results/swe_phase2_multistep")
    parser.add_argument("--iterations", type=int, default=500)
    parser.add_argument("--horizon", type=int, default=4)
    args = parser.parse_args()
    result = train_multistep(
        args.seed,
        args.pretrained_dir,
        args.output_dir,
        args.iterations,
        args.horizon,
    )
    print(result.to_string(index=False))
