"""Standalone IR-MLP comparison on the SCaSML LCD, VB, and DR benchmarks.

Main track:
- full-history MLP, n=2, M=10
- 1000 interior + 200 boundary test points
- standard Elworthy--Bismut--Li terminal-gradient normalization
- projection happens inside recursion

Use --include-upstream-normalization only as a diagnostic. The public
MLP_full_history.py divides the terminal z estimator by T-t; the main track
uses sqrt(T-t).

The defaults follow the executable public SCaSML repo where it conflicts with
paper prose: LCD evaluation uses the inherited [-0.5,0.5]^d test geometry and
the current Grad_Dependent_Nonlinear class uses sigma=0.25.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PAPER = {
    "LCD": {10: (.227, .0274), 20: (.235, .0472), 30: (.238, .0972), 60: (.239, .132)},
    "VB": {20: (.0836, .00403), 40: (.104, .0292), 60: (.117, .0288), 80: (.119, .0564)},
    "DR": {100: (.0899, .0111), 120: (.0913, .0103), 140: (.0897, .0300), 160: (.0900, .0322)},
}


def sigmoid(a):
    out = np.empty_like(a, dtype=np.float64)
    pos = a >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-a[pos]))
    ea = np.exp(a[~pos])
    out[~pos] = ea / (1.0 + ea)
    return out


@dataclass
class LCD:
    d: int
    T: float = .5

    def __post_init__(self):
        self.mu = -1.0 / self.d
        self.sigma = math.sqrt(2.0)
        self.norm = .5 * (self.d + 1)

    def g(self, xt):
        return (xt[:, :-1].sum(1) + self.T)[:, None]

    def f(self, xt, u, z):
        return np.zeros_like(u)

    def exact(self, xt):
        return xt[:, :-1].sum(1) + xt[:, -1]

    def project(self, out, xt):
        # On the executable test domain x in [-.5,.5]^d, t in [0,.5].
        u = np.clip(out[:, 0], -.5 * self.d, .5 * self.d + self.T)
        # Exact affine solution: grad u=1, hence z=sigma*1.
        z = np.full_like(out[:, 1:], self.sigma)
        return np.c_[u, z]


@dataclass
class VB:
    d: int
    sigma: float = .25
    T: float = .5

    def __post_init__(self):
        self.mu = -(1.0 / self.d + self.sigma**2 / 2.0)
        self.norm = 1.0
        self.z_radius = self.sigma * math.sqrt(self.d) / 4.0

    def g(self, xt):
        return sigmoid(self.T + xt[:, :-1].sum(1))[:, None]

    def f(self, xt, u, z):
        return self.sigma * u * np.sum(z, axis=1, keepdims=True)

    def exact(self, xt):
        return sigmoid(xt[:, -1] + xt[:, :-1].sum(1))

    def project(self, out, xt):
        # u=sigmoid(.) => 0<=u<=1, each grad coordinate in [0,1/4].
        u = np.clip(out[:, 0], 0.0, 1.0)
        z = np.maximum(out[:, 1:], 0.0)
        zn = np.linalg.norm(z, axis=1)
        scale = np.minimum(1.0, self.z_radius / np.maximum(zn, 1e-12))
        return np.c_[u, z * scale[:, None]]


@dataclass
class DR:
    d: int
    T: float = 1.0
    kappa: float = 1.6
    lam: float = .1

    def __post_init__(self):
        self.mu = 0.0
        self.sigma = 1.0
        self.norm = 10.0

    def amp(self, t):
        return np.exp(self.lam**2 * self.d * (t - self.T) / 2.0)

    def exact(self, xt):
        x, t = xt[:, :-1], xt[:, -1]
        return self.kappa + np.sin(self.lam * x.sum(1)) * self.amp(t)

    def g(self, xt):
        return self.exact(xt)[:, None]

    def f(self, xt, u, z):
        e = u[:, 0] - self.exact(xt)
        return np.minimum(1.0, e * e)[:, None]

    def project(self, out, xt):
        a = self.amp(xt[:, -1])
        u = np.minimum(np.maximum(out[:, 0], self.kappa - a), self.kappa + a)
        z = out[:, 1:].copy()
        radius = self.lam * math.sqrt(self.d) * a
        zn = np.linalg.norm(z, axis=1)
        scale = np.minimum(1.0, radius / np.maximum(zn, 1e-12))
        return np.c_[u, z * scale[:, None]]


class FullHistoryMLP:
    def __init__(self, eq, hard=False, terminal_norm="corrected", seed=0):
        self.eq = eq
        self.hard = hard
        self.terminal_norm = terminal_norm
        self.seed = int(seed)

    def solve(self, n, xt, M=10):
        eq = self.eq
        T, d = eq.T, eq.d
        b = xt.shape[0]
        x, t = xt[:, :-1], xt[:, -1]
        delta = T - t

        rng = np.random.default_rng(self.seed + 100003 * n + 17 * b)
        mcg = M**n
        std = rng.standard_normal((b, mcg, d))
        dW = np.sqrt(delta)[:, None, None] * std
        X = x[:, None, :] + eq.mu * delta[:, None, None] + eq.sigma * dW
        terminal = np.concatenate(
            [X, np.full((b, mcg, 1), T)], axis=2
        ).reshape(-1, d + 1)
        gout = eq.g(terminal).reshape(b, mcg, 1)

        u = np.mean(gout, axis=1)
        if self.terminal_norm == "upstream":
            denom = np.maximum(delta[:, None], 1e-8)
        else:
            denom = np.maximum(np.sqrt(delta)[:, None], 1e-8)
        z = np.mean(gout * std, axis=1) / denom
        state = np.c_[u, z]
        if n == 0:
            return np.zeros_like(state)

        for level in range(n):
            mcf = M ** (n - level)
            rng_l = np.random.default_rng(
                self.seed + 900001 + 1009 * n + 97 * level + 13 * b
            )
            tau = rng_l.random((b, mcf))
            steps = tau * delta[:, None]
            std2 = rng_l.standard_normal((b, mcf, d))
            X2 = (
                x[:, None, :]
                + eq.mu * steps[:, :, None]
                + eq.sigma * np.sqrt(steps)[:, :, None] * std2
            )
            inter = np.concatenate(
                [X2, (t[:, None] + steps)[:, :, None]], axis=2
            )
            flat = inter.reshape(-1, d + 1)

            sim = self.solve(level, flat, M).reshape(b, mcf, d + 1)
            y = eq.f(
                flat,
                sim[:, :, 0].reshape(-1, 1),
                sim[:, :, 1:].reshape(-1, d),
            ).reshape(b, mcf, 1)
            u += delta[:, None] * np.mean(y, axis=1)
            ds = np.sqrt(steps + 1e-6)
            z += delta[:, None] * np.mean(
                y * std2 / ds[:, :, None], axis=1
            )

            if level:
                prev = self.solve(level - 1, flat, M).reshape(b, mcf, d + 1)
                yprev = eq.f(
                    flat,
                    prev[:, :, 0].reshape(-1, 1),
                    prev[:, :, 1:].reshape(-1, d),
                ).reshape(b, mcf, 1)
                u -= delta[:, None] * np.mean(yprev, axis=1)
                z -= delta[:, None] * np.mean(
                    yprev * std2 / ds[:, :, None], axis=1
                )

        out = np.c_[u, z]
        if self.hard:
            return eq.project(out, xt)
        return np.clip(out, -eq.norm, eq.norm)


def hypercube(rng, n, d, lo, hi, T, boundary=False):
    x = rng.uniform(lo, hi, (n, d))
    t = rng.uniform(0.0, T, n)
    if boundary:
        idx = rng.integers(0, d, n)
        side = rng.integers(0, 2, n)
        x[np.arange(n), idx] = np.where(side == 0, lo, hi)
    return np.c_[x, t]


def ball(rng, n, d, T, boundary=False):
    z = rng.standard_normal((n, d))
    z /= np.linalg.norm(z, axis=1, keepdims=True)
    r = np.ones(n) if boundary else rng.random(n) ** (1.0 / d)
    return np.c_[z * r[:, None], rng.uniform(0.0, T, n)]


def make_test(eq, seed, n_domain=1000, n_boundary=200):
    rng = np.random.default_rng(seed)
    if isinstance(eq, (LCD, VB)):
        a = hypercube(rng, n_domain, eq.d, -.5, .5, eq.T)
        b = hypercube(rng, n_boundary, eq.d, -.5, .5, eq.T, True)
    else:
        a = ball(rng, n_domain, eq.d, eq.T)
        b = ball(rng, n_boundary, eq.d, eq.T, True)
    return np.r_[a, b]


def rel_l2(pred, truth):
    return float(np.linalg.norm(pred - truth) / np.linalg.norm(truth))


def evaluate(eq, reps=3, batch=100, include_upstream=False):
    raw = []
    for rep in range(reps):
        xt = make_test(eq, 202600 + rep)
        truth = eq.exact(xt)
        row = {}
        variants = [
            ("baseline", False, "corrected"),
            ("hard", True, "corrected"),
        ]
        if include_upstream:
            variants += [
                ("upstream_baseline", False, "upstream"),
                ("upstream_hard", True, "upstream"),
            ]
        for label, hard, norm in variants:
            pred = []
            for i in range(0, len(xt), batch):
                solver = FullHistoryMLP(
                    eq, hard=hard, terminal_norm=norm, seed=7000 + rep
                )
                pred.append(solver.solve(2, xt[i:i + batch], 10)[:, 0])
            row[label] = rel_l2(np.concatenate(pred), truth)
        raw.append(row)

    summary = {}
    for key in raw[0]:
        vals = np.asarray([r[key] for r in raw])
        summary[key] = float(vals.mean())
        summary[key + "_std"] = (
            float(vals.std(ddof=1)) if reps > 1 else 0.0
        )
    summary["hard_reduction_pct"] = (
        100.0 * (summary["baseline"] - summary["hard"])
        / summary["baseline"]
    )
    return summary, raw


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--batch", type=int, default=100)
    p.add_argument("--output", type=Path,
                   default=Path("scasml_benchmark_completion.json"))
    p.add_argument("--include-upstream-normalization", action="store_true")
    p.add_argument("--vb-sigma-paper", action="store_true")
    args = p.parse_args()

    vb_sigma = math.sqrt(2.0) if args.vb_sigma_paper else .25
    suites = [
        ("LCD", [10, 20, 30, 60], lambda d: LCD(d)),
        ("VB", [20, 40, 60, 80], lambda d: VB(d, sigma=vb_sigma)),
        ("DR", [100, 120, 140, 160], lambda d: DR(d)),
    ]
    result = {
        "protocol": {
            "n": 2,
            "M": 10,
            "interior_points": 1000,
            "boundary_points": 200,
            "repetitions": args.reps,
            "terminal_gradient_normalization": "standard EBL",
            "vb_sigma": vb_sigma,
        },
        "benchmarks": {},
    }

    for name, dims, ctor in suites:
        result["benchmarks"][name] = {}
        for d in dims:
            print("RUN", name, d, flush=True)
            summary, raw = evaluate(
                ctor(d), args.reps, args.batch,
                args.include_upstream_normalization
            )
            summary["paper_mlp_rel_l2"], summary["paper_scasml_rel_l2"] = \
                PAPER[name][d]
            summary["raw_repetitions"] = raw
            result["benchmarks"][name][str(d)] = summary
            print(json.dumps(summary), flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("wrote", args.output)


if __name__ == "__main__":
    main()
