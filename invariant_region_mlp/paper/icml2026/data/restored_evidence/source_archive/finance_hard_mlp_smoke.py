"""
Standalone smoke test for hard-constrained MLP on the 100D nonlinear
option-pricing PDE with different borrowing/lending rates.

This script follows the full-history gradient-dependent MLP structure used in:
  E, Hutzenthaler, Jentzen, Kruse (2019)
  "On multilevel Picard numerical approximations..."
and the public philipschw/mlpSolver implementation.

It is intentionally lightweight: NumPy only, single-process, low-budget.
It is NOT a drop-in patch of the public repository.

Main benchmark parameters:
  d=100, T=0.5, sigma=0.2, R_l=0.04, R_b=0.06, mu=0.06
  x0=(100,...,100)
  g(x) = [max(x)-120]_+ - 2[max(x)-150]_+
Paper reference at x0: approximately 21.299.

Hard constraint:
  Delta = grad_x u
  ||Delta(t,x)||_2 <= exp(0.5*sigma^2*(T-t))
because g is 1-Lipschitz in Euclidean norm and the controlled
Black-Scholes HJB preserves a corresponding Lipschitz envelope.

Since the MLP gradient state is z_i = sigma*x_i*Delta_i, enforce
  sum_i (z_i/(sigma*x_i))^2 <= exp(sigma^2*(T-t)).
"""

import math
import numpy as np


class FundingMaxSpread:
    def __init__(
        self,
        d=100,
        T=0.5,
        sigma=0.2,
        mu=0.06,
        Rl=0.04,
        Rb=0.06,
    ):
        self.d = d
        self.T = T
        self.sigma = sigma
        self.mu = mu
        self.Rl = Rl
        self.Rb = Rb

    def g(self, x):
        m = float(np.max(x))
        return max(m - 120.0, 0.0) - 2.0 * max(m - 150.0, 0.0)

    def f(self, t, x, v):
        y = float(v[0])
        z = np.asarray(v[1:])
        sum_z = float(np.sum(z))
        return (
            -self.Rl * y
            - ((self.mu - self.Rl) / self.sigma) * sum_z
            + (self.Rb - self.Rl) * max(sum_z / self.sigma - y, 0.0)
        )

    def sample(self, t, s, x, rng):
        dt = s - t
        W = rng.normal(0.0, math.sqrt(dt), size=self.d)
        X = x * np.exp(
            (self.mu - 0.5 * self.sigma**2) * dt + self.sigma * W
        )
        dI = np.concatenate(([1.0], W / dt))
        return X, dI

    def delta_radius(self, t):
        # L_g = 1 for the max-call-spread payoff.
        return math.exp(0.5 * self.sigma**2 * (self.T - t))


def seeded_rng(seed, path):
    vals = [int(seed)] + [int(v) & 0xFFFFFFFF for v in path]
    return np.random.default_rng(np.random.SeedSequence(vals))


class FullHistoryMLP:
    def __init__(
        self,
        equation,
        M=3,
        time_dist_exponent=0.5,
        hard_constraint=False,
        radius_factor=1.0,
        seed=0,
    ):
        self.eq = equation
        self.M = M
        self.alpha = time_dist_exponent
        self.hard_constraint = hard_constraint
        self.radius_factor = radius_factor
        self.seed = seed

    def project(self, v, x, t):
        if not self.hard_constraint:
            return v

        y = float(v[0])
        z = np.asarray(v[1:], dtype=float).copy()

        # MLP state: z_i = sigma * x_i * Delta_i.
        delta = z / (self.eq.sigma * np.maximum(x, 1e-12))

        radius = self.radius_factor * self.eq.delta_radius(t)
        norm = np.linalg.norm(delta)

        if norm > radius:
            delta *= radius / norm

        z = self.eq.sigma * x * delta
        return np.concatenate(([y], z))

    def solve(self, n, t, x, path=(0,)):
        d = self.eq.d
        if n == 0:
            return np.zeros(d + 1)

        M = self.M
        T = self.eq.T
        a = self.alpha

        gx = self.eq.g(x)

        # Terminal contribution.
        rhs_g = np.zeros(d + 1)
        for i in range(M**n):
            rng = seeded_rng(self.seed, path + (11, n, i))
            XT, dI = self.eq.sample(t, T, x, rng)
            rhs_g += (self.eq.g(XT) - gx) * dI

        rhs = np.concatenate(([gx], np.zeros(d))) + rhs_g / (M**n)

        # Multilevel nonlinear contributions.
        for level in range(n):
            N = M ** (n - level)
            rhs_f = np.zeros(d + 1)

            for i in range(N):
                rng = seeded_rng(self.seed, path + (21, n, level, i))
                r = rng.power(a)
                R = t + (T - t) * r

                XR, dI = self.eq.sample(t, R, x, rng)

                v_level = self.solve(
                    level,
                    R,
                    XR,
                    path + (100 + level, i, 1),
                )
                f_level = self.eq.f(R, XR, v_level)

                if level == 0:
                    f_previous = 0.0
                else:
                    v_previous = self.solve(
                        level - 1,
                        R,
                        XR,
                        path + (100 + level, i, 2),
                    )
                    f_previous = self.eq.f(R, XR, v_previous)

                rhs_f += (
                    r ** (1.0 - a)
                    * (f_level - f_previous)
                    * dI
                )

            rhs += (T - t) * rhs_f / (a * N)

        return self.project(rhs, x, t)


def summarize(values, reference):
    values = np.asarray(values)
    errors = np.abs(values - reference)
    return {
        "mean": float(np.mean(values)),
        "mae": float(np.mean(errors)),
        "median_abs_error": float(np.median(errors)),
        "std": float(np.std(values)),
    }


def run(
    d=100,
    n=3,
    M=3,
    seeds=50,
    reference=21.299,
    radius_factor=1.0,
):
    eq = FundingMaxSpread(d=d)
    x0 = np.full(d, 100.0)

    results = {}

    for name, hard in [("baseline", False), ("hard", True)]:
        values = []
        delta_norms = []

        for seed in range(seeds):
            solver = FullHistoryMLP(
                eq,
                M=M,
                hard_constraint=hard,
                radius_factor=radius_factor,
                seed=seed,
            )
            out = solver.solve(n, 0.0, x0)
            values.append(out[0])

            delta = out[1:] / (eq.sigma * x0)
            delta_norms.append(np.linalg.norm(delta))

        stats = summarize(values, reference)
        stats["median_delta_l2"] = float(np.median(delta_norms))
        stats["p90_delta_l2"] = float(np.quantile(delta_norms, 0.90))
        results[name] = stats

    print(f"d={d}, n={n}, M={M}, seeds={seeds}")
    print(f"paper reference ~ {reference:.3f}")
    print(f"certified delta radius at t=0: {eq.delta_radius(0.0):.6f}")
    print()

    for name, stats in results.items():
        print(name)
        for key, val in stats.items():
            print(f"  {key}: {val:.6f}")
        print()

    base = results["baseline"]
    hard = results["hard"]
    print(
        "MAE reduction: "
        f"{100.0 * (base['mae'] - hard['mae']) / base['mae']:.2f}%"
    )
    print(
        "Std reduction: "
        f"{100.0 * (base['std'] - hard['std']) / base['std']:.2f}%"
    )


if __name__ == "__main__":
    run()
