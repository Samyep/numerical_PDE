"""Extension of screen_effdim.py to larger effective dimension k, plus a seed-stability check.

Families: rand (K = 2k random directions), rand_min (K = k + 1 random directions: the minimum
for rank k), simplex (K = k + 1). Grid: k in {10, 15, 20, 30, 50}, scale s in {2, 3, 4, 5},
T in {0.05, 0.1, 0.25}. Same diagnostics and the same pass rule as screen_effdim.py
(G2 >= 0.15 and LT <= 1.5). Configurations that pass are re-screened with 4 further
direction seeds (seed stability).
"""
import json
import math
import sys

import numpy as np

import screen_effdim as base

base.S = 3000
_orig_family = base.family


def family(name, k, rng):
    if name == "rand_min":
        B = rng.standard_normal((k + 1, k)); B /= np.linalg.norm(B, axis=1, keepdims=True)
        return B * rng.uniform(0.5, 1.5, (k + 1, 1))
    return _orig_family(name, k, rng)


base.family = family


if __name__ == "__main__":
    out = sys.argv[1]
    rows, passed = [], []
    for name in ("rand", "rand_min", "simplex"):
        for k in (10, 15, 20, 30, 50):
            for s in (2.0, 3.0, 4.0, 5.0):
                for T in (0.05, 0.1, 0.25):
                    r = base.screen(name, k, s, T); r["seed"] = 0; rows.append(r)
                    ok = r["G2"] >= 0.15 and r["LT"] <= 1.5
                    if ok:
                        passed.append((name, k, s, T))
                    print(f"{name:8s} k={k:2d} s={s:3.1f} T={T:4.2f}  S_NL={r['S_NL']:.3f}  G2={r['G2']:.3f}  "
                          f"CV={r['CV_gen']:.2f}  LT={r['LT']:.2f}  {'PASS' if ok else ''}", flush=True)
    print("\n== seed stability of passing configurations (seeds 1-4)")
    for name, k, s, T in passed:
        g2 = []
        for seed in (1, 2, 3, 4):
            r = base.screen(name, k, s, T, seed=seed); r["seed"] = seed; rows.append(r); g2.append(r["G2"])
        print(f"{name:8s} k={k:2d} s={s:3.1f} T={T:4.2f}  G2 over seeds: {np.round(g2, 3).tolist()}  min={min(g2):.3f}", flush=True)
    json.dump(rows, open(out, "w"), indent=1)
