"""Analytic screen for benchmarks with effective dimension k > 1 (no MLP).

PDE (sigma = sqrt2, mu = 0):  u_t + Lap u - |grad u|^2 = 0,  u(T, x) = g(x),
    g(x) = -log sum_j c_j exp(a_j . x),   a_j = Q b_j,  Q in R^{d x k} orthonormal, b_j in R^k.
Hopf-Cole closed form:  u(t, x) = -log sum_j c_j exp(a_j . x + |a_j|^2 (T - t)),
    grad u = -sum_j pi_j a_j   (softmax weights), so z = sqrt2 grad u lies in -sqrt2 conv{a_j}:
    a PDE-derived k-dimensional polytope certificate. The solution is not separable and has no
    symmetry that reduces it below k dimensions for generic b_j.

Families of {b_j} in R^k (K vectors):
  rand  : K = 2k Gaussian directions, norms uniform in [0.5, 1.5] (times scale s)
  cross : the 2k vectors +-e_i (times s)  -- sum of exponentials, not a product: non-separable
  simplex: K = k+1 centred simplex vertices (times s)
Diagnostics on 800 test points t ~ U[0,T), x ~ U[-1,1]^d (d = 100):
  S_NL  = skill of the f=0 solution E g(x + sqrt2 W_h)       (size of the nonlinearity)
  G2    = skill after the best time-only source a h + b h^2   (spatial dependence; gate >= 0.15)
  CV    = coefficient of variation of the exact generator |grad u|^2 over test points
  LT    = sqrt2 max_j |a_j| * T                               (Picard-difficulty proxy; guide <~ 1.5)
Skill = RMSE / std(u*).
"""
import json
import math
import sys

import numpy as np

D = 100
J = 800
S = 6000


def family(name, k, rng):
    if name == "rand":
        B = rng.standard_normal((2 * k, k)); B /= np.linalg.norm(B, axis=1, keepdims=True)
        return B * rng.uniform(0.5, 1.5, (2 * k, 1))
    if name == "cross":
        return np.concatenate([np.eye(k), -np.eye(k)])
    if name == "simplex":
        V = np.eye(k + 1) - 1.0 / (k + 1)                     # k+1 points in a k-dim hyperplane of R^{k+1}
        U, _, _ = np.linalg.svd(V.T, full_matrices=False)     # orthonormal coordinates of that hyperplane
        P = V @ U[:, :k]
        return P / np.linalg.norm(P, axis=1, keepdims=True)
    raise ValueError(name)


def lse(Z):
    m = Z.max(-1, keepdims=True)
    return m[..., 0] + np.log(np.exp(Z - m).sum(-1))


def screen(name, k, s, T, seed=0):
    rng = np.random.default_rng(np.random.SeedSequence([seed, k, {"rand": 1, "cross": 2, "simplex": 3, "rand_min": 4}[name]]))
    B = s * family(name, k, rng)                              # (K, k)
    K = len(B); logc = -math.log(K) * np.ones(K)
    Q, _ = np.linalg.qr(rng.standard_normal((D, k)))
    prng = np.random.default_rng(7)
    x = prng.uniform(-1, 1, (J, D)); t = prng.uniform(0, T, J); h = T - t
    y = x @ Q                                                 # coordinates in the active subspace
    n2 = (B * B).sum(1)
    logits = logc + y @ B.T + h[:, None] * n2
    u = -lse(logits)
    pi = np.exp(logits - lse(logits)[:, None])
    grad_sub = -(pi @ B)                                      # grad u in subspace coordinates
    gen = (grad_sub ** 2).sum(1)                              # |grad u|^2 (= -f with z = sqrt2 grad u)
    # f = 0 solution: E g(x + sqrt2 W_h); only the k subspace coordinates matter
    G = np.random.default_rng(11).standard_normal((S, k))
    f0 = np.empty(J)
    for i in range(J):
        Y = y[i] + math.sqrt(2 * h[i]) * G
        f0[i] = (-lse(logc + Y @ B.T)).mean()
    sd = u.std(); sk = lambda p: float(np.sqrt(np.mean((p - u) ** 2)) / sd)
    A = np.stack([h, h ** 2], 1)
    lin = f0 + A @ np.linalg.lstsq(A, u - f0, rcond=None)[0]
    return dict(family=name, k=k, scale=s, T=T, K=K, S_NL=sk(f0), G2=sk(lin),
                CV_gen=float(gen.std() / gen.mean()), LT=float(math.sqrt(2) * math.sqrt(n2.max()) * T))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "effdim_screen.json"
    rows = []
    for name in ("rand", "cross", "simplex"):
        for k in (3, 5, 8, 10):
            for s in (0.5, 1.0, 1.5, 2.0, 3.0):
                for T in (0.1, 0.25, 0.5):
                    r = screen(name, k, s, T); rows.append(r)
                    ok = "PASS" if (r["G2"] >= 0.15 and r["LT"] <= 1.5) else ""
                    print(f"{name:7s} k={k:2d} s={s:3.1f} T={T:4.2f}  S_NL={r['S_NL']:.3f}  G2={r['G2']:.3f}  "
                          f"CV={r['CV_gen']:.2f}  LT={r['LT']:.2f}  {ok}", flush=True)
    json.dump(rows, open(out, "w"), indent=1)
