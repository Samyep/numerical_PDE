"""Probe: surrogate control variates (Stein / first-order) inside the MLP gradient estimators.

Question: if a network supplied an approximate gradient field z_th ~ sigma grad u and f-field
F_th ~ f(u, z) with gradient, how much of the noise-into-f error would disappear, as a function of
the surrogate's relative error eps?  No network is trained here: the "surrogate" is the exact
solution plus a smooth random perturbation of controlled relative size eps (per point).

Estimators (Gaussian integration by parts: E[(a . G) G] = a, so every control variate has known mean
and the estimator stays unbiased for ANY surrogate):
  terminal : z_hat = mean_k (Delta_k / sqrt h - z_th . G_k) G_k + z_th
  level 1  : z_inc = mean_k I_k (D_k - F_th - sigma sqrt(dt_k) gF_th . G_k) G_k / sqrt(dt_k) + mean_k I_k sigma gF_th
             (D_k = f(U_1)(R_k, X_k); F_th, gF_th evaluated at (R_k, x); levels > 1 unchanged)
Value estimates are left unchanged so that only the z channel is affected.
Surrogates: cv_eX (oracle + perturbation of relative size X), cv_gradg (data-free: the terminal
condition used as surrogate, i.e. the closed form evaluated with T - t = 0).
Usage: python cv_probe.py OUT.jsonl [workers]
"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import effdim_mlp as E

MM = E.MM
SQ2 = math.sqrt(2.0)


def _fields(eq, t, x, at_terminal=False):
    """z*, F* = f(z*), grad_x F* for the multi-ridge HJB (closed form)."""
    h = np.zeros_like(t) if at_terminal else (eq.T - t)
    lg = eq.logc + x @ eq.A.T + h[..., None] * eq.n2
    pi = np.exp(lg - E._lse(lg)[..., None])
    z = -SQ2 * (pi @ eq.A)
    F = -0.5 * np.sum(z * z, -1)
    Az = z @ eq.A.T
    w = pi * Az - pi * np.sum(pi * Az, -1, keepdims=True)
    gF = SQ2 * (w @ eq.A)
    return z, F, gF


def _pert(eq, x, which):
    rng = np.random.default_rng([777, eq.d, which])
    R = rng.standard_normal((eq.d, 64)) * (3.0 / math.sqrt(eq.d))
    P = rng.standard_normal((64, eq.d))
    N = np.tanh(x @ R + rng.uniform(-1, 1, 64)) @ P
    return N / np.maximum(np.linalg.norm(N, axis=-1, keepdims=True), 1e-300)


def surrogate(eq, t, x, mode):
    if mode == "gradg":
        return _fields(eq, t, x, at_terminal=True)
    eps = float(mode)
    z, F, gF = _fields(eq, t, x)
    if eps == 0.0:
        return z, F, gF
    z = z + eps * np.linalg.norm(z, axis=-1, keepdims=True) * _pert(eq, x, 1)
    gF = gF + eps * np.linalg.norm(gF, axis=-1, keepdims=True) * _pert(eq, x, 2)
    F = F * (1.0 + eps * np.tanh(x[..., 0]))
    return z, F, gF


MODE = {}


class CVMLP(MM.MechanismMLP):
    @property
    def _mode(self):
        return MODE.get(self.method.name)

    def _terminal_estimate(self, n, t, x):
        if self._mode is None:
            return super()._terminal_estimate(n, t, x)
        eq = self.equation
        b, d = len(t), eq.d
        h = np.maximum(eq.T - t, 1e-300)
        normal = self._normal((b, self.M ** n, d))
        gx = eq.terminal(x)
        xt = x[:, None, :] + eq.sigma * np.sqrt(h)[:, None, None] * normal
        diff = eq.terminal(xt) - gx[:, None]
        self.stats.terminal_g_evals += b * (1 + self.M ** n)
        self.stats.terminal_samples += b * self.M ** n
        out = np.zeros((b, d + 1))
        out[:, 0] = gx + diff.mean(1)
        term, _ = self._mode
        if term == "path":   # Stein with phi = g itself: sigma * mean grad g(X_k)  (pathwise estimator)
            zX, _, _ = _fields(eq, np.full(xt.shape[:2], eq.T), xt, at_terminal=True)
            out[:, 1:] = zX.mean(1)
        elif term == "bismut":
            out[:, 1:] = np.mean((diff / np.sqrt(h)[:, None])[:, :, None] * normal, 1)
        else:
            raise ValueError(term)
        return out

    def solve(self, n, t, x, *, collect_root=False):
        if self._mode is None:
            return super().solve(n, t, x, collect_root=collect_root)
        t = np.asarray(t, dtype=np.float64).reshape(-1)
        eq = self.equation
        x = np.asarray(x, dtype=np.float64).reshape(len(t), eq.d)
        b = len(t)
        self.stats.recursively_evaluated_states += b
        if n <= 0:
            return np.zeros((b, eq.d + 1))
        output = self._terminal_estimate(n, t, x)
        remaining = eq.T - t
        alpha = self.time_beta_alpha
        for level in range(1, n):
            sib = self.M ** (n - level)
            r = self._power_time((b, sib))
            dt = remaining[:, None] * r
            normal = self._normal((b, sib, eq.d))
            xr = x[:, None, :] + eq.sigma * np.sqrt(dt)[:, :, None] * normal
            tr = t[:, None] + dt
            high = self.solve(level, tr.reshape(-1), xr.reshape(-1, eq.d)).reshape(b, sib, eq.d + 1)
            diff = self._generator(high, tr, xr)
            if level > 1:
                low = self.solve(level - 1, tr.reshape(-1), xr.reshape(-1, eq.d)).reshape(b, sib, eq.d + 1)
                diff = diff - self._generator(low, tr, xr)
            imp = remaining[:, None] * r ** (1.0 - alpha) / alpha
            sdt = np.sqrt(np.maximum(dt, 1e-300))
            output[:, 0] += np.mean(imp * diff, 1)
            sur = self._mode[1]
            if level == 1 and sur is not None:
                # full Stein control variate phi(G) = F_th(R, x + sigma sqrt(dt) G): mean of phi(G) G / sqrt(dt)
                # equals sigma * E grad F_th(R, X) by Gaussian integration by parts
                _, Fs, gFs = surrogate(eq, tr, xr, sur)
                output[:, 1:] += np.mean((imp * (diff - Fs) / sdt)[:, :, None] * normal, 1) \
                    + np.mean(imp[:, :, None] * eq.sigma * gFs, 1)
            else:
                output[:, 1:] += np.mean((imp * diff / sdt)[:, :, None] * normal, 1)
        return output


MM.MechanismMLP = CVMLP

MODE.update({"path": ("path", None), "bismut_cv_0": ("bismut", "0"), "bismut_cv_0.3": ("bismut", "0.3")})
for m in ["0", "0.1", "0.3", "1", "gradg"]:
    MODE[f"path_cv_{m}"] = ("path", m)
METHODS = [MM.MethodSpec("raw", "raw", 1.0), MM.MethodSpec("oracle_state", "oracle_state", 1.0),
           MM.MethodSpec("centre", "centre", 1.0)] + [MM.MethodSpec(m, "raw", 1.0) for m in MODE]


def task(args):
    st, d, n, M, meth, rep = args
    eq = E.MultiRidgeHJB(d=d, **st)
    rng = np.random.default_rng(np.random.SeedSequence([20261008, d, st["k"], E.NPTS]))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), E.NPTS); x = rng.uniform(-1, 1, (E.NPTS, d))
    res = MM.run_single_repetition(pde_id=f"EK_k{st['k']}", equation=eq, method=meth, n=n, M=M, repetition=rep,
                                   t=t, x=x, is_validation=np.zeros(E.NPTS, bool), base_seed=20261009, chunk_size=4)
    md = res["metadata"]; g = md["work"].get("generator") or {}
    return dict(k=st["k"], d=d, n=n, M=M, method=meth.name, rep=rep, skill=md["metrics"]["all"]["skill"],
                generator_bias=g.get("bias"), seconds=md["wall_clock_seconds"], f_calls=md["work"].get("f_evals"))


if __name__ == "__main__":
    settings = [dict(k=10, scale=4.0, T=0.1), dict(k=50, scale=5.0, T=0.1)]
    jobs = [(st, d, n, M, m, r) for st in settings for d in (100, 400) for (n, M) in [(2, 32), (3, 6)]
            for m in METHODS for r in range(2)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(sys.argv[2]) if len(sys.argv) > 2 else 2) as pool, \
            open(sys.argv[1], "w") as fh:
        for row in pool.map(task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print(f"done {len(jobs)} in {time.time() - t0:.0f}s")
