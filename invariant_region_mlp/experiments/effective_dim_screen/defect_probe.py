"""Probe: SCaSML-style defect correction on the multi-ridge HJB with a controlled surrogate.

Surrogate u_th = u* + eps * psi(x), psi a smooth random tanh field (32 ridge units in random ambient
directions, so its gradient error is NOT confined to the active subspace, as for a trained network),
normalised so that RMS sigma|grad psi| = RMS |z*| on the test points: eps = relative gradient error.
Defect e = u - u_th solves  e_t + Lap e + f_e(t, x, z_e) = 0,  e(T) = g - u_th(T) = -eps psi,
    f_e = f(z_th + z_e) - f(z_th) + r_th,   r_th = d_t u_th + Lap u_th + f(z_th) = eps Lap psi - f(z*) + f(z_th)
(r_th is the surrogate's PDE residual; SCaSML computes it by autodiff, here it is closed form).
The unchanged MLP recursion is run on the defect equation; prediction u = u_th + e_hat.
Terminal gradient: Bismut (as in the base code) or pathwise sigma grad e(T, X).
Skill is reported in units of std(u*) so it is comparable with the plain-MLP rows.
Usage: python defect_probe.py OUT.jsonl [workers]
"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

import cv_probe as C

E, MM = C.E, C.MM
SQ2 = math.sqrt(2.0)


@dataclass(frozen=True)
class DefectEq(E.MultiRidgeHJB):
    eps: float = 0.1
    family: str = "defect"
    name: str = "EK_multiridge_defect"

    def __post_init__(self):
        super().__post_init__()
        rng = np.random.default_rng([4242, self.d, self.k])
        R = rng.standard_normal((32, self.d)); R *= 1.5 / np.linalg.norm(R, axis=1, keepdims=True)
        bb = rng.uniform(-1, 1, 32); c = rng.standard_normal(32)
        object.__setattr__(self, "R", R); object.__setattr__(self, "bb", bb); object.__setattr__(self, "c", c)
        xs = np.random.default_rng(5).uniform(-1, 1, (2000, self.d))
        ts = np.random.default_rng(6).uniform(0, self.T, 2000)
        zs = super().exact_z(ts, xs)
        gp = self.sigma * self._gpsi(xs)
        object.__setattr__(self, "c", c * np.sqrt(np.mean(np.sum(zs ** 2, -1)) / np.mean(np.sum(gp ** 2, -1))))

    def _psi(self, x):
        return np.tanh(x @ self.R.T + self.bb) @ self.c

    def _gpsi(self, x):
        s = 1.0 / np.cosh(x @ self.R.T + self.bb) ** 2
        return (s * self.c) @ self.R

    def _lap_psi(self, x):
        a = x @ self.R.T + self.bb
        th = np.tanh(a); s = 1.0 - th ** 2
        return (-2.0 * th * s * self.c * np.sum(self.R ** 2, 1)) @ np.ones(len(self.c))

    # defect solution and data
    def exact_u(self, t, x):
        _, xx = E._as_points(t, x, self.d)
        return -self.eps * self._psi(xx)

    def exact_z(self, t, x):
        _, xx = E._as_points(t, x, self.d)
        return -self.eps * self.sigma * self._gpsi(xx)

    def terminal(self, x):
        return -self.eps * self._psi(np.asarray(x))

    def sigma_grad_terminal(self, x):
        return -self.eps * self.sigma * self._gpsi(x)

    def u_star(self, t, x):
        return super().exact_u(t, x)

    def z_star(self, t, x):
        return super().exact_z(t, x)

    def defect_generator(self, t, x, z_e):
        zs = self.z_star(t, x)
        zth = zs + self.eps * self.sigma * self._gpsi(x)
        f = lambda z: -0.5 * np.sum(z * z, -1)
        return f(zth + z_e) + self.eps * self._lap_psi(x) - f(zs)


class DefectMLP(C.CVMLP):
    def _terminal_estimate(self, n, t, x):
        eq = self.equation
        if not isinstance(eq, DefectEq) or self._mode is None or self._mode[0] != "path":
            return super()._terminal_estimate(n, t, x)
        b, d = len(t), eq.d
        h = np.maximum(eq.T - t, 1e-300)
        normal = self._normal((b, self.M ** n, d))
        xt = x[:, None, :] + eq.sigma * np.sqrt(h)[:, None, None] * normal
        out = np.zeros((b, d + 1))
        out[:, 0] = eq.terminal(xt).mean(1)
        out[:, 1:] = eq.sigma_grad_terminal(xt).mean(1)
        self.stats.terminal_samples += b * self.M ** n
        return out

    def _generator(self, state, t, x):
        eq = self.equation
        if not isinstance(eq, DefectEq):
            return super()._generator(state, t, x)
        exact = eq.exact_state(t, x)
        corrected = self._transform(state, exact)
        estimate = eq.defect_generator(t, x, corrected[..., 1:])
        truth = eq.eps * eq._lap_psi(x)
        count = state.shape[0] * state.shape[1]
        self.stats.f_evals += count
        self.stats.generator_count += count
        self.stats.generator_error_sum += float(np.sum(estimate - truth))
        return estimate


MM.MechanismMLP = DefectMLP
C.MODE.update({"defect_bismut": ("bismut", None), "defect_path": ("path", None)})
DEFECT_METHODS = [MM.MethodSpec("defect_bismut", "raw", 1.0), MM.MethodSpec("defect_path", "raw", 1.0)]
EPS = [0.03, 0.1, 0.3]


def task(args):
    st, d, n, M, meth, rep, eps = args
    if eps is None:
        return C.task((st, d, n, M, meth, rep))
    eq = DefectEq(d=d, eps=eps, **st)
    rng = np.random.default_rng(np.random.SeedSequence([20261008, d, st["k"], E.NPTS]))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), E.NPTS); x = rng.uniform(-1, 1, (E.NPTS, d))
    res = MM.run_single_repetition(pde_id=f"EK_k{st['k']}", equation=eq, method=meth, n=n, M=M, repetition=rep,
                                   t=t, x=x, is_validation=np.zeros(E.NPTS, bool), base_seed=20261009, chunk_size=4)
    pred = res["prediction"] if "prediction" in res else None
    md = res["metadata"]
    sk_e = md["metrics"]["all"]["skill"]
    sd_e, sd_u = np.std(eq.exact_u(t, x)), np.std(eq.u_star(t, x))
    surrogate_skill = float(np.sqrt(np.mean((eq.eps * eq._psi(x)) ** 2)) / sd_u)
    return dict(k=st["k"], d=d, n=n, M=M, method=f"{meth.name}_e{eps}", rep=rep, skill=sk_e * sd_e / sd_u,
                surrogate_skill=surrogate_skill, seconds=md["wall_clock_seconds"], f_calls=md["work"].get("f_evals"))


if __name__ == "__main__":
    settings = [dict(k=10, scale=4.0, T=0.1), dict(k=50, scale=5.0, T=0.1)]
    base = [m for m in C.METHODS if m.name in ("raw", "oracle_state", "centre", "path", "path_cv_0", "path_cv_gradg")]
    jobs = []
    for st in settings:
        for d in (100, 400):
            for (n, M) in [(2, 32), (3, 6)]:
                for r in range(2):
                    jobs += [(st, d, n, M, m, r, None) for m in base]
                    jobs += [(st, d, n, M, m, r, e) for m in DEFECT_METHODS for e in EPS]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(sys.argv[2]) if len(sys.argv) > 2 else 2) as pool, \
            open(sys.argv[1], "w") as fh:
        for row in pool.map(task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print(f"done {len(jobs)} in {time.time() - t0:.0f}s")
