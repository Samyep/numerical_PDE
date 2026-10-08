"""MLP screen for effective dimension k > 1 (multi-direction log-sum-exp quadratic HJB).

PDE (sigma = sqrt2, mu = 0): u_t + Lap u - |grad u|^2 = 0, i.e. f(z) = -|z|^2 / 2, z = sqrt2 grad u.
Terminal g(x) = -log mean_j exp(a_j . x), a_j = Q b_j, Q in R^{d x k} orthonormal (random), b_j in R^k.
Closed form u(t,x) = -log mean_j exp(a_j . x + |a_j|^2 (T-t)), z = -sqrt2 sum_j pi_j a_j.

Certificates (all PDE-derived; Hopf-Cole tilted expectation => grad u in conv{-a_j}):
  sub_box  : z in span(Q_hat) and, in Q_hat coordinates, inside the bounding box of {-sqrt2 b_j}
             (Q_hat is ESTIMATED from samples of grad g only; principal angles to Q are reported)
  subspace : z in span(Q_hat) only (no magnitude bound)
  box      : ambient coordinatewise hull of {-sqrt2 a_j}
  ball     : |z| <= sqrt2 max_j |a_j|
Controls: raw, oracle_state, centre (box-hull centre, data-free), f_zero.
The recursion is the unchanged mechanism-suite MechanismMLP; only _project_segment / _project_span
are redirected for this family ("segment" -> sub_box, "span_only" -> subspace).
Usage: python effdim_mlp.py OUT.jsonl [workers]
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from mechanism_suite import mechanism_mlp as MM  # noqa: E402
from mechanism_suite.equations import ExactEquation, _as_points  # noqa: E402

SQ2 = math.sqrt(2.0)


def _lse(z):
    m = z.max(-1, keepdims=True)
    return m[..., 0] + np.log(np.exp(z - m).sum(-1))


@dataclass(frozen=True)
class MultiRidgeHJB(ExactEquation):
    d: int
    k: int
    scale: float
    T: float
    seed: int = 0
    n_grad_samples: int = 4000
    sigma: float = SQ2
    name: str = "EK_multiridge_lse_hjb"
    family: str = "multiridge"

    def __post_init__(self):
        rng = np.random.default_rng(np.random.SeedSequence([20261008, self.k, self.d, self.seed]))
        B = rng.standard_normal((2 * self.k, self.k)); B /= np.linalg.norm(B, axis=1, keepdims=True)
        B *= self.scale * rng.uniform(0.5, 1.5, (2 * self.k, 1))
        Q, _ = np.linalg.qr(rng.standard_normal((self.d, self.k)))
        A = B @ Q.T                                                   # rows a_j in R^d
        object.__setattr__(self, "B", B); object.__setattr__(self, "Q", Q); object.__setattr__(self, "A", A)
        object.__setattr__(self, "mu", 0.0)
        object.__setattr__(self, "logc", -math.log(len(B)) * np.ones(len(B)))
        object.__setattr__(self, "n2", (A * A).sum(1))
        # --- subspace estimated from samples of grad g only (no solution used) ---
        xs = np.random.default_rng(99).uniform(-1, 1, (self.n_grad_samples, self.d))
        lg = self.logc + xs @ A.T
        pi = np.exp(lg - _lse(lg)[:, None])
        G = -(pi @ A)                                                 # grad g samples
        U, S, _ = np.linalg.svd(G.T, full_matrices=False)
        Qh = U[:, :self.k]
        cosines = np.linalg.svd(Q.T @ Qh, compute_uv=False)
        object.__setattr__(self, "Qh", Qh)
        object.__setattr__(self, "max_principal_angle", float(np.arccos(np.clip(cosines.min(), -1, 1))))
        object.__setattr__(self, "sv_gap", float(S[self.k - 1] / max(S[self.k], 1e-300)) if len(S) > self.k else float("inf"))
        V = -SQ2 * B @ (Q.T @ Qh)                                     # vertices -sqrt2 b_j in Q_hat coordinates
        object.__setattr__(self, "sub_lo", V.min(0)); object.__setattr__(self, "sub_hi", V.max(0))
        verts = -SQ2 * A
        object.__setattr__(self, "_box_lo", verts.min(0)); object.__setattr__(self, "_box_hi", verts.max(0))

    def _logits(self, t, x):
        time_, state = _as_points(t, x, self.d)
        return self.logc + state @ self.A.T + (self.T - time_)[..., None] * self.n2

    def exact_u(self, t, x):
        return -_lse(self._logits(t, x))

    def exact_z(self, t, x):
        lg = self._logits(t, x)
        pi = np.exp(lg - _lse(lg)[..., None])
        return -SQ2 * (pi @ self.A)

    def terminal(self, x):
        return -_lse(self.logc + np.asarray(x) @ self.A.T)

    def generator(self, u, z):
        del u
        z = np.asarray(z, dtype=np.float64)
        return -0.5 * np.sum(z * z, axis=-1)

    @property
    def box_low(self):
        return self._box_lo

    @property
    def box_high(self):
        return self._box_hi

    @property
    def z_ball_radius(self):
        return SQ2 * float(np.sqrt(self.n2.max()))

    @property
    def z_center(self):
        return 0.5 * (self._box_lo + self._box_hi)

    @property
    def dose_scale(self):
        return float(np.mean(0.5 * (self._box_hi - self._box_lo)))

    @property
    def picard_proxy(self):
        return self.z_ball_radius * self.T

    def certificate_labels(self):
        return {"sub_box": "PDE-derived Hopf-Cole hull, subspace estimated from grad g samples",
                "subspace": "span of grad g samples (translation invariance)",
                "box": "ambient coordinatewise hull of -sqrt2 a_j", "ball": "sqrt2 max |a_j|"}

    def certificate_parameters(self):
        return {"k": self.k, "scale": self.scale, "seed": self.seed,
                "max_principal_angle": self.max_principal_angle, "sv_gap": self.sv_gap}


_orig_segment, _orig_span = MM._project_segment, MM._project_span


def _segment(eq, state):
    if eq.family != "multiridge":
        return _orig_segment(eq, state)
    out = np.asarray(state, dtype=np.float64).copy()
    c = np.clip(out[..., 1:] @ eq.Qh, eq.sub_lo, eq.sub_hi)
    out[..., 1:] = c @ eq.Qh.T
    return out


def _span(eq, state):
    if eq.family != "multiridge":
        return _orig_span(eq, state)
    out = np.asarray(state, dtype=np.float64).copy()
    out[..., 1:] = (out[..., 1:] @ eq.Qh) @ eq.Qh.T
    return out


MM._project_segment, MM._project_span = _segment, _span

METHODS = [MM.MethodSpec("raw", "raw", 1.0), MM.MethodSpec("sub_box", "segment", 1.0),
           MM.MethodSpec("subspace", "span_only", 1.0), MM.MethodSpec("box", "box", 1.0),
           MM.MethodSpec("ball", "ball", 1.0), MM.MethodSpec("oracle_state", "oracle_state", 1.0),
           MM.MethodSpec("centre", "centre", 1.0), MM.F_ZERO]
SETTINGS = [dict(k=10, scale=4.0, T=0.1), dict(k=50, scale=5.0, T=0.1)]
DIMS = [100, 400]
CELLS = [(2, 32), (3, 6), (4, 6)]
NPTS = 240


def task(args):
    st, d, n, M, meth, rep = args
    eq = MultiRidgeHJB(d=d, **st)
    rng = np.random.default_rng(np.random.SeedSequence([20261008, d, st["k"], NPTS]))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), NPTS); x = rng.uniform(-1, 1, (NPTS, d))
    res = MM.run_single_repetition(pde_id=f"EK_k{st['k']}", equation=eq, method=meth, n=n, M=M, repetition=rep,
                                   t=t, x=x, is_validation=np.zeros(NPTS, bool), base_seed=20261009, chunk_size=4)
    md = res["metadata"]; w = md["work"]; g = w.get("generator") or {}
    return dict(k=st["k"], scale=st["scale"], T=st["T"], d=d, n=n, M=M, method=meth.name, rep=rep,
                skill=md["metrics"]["all"]["skill"], seconds=md["wall_clock_seconds"], f_calls=w.get("f_evals"),
                generator_bias=g.get("bias"), angle=eq.max_principal_angle,
                violation=w.get("pre_box_violation_rate"))


def main(out, workers):
    jobs = [(st, d, n, M, m, r) for st in SETTINGS for d in DIMS if st["k"] < d
            for (n, M) in CELLS for m in METHODS for r in range(2)]
    jobs.sort(key=lambda j: (j[1], j[2] * 100 + j[3]))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool, open(out, "a") as fh:
        for row in pool.map(task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
    print(f"done {len(jobs)} jobs in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else max(1, (os.cpu_count() or 2)))
