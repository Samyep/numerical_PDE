"""Probe: data-driven spectral projection of Monte Carlo gradient estimates (no PDE certificate).

Inside each solver instance, every gradient estimate handed to f is accumulated into an uncentred
second-moment matrix S. At each generator call the signal subspace is the span of eigenvectors of
S/n whose eigenvalues exceed tau * (1 + sqrt(d/n))^2 * median(eigenvalues) (Marchenko-Pastur edge with
the noise level estimated by the median eigenvalue; tau = 1.5 safety factor). Each estimate is then
orthogonally projected onto that subspace. No knowledge of k, of the PDE, or of g is used.
Usage: python spectral_probe.py OUT.jsonl [workers]
"""
import json, sys, time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import numpy as np
import effdim_mlp as E
MM = E.MM


@dataclass(frozen=True)
class SpectralEq(E.MultiRidgeHJB):
    family: str = "multiridge_spectral"


_prev_span = MM._project_span
_STATE = {}


def _spectral(eq, state):
    if eq.family != "multiridge_spectral":
        return _prev_span(eq, state)
    out = np.asarray(state, dtype=np.float64).copy()
    Z = out[..., 1:].reshape(-1, eq.d)
    key = id(_STATE.get("solver"))
    acc = _STATE.setdefault(("acc", key), [np.zeros((eq.d, eq.d)), 0])
    acc[0] += Z.T @ Z; acc[1] += len(Z)
    S, n = acc[0] / acc[1], acc[1]
    lam, V = np.linalg.eigh(S)
    edge = 1.5 * (1 + np.sqrt(eq.d / n)) ** 2 * np.median(lam)
    U = V[:, lam > edge]
    _STATE.setdefault(("kept", key), []).append(U.shape[1])
    out[..., 1:] = (out[..., 1:] @ U) @ U.T
    return out


MM._project_span = _spectral
_orig_init = MM.MechanismMLP.__init__


def _init(self, *a, **k):
    _orig_init(self, *a, **k)
    _STATE["solver"] = self


MM.MechanismMLP.__init__ = _init

METHODS = [MM.MethodSpec("raw", "raw", 1.0), MM.MethodSpec("spectral", "span_only", 1.0),
           MM.MethodSpec("oracle_state", "oracle_state", 1.0), MM.MethodSpec("centre", "centre", 1.0)]


def task(args):
    st, d, n, M, meth, rep = args
    eq = SpectralEq(d=d, **st)
    rng = np.random.default_rng(np.random.SeedSequence([20261008, d, st["k"], E.NPTS]))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), E.NPTS); x = rng.uniform(-1, 1, (E.NPTS, d))
    _STATE.clear()
    res = MM.run_single_repetition(pde_id=f"EK_k{st['k']}", equation=eq, method=meth, n=n, M=M, repetition=rep,
                                   t=t, x=x, is_validation=np.zeros(E.NPTS, bool), base_seed=20261009, chunk_size=4)
    kept = [v for key, val in _STATE.items() if key[0] == "kept" for v in val]
    md = res["metadata"]; g = md["work"].get("generator") or {}
    return dict(k=st["k"], d=d, n=n, M=M, method=meth.name, rep=rep, skill=md["metrics"]["all"]["skill"],
                generator_bias=g.get("bias"), kept_median=float(np.median(kept)) if kept else None,
                kept_min=int(min(kept)) if kept else None, kept_max=int(max(kept)) if kept else None)


if __name__ == "__main__":
    jobs = [(st, 100, n, M, m, r) for st in E.SETTINGS[:1] + [dict(k=20, scale=4.0, T=0.1)]
            for (n, M) in [(2, 32), (3, 6)] for m in METHODS for r in range(2)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=int(sys.argv[2]) if len(sys.argv) > 2 else 2) as pool, open(sys.argv[1], "w") as fh:
        for row in pool.map(task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print(f"done {len(jobs)} in {time.time() - t0:.0f}s")
