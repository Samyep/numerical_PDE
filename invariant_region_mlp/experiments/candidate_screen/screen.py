"""MLP-level screen of candidate PDEs (240 points, 2 paired repetitions).

Methods: raw, segment (certified), box (certified hull), oracle_state, centre
(data-free certificate centre), f_zero.  The recursion is the unchanged
mechanism-suite `MechanismMLP`; only `_project_segment` is generalised to the
new families (z = sigma c w with c clipped to `segment_coefficients`).
Usage: python screen.py OUT.jsonl [workers]
"""
from __future__ import annotations

import json
import warnings
warnings.filterwarnings("ignore")
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))

import candidates as C  # noqa: E402
from mechanism_suite import mechanism_mlp as MM  # noqa: E402
from mechanism_suite.equations import RidgeLSEHJB  # noqa: E402
from mechanism_suite.norm_hjb import NormDriverHJB  # noqa: E402

_orig_segment = MM._project_segment


def _generic_segment(equation, state):
    if equation.family in ("ridge_lse", "norm_hjb"):
        return _orig_segment(equation, state)
    out = np.asarray(state, dtype=np.float64).copy()
    w = np.asarray(equation.w)
    lo, hi = equation.segment_coefficients
    c = np.clip(np.einsum("...d,d->...", out[..., 1:], w) / equation.sigma, lo, hi)
    out[..., 1:] = equation.sigma * c[..., None] * w
    return out


MM._project_segment = _generic_segment

METHODS = [MM.MethodSpec("raw", "raw", 1.0), MM.MethodSpec("segment", "segment", 1.0),
           MM.MethodSpec("box", "box", 1.0), MM.MethodSpec("oracle_state", "oracle_state", 1.0),
           MM.MethodSpec("centre", "centre", 1.0), MM.F_ZERO]
CELLS = [(2, 32), (3, 6), (4, 6)]
DIMS = [20, 100]


def make(cfg, d):
    kind = cfg["kind"]
    if kind == "P1":
        return RidgeLSEHJB(d=d)
    if kind == "P4":
        return NormDriverHJB(d=d, reference_path=str(C.REF_P4), T=0.5)
    if kind == "C1":
        return C.L1ControlHJB(d=d)
    if kind == "C2":
        return C.LQGameHJB(d=d, a=cfg["a"], b=cfg["b"], frac_A=cfg["frac_A"])
    if kind == "C3":
        return C.CubicHJB(d=d, kappa=0.5, A=2.0, T=0.25)
    raise ValueError(kind)


CONFIGS = [
    dict(label="P1_anchor", kind="P1"),
    dict(label="P4_anchor", kind="P4"),
    dict(label="C1_l1_control", kind="C1"),
    dict(label="C2_game_convex(a=2,b=0,dA=d/2)", kind="C2", a=2.0, b=0.0, frac_A=0.5),
    dict(label="C2_game_cancel(a=3,b=1,dA=d/4)", kind="C2", a=3.0, b=1.0, frac_A=0.25),
    dict(label="C2_game_flip(a=3,b=1,dA=d/8)", kind="C2", a=3.0, b=1.0, frac_A=0.125),
    dict(label="C3_cubic(k=0.5,A=2)", kind="C3"),
]


def task(args):
    cfg, d, n, M, meth, rep = args
    eq = make(cfg, d)
    rng = np.random.default_rng(np.random.SeedSequence([20261007, d, 240]))
    t = rng.uniform(0.0, np.nextafter(eq.T, 0.0), 240); x = rng.uniform(-1, 1, (240, d))
    res = MM.run_single_repetition(pde_id=cfg["label"], equation=eq, method=meth, n=n, M=M, repetition=rep,
                                   t=t, x=x, is_validation=np.zeros(240, bool), base_seed=20261008,
                                   chunk_size=8 if d <= 20 else 4)
    md = res["metadata"]; work = md["work"]
    out = dict(label=cfg["label"], d=d, n=n, M=M, method=meth.name, rep=rep,
               skill=md["metrics"]["all"]["skill"], seconds=md["wall_clock_seconds"],
               f_calls=work.get("f_evals"),
               generator_bias=(work.get("generator") or {}).get("bias"),
               generator_mae=(work.get("generator") or {}).get("mae"))
    if cfg["kind"] == "C2":
        out["predicted_bias_factor"] = eq.predicted_orthogonal_bias_factor()
    return out


def main(out_path, workers):
    jobs = [(cfg, d, n, M, m, r) for cfg in CONFIGS for d in DIMS for (n, M) in CELLS for m in METHODS for r in range(2)]
    jobs.sort(key=lambda j: (j[1], j[2] * 10 + j[3]))
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool, open(out_path, "a") as fh:
        for row in pool.map(task, jobs, chunksize=1):
            fh.write(json.dumps(row) + "\n"); fh.flush()
    print(f"done {len(jobs)} jobs in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else max(1, os.cpu_count() - 1))
