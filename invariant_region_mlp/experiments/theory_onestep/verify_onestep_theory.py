"""Monte Carlo check of the one-step gradient-noise theory
(docs/theory/onestep_gradient_noise.tex).

Setting: ridge terminal g(x) = phi(w.x), |w| = 1, forward X_T = x + sigma*sqrt(h)*G (no drift here),
centred terminal-block gradient estimator
    z_hat = (1/M) sum_k [g(x + sigma sqrt(h) G_k) - g(x)] G_k / sqrt(h).
Checked statements (brute force in R^d, no use of the lemma's conditional law):
  L1(iv)  E||z_perp||^2 = (d-1) c_h / M,                c_h = E[Delta^2]/h
  P2      E f(z_hat) - f(E z_hat) = -(lam/2)(Var z_par + (d-1) c_h / M)      for f = -(lam/2)||z||^2
  P3      lam*E sqrtQ*(d-1)/sqrt((d+1)M) <= lam*E||z_hat|| <= lam*sqrt(E||z_hat||^2)   for f = -lam||z||
  P4      law of f(Pi_C z_hat) does not depend on d when C is inside span(w)

Problems: P1 (ridge LSE, lam=(0.5,1.5,3), sigma=sqrt2, f=-|z|^2/2, certificate segment)
          P4 (log cosh(beta s)/beta, beta=2, sigma=sqrt2, f=-(1/sqrt2)|z|, certificate segment |c|<=1).
Usage: python verify_onestep_theory.py [out_dir]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

SQ2 = math.sqrt(2.0)
XG, WG = np.polynomial.hermite_e.hermegauss(120)
WG = WG / WG.sum()


def phi_p1(s):
    lam = np.array([0.5, 1.5, 3.0])
    a = np.log(1.0 / 3.0) + lam * np.asarray(s)[..., None]
    m = a.max(-1, keepdims=True)
    return -(m[..., 0] + np.log(np.exp(a - m).sum(-1)))


def phi_p4(s, beta=2.0):
    s = np.asarray(s)
    return (np.abs(beta * s) + np.log1p(np.exp(-2 * np.abs(beta * s))) - math.log(2.0)) / beta


PROBLEMS = {
    # name: (phi, sigma, f(norm^2, norm), Lipschitz-in-c generator on the certificate, certificate interval for c in z = c w)
    'P1_quadratic': dict(phi=phi_p1, sigma=SQ2, kind='quadratic', lam=1.0, interval=(-SQ2 * 3.0, -SQ2 * 0.5)),
    'P4_norm': dict(phi=phi_p4, sigma=SQ2, kind='norm', lam=1.0 / SQ2, interval=(-SQ2, SQ2)),
}


def f_of(kind, lam, z):
    n2 = np.sum(z * z, axis=-1)
    return -0.5 * lam * n2 if kind == 'quadratic' else -lam * np.sqrt(n2)


def run_case(name, spec, d, M, s0=0.3, h=0.2, reps=20000, chunk=2000, seed=0):
    rng = np.random.default_rng(np.random.SeedSequence([20261007, d, M, 1 if name.startswith('P1') else 4]))
    w = rng.normal(size=d); w /= np.linalg.norm(w)
    r0 = 0.5 * rng.uniform(-1, 1, d)
    x = s0 * w + (r0 - (w @ r0) * w)               # w.x = s0 exactly; orthogonal part arbitrary
    s = float(w @ x)
    phi, sig, kind, lam = spec['phi'], spec['sigma'], spec['kind'], spec['lam']
    lo, hi = spec['interval']
    # exact 1-D quantities by quadrature
    delta_q = phi(s + sig * math.sqrt(h) * XG) - phi(s)
    c_h = float((WG * delta_q**2).sum() / h)
    zbar_par = float((WG * delta_q * XG).sum() / math.sqrt(h))   # E z_par (gradient of the f=0 solution)
    Z, ZPAR, ZPERP2, SQRTQ = [], [], [], []
    for b in range(0, reps, chunk):
        r = min(chunk, reps - b)
        G = rng.standard_normal((r, M, d))
        delta = phi((x + sig * math.sqrt(h) * G) @ w) - phi(s)          # (r, M)
        zh = np.einsum('rm,rmd->rd', delta, G) / (M * math.sqrt(h))       # brute force in R^d
        zp = zh @ w
        Z.append(f_of(kind, lam, zh)); ZPAR.append(zp)
        ZPERP2.append(np.sum(zh * zh, axis=1) - zp**2)
        SQRTQ.append(np.sqrt(np.mean(delta**2, axis=1) / h))
    fz = np.concatenate(Z); zpar = np.concatenate(ZPAR); zperp2 = np.concatenate(ZPERP2); sq = np.concatenate(SQRTQ)
    out = dict(problem=name, d=d, M=M, s=s, h=h, reps=reps, c_h=c_h,
               E_zperp2_mc=float(zperp2.mean()), E_zperp2_theory=(d - 1) * c_h / M,
               E_zpar_mc=float(zpar.mean()), E_zpar_theory=zbar_par)
    zbar = np.array([zpar.mean()])
    if kind == 'quadratic':
        gap_mc = float(fz.mean() - (-0.5 * lam * zpar.mean()**2))
        gap_th = float(-0.5 * lam * (zpar.var() + (d - 1) * c_h / M))
        out.update(jensen_gap_mc=gap_mc, jensen_gap_theory=gap_th)
    else:
        En = float((-fz / lam).mean())
        lower = float(sq.mean() * (d - 1) / math.sqrt((d + 1) * M))
        upper = float(math.sqrt(zpar.mean()**2 + zpar.var() + (d - 1) * c_h / M))
        out.update(E_norm_mc=En, E_norm_lower=lower, E_norm_upper=upper)
    # certified projection onto C = {c w : c in [lo, hi]}: depends on z_par only
    cproj = np.clip(zpar, lo, hi)
    f_proj = -0.5 * lam * cproj**2 if kind == 'quadratic' else -lam * np.abs(cproj)
    ref_c = float(np.clip(zbar_par, lo, hi))       # f=0 gradient lies in C here; used as a common reference
    f_ref = -0.5 * lam * ref_c**2 if kind == 'quadratic' else -lam * abs(ref_c)
    out.update(raw_generator_bias_vs_ref=float(fz.mean() - f_ref),
               proj_generator_bias_vs_ref=float(f_proj.mean() - f_ref),
               proj_generator_q10_q50_q90=[float(v) for v in np.quantile(f_proj - f_ref, [.1, .5, .9])],
               raw_generator_rmse_vs_ref=float(np.sqrt(np.mean((fz - f_ref)**2))),
               proj_generator_rmse_vs_ref=float(np.sqrt(np.mean((f_proj - f_ref)**2))))
    return out


def main(out_dir='results/theory_onestep'):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, spec in PROBLEMS.items():
        for M in (4, 16):
            for d in (10, 50, 200, 1000):
                reps = 20000 if d <= 200 else 6000
                rows.append(run_case(name, spec, d, M, reps=reps, chunk=2000 if d <= 200 else 500))
                print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}), flush=True)
    (out_dir / 'onestep_theory_check.json').write_text(json.dumps(rows, indent=1))
    lines = ['| problem | d | M | E‖z⊥‖² MC | theory | Jensen gap / E‖ẑ‖ MC | theory or [lower, upper] | raw gen. bias | certified gen. bias | certified q10/q50/q90 |',
             '|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        if 'jensen_gap_mc' in r:
            mid = f"{r['jensen_gap_mc']:.4f} | {r['jensen_gap_theory']:.4f}"
        else:
            mid = f"{r['E_norm_mc']:.4f} | [{r['E_norm_lower']:.4f}, {r['E_norm_upper']:.4f}]"
        q = r['proj_generator_q10_q50_q90']
        lines.append(f"| {r['problem']} | {r['d']} | {r['M']} | {r['E_zperp2_mc']:.4f} | {r['E_zperp2_theory']:.4f} | {mid} | "
                     f"{r['raw_generator_bias_vs_ref']:.4f} | {r['proj_generator_bias_vs_ref']:.4f} | {q[0]:.3f} / {q[1]:.3f} / {q[2]:.3f} |")
    (out_dir / 'onestep_theory_check.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:])
