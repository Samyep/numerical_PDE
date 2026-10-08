"""Han-Jentzen-E (PNAS 2018) HJB: u_t + Lap u - lam |grad u|^2 = 0, g = log((1+|x|^2)/2), d=100, T=1.
Hopf-Cole: u = -(1/lam) log E exp(-lam g(x + sqrt2 W_{T-t})).  f=0 solution: E g(x + sqrt2 W_{T-t}).
|x + sqrt2 W_h|^2 ~ 2h * noncentral chi^2_d(|x|^2/(2h)) -> 1-D quadrature over the chi^2 law (exact)."""
import numpy as np, math
from scipy import stats, integrate
d, T = 100, 1.0
def expect(fun, r2, h, lam=None):
    nc = r2/(2*h) if h > 0 else 0.0
    dist = stats.ncx2(d, nc) if nc > 0 else stats.chi2(d)
    lo, hi = dist.ppf(1e-12), dist.ppf(1-1e-12)
    val, _ = integrate.quad(lambda s: fun(2*h*s)*dist.pdf(s), lo, hi, limit=400)
    return val
g = lambda R2: np.log((1+R2)/2)
def u_exact(r2, h, lam):
    if h == 0: return g(r2)
    # -(1/lam) log E exp(-lam g) = -(1/lam) log E ((1+R2)/2)^(-lam)
    return -math.log(expect(lambda R2: ((1+R2)/2)**(-lam), r2, h))/lam
def u_lin(r2, h): return g(r2) if h == 0 else expect(g, r2, h)
print('u(0,0) for lam = 1, 10, 50 (PNAS reports 4.5901 at lam=1):')
for lam in (1, 10, 50):
    print(f'  lam={lam}: u={u_exact(0.0, T, lam):.4f}   f=0: {u_lin(0.0, T):.4f}')
# test distribution: t ~ U[0,T), x ~ U[-1,1]^d  (|x|^2 ~ d/3 concentrated); also x ~ N(0, I) like the PNAS sample paths
rng = np.random.default_rng(0)
for name, X in [('x ~ U[-1,1]^d', rng.uniform(-1, 1, (300, d))), ('x = sqrt2 W_t (forward paths)', None)]:
    t = rng.uniform(0, T, 300)
    if X is None: X = math.sqrt(2) * np.sqrt(t)[:, None] * rng.standard_normal((300, d))
    r2 = (X**2).sum(1); h = T - t
    for lam in (1, 10):
        U = np.array([u_exact(a, b, lam) for a, b in zip(r2, h)]); L = np.array([u_lin(a, b) for a, b in zip(r2, h)])
        A = np.stack([h, h**2], 1); Lt = L + A @ np.linalg.lstsq(A, U-L, rcond=None)[0]
        sd = U.std(); sk = lambda p: np.sqrt(np.mean((p-U)**2))/sd
        print(f'{name:30s} lam={lam:2d}: nonlinear share ||u-u_lin||/||u|| = {np.linalg.norm(U-L)/np.linalg.norm(U):.4f}   S_NL={sk(L):.3f}   G2={sk(Lt):.3f}')
