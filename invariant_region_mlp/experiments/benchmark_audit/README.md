# Audit of classic high-dimensional benchmarks (exact, no MLP)

`hje_hjb_audit.py`: Han-Jentzen-E (PNAS 2018) HJB, d=100, T=1, g=log((1+|x|^2)/2), via Hopf-Cole and
one-dimensional noncentral chi-square quadrature (exact up to quadrature error).

Output (2026-10-08):
- u(0,0): lam=1 -> 4.5902 (published reference 4.5901); f=0 (linear) solution -> 4.6002 (0.22% difference).
- On x ~ U[-1,1]^100 and on forward paths: nonlinear share ||u-u_lin||/||u|| = 0.16-0.21% at lam=1;
  best time-only correction leaves skill 0.003 (gate G2), i.e. the nonlinear effect is a spatially constant drift.
- lam=10: share 1.7-2.3%, G2 = 0.03.

For comparison, Allen-Cahn (PNAS 2018, d=100, T=0.3, u-dependent nonlinearity): f=0 value 0.0391 vs reference
0.0528 (26% difference): the u-channel nonlinearity there is real.
Earlier audits: Rosenbrock HJB (SCaSML) 0.26%; 100D funding / different interest rates: the solution lies
between the two constant-rate Black-Scholes prices, nonlinear premium 0.94%.
