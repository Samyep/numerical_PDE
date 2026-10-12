# HN 2025: a minimal failure-and-repair introduction figure

**Status:** Executed on CPU with the full-history Eq. (6) structure and Beta(1/2,1) random-time sampling.

**Source:** Hutzenthaler and Nguyen, arXiv:2506.23969 (2025), Eq. (6), https://arxiv.org/pdf/2506.23969 .

**PDE:** $u_t+\frac12\Delta u+\|\nabla_{2:d}u\|_2=0$, $u(1,x)=|x_1|$. Exact $u(0,0)=\sqrt{2/\pi}$, gradient zero.

**Protocol:** $d=50,100,200,400$; $n=M=3$, 64 independent paired repetitions per dimension; $t=0$, $x=0$; float64. No NN training, actor, double, learned baseline or projection. Modified method changes only the terminal gradient from centered BEL to the a.e. pathwise derivative.

## Results

| d | Vanilla joint RMSE | Pathwise joint RMSE | Vanilla value RMSE | Pathwise value RMSE |
|---:|---:|---:|---:|---:|
| 50 | 171.57 | 0.22055 | 19.343 | 0.11467 |
| 100 | 492.43 | 0.2419 | 39.234 | 0.12161 |
| 200 | 1393.2 | 0.2232 | 77.403 | 0.11457 |
| 400 | 4111.3 | 0.20458 | 151.8 | 0.1104 |

For all 95% bootstrap intervals (5,000 resamples), measured biases, gradient errors, per-run work and run time, consult `summary.csv` and `results.csv`.

**Independent theoretical pathwise sanity check:** with $N=M^n=27$ root terminal samples at $(0,0)$, the true pathwise value RMSE is $\sqrt{(1-2/\pi)/27}\approx0.11601$ and full-state RMSE is $\sqrt{(2-2/\pi)/27}\approx0.22471$ (expectation over algorithm randomness). Our finite-64-repetition estimates are compatible with these exact scales.

**Observed growth exponents (four points, descriptive only):** the vanilla value RMSE scales approximately as $d^{0.99}$ and vanilla joint-state RMSE as $d^{1.52}$; these are finite-budget fits, not theoretical exponents.

## Internal validation

- Original HN convention retained: diffusion coefficient 1, terminal T=1, Gaussian increment sqrt(h) G, gradient weight G/sqrt(h).
- Terminal is `(g(x),0)+MC[g(X)-g(x)]*(1,G/sqrt(h))`, i.e. **centered BEL as in HN Eq. (6)**.
- Source random time R=t+h U^2, U Uniform(0,1); exact inverse density 2h U; weighted derivative kernel simplifies to 2 sqrt(h) G.
- Full-history telescoping: sum of level corrections for ell=1,...,n-1, using independent recursive calls for high and low components. The ell=0 term is identically zero since f(U0)=0.
- Both estimators consume the **same random stream per paired seed**, including the recursive calls. All sample counts were checked to agree between paired methods.
- Pathwise terminal gradients have zero inactive coordinates pathwise and all recursive driver evaluations vanished exactly (verified at n<=4 on small d and during production).
- The pathwise value output equals the terminal Monte Carlo sample mean (all levels) to floating-point arithmetic precision.
- The source theorem establishes an asymptotic lower bound on joint state error for a designated recursion; this finite (n,M) sweep does **not** prove the source asymptotic claim.

## Interpretation

The exact solution is one-dimensional and its Hamiltonian vanishes. The source full-history MLP applies a nonlinear norm to the spurious sampled gradient in d-1 irrelevant directions, so deeper nonlinear source corrections can grow with dimension. The pathwise terminal estimator keeps these irrelevant coordinates identically zero and therefore makes all source terms vanish recursively. **This is a diagnostic explanation of a special counterexample, not evidence that the learned actor is necessary or that general nonlinear PDEs are solved.**

## Costs and limitations

- Wall time to execute both methods and produce files: 1.67 seconds (CPU environment; timings include Python overhead).
- Sample counts per root: terminal 117, transition 21, Gaussian vectors 138, generator evaluations 24.
- The run evaluates only (t,x)=(0,0), not a spatial test grid.
- The pathwise estimate includes stochastic terminal value and active-gradient noise; it removes only the irrelevant-gradient nonlinear feedback for this specially structured source.
- No training cost or GPU execution occurs; this is intentionally a small introduction-only experiment.

## Reproduce

From this directory: `python run_hn_intro.py --reps 64 --dims 50 100 200 400 --n 3 --m 3`

Outputs: `results.csv`, `summary.csv`, `figure1_hn_failure_repair.pdf`, `figure1_hn_failure_repair.png`, `config.json`.