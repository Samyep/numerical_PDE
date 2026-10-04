# Invariant-Region Multilevel Picard (IR-MLP)

**Canonical development location.** This directory is now the maintained IR-MLP repository. The older standalone `Samyep/invariant-region-mlp` repository is archival.

IR-MLP inserts PDE-certified projection inside stochastic Picard recursion, before noisy intermediate states are reused by nonlinear generators.

## Current central result: IR-MLP breaks a known MLP counterexample

For the Hutzenthaler--Nguyen (2025) HJB counterexample
[
u_t + \tfrac12\Delta u + (\sum_{j=2}^d |\partial_{x_j}u|^2)^{1/2}=0,\qquad u(1,x)=|x_1|,
]
the true gradient lies in the certified subspace `span(e1)`. Projecting each recursive gradient state to this subspace makes the nonlinear driver vanish pathwise, removes the spurious inactive-coordinate feedback, and reduces the scheme to one-dimensional terminal Monte Carlo.

At `(t,x)=(0,0)`, with `N=m^n` terminal samples,
[
\mathbb E\|\widehat Y^{IR}-Y^*\|_2^2=(4-2/\pi)/N,
]
which is independent of the ambient dimension.

At `n=m=3`, 512 paired replicas give full-state RMSE about 0.356 for recursive subspace IR from d=2 through d=1000, while the original MLP grows from about 1.20 (d=2) to about 1.70e4 (d=1000). Final-output-only projection does not repair the value error; projection must occur before nonlinear feedback.

See:
- `experiments/counterexample_rescue/` for the simulator and validation code.
- `results/counterexample_rescue/` for aggregate and raw replica results.
- `docs/counterexample_rescue_theory.md` for the proof notes.
- `docs/COUNTEREXAMPLE_REPORT.md` for the audited experiment report.

## Broader project

Other current evidence includes high-dimensional quadratic HJB, nonlinear finance, Neufeld--Wu 100--300D gradient-dependent MLP benchmarks, mechanism diagnostics based on invariant-region violations, and negative controls. Older SCaSML integration material is retained for reference but is no longer the center of the project.

The maintained research question is:

> Can PDE-certified invariant regions remove statistically spurious directions in high-dimensional stochastic fixed-point solvers before nonlinearities amplify them?
