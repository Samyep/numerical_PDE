# Joint nonseparable PDE extension to the ICE single-point benchmark

Repository `Samyep/numerical_PDE`, branch
`ir-mlp-icml-joint-pde`, based on completed ICE suite commit `32ef8b074`.
This document must be committed unchanged before implementation, reference
construction, or solver execution.  The completed four-PDE suite, its frozen
protocol, code, and results remain unchanged.

## Motivation and prior information

The completed suite already contains a joint driver in `dense_fourier`,

\[
f(y,z)=\tfrac12\tanh(y)\lvert z\rvert^2.
\]

This extension therefore does **not** claim that joint value-gradient
dependence is absent from the original suite.  It asks a narrower question:
do the three methods remain competitive for a genuinely nonseparable driver
with a nonzero mixed derivative, rather than a scalar value-dependent
coefficient multiplying \(\lvert z\rvert^2\)?

The equation family and Tangent-Actor-only behavior were studied previously
in `convex_concave_semilinear_extended`.  In particular, a 100D three-seed
pilot reported a final Tangent-Actor RMSE near `0.0047`.  No FH-MLP or Deep
BSDE comparison on the query banks below has been observed.  Consequently,
this is a prospective strong-baseline extension, not a novel-PDE discovery,
and all claims must disclose this prior information.

## PDE

For \(d\in\{20,50,100,200,400\}\), solve

\[
\partial_t u+\Delta u+f(u,z)=0,\qquad z=\sqrt2\nabla u,
\qquad T=0.5,
\]

with

\[
f(y,z)=y-y^3-\tfrac12\lvert z+yv\rvert^2,
\]

and terminal condition

\[
g(x)=0.5+0.15\sin(\beta_1^\top x)
          +0.15\cos(\beta_2^\top x).
\]

The columns \(\beta_1,\beta_2\) are the sign-canonicalized thin-QR
orthonormalization of a `d x 2` standard-normal matrix drawn from
`default_rng(SeedSequence([20261232, d, 17]))`; \(v=\beta_1\).  The same
matrix is used by every method and is frozen before solver execution.

This driver is jointly concave on \(y\in[0,1]\), has

\[
f_y=1-3y^2-(z+yv)\cdot v,\qquad f_z=-(z+yv),
\qquad \nabla_z f_y=-v\ne0,
\]

so it is genuinely nonseparable in \((y,z)\).  The terminal data lie in
`[0.2, 0.8]`; the invariant-domain check uses `[0,1]` but no method may use
the numerical reference for projection, labels, training, stopping, or model
selection.

## Independent reference

Because the data depend only on
`s=(beta_1^T x, beta_2^T x)`, the exact ambient solution reduces to a
two-dimensional periodic PDE.  References are constructed by the existing
independent Fourier solver, which is not imported by any learning method.

- primary grid: `256 x 256`, internal step at most `1.25e-4`;
- audit grid: `128 x 128`, internal step at most `2.5e-4`;
- periodic domain and interpolation follow the independently validated
  `convex_concave_semilinear_extended` solver;
- primary/audit RMSE must be at most `5e-5`, maximum difference at most
  `2e-4`, PDE-residual RMSE at most `2e-6`, and terminal interpolation maximum
  error at most `1e-4`;
- the numerical solution must remain in `[0,1]` up to `1e-8`.

Failure of any reference gate stops all GPU solver runs and is reported as
`BLOCKED_REFERENCE_VALIDATION`.

## Fresh query banks and randomness

This extension uses no query point or stochastic solver draw from the
completed suite or the earlier nonseparable pilot.

- 32 time-zero queries per dimension;
- row zero is exactly the origin;
- the remaining 31 rows are iid `N(0,I_d)` from
  `default_rng(SeedSequence([20261232, d, 32]))`;
- all stochastic method seeds are derived from the revision string
  `single_point_joint_nonseparable_v1`, method, dimension, query index,
  configuration, and repeat using SHA-256 as in the completed suite;
- Tangent-Actor uses training seed `0`; Deep BSDE uses training seeds
  `{0,1,2}`; FH-MLP uses repeats `0,...,9`.  The revision namespace makes
  every underlying draw fresh even when an integer index matches an older
  study.

The frozen parameter and query files, references, and their SHA-256 hashes
must be committed before the first formal solver job.

## Methods and fixed budgets

All implementations import the completed suite or its already audited source
modules; the extension adds adapters and dispatch only.  No existing solver
module is edited.

1. **Tangent-Actor.** The completed suite's joint terminal-aware Actor:
   8,192 training states, 512 validation states, 128 paths per label, 32 label
   steps, 4x384 SiLU, 5,000 AdamW steps per update, four updates, joint
   value-gradient MCV labels, the frozen reference-free selector, and seed 0.
   Final evaluation uses 1,024 paths, 64 Euler steps, and four independent
   inference replicates.
2. **Standard Full-History MLP.** The completed generic float64 recursion,
   pathwise terminal gradient, untruncated Bismut level gradient, and no
   truth-based selection.  Run `(n,M)` in
   `{(2,16),(2,32),(3,8),(3,12),(3,16),(3,24)}`, ten repeats each.  The
   preregistered main cell is `(3,16)`; all cells are reported.
3. **Deep BSDE.** The completed single-initial-state implementation, 20 time
   intervals, batch size 64, validation size 256, float64, three independent
   trainings at every query, and the same fixed training schedule used by the
   non-Allen equations in the completed suite.  Truth is never used for
   optimization or checkpoint selection.

No hyperparameter search is permitted.  Infrastructure failures may be
resumed from exact checkpoints; a mathematical or numerical failure is
reported, not replaced by another configuration.

## Validation before formal runs

The CPU/CUDA preflight must verify:

- terminal identity and terminal gradient against central differences;
- analytic tangent coefficients against autograd to `1e-11` in float64;
- a nonzero mixed derivative and joint-Hessian negative semidefiniteness on
  at least 1,024 frozen states;
- the two-dimensional reduction against the ambient driver and gradient;
- reference convergence and residual gates above;
- one draw-tape replay for FH-MLP and one checkpoint/resume replay for Deep
  BSDE;
- absence of nonterminal reference access in all three method call graphs.

Run a single 100D end-to-end smoke cell before formal dispatch.  The smoke is
an infrastructure gate only: finite outputs and exact artifact replay are
required, but observed accuracy cannot change any budget, threshold, or
method.

## Metrics and preregistered judgments

Primary metric: RMSE over the 32 points and formal repeats, with nonfinite
predictions making the corresponding method-cell outcome a numerical failure.
Also report relative L2, bias, repeat dispersion, nonfinite counts, synchronized
GPU time, peak memory, and reference uncertainty.

- **J-1 (best-method coverage):** Tangent-Actor has the lowest primary RMSE
  among the three preregistered methods in at least 4 of 5 dimensions.
- **J-2 (versus standard FH-MLP):** Tangent-Actor RMSE is at most `0.5` times
  the `(n=3,M=16)` FH-MLP RMSE in at least 4 of 5 dimensions.
- **J-3 (versus Deep BSDE):** Tangent-Actor RMSE is at most `0.5` times the
  Deep-BSDE RMSE in at least 3 of 5 dimensions.
- **J-4 (dimension robustness):** Tangent-Actor's 400D/20D RMSE ratio is at
  most `3.0`.

These are extension-specific judgments and do not alter, overwrite, or
relabel any conclusion of the completed four-PDE study.  Every failure is
reported verbatim.  No method is selected using truth, and descriptive best
FH cells cannot replace the preregistered main cell in J-2.

## Execution order and compute ceiling

1. Freeze this document.
2. Implement adapters/tests and commit them.
3. Freeze parameters, query banks, references, and hashes; run preflight.
4. Run the 100D smoke.
5. Run formal dimensions in order `20, 50, 100, 200, 400`.
6. Audit, aggregate, report, commit, and push.

Use PACE ICE H100 GPUs.  The extension has a ceiling of 12 allocated H100
GPU-hours.  A formal allocation has an 8-hour wall limit.  If the ceiling is
reached, completed cells remain valid and missing judgments are marked
`NOT EVALUATED (COMPUTE)`; dimensions, repeats, or budgets may not be silently
reduced.

## Outputs

New code and artifacts are isolated from the frozen suite:

```text
invariant_region_mlp/experiments/single_point_joint_extension/
invariant_region_mlp/results/single_point_joint_extension/
invariant_region_mlp/docs/ICML_SINGLE_POINT_JOINT_EXTENSION_REPORT.md
```

The compact committed outputs include frozen manifests, validation audits,
method/result CSV files, the five-dimension RMSE table, resource accounting,
an artifact audit, and publication figures.  Checkpoints, label caches,
per-query predictions, raw traces, and Slurm logs remain on ICE and are bound
to the compact artifacts by hashes and row counts.
