# Frozen Protocol: 32-Point Single-State Solver Comparison on PACE ICE

Status: **FROZEN BEFORE INPUT GENERATION, PILOT EXECUTION, OR FORMAL TRAINING**  
Revision: `single_point_32_ice_v1`  
Freeze date: 2026-10-23  
Branch: `ir-mlp-icml-single-point-32-ice`  
Base commit: `027240a72fe3028b4308275e5ec8f06722b79c33`

This document preregisters the formal comparison of Tangent-Actor + MCV,
standard Full-History Multilevel Picard (FH-MLP), and the audited Deep BSDE
PyTorch port. It is frozen before generating the 32-point banks and before any
pilot or formal run. Exact solution values are evaluation-only and may not be
used for training, checkpoint selection, solver-budget selection, early
stopping, or failure recovery.

## 1. Scientific question and scope

The primary question is whether one globally trained Tangent-Actor can return
accurate values at multiple initial states as ambient dimension grows, when it
is compared fairly with two traditional single-point solvers that are rerun at
each requested state.

Only these methods are in scope:

1. `tangent_actor`: one freshly trained global Actor per PDE, dimension, and
   training seed, followed by MCV inference at all 32 frozen points.
2. `fh_mlp_standard`: an untrained recursive Monte Carlo estimator rerun at
   each frozen point and replicate.
3. `deep_bsde`: one independently trained scalar-initial-value Deep BSDE model
   per frozen point and training seed.

DPI, PINN-HTE, SDGD-PINN, amortized Deep BSDE, oracle features, and reuse of a
previous Actor as a formal result are excluded. The two point solvers must not
be converted into global function learners.

## 2. Frozen PDEs and dimensions

Dimensions are exactly `20, 50, 100, 200, 400`. Query time is always `t0=0`.

### P1: `hjb_log`

\[
u_t+\Delta u-\lVert\nabla u\rVert^2=0,\qquad
u(1,x)=\log((1+\lVert x\rVert^2)/2).
\]

The reference is the existing Hopf-Cole radial reference, independently
recomputed and convergence-audited for every frozen radius.

### P2: `allen_cahn`

\[
u_t+\Delta u+u-u^3=0,\qquad
u(0.3,x)=1/(2+0.4\lVert x\rVert^2).
\]

The reference is the existing radial PDE machinery, rerun on a domain covering
the newly frozen `t=0` spatial radii. These points are not Brownian-marginal
samples. A nested domain/time/space refinement audit is mandatory.

### P3: `b6_fullrank`

\[
u_t+\Delta u-\tfrac12\lVert z\rVert^2=0,\quad
z=\sqrt2\nabla u,\quad T=0.25,
\]

\[
g(x)=-\log\left((2d)^{-1}\sum_{m=1}^{2d}e^{b_m^\top x}\right).
\]

The exact value is

\[
u^*(t,x)=-\log\left((2d)^{-1}\sum_{m=1}^{2d}
e^{b_m^\top x+(T-t)\lVert b_m\rVert^2}\right).
\]

The parameter arrays must be byte-identical to the completed ICE expansion.
The frozen `(beta, gamma)` pairs are `(1.0,1.0)`, `(1.5,1.0)`,
`(1.5,1.0)`, `(2.0,1.0)`, `(2.0,1.0)` for dimensions
`20,50,100,200,400`, respectively. Recalibration is forbidden.

### P4: `dense_fourier`

\[
u_t+\Delta u+\tfrac12\tanh(u)\lVert z\rVert^2=0,\quad
z=\sqrt2\nabla u,\quad T=0.5,
\]

\[
g(x)=\operatorname{arsinh}\left[\tfrac12+
\frac2{\sqrt{2d}}\sum_{j=1}^{2d}a_j
\sin(\omega_j^\top x+\phi_j)\right],
\]

\[
u^*(t,x)=\operatorname{arsinh}\left[\tfrac12+
\frac2{\sqrt{2d}}\sum_{j=1}^{2d}a_je^{-(T-t)}
\sin(\omega_j^\top x+\phi_j)\right].
\]

The Fourier directions, signs, and phases must be byte-identical to the
completed ICE expansion and each direction must have unit norm. No new
frequency is permitted. This PDE uses the complete joint value-gradient Actor.

## 3. Frozen query banks

Each PDE-dimension cell has 32 float64 states. Row zero is exactly the zero
vector. Rows 1--31 are drawn once with NumPy `Generator(PCG64(seed))`:

- `hjb_log`, `allen_cahn`: independent `N(0,I_d)` entries.
- `b6_fullrank`, `dense_fourier`: independent `U[-1,1]` entries.

Seeds are the unsigned little-endian integer represented by the first four
bytes of
`SHA256("single_point_32_ice_v1::query_bank::<pde>::<dimension>")`:

| PDE | 20 | 50 | 100 | 200 | 400 |
|---|---:|---:|---:|---:|---:|
| hjb_log | 1014295083 | 1870289811 | 1898610954 | 2768728099 | 1832156210 |
| allen_cahn | 1716816384 | 259137838 | 4074736349 | 1881011736 | 255990953 |
| b6_fullrank | 214171249 | 697137073 | 1737294624 | 2575111948 | 3771049365 |
| dense_fourier | 3392201883 | 1193042412 | 2884838781 | 2659431886 | 2229501467 |

For every bank, the preparation stage must save `t.npy`, `x.npy`, and
`reference.npy`, plus a JSON record of the seed, distribution, array shape,
dtype, per-file SHA-256, canonical concatenated-array SHA-256, parameter hash,
reference algorithm, and numerical uncertainty. No point may be replaced after
seeing any method's output.

The Actor training and validation generators use disjoint namespaces. The
preflight must assert that no query row appears byte-identically in a training
or validation bank.

## 4. Reference validation

- P1 uses two independent numerical implementations of the Hopf-Cole radial
  expectation (the mature repository implementation and a higher-precision
  one-dimensional quadrature). The required maximum disagreement is `5e-8`.
- P2 uses nested radial solves. The outer radius must exceed the largest query
  radius by the larger of `8` and four diffusion standard deviations. Increase
  the outer radius by 25%, halve the radial mesh width, and halve the time step.
  The maximum change at the 32 query radii must be at most `2e-5`; otherwise
  refine again. The reported per-point uncertainty is the maximum observed
  change across the last two accepted comparisons.
- P3 and P4 use their analytic formulas in float64 with stabilized log-sum-exp
  or pairwise summation. A second long-double/independent-order computation is
  required where supported. The uncertainty record includes the discrepancy
  and a roundoff bound.

A reference failure blocks the affected cell; it is never replaced by a model
estimate.

## 5. Tangent-Actor protocol

Formal training seeds are exactly `0,1,2`; each maps through the stable SHA-256
seed namespace `single_point_32_ice_v1::ta::<pde>::<d>::<seed>`.

Every run is fresh and uses:

- 8,192 independent training states and 512 independent validation states;
- 128 paths per tangent label and 32 label time steps;
- four hidden layers, width 384, SiLU;
- Adam learning rate `1.5e-4`, batch size 512, gradient clip 10;
- 5,000 optimizer steps per Actor update and four Actor updates;
- MCV labels and MCV inference;
- parent batch 8 and path chunk 32 (these may be lowered only as
  memory-layout parameters while preserving the total mathematical sample
  count);
- checkpoints every 250 optimizer steps;
- terminal-aware parameterization and no exact nonterminal feature.

P1 and P3 use the mature gradient Actor, P2 the mature value Actor, and P4 the
mature joint value-gradient Actor including tangent reaction weight, drift,
and source integral.

Selection among `a0,...,a4` is reference-free. The primary key is the frozen
held-out tangent-label/fixed-point validation score; numerical-validity gates
are applied first, and exact ties choose the earlier update. Query truth and
query predictions are unavailable to selection code.

Each selected Actor is evaluated on all 32 points using 1,024 MCV paths,
64 Euler steps, and four independent inference replicates. Training and
inference seeds are disjoint. Euler sensitivity at 128 and 256 steps is a
secondary audit for seed 0 at 100D and 400D and cannot change the 64-step main
result. Raw prediction for every point, training seed, and inference replicate
is retained.

## 6. Deep BSDE protocol

Deep BSDE is independently trained at every `(pde,d,point,seed)` with formal
training seeds `0,1,2`: `20*32*3 = 1,920` models. The forward state is initialized
at the frozen `x_i`, never at a hard-coded zero. A trained scalar `Y0` is used
only for that same point.

The audited architecture and numerical scheme are frozen:

- scalar trainable `Y0`, free `Z0`, and 19 time-specific nonshared Z networks;
- two hidden layers of width `d+10` per Z network, with batch normalization;
- 20 Euler intervals, `sigma=sqrt(2) I`, float64;
- `Y[k+1] = Y[k] - f(Y[k],Z[k])*dt + Z[k]^T dW[k]`;
- source-compatible clipped terminal loss with Huber transition at absolute
  terminal residual 50;
- Adam epsilon `1e-8`, batch size 64, validation batch 256;
- P1/P3/P4: 2,000 iterations, learning rate `1e-2`, `Y0~U[0,1]`;
- P2: 4,000 iterations, learning rate `5e-4`, `Y0~U[0.3,0.5]`;
- fixed validation Brownian increments and logging/checkpointing every 100
  iterations; the final scheduled iterate is the formal estimate. Validation
  loss is diagnostic and cannot select using truth.

P1/P2 are `native_source_family` after the 100D zero-state fidelity anchors
pass. P3/P4 are `faithful_transfer`: only the PDE driver, horizon, terminal
function, frozen parameters, dimension, and fixed initial state change. If a
different initialization, clipping rule, architecture, optimizer, or schedule
is later required, that run is separately named `adapted_transfer`, is fully
audited, and cannot silently replace the faithful-transfer main result.

Exact resume checkpoints must include model state, optimizer state, completed
iteration, all CPU/CUDA/NumPy/Python RNG states, training-generator state,
fixed validation increments (or their content hash and exact seed), elapsed
timing, and the best diagnostic record. A resume equivalence test must show
bitwise-identical CPU continuation and tolerance-equivalent CUDA continuation.

## 7. Full-History MLP protocol

The standard recursive estimator is run at every point for all six frozen
configurations

`(n,M) in {(2,16),(2,32),(3,8),(3,12),(3,16),(3,24)}`

with ten independent MC replicates per configuration. The preregistered main
configuration is `(n=3,M=16)` for every PDE and dimension. No test truth may
change this choice. All six configurations remain in the diagnostic output.

The audit must verify the terminal/value/gradient estimators, multilevel
telescoping signs, independent random streams across levels and branches, time
randomization densities and weights, diffusion scaling, and PDE driver sign.
No exact nonterminal value, Actor feature, or test-truth tuning is allowed.

For P4, a joint value-gradient recursion must be derived and tested directly;
the quadratic-driver double estimator must not be reused without a valid joint
nonlinear derivation. The previously validated repaired B6 estimator may be
reported as `fh_mlp_b6_repaired` only as a separate variant. Standard FH-MLP is
always retained, including nonfinite or unstable outcomes. A nonlinear plug-in
estimate is not described as unbiased unless proved for that exact estimator.

## 8. Stage gates

### Stage A: mathematical preflight

1. Generate and hash the 20 query banks and copy/verify frozen P3/P4 parameters.
2. Complete P1/P2 reference convergence audits.
3. Generalize Deep BSDE to arbitrary fixed `X0` and pass native 100D anchors.
4. Derive/test the P4 FH-MLP joint recursion and audit all four drivers.
5. Assert Actor query isolation and absence of oracle features.
6. Pass unit tests, CPU mathematical tests, H100 CUDA smoke tests, exact-resume
   tests, duplicate-writer tests, and artifact-schema tests.
7. Commit the frozen input/parameter manifest and Stage-A audit before Stage B.

Any mathematical failure blocks downstream execution.

### Stage B: nonformal 100D representative pilot

Use the first four frozen points (origin plus the first three random points),
seed 0, every method, and the small FH-MLP grid. Pilot outputs are namespaced
`pilot` and are excluded from all formal tables. The gate checks end-to-end
correctness, finite values, H100 memory, wall time, signs, reference lookup, and
checkpoint continuation.

### Stage C: formal dimension scaling

Release dimensions in order `20 -> 50 -> 100 -> 200 -> 400`. Independent
single-GPU jobs may run concurrently, with a default throttle of four H100s.
Every task has one deterministic ID, one atomic status file, and at most one
writer lease. Retrying an incomplete task resumes its exact checkpoint; a
complete task is immutable. Reducing microbatch or path chunk is permitted for
OOM only when total samples and equations remain unchanged and the change is
recorded.

### Stage D: independent audit

The auditor loads raw prediction artifacts rather than aggregate tables,
verifies all bank/parameter/algorithm hashes and seeds, and recomputes every
metric. A PDE-dimension cell is `COMPLETE` only after all formal repeats of all
three main methods exist and pass integrity checks. Otherwise it is `PARTIAL`
and excluded from complete-cell rankings.

### Stage E: aggregation and publication artifacts

Generate the three specified figures, full tables, LaTeX tables, resource
accounting, artifact audit, and the report
`invariant_region_mlp/docs/ICML_SINGLE_POINT_32_ICE_REPORT.md`. Commit and push
the complete results and report, and record the final Git SHA.

## 9. Metrics and uncertainty

For predictions `p[i,r]` and frozen references `u[i]`, the primary statistic is

\[
RMSE_{32,R}=\sqrt{(32R)^{-1}\sum_{i=1}^{32}\sum_{r=1}^R
(p_{i,r}-u_i)^2}.
\]

Report RMSE, MAE, relative L2, signed bias, maximum absolute error, RMSE on the
31 random points, absolute error at the origin, nonfinite count, and all raw
repeat predictions. Relative L2 is
`sqrt(sum error^2 / sum reference^2)` over the same finite repeat-expanded
sample; a zero denominator is reported rather than regularized silently.

Between-training-seed variation and MC-inference variation are separate for
Tangent-Actor. Deep BSDE reports training-seed variation. FH-MLP reports
recursive-MC variation. Numerical-reference uncertainty is stored separately
and is never folded invisibly into method error. Pointwise bootstrap summaries
must state that the bank includes one deterministic origin and that training
randomness is not query-sampling randomness. No best-seed reporting is allowed.

## 10. Compute and failure policy

The experiment uses the existing PACE ICE H100 environment and account limits.
Jobs resolve the repository root from the submitting checkout and use unique
job names and result directories. Existing whole-function results, checkpoints,
protocols, and Slurm ledgers are read-only.

Record wall time, GPU allocation time, peak GPU memory, CPU/RAM request, Slurm
job/array/task IDs, host/GPU/environment metadata, restart count, and checkpoint
hash. Scheduler, quota, OOM, timeout, nonfinite, and mathematical failures have
distinct statuses. Results are never fabricated, and account limits are not
increased automatically.

## 11. Required artifacts

Code lives in `invariant_region_mlp/experiments/single_point_32_ice/`; results
live in `invariant_region_mlp/results/single_point_32_ice/`. Required outputs
are:

- `FROZEN_PROTOCOL.json`, `frozen_query_banks/`, `parameter_manifest.json`,
  `reference_validation.json`, `training_manifest.csv`,
  `method_fidelity_audit.md`;
- `ta_results.csv`, `fh_mlp_results.csv`, `deep_bsde_results.csv`,
  `main_pointwise_results.csv`, `main_dimension_scaling.csv`,
  `randomness_decomposition.csv`, `resource_accounting.csv`;
- `execution_status.json`, `artifact_audit.json`, `final_summary.json`,
  `checkpoints/`, `predictions/`, `figures/`, and `tables/`.

Figures are (1) log-scale RMSE versus dimension with one panel per PDE,
(2) pointwise absolute-error distributions at 100D and 400D, and
(3) method-specific numerical-stability diagnostics. No new formal
accuracy-compute frontier is part of this experiment.

## 12. Source provenance and amendment rule

Primary source commits are:

- completed B6/Dense ICE expansion and branch base:
  `027240a72fe3028b4308275e5ec8f06722b79c33`;
- mature whole-function suite:
  `136e892ba6ffd5fcf33d5d3dc618828584675f44`;
- validated B6 source:
  `281977883e0a5ce4f05895878680b582620461be`;
- Dense Fourier pilot:
  `3d3d4c80f0b29af39e7d509a143812a1f5a702f0`;
- strong-baseline source:
  `66495efb3510ab99514942ed7f9ee62259b2d21d`.

Implementation commits after this freeze may fix or add code needed to execute
the protocol, but may not change banks, formal budgets, main estimators,
metrics, or selection rules. A scientifically material deviation requires an
explicit versioned amendment committed before affected runs, separate artifact
names, and disclosure in the final report. It may not overwrite v1 results.

