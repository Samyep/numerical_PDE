# Three-seed, two-PDE hyperparameter-specific accuracy--compute scaling

This document is frozen before any new timing calibration, formal training, or
test evaluation for this study.  The study is deliberately separate from
`whole_function_compute_frontier`: its outputs do not overwrite the previous
single-seed experiment or its failed fidelity check.

## Scientific question and terminology

For B6 Full-Rank HJB and Dense Fourier in 100 dimensions, where does each
method benefit from spending additional computation: additional nonlinear
iterations, a more accurate stochastic estimator, or more network
optimization/sampling per nonlinear iterate?

The result is called **hyperparameter-specific accuracy--compute scaling**.
It is not called an optimal frontier: only three one-parameter sweeps per
method are examined, no joint hyperparameter search is performed, and no
test-selected lower envelope is reported.

## Frozen scope

- PDEs: B6 Full-Rank HJB and Dense Fourier, dimension 100.
- Methods: Tangent-Actor (TA), DPI, and PINN-HTE.
- Sweeps: C0, C1, and C2 as defined below.
- Independent training seeds: exactly `0`, `1`, and `2`.
- Hardware: NVIDIA H100 80GB for calibration, training, and inference.
- Test data: the existing frozen 960-state whole-function bank for each PDE.
- Formal allocation cap: 90 H100 GPU-hours, including calibration, failed
  attempts, retries, evaluation, and aggregation.
- No new SDGD-PINN training is part of this study.

The expected frozen input hashes are inherited from the audited preceding
study:

| PDE | test-input SHA-256 | truth-value SHA-256 |
|---|---|---|
| B6 Full-Rank HJB | `0ee05fff7ce16f30d0d31563f60df7db3eb23a9c46bbe7bc9be8eb3dbee1311b` | `6262db7f6b7ee2c41cb6011dadb7288226157533488f708be7d6ea6773f2e13a` |
| Dense Fourier | `d409296128f48c5ae00c5de2c6ddd998768c01dff18a89c385245a50ff581b70` | `3b33ea36e21595b5c20e5d91097d1fadfcbf3f7884c0328f20579e69995fc816` |

## What each sweep changes

Exactly one named hyperparameter varies along each line.  All other
method/PDE choices, including network, objective, optimizer, learning-rate
schedule, state law, terminal loss, and discretization, stay fixed at the
audited mature configuration.

| Method | C0: anytime trajectory | C1: stochastic-estimator accuracy | C2: optimization/sampling amount |
|---|---|---|---|
| TA | Actor updates | Monte Carlo paths per label | Optimizer steps per update |
| DPI | Picard iterations | Integral Monte Carlo samples | Epochs per Picard iteration |
| PINN-HTE | Training epochs | Hutchinson probes | Batch size |

The fixed baseline tuple is:

- TA: 4 Actor updates, 128 paths per label, 5,000 optimizer steps per update.
- DPI: 80 Picard iterations, 4,096 integral samples, 16 epochs per Picard.
- PINN-HTE: 30,000 epochs, 16 Hutchinson probes, batch size 512.

C0 is one continuous run per `(PDE, method, seed)`, observed at genuine
algorithmic checkpoints.  C1 and C2 are independent fresh training runs for
every `(PDE, method, sweep, parameter, seed)`.  Even when a C1/C2 parameter
equals the baseline tuple, its run is not reused from C0 or the other sweep.

No RMSE is replaced by a running minimum, smoothed to enforce monotonicity, or
formed from a test-selected envelope.  Nonmonotone curves remain nonmonotone.

## Timing-only calibration and deterministic grid freeze

The candidate ladders are:

| Method/sweep | ordered candidates |
|---|---|
| TA C0 updates | `1,2,3,4,5,6,8,10,12,16` |
| TA C1 paths | `32,64,128,192,256,384,512,768,1024` |
| TA C2 fit steps | `1250,2500,5000,7500,10000,15000,20000,30000,40000` |
| DPI C0 Picard | `20,40,80,160,320,640,1280` |
| DPI C1 MC samples | `512,1024,2048,4096,8192,16384,32768,65536` |
| DPI C2 epochs/Picard | `4,8,16,32,64,128,256` |
| PINN C0 epochs | `10000,30000,60000,120000,240000,480000` |
| PINN C1 probes | `4,8,16,32,64,128,256` |
| PINN C2 batch | `128,256,512,1024,2048,4096,8192,16384` |

Before any formal run, a timing-only calibration uses calibration seed
`314159`.  Calibration may inspect synchronized elapsed GPU time, peak memory,
finite/nonfinite status, and OOM status only.  It must not load the frozen
nonterminal truth, compute a test error, or retain a calibration model as a
formal result.  Existing timings may guide which anchors are cheap to measure,
but every method/PDE receives a fresh H100 timing anchor.

For each of the 18 `(PDE, method, sweep)` lines, measured anchors define a
monotone piecewise log-linear time model after replacing a downward timing
fluctuation by the preceding cumulative maximum.  Candidate values nearest
the target elapsed times `5, 10, 15, 20, 30` minutes are selected and
deduplicated.  At least three distinct values are required.  The final value
must be the candidate whose predicted time is closest to 30 minutes among
values in `[25,35]` minutes.  If the initial ladder remains below 25 minutes,
the varied parameter is repeatedly doubled; this is the predeclared response
to slow time growth.  If adjacent candidates jump from below 25 to above 35
minutes, one rounded geometric-intermediate parameter is inserted and timed.
No point is chosen or removed using accuracy.

Calibration produces `CALIBRATION_MANIFEST.json`, containing raw anchors,
the fitted time model, memory observations, selected ladders, projected GPU
cost, code revision, and hashes.  It is committed before formal training.

The projected formal allocation may use at most 82 GPU-hours, reserving eight
hours for calibration, inference, failures, and audits.  If the five-target
grids exceed 82 projected hours, every C1/C2 line first keeps the candidates
nearest 5, 15, and 30 minutes; remaining candidates nearest 10 and 20 minutes
are then added round-robin in lexical `(PDE, method, sweep)` order while the
projection remains within 82 hours.  C0 trajectories are never pruned.  Thus
each sensitivity line has three to five points, and different lines need not
have equal point counts.

Once the calibration manifest is committed, formal grids are immutable.
Actual formal times may fall outside 25--35 minutes and are reported as
measured; points and lines are never extended, repeated, or interpolated to
manufacture a 30-minute endpoint.

## Randomness, training, and inference

Formal training seeds are exactly `0,1,2`.  Configuration identity, seed,
initialization hash, final/checkpoint hash, and the relevant random-state
metadata are recorded.  Formal training contains no exact nonterminal-solution
callback.

The preceding study's Dense Fourier DPI C0 mismatch remains recorded as a
failed historical-consistency audit: historical RMSE `0.0646384` versus the
truth-gated run `0.1151372`.  Its post-freeze RNG-only diagnostic obtained
`0.0662925`, a 2.56% difference from the historical RMSE and a 96.72%
reduction of the discrepancy.  This strongly implicates the historical
`EVAL.FREQ=8` exact-solution callback, which advanced the global Torch RNG.
It was not a bitwise reproduction and does not turn the earlier failure into a
pass.  In the new study `EVAL.FREQ=null`; exact callbacks are absent, and
their time is not mixed into the compute axis.

TA inference is fixed at 1,024 paths, 64 Euler steps, the mature martingale
control variate, and three independent inference replicates per trained
checkpoint.  No prediction ensemble is formed.  For one training seed, the
reported RMSE is the arithmetic mean of its three replicate RMSEs.  The
within-seed replicate standard deviation is stored separately.  DPI and
PINN-HTE use deterministic checkpoint inference under the frozen test bank;
any remaining evaluation randomness is explicitly seeded and recorded.

When a large stochastic configuration does not fit in memory, only the
following estimator-preserving fallbacks are allowed:

- DPI integral paths may be evaluated in chunks and accumulated before the
  same sample mean/target is formed.
- TA label paths may be chunked and accumulated before the same label mean is
  formed.
- A PINN batch may use exact loss-sum gradient accumulation over microbatches,
  followed by one optimizer step, provided the network has no batch-dependent
  layers.  The full batch's fixed Hutchinson draws are preserved.

Fallbacks, OOM attempts, and their allocated GPU time are reported.  No
smaller completed subset may be assigned the nominal larger parameter value.

## Truth gate and metrics

Calibration and all formal training finish before the nonterminal frozen
truth gate opens.  Test truth cannot affect calibration, grid choice,
training, checkpoint choice, stopping, retries, or extensions.  C0 reports the
model at each predeclared algorithmic checkpoint; C1/C2 report the terminal
model of the fixed run.  Reference-free diagnostics may be recorded but do
not select another checkpoint.

For every seed and inference replicate, report:

- whole-function value RMSE (primary);
- relative L2, RelMAE, bias, nonfinite count;
- synchronized active training seconds;
- inference seconds, peak GPU memory, job wall time, and allocated GPU time.

For a fixed hyperparameter value `lambda`, the plotted point is

`(mean_seed training_time(seed,lambda), mean_seed RMSE(seed,lambda))`.

The vertical error bar is the sample standard deviation (`ddof=1`) of the
three per-training-seed RMSE values.  It is not the uncertainty of an ensemble
prediction.  Timing variability and TA inference-replicate variability are
reported separately.

## Figure and reporting

The main figure has B6 and Dense Fourier panels.  Method colors are fixed
across panels; C0/C1/C2 use solid/dashed/dotted line styles.  Points use their
measured mean training time, log-scale RMSE, and training-seed standard-
deviation bars.  Parameter values are available in the table/metadata and may
be annotated sparingly.  There is no running-minimum curve, interpolation,
test-selected envelope, or claim of an optimal Pareto frontier.

The report begins with a completeness/audit table, then the main mean±std
results, timing calibration, failures, Dense Fourier DPI provenance, resource
accounting, and exact commit/input/artifact hashes.  Negative results and
nonmonotonicity are reported without threshold changes.

Required outputs are placed only in the new study directories:

```
invariant_region_mlp/results/whole_function_compute_scaling_3seed/
  CALIBRATION_MANIFEST.json
  run_manifest.json
  seed_metrics.csv
  scaling_summary.csv
  resource_accounting.csv
  artifact_audit.json
  figures/two_pde_three_seed_scaling.{png,pdf}
invariant_region_mlp/docs/WHOLE_FUNCTION_COMPUTE_SCALING_3SEED_REPORT.md
```

The original frontier results and manuscript are not modified.
