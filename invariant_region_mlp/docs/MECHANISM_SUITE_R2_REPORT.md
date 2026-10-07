# Mechanism suite, round 2: confirmatory re-run

Round 2 used criteria written after the round-1 results were seen. It is therefore a fresh-randomness confirmatory study, not a replacement for the original pre-registered verdicts. Round-1 results and files were left unchanged.

## Round-1 outcome (unchanged)

| PDE | C1 | C2 | C3 | C3 joint-cell fraction | certificate-prior dominated | centre better fraction | shrink-centre better fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | PASS | PASS | FAIL | 0.5000 | no | 0.0000 | 0.1250 |
| P2_a4 | PASS | FAIL | FAIL | 0.5000 | yes | 0.0000 | 0.3750 |
| P2_a8 | PASS | FAIL | PASS | 0.7500 | no | 0.0000 | 0.2500 |
| P3_rho1 | PASS | FAIL | FAIL | 0.6250 | yes | 0.0000 | 0.3750 |
| P3_rho2 | PASS | FAIL | FAIL | 0.6250 | yes | 0.0000 | 0.3750 |
| P4 | PASS | PASS | FAIL | 0.0000 | yes | 0.0000 | 0.7500 |

## Round-2 outcome

| PDE | C2-deep | C2-shallow | C3-A | C3-B | C3-C | class-A dominated |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | PASS | PASS | PASS | PASS | PASS | no |
| P2_a4 | PASS | PASS | PASS | PASS | PASS | no |
| P2_a8 | PASS | PASS | PASS | PASS | PASS | no |
| P3_rho1 | PASS | PASS | PASS | PASS | PASS | no |
| P3_rho2 | PASS | PASS | PASS | PASS | PASS | no |
| P4 | PASS | PASS | FAIL | PASS | FAIL | yes |

The pooled C2-graded Spearman correlation is 0.8115 (PASS). C-tight-1=PASS, C-tight-2=FAIL, C-loc-1=PASS, C-loc-2=FAIL, and C-dim=PASS.

## Pre-registered criteria and mechanical verdicts

### G5-R3: PASS

> R3 proceeds only if the violation decreases under refinement and the finest-grid violation is below 1e-5.

decreasing=True; finest=7.63993e-06

### C2-deep: PASS

> C2-deep: primary has `Gc >= 0.5` in >= 75% of deep cells, for every PDE.

P1=1.000; P2_a4=1.000; P2_a8=1.000; P3_rho1=1.000; P3_rho2=1.000; P4=1.000

### C2-shallow: PASS

> C2-shallow (boundary prediction): for P2/P3, primary has `Gc < 0.5` in >= 50% of shallow cells. For P1 and P4, `Gc >= 0.5` in >= 75% of shallow cells.

P1=0.917 (Gc >= 0.5); P2_a4=1.000 (Gc < 0.5); P2_a8=1.000 (Gc < 0.5); P3_rho1=1.000 (Gc < 0.5); P3_rho2=1.000 (Gc < 0.5); P4=1.000 (Gc >= 0.5)

### C2-graded: PASS

> C2-graded: pooled over all P2/P3 cells, Spearman correlation between `Gc` and `log(raw/oracle_state)` is >= 0.5.

rho=0.811493

### C3-A: FAIL

> C3-A: primary wins >= 7/10 vs centre in >= 75% of cells.

P1=1.000; P2_a4=1.000; P2_a8=1.000; P3_rho1=1.000; P3_rho2=1.000; P4=0.000

### C3-B: PASS

> C3-B: primary wins >= 7/10 vs best class-B control in >= 75% of cells.

P1=1.000; P2_a4=0.875; P2_a8=0.875; P3_rho1=0.875; P3_rho2=1.000; P4=1.000

### C3-C: FAIL

> C3-C: primary wins >= 7/10 vs best class-C control in >= 75% of cells.

P1=0.750; P2_a4=1.000; P2_a8=1.000; P3_rho1=1.000; P3_rho2=1.000; P4=0.000

### P4-round1-segment-C: PASS

> P4 with the round-1 `segment` is predicted to fail C (round-1 replication).

prediction confirmed

### C-tight-1: PASS

> C-tight-1: test skill decreases monotonically in theta (mean over reps) in >= 75% of the tightness-family cells.

fraction=1.000

### C-tight-2: FAIL

> C-tight-2: tight_segment passes C3-A, C3-B and C3-C (R2 thresholds).

all P4 classes pass=False

### C-loc-1: PASS

> C-loc-1: for every dose >= 1, the central bin `[-0.25, 0.25)` has the largest mean |dose-induced generator error|.

dose 1: central=5.616, outer=5.227; dose 2: central=11.36, outer=10.95; dose 4: central=22.85, outer=22.42

### C-loc-2: FAIL

> C-loc-2: in raw E1 rows, the normalised bias is largest in the central bin.

central=0.4265, outer=0.4279

### R5-direction: PASS

> Directional prediction: soft has lower mean skill than hard in >= 50% of shallow cells, and hard has lower mean skill than soft in >= 75% of deep cells.

shallow=0.875; deep=1.000

### C-dim: PASS

> C-dim: `max_d(box/oracle) / min_d(box/oracle) <= 1.25`.

max/min=1.21218

### all-C3-predictions: FAIL

> Predictions: P1–P3 pass A, B, C. P4 with `tight_segment` passes A, B, C.

joint summary of the separately evaluated C3 classes

## All failures

- C3-A: P1=1.000; P2_a4=1.000; P2_a8=1.000; P3_rho1=1.000; P3_rho2=1.000; P4=0.000
- C3-C: P1=0.750; P2_a4=1.000; P2_a8=1.000; P3_rho1=1.000; P3_rho2=1.000; P4=0.000
- C-tight-2: all P4 classes pass=False
- C-loc-2: central=0.4265, outer=0.4279
- all-C3-predictions: joint summary of the separately evaluated C3 classes
- P4 C3-A fraction=0.0000.
- P4 C3-C fraction=0.0000.

## Protocol and provenance

Base seed: `20261107`. The fixed 1,200-point data sets, including the 20% validation split, were regenerated from this seed. Tree seeds are `SeedSequence([20261107,d,n,M,rep,chunk_index])`; methods within a study use identical chunk sizes and paired trees. Every numerical array is float64. Value accuracy alone determines verdicts; gradient diagnostics are retained but do not determine a claim.

The round-2 pre-registration was frozen at commit `f83ce81fbeaff6b7c12aeec90906fc6e6124a71e` before implementation or computation. The containment artifact records code commit `421c4dbd33a9dce9af2088f83db4fb05f3135953` on branch `ir-mlp-mechanism-suite-r2`. The unchanged FullHistoryMLP source has SHA-256 `f4581babdd5ceed8354337ee7783112eb6b2590facc2a0115e7c02ccd50990b8`.

The P1/P4 directions and all equation parameters are inherited unchanged from round 1. No round-1 result is reused as a round-2 random draw; only mathematically identical rows within round 2 are hard-linked and counted once in work accounting.

## R3: tight P4 certificate and containment gate

For `v=psi_s`, differentiation gives `v_t + v_ss + b v_s = 0` with `b=-lambda_f sign(v)` and `|b|<=lambda_f`. Monotonicity of `tanh` and constant-drift comparison give

```text
v_minus(t,s) <= psi_s(t,s) <= v_plus(t,s)
v_pm(t,s) = E tanh(beta*(s +/- lambda_f*(T-t) + sqrt(2*(T-t))*xi))
xi ~ N(0,1).
```

The bounds were computed with 80-node Gauss-Hermite quadrature. The recursive solver uses an audited cubic cache of those quadrature values; it clips exactly to the interpolated endpoints and applies no inward safety margin. The sharp-bound G5 audit used 100,000 fresh test-distribution points plus 10,000 points with `T-t in [0,0.1]` and `s in [-0.25,0.25]`.

The first implementation audited only the first four levels of the already-existing six-level P4 reference hierarchy and stopped at `3.023e-5`. That failed audit is preserved under `audit_history/`. Before any round-2 MLP run, the implementation was corrected to use all six pre-existing levels; the `1e-5` threshold and every scientific criterion remained unchanged.

| n_space | n_steps | ds | maximum_violation | positive_violation_count |
| --- | --- | --- | --- | --- |
| 2401 | 512 | 0.0050 | 2.851e-04 | 4191 |
| 4801 | 1024 | 0.0025 | 1.337e-04 | 3669 |
| 9601 | 2048 | 0.0013 | 6.324e-05 | 3270 |
| 19201 | 4096 | 6.250e-04 | 3.023e-05 | 2916 |
| 38401 | 8192 | 3.125e-04 | 1.484e-05 | 2607 |
| 76801 | 16384 | 1.563e-04 | 7.640e-06 | 2376 |

The accepted Richardson reference has maximum bound violation 1.466e-05. The quadrature-cache interpolation audit maximum is 5.447e-08.

### Tightness family

| d | n | M | monotone_nonincreasing | mean_skill_theta_0 | mean_skill_theta_0.25 | mean_skill_theta_0.5 | mean_skill_theta_0.75 | mean_skill_theta_1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 3 | 6 | PASS | 0.3371 | 0.3254 | 0.3111 | 0.2882 | 0.2458 |
| 20 | 4 | 6 | PASS | 0.2156 | 0.2096 | 0.2015 | 0.1875 | 0.1316 |
| 50 | 3 | 6 | PASS | 0.3922 | 0.3785 | 0.3613 | 0.3318 | 0.2812 |
| 50 | 4 | 6 | PASS | 0.2474 | 0.2399 | 0.2303 | 0.2140 | 0.1538 |

[Per-repetition R3 rows](../results/mechanism_suite_r2/r3_tightness.csv) · [tightness figure](../results/mechanism_suite_r2/figures/r3_tightness.png)

## R1: C2 regime hypothesis

| PDE | deep Gc>=0.5 fraction | deep verdict | shallow rule | shallow fraction | shallow verdict |
| --- | --- | --- | --- | --- | --- |
| P1 | 1.0000 | PASS | Gc >= 0.5 | 0.9167 | PASS |
| P2_a4 | 1.0000 | PASS | Gc < 0.5 | 1.0000 | PASS |
| P2_a8 | 1.0000 | PASS | Gc < 0.5 | 1.0000 | PASS |
| P3_rho1 | 1.0000 | PASS | Gc < 0.5 | 1.0000 | PASS |
| P3_rho2 | 1.0000 | PASS | Gc < 0.5 | 1.0000 | PASS |
| P4 | 1.0000 | PASS | Gc >= 0.5 | 1.0000 | PASS |

Pooled P2/P3 Spearman rho=0.8115; the threshold was 0.5.

[All R1 rows](../results/mechanism_suite_r2/r1_regime.csv) · [regime figure](../results/mechanism_suite_r2/figures/r1_gc_regime.png)

## R2: corrected C3 comparator classes

| PDE | C3-A fraction | C3-B fraction | C3-C fraction | centre better fraction | class-A dominated |
| --- | --- | --- | --- | --- | --- |
| P1 | 1.0000 | 1.0000 | 0.7500 | 0.0000 | no |
| P2_a4 | 1.0000 | 0.8750 | 1.0000 | 0.0000 | no |
| P2_a8 | 1.0000 | 0.8750 | 1.0000 | 0.0000 | no |
| P3_rho1 | 1.0000 | 0.8750 | 1.0000 | 0.0000 | no |
| P3_rho2 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | no |
| P4 | 0.0000 | 1.0000 | 0.0000 | 1.0000 | yes |

Class B and C choices were selected only by mean validation skill and then frozen for the test split. The full choices, including every candidate validation score, are in `r2_tuning_choices.json`.

[All R2 rows](../results/mechanism_suite_r2/r2_controls.csv) · [class heatmap](../results/mechanism_suite_r2/figures/r2_comparator_classes.png)

## R4: homoscedastic rectification localization

### Dose-induced generator error

| dose | bin_index | bin_label | is_central_bin | generator_count | mean_absolute_dose_generator_error | mean_dose_generator_error |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5000 | 0 | [-inf,-2) | FAIL | 39715 | 1.9661 | -1.9661 |
| 0.5000 | 1 | [-2,-1) | FAIL | 513081 | 2.0650 | -2.0650 |
| 0.5000 | 2 | [-1,-0.25) | FAIL | 1649956 | 2.3764 | -2.3764 |
| 0.5000 | 3 | [-0.25,0.25) | PASS | 1525658 | 2.7466 | -2.7466 |
| 0.5000 | 4 | [0.25,1) | FAIL | 1593935 | 2.3874 | -2.3874 |
| 0.5000 | 5 | [1,2) | FAIL | 476302 | 2.1060 | -2.1060 |
| 0.5000 | 6 | [2,inf) | FAIL | 33353 | 2.0159 | -2.0159 |
| 1.0000 | 0 | [-inf,-2) | FAIL | 39715 | 4.6635 | -4.6635 |
| 1.0000 | 1 | [-2,-1) | FAIL | 513081 | 4.8143 | -4.8143 |
| 1.0000 | 2 | [-1,-0.25) | FAIL | 1649956 | 5.2094 | -5.2094 |
| 1.0000 | 3 | [-0.25,0.25) | PASS | 1525658 | 5.6162 | -5.6162 |
| 1.0000 | 4 | [0.25,1) | FAIL | 1593935 | 5.2272 | -5.2272 |
| 1.0000 | 5 | [1,2) | FAIL | 476302 | 4.8980 | -4.8980 |
| 1.0000 | 6 | [2,inf) | FAIL | 33353 | 4.7653 | -4.7653 |
| 2.0000 | 0 | [-inf,-2) | FAIL | 39715 | 10.1859 | -10.1859 |
| 2.0000 | 1 | [-2,-1) | FAIL | 513081 | 10.4197 | -10.4197 |
| 2.0000 | 2 | [-1,-0.25) | FAIL | 1649956 | 10.9215 | -10.9215 |
| 2.0000 | 3 | [-0.25,0.25) | PASS | 1525658 | 11.3589 | -11.3589 |
| 2.0000 | 4 | [0.25,1) | FAIL | 1593935 | 10.9517 | -10.9517 |
| 2.0000 | 5 | [1,2) | FAIL | 476302 | 10.5867 | -10.5867 |
| 2.0000 | 6 | [2,inf) | FAIL | 33353 | 10.3901 | -10.3901 |
| 4.0000 | 0 | [-inf,-2) | FAIL | 39715 | 21.2985 | -21.2985 |
| 4.0000 | 1 | [-2,-1) | FAIL | 513081 | 21.6867 | -21.6867 |
| 4.0000 | 2 | [-1,-0.25) | FAIL | 1649956 | 22.3692 | -22.3692 |
| 4.0000 | 3 | [-0.25,0.25) | PASS | 1525658 | 22.8461 | -22.8461 |
| 4.0000 | 4 | [0.25,1) | FAIL | 1593935 | 22.4237 | -22.4237 |
| 4.0000 | 5 | [1,2) | FAIL | 476302 | 22.0192 | -22.0192 |
| 4.0000 | 6 | [2,inf) | FAIL | 33353 | 21.7062 | -21.7062 |

### Raw R1 normalized bias

| bin_index | bin_label | is_central_bin | generator_count | generator_bias | z_error_rms | normalized_absolute_bias |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | [-inf,-2) | FAIL | 39715 | -2.0907 | 4.9401 | 0.4232 |
| 1 | [-2,-1) | FAIL | 513081 | -2.1178 | 5.0035 | 0.4233 |
| 2 | [-1,-0.25) | FAIL | 1649956 | -1.7951 | 4.2779 | 0.4196 |
| 3 | [-0.25,0.25) | PASS | 1525658 | -1.5420 | 3.6155 | 0.4265 |
| 4 | [0.25,1) | FAIL | 1593935 | -1.8168 | 4.3196 | 0.4206 |
| 5 | [1,2) | FAIL | 476302 | -2.1607 | 5.0861 | 0.4248 |
| 6 | [2,inf) | FAIL | 33353 | -2.1032 | 4.9151 | 0.4279 |

[Expanded bin rows](../results/mechanism_suite_r2/r4_localization.csv) · [localization figure](../results/mechanism_suite_r2/figures/r4_localization.png)

## R5: hard projection versus soft contraction

Across shallow cells, soft contraction has lower mean skill in 0.8750 of cells (threshold 0.5). Across deep cells, hard projection has lower mean skill in 1.0000 of cells (threshold 0.75).

| PDE | soft better shallow fraction | hard better deep fraction |
| --- | --- | --- |
| P2_a4 | 1.0000 | 1.0000 |
| P2_a8 | 0.7500 | 1.0000 |
| P3_rho1 | 0.8750 | 1.0000 |
| P3_rho2 | 0.8750 | 1.0000 |

Each cell's mean paired difference and deterministic 10,000-draw paired bootstrap interval is in the CSV. This study is descriptive with a pre-registered direction and does not determine a C-claim.

[R5 cell rows](../results/mechanism_suite_r2/r5_soft_vs_hard.csv) · [paired intervals](../results/mechanism_suite_r2/figures/r5_soft_vs_hard.png)

## R6: dimension stability

| d | box_mean_skill | oracle_state_mean_skill | box_over_oracle_state |
| --- | --- | --- | --- |
| 20 | 0.1399 | 0.0287 | 4.8665 |
| 50 | 0.1377 | 0.0292 | 4.7184 |
| 100 | 0.1431 | 0.0291 | 4.9112 |
| 200 | 0.1569 | 0.0299 | 5.2569 |
| 400 | 0.1639 | 0.0287 | 5.7196 |

`max/min=1.2122` against the pre-registered threshold 1.25: **PASS**.

[All R6 rows](../results/mechanism_suite_r2/r6_dimension.csv) · [dimension figure](../results/mechanism_suite_r2/figures/r6_dimension.png)

## R7: what remains between box and oracle

| pde | mean_per_rep_skill_gap | mean_repetition_averaged_skill_gap | mean_box_generator_signed_error | mean_oracle_generator_signed_error |
| --- | --- | --- | --- | --- |
| P1 | 0.1108 | 0.1166 | 0.1607 | 0.0000 |
| P2_a4 | 0.1723 | 0.2088 | -0.1954 | 0.0000 |
| P2_a8 | 0.2060 | 0.2458 | -0.1139 | 0.0000 |
| P3_rho1 | 0.1709 | 0.2119 | -0.2035 | 0.0000 |
| P3_rho2 | 0.1740 | 0.2170 | -0.2170 | 0.0000 |
| P4 | 0.2609 | 0.3343 | -0.2016 | 0.0000 |

The mean per-repetition box-minus-oracle skill gap is 0.177; the mean gap after averaging the ten predictions is 0.2143. The repetition-averaged gap is 1.211 times the mean per-repetition gap, so averaging does not close the gap; it enlarges it by 21.1%. Because normalised RMSE is nonlinear, these two gaps are not an additive bias--variance decomposition. Their persistence under averaging nevertheless supports a mainly systematic, bias-like residual rather than one dominated by Monte Carlo variation. The mean signed child-generator errors are -0.1062 for box and 0 for oracle_state. These numbers answer the question descriptively without introducing a post-hoc pass threshold.

[Paired R7 rows](../results/mechanism_suite_r2/r7_gap.csv) · [gap figure](../results/mechanism_suite_r2/figures/r7_gap.png)

## Work and integrity accounting

The six computational outputs contain 13,250 formal method-repetition rows, including 1,400 exact within-round hard-link reuses and 11,850 newly computed rows. Non-duplicated work totals 4,919,004,000 generator calls, 4,933,224,000 recursive states, 68,248,752,000 stochastic samples, and 27.69 summed worker-hours. Rows with any non-finite state/generator count: 0.

The manuscript, FullHistoryMLP source, and all round-1 results were not modified.
