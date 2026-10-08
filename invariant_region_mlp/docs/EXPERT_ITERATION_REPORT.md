# Noise-aware defect correction and expert iteration: decision-subset report

> Status: Study A decision subset only. The user explicitly limited this run to d={20,100,160}; Study B and omitted Study-A dimensions/cells are not evaluated, not negative results.

## Outcome table

| Gate / criterion | Verdict |
|---|---|
| A-G1 | PASS |
| A-G0 | REPORTED |
| A-1 | NOT EVALUATED (USER-DIRECTED SUBSET) |
| A-2 | NOT EVALUATED (USER-DIRECTED SUBSET) |
| A-3 | PASS |
| A-4 | NOT EVALUATED (USER-DIRECTED SUBSET) |
| A-5 | EXPLORATORY |
| B-pilot | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |
| B-1 | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |
| B-2 | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |
| B-3 | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |
| B-4 | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |
| B-5 | NOT EVALUATED (USER-DIRECTED DECISION SUBSET) |

## Decision-point result

Median test relative L2 error after first taking the median over five MC repetitions within each network seed, then the median over three network seeds:

| method | d=20 | d=100 | d=160 |
|---|---:|---:|---:|
| surrogate | 0.777174 | 0.845138 | 0.857251 |
| f_zero | 0.0232327 | 0.00299977 | 0.00171256 |
| mlp | 0.548921 | 2.19655 | 3.38165 |
| mlp_clip | 0.548921 | 2.19433 | 2.95152 |
| scasml | 0.741433 | 0.82176 | 0.83605 |
| scasml_noclip | 0.552094 | 2.1976 | 3.38264 |
| path | 0.0151147 | 0.00392384 | 0.00280498 |
| path_clip | 0.741433 | 0.82176 | 0.83605 |
| oracle_state | 0.0152096 | 0.00389162 | 0.00280243 |

Decision-point reading:

- Pathwise correction improves over the surrogate by 51.4x, 215.4x, and 305.6x at d=20,100,160, respectively.
- It improves over clipped SCaSML by 49.1x, 209.4x, and 298.1x, and is essentially coincident with the oracle-state lower bound.
- The benchmark itself becomes nearly linear with dimension: f=0 is 2.32% at d=20, 0.300% at d=100, and 0.171% at d=160. Pathwise correction is better than f=0 at d=20 but worse than f=0 at d=100 and d=160.
- Thus the requested dramatic repair is reproduced, while the nonlinear benchmark is simultaneously exposed as weak in the high-dimensional regime. Both facts are retained; neither is treated as cancelling the other.

## Gates and criteria (verbatim)

- **A-G0:** Report the relative L2 error of the f=0 solution. Result: `{20: 0.023232732569262535, 100: 0.002999770977715469, 160: 0.0017125607857761578}`.
- **A-G1:** MC estimated SE of u <=1e-4 at every test point. Verdict: **PASS**. Details: `{20: True, 100: True, 160: True}`.

### A-1

`mlp` relative L2 >= 1 at d >= 100, and `scasml` / `surrogate` error ratio >= 0.5 at d = 160 in the primary cell. (If not reproduced, report and continue.)

Verdict: **NOT EVALUATED (USER-DIRECTED SUBSET)**.

Mechanical details: `{"missing_dimensions": [120, 140], "mlp_median_errors": {"100": 2.196547969964566, "160": 3.381651206214368, "20": 0.5489206319579321}, "restricted_subset_pass": true, "scasml_over_surrogate_d160": 0.9752685490285039, "verdict": "NOT EVALUATED (USER-DIRECTED SUBSET)"}`

### A-2

In the primary cell, `path` <= 0.2 x `surrogate` and `path` <= 0.5 x `scasml` at every d in {100, ..., 160}.

Verdict: **NOT EVALUATED (USER-DIRECTED SUBSET)**.

Mechanical details: `{"missing_dimensions": [120, 140], "path_over_scasml": {"100": 0.004772325158735859, "160": 0.0033573596478880543}, "path_over_surrogate": {"100": 0.004640377455920984, "160": 0.003274274885338305}, "restricted_subset_pass": true, "verdict": "NOT EVALUATED (USER-DIRECTED SUBSET)"}`

### A-3

The ratio `path / surrogate` at d = 160 is <= 1.3 x its value at d = 20, while the ratio `scasml_noclip / surrogate` at d = 160 is >= 1.5 x its value at d = 20.

Verdict: **PASS**.

Mechanical details: `{"path_ratio_growth_d160_over_d20": 0.16835831513707536, "scasml_noclip_ratio_growth_d160_over_d20": 5.547045998389812, "verdict": "PASS"}`

### A-4

Mean generator bias for `scasml_noclip` is negative at every d and its magnitude at d = 160 is >= 3 x that at d = 20; for `path` the magnitude at d = 160 is <= 1.5 x that at d = 20.

Verdict: **NOT EVALUATED (USER-DIRECTED SUBSET)**.

Mechanical details: `{"median_generator_bias": {"100": {"path": -0.01089984347861444, "scasml_noclip": -100.93433039202368}, "160": {"path": -0.00806906681046694, "scasml_noclip": -201.05102901354684}, "20": {"path": -0.030164641876552302, "scasml_noclip": -9.637495860290818}}, "missing_dimensions": [50, 120, 140], "path_bias_growth": 0.26750083238147837, "restricted_subset_pass": true, "scasml_noclip_bias_growth": 20.86133492849874, "verdict": "NOT EVALUATED (USER-DIRECTED SUBSET)"}`

### A-5

Plot `path` and `scasml_noclip` error against surrogate error across the three checkpoints; report whether `scasml_noclip` is ever worse than the surrogate alone (exploratory).

Verdict: **EXPLORATORY**.

Mechanical details: `{"comparisons": 27, "ever_worse": true, "scasml_noclip_worse_count": 18, "verdict": "EXPLORATORY"}`

### B-pilot

Proceed only if `EI-path` test error after round 3 is <= 0.5 x after round 0 on B2.

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

### B-1

`EI-path` error after round 3 <= 0.2 x its error after round 0, on B1, B2, B3.

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

### B-2

On B4, `EI-bismut` error after round 3 >= its error after round 1, while `EI-path` satisfies B-1 on B4.

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

### B-3

At the wall-clock of `EI-path` round 6, `EI-path` error <= 0.5 x the best DPI error within the same wall-clock, on at least 2 of B1-B3.

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

### B-4

Report whether DPI at 4x the wall-clock reaches `EI-path`'s round-6 error (no verdict).

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

### B-5

For the final `EI-path` network, corrected error <= 0.3 x network error on B1-B3.

Verdict: **NOT EVALUATED (USER-DIRECTED DECISION SUBSET)**.

Mechanical details: `{"details": "The user explicitly requested stopping after the inexpensive Study-A decision point.", "verdict": "NOT EVALUATED (USER-DIRECTED DECISION SUBSET)"}`

## Exact-Laplacian control at d=100

`{"path": {"exact": 0.0039893386074845, "hutchinson": 0.0039896849180583, "relative_difference": 8.68090197083018e-05}, "scasml_noclip": {"exact": 2.1975964925264533, "hutchinson": 2.1975984797469112, "relative_difference": 9.042699443343632e-07}, "status": "REPORTED"}`

## Failures, limitations, and deviations

- The run is intentionally the user-requested decision subset d={20,100,160}. Dimensions 50, 120, and 140, secondary cells, and all of Study B are marked not evaluated.
- The frozen prompt's interpretation is lambda=1, zero drift, T=0.5, and unit-ball test points. The audited public SCaSML LQG markdown says T=1; this run follows the frozen prompt.
- The public SCaSML training code samples d/4 coordinate Hessian entries and also subsamples the gradient norm. The frozen prompt specifically says Hutchinson d/4 for the Laplacian; this implementation uses Rademacher Hutchinson for the Laplacian and the full gradient norm.
- Reported test errors use the required adaptive antithetic MC reference. A scaled Gauss--Laguerre Hopf--Cole integral supplies truth only for recursive child-state diagnostics/oracle calls and is audited against the MC reference.
- Antithetic reference points use common random numbers within each eight-point chunk. Each point retains the correct marginal iid Gaussian sample and its own SE; correlations across test-point errors do not enter any registered threshold.
- The pathwise change is restricted to the terminal gradient estimator; level terms retain the unchanged Bismut estimator.

## Provenance and totals

- Frozen preregistration commit: `85eaac9f7bfc6ba7c0ed291fd258cb3be86759c6`
- Code commit at analysis: `4b057e24bbeb8886c1b70abb3e446b65567127c1`
- Rows: 685; aggregate inference wall-clock: 1224.9 s
- Environment: `{"code_commit": "4b057e24bbeb8886c1b70abb3e446b65567127c1", "cuda_available": true, "cuda_runtime": "12.8", "device": "NVIDIA GeForce RTX 5070 Ti", "logical_cpus": 16, "numpy": "2.3.4", "platform": "Windows-11-10.0.26200-SP0", "python": "3.12.11", "torch": "2.11.0+cu128"}`
- Reference summary: `{"100": {"dimension": 100, "f_zero_relative_l2_test": 0.002999770977715469, "gate_A_G1": true, "max_antithetic_pairs": 2588672, "max_u_standard_error_all": 9.99956505155432e-05, "max_u_standard_error_test": 9.99956505155432e-05, "mc_vs_quadrature_u_relative_l2_test": 2.3013319211974058e-05, "mc_vs_quadrature_z_relative_l2_test": 0.00018424117307634912, "min_antithetic_pairs": 2506752, "quadrature_role": "independent audit and recursive-child truth; reported errors use MC", "reference_file": "results/expert_iteration/artifacts/references/lqg_d100.npz", "seed_scheme": "SeedSequence([20261101,1,1,d,0,20,point_chunk])", "terminal_samples_total": 6153043968, "tolerance": 0.0001, "wall_clock_seconds": 852.1834764999803}, "160": {"dimension": 160, "f_zero_relative_l2_test": 0.0017125607857761578, "gate_A_G1": true, "max_antithetic_pairs": 1622016, "max_u_standard_error_all": 9.999607638469564e-05, "max_u_standard_error_test": 9.999520598739125e-05, "mc_vs_quadrature_u_relative_l2_test": 2.0749524276421675e-05, "mc_vs_quadrature_z_relative_l2_test": 0.00018580624704691762, "min_antithetic_pairs": 1589248, "quadrature_role": "independent audit and recursive-child truth; reported errors use MC", "reference_file": "results/expert_iteration/artifacts/references/lqg_d160.npz", "seed_scheme": "SeedSequence([20261101,1,1,d,0,20,point_chunk])", "terminal_samples_total": 3861905408, "tolerance": 0.0001, "wall_clock_seconds": 820.9502270000521}, "20": {"dimension": 20, "f_zero_relative_l2_test": 0.023232732569262535, "gate_A_G1": true, "max_antithetic_pairs": 14336000, "max_u_standard_error_all": 9.999969285980934e-05, "max_u_standard_error_test": 9.999969285980934e-05, "mc_vs_quadrature_u_relative_l2_test": 2.936066874793782e-05, "mc_vs_quadrature_z_relative_l2_test": 0.000182110495792915, "min_antithetic_pairs": 12255232, "quadrature_role": "independent audit and recursive-child truth; reported errors use MC", "reference_file": "results/expert_iteration/artifacts/references/lqg_d20.npz", "seed_scheme": "SeedSequence([20261101,1,1,d,0,20,point_chunk])", "terminal_samples_total": 32313180160, "tolerance": 0.0001, "wall_clock_seconds": 1460.2730165999383}}`

## Figures

- `results/expert_iteration/figures/A_error_vs_d.png`
- `results/expert_iteration/figures/A_generator_bias_vs_d.png`
- `results/expert_iteration/figures/A_error_vs_surrogate.png`
