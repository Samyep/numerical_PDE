# Non-ridge spectral HJB benchmark: pre-registered report

## Outcome table

| Block | N-1 | N-2 | N-3 | N-4 | N-5 |
|---|---|---|---|---|---|
| Stage A: S0 + S05 | PASS | PASS | PASS | PASS | PASS |
| Stage B: S1 | PASS | PASS | PASS | PASS | PASS |

## Stage A criteria (verbatim)

### N-1

> **N-1 (pathwise matters off-ridge):** `path` skill <= 0.5 x `raw` skill at (2,64) and at (3,16), in all 8 (family, d, cell) combinations. A non-finite `raw` mean counts as satisfied.

**Verdict: PASS (8/8; required 8).**

### N-2

> **N-2 (baseline off-ridge):** `path_base` skill <= 0.5 x `path` skill at (3,16), in all 4 (family, d) combinations.

**Verdict: PASS (4/4; required 4).**

### N-3

> **N-3 (double complements baseline off-ridge):** `double_path_base` skill <= `double_path` skill at (3,16) and (3,24), in at least 7 of the 8 (family, d, M) combinations.

**Verdict: PASS (8/8; required 7).**

### N-4

> **N-4 (n = 3 below the n = 2 floor at d = 400):** for each family at d = 400, the minimum over the n = 3 cells and methods {`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over the n = 2 `path` cells (M = 16, 64, 256). PASS if it holds for both families. Report the same quantity at d = 100 (no verdict).

**Verdict: PASS (2/2; required 2).**

### N-5

> **N-5 (dimension robustness):** `double_path_base` skill ratio d = 400 / d = 100 <= 1.3 at (3,16) and (3,24), for both families (4 comparisons; PASS only if all hold).

**Verdict: PASS (4/4; required 4).**

## Stage A mechanical details

### N-1 comparisons

| family | d | cell/M | ratio | met |
|---|---:|---|---:|---|
| S0 | 100 | (2,64) | 0.0278853 | yes |
| S0 | 100 | (3,16) | 0.00487061 | yes |
| S0 | 400 | (2,64) | 0.00276887 | yes |
| S0 | 400 | (3,16) | 5.95424e-05 | yes |
| S05 | 100 | (2,64) | 0.0406661 | yes |
| S05 | 100 | (3,16) | 0.00683939 | yes |
| S05 | 400 | (2,64) | 0.00476087 | yes |
| S05 | 400 | (3,16) | 0.000139425 | yes |

### N-2 comparisons

| family | d | cell/M | ratio | met |
|---|---:|---|---:|---|
| S0 | 100 | (3,16) | 0.0124585 | yes |
| S0 | 400 | (3,16) | 0.00308496 | yes |
| S05 | 100 | (3,16) | 0.0129164 | yes |
| S05 | 400 | (3,16) | 0.00467752 | yes |

### N-3 comparisons

| family | d | cell/M | ratio | met |
|---|---:|---|---:|---|
| S0 | 100 | (3,16) | 0.351773 | yes |
| S0 | 100 | (3,24) | 0.362526 | yes |
| S0 | 400 | (3,16) | 0.206381 | yes |
| S0 | 400 | (3,24) | 0.223734 | yes |
| S05 | 100 | (3,16) | 0.331338 | yes |
| S05 | 100 | (3,24) | 0.364104 | yes |
| S05 | 400 | (3,16) | 0.203319 | yes |
| S05 | 400 | (3,24) | 0.229196 | yes |

### N-5 comparisons

| family | d | cell/M | ratio | met |
|---|---:|---|---:|---|
| S0 | 400/100 | (3,16) | 0.76306 | yes |
| S0 | 400/100 | (3,24) | 0.768212 | yes |
| S05 | 400/100 | (3,16) | 0.785115 | yes |
| S05 | 400/100 | (3,24) | 0.772516 | yes |

### N-4 floor comparison

| family | d | n=2 floor | best n=3 | ratio | verdict scope | met |
|---|---:|---|---|---:|---|---|
| S0 | 100 | M=256; 0.0168812 | double_path_base M=24; 0.0105186 | 0.623096 | reported only | yes |
| S0 | 400 | M=256; 0.0112016 | double_path_base M=24; 0.00808051 | 0.721371 | yes | yes |
| S05 | 100 | M=256; 0.0285142 | double_path_base M=24; 0.0121543 | 0.426253 | reported only | yes |
| S05 | 400 | M=256; 0.0180238 | double_path_base M=24; 0.00938937 | 0.520944 | yes | yes |

## Stage B: S1 pre-registered block

The same N-1 through N-5 inequalities and thresholds are applied to S1 alone; N-3 requires all 4 S1 comparisons.

- **N-1: PASS** (4/4; required 4).
- **N-2: PASS** (2/2; required 2).
- **N-3: PASS** (4/4; required 4).
- **N-4: PASS** (1/1; required 1).
- **N-5: PASS** (2/2; required 2).

## Extended M=32 / M=512 floor re-evaluation

| family | d | initial ratio | extended ratio | extended attaining cells | <= 0.9 |
|---|---:|---:|---:|---|---|
| S0 | 100 | 0.623096 | 0.555441 | n2 M=512; n3 double_path_base M=32 | yes |
| S0 | 400 | 0.721371 | 0.702978 | n2 M=512; n3 double_path_base M=32 | yes |
| S05 | 100 | 0.426253 | 0.323773 | n2 M=512; n3 double_path_base M=32 | yes |
| S05 | 400 | 0.520944 | 0.425758 | n2 M=512; n3 double_path_base M=32 | yes |
| S1 | 100 | 0.241169 | 0.192854 | n2 M=512; n3 double_path_base M=32 | yes |
| S1 | 400 | 0.240727 | 0.197262 | n2 M=512; n3 double_path_base M=32 | yes |

## Calibration (reported, no verdict)

| PDE | d | std u* | mean f(z*) | std f(z*) | RMS norm(z*) | norm(mean z*) | heat skill | z var outside top-1 | outside top-10 | Hessian Fro | Hessian eff. rank |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S0 | 100 | 1.32382 | -3.23288 | 0.180904 | 2.54279 | 2.45609 | 0.360544 | 0.982471 | 0.840724 | 0.891843 | 86.0472 |
| S0 | 400 | 1.71717 | -3.17431 | 0.0873085 | 2.51965 | 2.4338 | 0.278514 | 0.993239 | 0.935977 | 0.892033 | 344.281 |
| S05 | 100 | 1.22642 | -3.39732 | 0.45763 | 2.60665 | 2.51896 | 0.41161 | 0.790364 | 0.408895 | 0.891044 | 16.0789 |
| S05 | 400 | 1.49684 | -3.28196 | 0.353175 | 2.56202 | 2.47378 | 0.321923 | 0.825105 | 0.524728 | 0.889541 | 26.9013 |
| S1 | 100 | 1.27805 | -3.84981 | 1.33407 | 2.77482 | 2.67751 | 0.434908 | 0.333672 | 0.0414929 | 0.877263 | 2.82224 |
| S1 | 400 | 1.31093 | -3.92665 | 1.33002 | 2.80237 | 2.70901 | 0.445532 | 0.340141 | 0.0456008 | 0.880521 | 2.81546 |
| P1 | 100 | 1.33351 | -4.77001 | 2.15948 | 3.08869 | 2.99384 | 0.505741 | 1.52623e-29 | 1.99579e-32 | 0.862765 | 1 |
| P1 | 400 | 1.3318 | -4.64719 | 2.16906 | 3.04867 | 2.95181 | 0.507433 | 8.71473e-30 | 7.8447e-32 | 0.861822 | 1 |

## Equation validation

- **V-1: PASS**; checks 6, failures 0.
- **V-2: PASS**; checks 2, failures 0.
- **V-3: PASS**; checks 7, failures 0.
- **V-4: PASS**; checks 6, failures 0.
- **V-5: PASS**; checks 6, failures 0.
- **V-6: PASS**; checks 8, failures 0.

The GPU rotation uses the Makhoul FFT DCT. Timings are for one 1,024-point forward DCT:

| family | d | milliseconds/call | max equation relative error |
|---|---:|---:|---:|
| S0 | 100 | 0.213095 | 4.7933e-16 |
| S0 | 400 | 0.203745 | 3.11144e-16 |
| S05 | 100 | 0.21195 | 3.8201e-16 |
| S05 | 400 | 0.204755 | 3.24989e-16 |
| S1 | 100 | 0.208335 | 3.42408e-16 |
| S1 | 400 | 0.207975 | 3.60948e-16 |

## Ratio-level comparison with Study F (reported, unpaired)

| family | d | source | path_base/path M16 | double_base/double M16 | best n3/n2 floor |
|---|---:|---|---:|---:|---:|
| S0 | 100 | nonridge GPU | 0.0124585 | 0.351773 | 0.623096 |
| S0 | 400 | nonridge GPU | 0.00308496 | 0.206381 | 0.721371 |
| S05 | 100 | nonridge GPU | 0.0129164 | 0.331338 | 0.426253 |
| S05 | 400 | nonridge GPU | 0.00467752 | 0.203319 | 0.520944 |
| S1 | 100 | nonridge GPU | 0.0275018 | 0.357521 | 0.241169 |
| S1 | 400 | nonridge GPU | 0.0243345 | 0.194788 | 0.240727 |
| P1 | 100 | Study F CPU; ratios only, unpaired | 0.0442866 | 0.376181 | 0.2125 |
| P1 | 400 | Study F CPU; ratios only, unpaired | 0.0438045 | 0.221509 | 0.215155 |
| MR | 100 | Study F CPU; ratios only, unpaired | 0.0873141 | 0.323269 | 0.243863 |
| MR | 400 | Study F CPU; ratios only, unpaired | 0.0722564 | 0.215392 | 0.258066 |

## V1-collapse fit (reported, no verdict)

- S0: slope 0.7722 from 4 available points.
- S05: slope 0.917266 from 4 available points.
- S1: slope 1.06604 from 4 available points.

## Exploratory n=4 results (no verdict)

| family | d | M | skill | generator calls | Gaussian coordinates | GPU seconds |
|---|---:|---:|---:|---:|---:|---:|
| S0 | 400 | 6 | 0.0184328 | 2.982e+06 | 1.90627e+10 | 45.0919 |
| S0 | 400 | 10 | 0.00740477 | 1.25532e+07 | 1.39925e+11 | 332.102 |
| S05 | 400 | 6 | 0.0214743 | 2.982e+06 | 1.90627e+10 | 45.0813 |
| S05 | 400 | 10 | 0.00901484 | 1.25532e+07 | 1.39925e+11 | 332.138 |
| S1 | 400 | 6 | 0.057643 | 2.982e+06 | 1.90627e+10 | 45.0852 |
| S1 | 400 | 10 | 0.0246092 | 1.25532e+07 | 1.39925e+11 | 332.121 |

## Equal-cost frontiers (reported, no verdict)

Attaining cells for both generator-call and Gaussian-coordinate lower envelopes are recorded in `results/nonridge/frontier.csv` and annotated in the two frontier figures.

## Failures and deviations

- Mechanical criterion failures: none.
- No post-hoc thresholds were introduced; non-finite raw means are handled only by the pre-registered N-1 rule.
- The registered Stage-A grid contains `path` and `path_base` at M=16 but not M=24, although the reported V1-collapse sentence names both M=16 and M=24. No unregistered M=24 rows were added; the fit uses the available registered path/baseline points and this mismatch is retained explicitly.
- Stage-A >30-minute curtailments: 0 cells; details are in `results/nonridge/curtailed_cells.json` if nonzero.

## Provenance

- Frozen preregistration commit: `43cda1b3b5bbbb1dff8d6ddb7dbe47e76cc38f88`.
- Analysis code commit: `9b5c86f9d1bb7b83e2421c181758f4130e8f1a5f`.
- Solver/run commits represented in rows: `9b5c86f9d1bb7b83e2421c181758f4130e8f1a5f`, `ad4121678bbfd1929ce7ca07d43bd558f22d3c77`, `dbcb22b9a598d4f2aa901bde83fd24b2c1f39aac`.
- GPU: NVIDIA GeForce RTX 5070 Ti; driver 610.47; CUDA 12.8; PyTorch 2.11.0+cu128; DCT: Makhoul FFT type-II DCT/IDCT in PyTorch.
- Completed rows: 900; Stage A 440/440; S1 220/220; extension 180/180; n=4 60/60.
- Non-finite state values: 0; non-finite generator values: 0.
- GPU wall time represented by rows: 81193.5 seconds.
- Validation: V-1=PASS, V-2=PASS, V-3=PASS, V-4=PASS, V-5=PASS, V-6=PASS.
- Arithmetic is float64; verdicts use the frozen 960-point test split and 10 repetitions unless the pre-registered >30-minute rule is triggered.

## Figures

- `results/nonridge/figures/skill_vs_M.png`
- `results/nonridge/figures/v1_collapse.png`
- `results/nonridge/figures/skill_vs_dimension.png`
- `results/nonridge/figures/frontier_generator_calls.png`
- `results/nonridge/figures/frontier_gaussian_coordinates.png`
