# High-effective-dimension, symmetry-free benchmark

## Outcome table

| Criterion | Verdict | Wins / evaluated | Required |
|---|---|---:|---:|
| H1 | **NOT EVALUATED (COMPUTE)** | 0 / 3 | 4 |
| H2 | **NOT EVALUATED (COMPUTE)** | 0 / 0 | 2 |
| H3 | **NOT EVALUATED (COMPUTE)** | 3 / 3 | 4 |
| H4 | **NOT EVALUATED (COMPUTE)** | 2 / 2 | 3 |
| H5 | **PASS** | 2 / 2 | 2 |
| H6 | **NOT EVALUATED (COMPUTE)** | 0 / 1 | 2 |

## Pre-registered criteria (verbatim)

### H1

> RMSE(TA) <= 0.5 x RMSE(TA_a0) in all 4 (PDE, d) cells.

**Verdict: NOT EVALUATED (COMPUTE).**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| S05 | 100 | 0.00183597 | 0.00258117 | 0.711292 | no |
| B6 | 100 | 0.00133404 | 0.00164973 | 0.808644 | no |
| S05 | 400 | 0.00196131 | 0.00234057 | 0.837961 | no |
| B6 | 400 | NA | NA | NA | not evaluated |

### H2

> RMSE(TA) <= 0.8 x RMSE(a1), both under the same post-freeze common-random-number diagnostic, in at least 2 of the 4 cells.

**Verdict: NOT EVALUATED (COMPUTE).**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| S05 | 100 | 0.00183597 | NA | NA | not evaluated |
| B6 | 100 | 0.00133404 | NA | NA | not evaluated |
| S05 | 400 | 0.00196131 | NA | NA | not evaluated |
| B6 | 400 | NA | NA | NA | not evaluated |

### H3

> RMSE(TA) <= 0.5 x the minimum RMSE over the entire frozen FH_raw grid, in all 4 cells.

**Verdict: NOT EVALUATED (COMPUTE).**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| S05 | 100 | 0.00183597 | 0.277487 | 0.00661642 | yes |
| B6 | 100 | 0.00133404 | 0.0865748 | 0.0154091 | yes |
| S05 | 400 | 0.00196131 | 2.21073 | 0.000887174 | yes |
| B6 | 400 | NA | NA | NA | not evaluated |

### H4

> RMSE(TA) <= the minimum RMSE over all frozen FH_raw, FH_path and FH_rep cells (any cost), in at least 3 of the 4 cells.

**Verdict: NOT EVALUATED (COMPUTE).**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| S05 | 100 | 0.00183597 | 0.00705174 | 0.260357 | yes |
| B6 | 100 | 0.00133404 | 0.00297944 | 0.44775 | yes |
| S05 | 400 | 0.00196131 | NA | NA | not evaluated |
| B6 | 400 | NA | NA | NA | not evaluated |

### H5

> In at least 2 of the 4 cells, the best Full-History family does not reach 2 x RMSE(TA) within its frozen grid, or needs more GPU seconds than TA's total to do so.

**Verdict: PASS.**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| S05 | 100 | 0.00183597 | NA | NA | yes |
| B6 | 100 | 0.00133404 | NA | NA | yes |
| S05 | 400 | 0.00196131 | NA | NA | not evaluated |
| B6 | 400 | NA | NA | NA | not evaluated |

### H6

> For B6, r_glob >= d/4 at the frozen calibration, and H1 holds on B6 at both d.

**Verdict: NOT EVALUATED (COMPUTE).**

| PDE | d | TA RMSE / r_glob | comparison / threshold | ratio / H1 | met |
|---|---:|---:|---:|---:|---|
| B6 | 100 | 28.9865 | 25 | 0 | no |
| B6 | 400 | 77.1437 | 100 | NA | not evaluated |

## B6 calibration (reference-only, frozen before solver runs)

| d | beta | gamma | ESS med | z RMS | f CV | heat skill | r glob | outside top 10 | r loc | qualifies | selected |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 100 | 1 | 0 | 163.347 | 0.176325 | 0.271794 | 0.0297478 | 31.1258 | 0.524869 | 45.3792 | no |  |
| 100 | 1 | 1 | 163.306 | 1.42966 | 0.0196103 | 0.251741 | 31.1956 | 0.524584 | 45.3512 | no |  |
| 100 | 1.5 | 0 | 119.533 | 0.423287 | 0.617221 | 0.102136 | 28.6678 | 0.498974 | 30.7498 | no |  |
| 100 | 1.5 | 1 | 119.411 | 1.48893 | 0.071921 | 0.262954 | 28.9865 | 0.500913 | 30.5981 | yes | yes |
| 100 | 2 | 0 | 72.9563 | 0.910479 | 0.930358 | 0.243728 | 26.119 | 0.467348 | 18.7405 | yes |  |
| 100 | 2 | 1 | 72.6708 | 1.70465 | 0.274645 | 0.327916 | 26.6028 | 0.472447 | 18.6527 | yes |  |
| 100 | 3 | 0 | 20.7031 | 2.59038 | 0.880357 | 0.568573 | 24.3863 | 0.443397 | 7.63918 | yes |  |
| 100 | 3 | 1 | 20.7725 | 2.98126 | 0.667036 | 0.595925 | 24.6522 | 0.44627 | 7.67824 | yes |  |
| 100 | 4 | 0 | 7.78974 | 4.65262 | 0.692437 | 0.82716 | 23.9168 | 0.437991 | 4.28268 | no |  |
| 100 | 4 | 1 | 7.80671 | 4.89105 | 0.632657 | 0.843421 | 24.0379 | 0.438937 | 4.2774 | no |  |
| 400 | 1 | 0 | 650.715 | 0.0897055 | 0.168025 | 0.00516807 | 112.115 | 0.821044 | 176.582 | no |  |
| 400 | 1 | 1 | 650.811 | 1.41459 | 0.00459232 | 0.235134 | 112.188 | 0.821084 | 176.701 | no |  |
| 400 | 1.5 | 0 | 464.863 | 0.225977 | 0.572034 | 0.0310837 | 98.0497 | 0.792179 | 105.014 | no |  |
| 400 | 1.5 | 1 | 465.205 | 1.42753 | 0.0190336 | 0.225895 | 98.372 | 0.79298 | 105.241 | no |  |
| 400 | 2 | 0 | 263.721 | 0.553283 | 1.2238 | 0.119779 | 76.4604 | 0.739918 | 51.7587 | no |  |
| 400 | 2 | 1 | 263.652 | 1.51068 | 0.163176 | 0.241803 | 77.1437 | 0.742048 | 51.776 | yes | yes |
| 400 | 3 | 0 | 53.9525 | 2.01781 | 1.14344 | 0.412185 | 73.2173 | 0.739778 | 14.1266 | yes |  |
| 400 | 3 | 1 | 54.0834 | 2.45171 | 0.767671 | 0.446771 | 73.9951 | 0.741434 | 14.2162 | yes |  |
| 400 | 4 | 0 | 14.5699 | 4.04307 | 0.839659 | 0.674798 | 76.8241 | 0.752987 | 6.35026 | no |  |
| 400 | 4 | 1 | 14.3704 | 4.27081 | 0.749533 | 0.689511 | 77.3591 | 0.753465 | 6.31839 | no |  |

S05 global gradient-covariance effective ranks (same definition as `r_glob`): d=100: 14.5467, d=400: 21.5715.

## Validation

- **autodiff-pde-residual: PASS** (2 checks, 0 failures).
- **calibration-freeze: PASS** (2 checks, 0 failures).
- **draw-tape-replay: PASS** (4 checks, 0 failures).
- **gpu-adapter: PASS** (2 checks, 0 failures).
- **hopf-cole-monte-carlo: PASS** (2 checks, 0 failures).
- **s05-bank-reuse: PASS** (2 checks, 0 failures).
- **terminal-identity: PASS** (2 checks, 0 failures).
- **z-finite-difference: PASS** (2 checks, 0 failures).

## Actor selection and iteration diagnostics

| PDE | d | seed | actor | resolved | steps | RMSE | MC contribution | GPU s |
|---|---:|---:|---|---|---:|---:|---:|---:|
| B6 | 100 | 0 | a0 | a0 | 64 | 0.00164973 | 0.00159905 | 27.2393 |
| B6 | 100 | 0 | selected | a4 | 64 | 0.00133404 | 0.00134975 | 27.1994 |
| S05 | 100 | 0 | a0 | a0 | 64 | 0.00258117 | 0.00200584 | 45.2258 |
| S05 | 100 | 0 | selected | a4 | 64 | 0.00183597 | 0.00183033 | 45.6979 |
| S05 | 400 | 0 | a0 | a0 | 64 | 0.00234057 | 0.00196193 | 44.902 |
| S05 | 400 | 0 | selected | a4 | 64 | 0.00196131 | 0.00187815 | 44.3855 |

The headline uses the arithmetic mean of the four independently replicated 256-path RMSE values; the repetitions are not pooled into a 1,024-path estimator. Truth was loaded only after reference-free checkpoint selection.

## Full-History grid

| PDE | d | method | n | M | reps | reused | RMSE | skill | GPU s | nonfinite |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B6 | 100 | FH_rep | 3 | 16 | 3 | 0 | 0.0106588 | 0.0163939 | 9.80927 | 0 |
| B6 | 100 | FH_rep | 3 | 24 | 3 | 0 | 0.00669937 | 0.0103041 | 32.9594 | 0 |
| B6 | 100 | FH_rep | 3 | 32 | 3 | 0 | 0.00470947 | 0.00724346 | 78.1395 | 0 |
| B6 | 100 | FH_rep | 4 | 6 | 3 | 0 | 0.0171328 | 0.0263513 | 11.1174 | 0 |
| B6 | 100 | FH_rep | 4 | 10 | 3 | 0 | 0.00669291 | 0.0102941 | 76.2747 | 0 |
| B6 | 100 | FH_rep | 4 | 16 | 3 | 0 | 0.00297944 | 0.00458256 | 482.29 | 0 |
| B6 | 100 | FH_path | 2 | 16 | 3 | 0 | 0.0426869 | 0.0656551 | 0.180162 | 0 |
| B6 | 100 | FH_path | 2 | 64 | 3 | 0 | 0.0164603 | 0.025317 | 2.09848 | 0 |
| B6 | 100 | FH_path | 2 | 256 | 3 | 0 | 0.00847328 | 0.0130324 | 33.604 | 0 |
| B6 | 100 | FH_path | 3 | 8 | 3 | 0 | 0.518291 | 0.797164 | 0.67549 | 0 |
| B6 | 100 | FH_path | 3 | 16 | 3 | 0 | 0.252203 | 0.387904 | 4.34311 | 0 |
| B6 | 100 | FH_path | 3 | 24 | 3 | 0 | 0.16726 | 0.257256 | 14.5248 | 0 |
| B6 | 100 | FH_path | 4 | 6 | 3 | 0 | 19.8313 | 30.5018 | 3.21938 | 0 |
| B6 | 100 | FH_path | 4 | 10 | 3 | 0 | 3.84832 | 5.91895 | 22.3193 | 0 |
| B6 | 100 | FH_raw | 2 | 16 | 3 | 0 | 1.45919 | 2.24432 | 0.194804 | 0 |
| B6 | 100 | FH_raw | 2 | 64 | 3 | 0 | 0.357972 | 0.550583 | 1.21793 | 0 |
| B6 | 100 | FH_raw | 2 | 256 | 3 | 0 | 0.0865748 | 0.133157 | 19.3622 | 0 |
| B6 | 100 | FH_raw | 3 | 8 | 3 | 0 | 221.05 | 339.989 | 0.436529 | 0 |
| B6 | 100 | FH_raw | 3 | 16 | 3 | 0 | 27.135 | 41.7353 | 2.53926 | 0 |
| B6 | 100 | FH_raw | 3 | 24 | 3 | 0 | 8.49411 | 13.0645 | 8.47555 | 0 |
| B6 | 100 | FH_raw | 4 | 6 | 3 | 0 | 1.55323e+07 | 2.38896e+07 | 1.97461 | 0 |
| B6 | 100 | FH_raw | 4 | 10 | 3 | 0 | 302251 | 464881 | 13.1797 | 0 |
| S05 | 100 | FH_rep | 3 | 16 | 3 | 3 | 0.0237027 | 0.0193267 | 10.1359 | 0 |
| S05 | 100 | FH_rep | 3 | 24 | 3 | 3 | 0.0153439 | 0.0125111 | 34.064 | 0 |
| S05 | 100 | FH_rep | 3 | 32 | 3 | 3 | 0.010683 | 0.00871071 | 80.5085 | 0 |
| S05 | 100 | FH_rep | 4 | 6 | 3 | 0 | 0.0343818 | 0.0280343 | 10.9356 | 0 |
| S05 | 100 | FH_rep | 4 | 10 | 3 | 0 | 0.0147717 | 0.0120446 | 77.2876 | 0 |
| S05 | 100 | FH_rep | 4 | 16 | 3 | 0 | 0.00705174 | 0.00574986 | 495.664 | 0 |
| S05 | 100 | FH_path | 2 | 16 | 3 | 3 | 0.0979944 | 0.0799028 | 0.176841 | 0 |
| S05 | 100 | FH_path | 2 | 64 | 3 | 3 | 0.0491456 | 0.0400724 | 2.21223 | 0 |
| S05 | 100 | FH_path | 2 | 256 | 3 | 3 | 0.0346006 | 0.0282127 | 35.1947 | 0 |
| S05 | 100 | FH_path | 3 | 8 | 3 | 0 | 3.99759 | 3.25956 | 0.753427 | 0 |
| S05 | 100 | FH_path | 3 | 16 | 3 | 3 | 1.99156 | 1.62388 | 4.5223 | 0 |
| S05 | 100 | FH_path | 3 | 24 | 3 | 0 | 1.31745 | 1.07423 | 15.259 | 0 |
| S05 | 100 | FH_path | 4 | 6 | 3 | 0 | 1065.18 | 868.526 | 3.39603 | 0 |
| S05 | 100 | FH_path | 4 | 10 | 3 | 0 | 216.848 | 176.813 | 22.9944 | 0 |
| S05 | 100 | FH_raw | 2 | 16 | 3 | 0 | 4.91355 | 4.00641 | 0.198554 | 0 |
| S05 | 100 | FH_raw | 2 | 64 | 3 | 3 | 1.19515 | 0.974503 | 1.15075 | 0 |
| S05 | 100 | FH_raw | 2 | 256 | 3 | 0 | 0.277487 | 0.226257 | 18.3219 | 0 |
| S05 | 100 | FH_raw | 3 | 8 | 3 | 0 | 2377.48 | 1938.55 | 0.515481 | 0 |
| S05 | 100 | FH_raw | 3 | 16 | 3 | 3 | 291.136 | 237.387 | 2.40227 | 0 |
| S05 | 100 | FH_raw | 3 | 24 | 3 | 0 | 89.883 | 73.2889 | 8.01356 | 0 |
| S05 | 100 | FH_raw | 4 | 6 | 3 | 0 | 1.85639e+09 | 1.51367e+09 | 1.91163 | 0 |
| S05 | 100 | FH_raw | 4 | 10 | 3 | 0 | 3.56104e+07 | 2.9036e+07 | 12.2673 | 0 |
| S05 | 400 | FH_rep | 3 | 16 | 3 | 3 | 0.022487 | 0.015023 | 43.8146 | 0 |
| S05 | 400 | FH_rep | 3 | 24 | 3 | 3 | 0.014079 | 0.0094058 | 147.676 | 0 |
| S05 | 400 | FH_rep | 3 | 32 | 3 | 3 | 0.0102832 | 0.00686995 | 352.12 | 0 |
| S05 | 400 | FH_rep | 4 | 6 | 3 | 3 | 0.03233 | 0.0215988 | 45.0893 | 0 |
| S05 | 400 | FH_rep | 4 | 10 | 3 | 3 | 0.0131393 | 0.00877801 | 332.197 | 0 |
| S05 | 400 | FH_rep | 4 | 16 | 1 | 0 | 0.00629363 | 0.00420461 | 2109.83 | 0 |
| S05 | 400 | FH_path | 2 | 16 | 3 | 3 | 0.0927653 | 0.0619742 | 0.718062 | 0 |
| S05 | 400 | FH_path | 2 | 64 | 3 | 3 | 0.0421755 | 0.0281764 | 9.51981 | 0 |
| S05 | 400 | FH_path | 2 | 256 | 3 | 3 | 0.0272684 | 0.0182173 | 150.715 | 0 |
| S05 | 400 | FH_path | 3 | 8 | 3 | 0 | 15.1065 | 10.0923 | 2.69361 | 0 |
| S05 | 400 | FH_path | 3 | 16 | 3 | 3 | 7.44757 | 4.97554 | 19.4894 | 0 |
| S05 | 400 | FH_path | 3 | 24 | 3 | 0 | 4.97737 | 3.32526 | 65.7476 | 0 |
| S05 | 400 | FH_path | 4 | 6 | 3 | 0 | 53800.6 | 35942.8 | 13.5599 | 0 |
| S05 | 400 | FH_path | 4 | 10 | 3 | 0 | 11387.8 | 7607.9 | 98.8178 | 0 |
| S05 | 400 | FH_raw | 2 | 16 | 3 | 0 | 35.8169 | 23.9284 | 0.553182 | 0 |
| S05 | 400 | FH_raw | 2 | 64 | 3 | 3 | 8.89762 | 5.94427 | 4.86999 | 0 |
| S05 | 400 | FH_raw | 2 | 256 | 3 | 0 | 2.21073 | 1.47694 | 76.2904 | 0 |
| S05 | 400 | FH_raw | 3 | 8 | 3 | 0 | 454065 | 303349 | 1.51136 | 0 |
| S05 | 400 | FH_raw | 3 | 16 | 3 | 3 | 53531.9 | 35763.3 | 10.0745 | 0 |
| S05 | 400 | FH_raw | 3 | 24 | 3 | 0 | 15632.4 | 10443.6 | 34.2148 | 0 |
| S05 | 400 | FH_raw | 4 | 6 | 3 | 0 | 2.60158e+14 | 1.73805e+14 | 7.3011 | 0 |
| S05 | 400 | FH_raw | 4 | 10 | 3 | 0 | 5.67604e+12 | 3.79202e+12 | 51.7859 | 0 |

## Failures and deviations

- Mechanical criterion failures: none among evaluated criteria.
- Frozen Full-History rows not yet available: 68; affected criteria are marked `NOT EVALUATED (COMPUTE)` rather than removed.
- Runtime hardware deviation: the requested GH200 was unavailable; all new timings are from an NVIDIA GeForce RTX 5070 Ti. Historical reused S05 rows were generated on the same GPU model.
- The two final-suite DPI unit tests that import the external `picard` package were not runnable because that optional package is absent; the 16 relevant highdim/nonridge tests passed.
- Reused S05 artifacts do not contain per-point predictions, so their MC/non-MC decomposition is unavailable; no proxy was invented.
- B6 d=400 selected r_glob is below d/4. The frozen selector was applied unchanged and the resulting H6 failure is retained.
- Compute-curtailed: S05, d=400, FH_rep, n=4, M=16, rep=1 after rep0 took 2124.8 GPU seconds.
- Compute-curtailed: S05, d=400, FH_rep, n=4, M=16, rep=2 after rep0 took 2124.8 GPU seconds.

## Provenance

- Preregistration commit: `d97ea678c0daddf1b74117289237591593e9aeed`.
- Calibration freeze commit: `a3d3887817b84d3f6b2aaa57ed787dab3825e2dc`.
- Analysis commit/worktree head: `dba32fdf363aa8bc12c9eb7ff6f6eeb02c17e223`.
- Imported nonridge commit: `d384f4e8`; byte-identity audit: PASS (11/11 blobs identical).
- Frozen final-suite protocol SHA-256: `80e81e5388e71aabbd5f5553939b84c4e265b5871a889402882288037491ee79`.
- Frozen review-gate SHA-256: `ee59264b2e8348be33a0a19c3dc97941772df9c6edfadfa57c8688759f707f19`.
- GPU/software: NVIDIA GeForce RTX 5070 Ti; Python 3.12.11, NumPy 2.3.4, PyTorch 2.11.0+cu128, CUDA 12.8.
- Validation: PASS.

## Figures

- `results/highdim_nosym/figures/rmse_vs_gpu_seconds.png`
- `results/highdim_nosym/figures/actor_iterations.png`
- `results/highdim_nosym/figures/b6_calibration_heatmap.png`
