# GPU MLP port: pre-registered report

## Outcome table

| Criterion | Outcome |
|---|---|
| G-0 | RECORDED |
| G-1 | PASS |
| G-2 | PASS |
| G-3 | PASS |
| G-4 | REPORTED |
| G-5 | REPORTED |

## Pre-registered criteria (verbatim)

### G-0

> **G-0 (environment):** record GPU model, driver, CUDA, framework versions, FP64 throughput class. If no GPU is available, implement and validate on the framework's CPU backend and mark G-4 and G-5 "not evaluated (no GPU)".

**Outcome: RECORDED.**

### G-1

> **G-1 (implementation):** PyTorch (or JAX), **float64 only**, breadth-first batched like the numpy code, supporting `raw`, `path`, `double`, `double_path`, `path_base`, `double_path_base` for P1 and MR; points processed in chunks to fit memory. Re-implement the equations; verify terminal, exact state and generator against the numpy equations to relative 1e-12 on random inputs. Seeding: a device generator seeded from `SeedSequence([20261201, d, n, M, rep, chunk_index]).generate_state(2)`, identical for all methods of a cell, so GPU methods are paired with each other (they are not paired with CPU rows; say so).

**Outcome: PASS.**

All 12 equation checks passed; maximum relative error 4.81387e-16 against threshold 1e-12.

### G-2

> **G-2 (exactness, verdict):** with a draw-tape interface, record every primary and auxiliary draw of the numpy solver for one chunk and replay it in the GPU solver. For P1 and MR at d=100, cells (2,8), (3,4), (3,6), all six methods: every output (u and z) agrees with numpy to relative 1e-10. PASS only if all agree.

**Outcome: PASS.**

All 36/36 checks were evaluated; maximum relative error 5.28786e-15 against threshold 1e-10.

### G-3

> **G-3 (statistical agreement, verdict):** with native GPU randomness, P1 and MR, d in {100, 400}, cells (2,32) and (3,10), methods `path`, `double_path`, `path_base`, `double_path_base`: 30 GPU repetitions versus 30 new CPU repetitions (stored in `results/gpu_mlp/`). PASS if |difference of mean skill| <= 3 x Welch standard error in every cell.

**Outcome: PASS.**

| PDE | d | cell | method | CPU mean | GPU mean | |difference| / Welch SE | pass |
|---|---:|---|---|---:|---:|---:|---|
| MR | 100 | (2,32) | double_path | 0.159632 | 0.16071 | 1.46912 | yes |
| MR | 100 | (2,32) | double_path_base | 0.159632 | 0.16071 | 1.46912 | yes |
| MR | 100 | (2,32) | path | 0.156412 | 0.157516 | 1.51157 | yes |
| MR | 100 | (2,32) | path_base | 0.156412 | 0.157516 | 1.51157 | yes |
| MR | 100 | (3,10) | double_path | 0.282561 | 0.276898 | 1.19387 | yes |
| MR | 100 | (3,10) | double_path_base | 0.0709358 | 0.0715917 | 0.964671 | yes |
| MR | 100 | (3,10) | path | 7.08151 | 7.02444 | 1.46918 | yes |
| MR | 100 | (3,10) | path_base | 0.681741 | 0.676284 | 1.12061 | yes |
| MR | 400 | (2,32) | double_path | 0.125207 | 0.125529 | 0.434283 | yes |
| MR | 400 | (2,32) | double_path_base | 0.125207 | 0.125529 | 0.434283 | yes |
| MR | 400 | (2,32) | path | 0.123201 | 0.123475 | 0.364343 | yes |
| MR | 400 | (2,32) | path_base | 0.123201 | 0.123475 | 0.364343 | yes |
| MR | 400 | (3,10) | double_path | 0.375602 | 0.363788 | 1.36884 | yes |
| MR | 400 | (3,10) | double_path_base | 0.0589902 | 0.0589346 | 0.0798141 | yes |
| MR | 400 | (3,10) | path | 19.8266 | 19.6302 | 1.85853 | yes |
| MR | 400 | (3,10) | path_base | 1.60296 | 1.59867 | 0.356433 | yes |
| P1 | 100 | (2,32) | double_path | 0.133849 | 0.13322 | 0.881232 | yes |
| P1 | 100 | (2,32) | double_path_base | 0.133849 | 0.13322 | 0.881232 | yes |
| P1 | 100 | (2,32) | path | 0.133559 | 0.13292 | 0.89864 | yes |
| P1 | 100 | (2,32) | path_base | 0.133559 | 0.13292 | 0.89864 | yes |
| P1 | 100 | (3,10) | double_path | 0.178028 | 0.175225 | 0.844155 | yes |
| P1 | 100 | (3,10) | double_path_base | 0.0502725 | 0.0504343 | 0.385517 | yes |
| P1 | 100 | (3,10) | path | 4.30093 | 4.26012 | 1.80319 | yes |
| P1 | 100 | (3,10) | path_base | 0.210914 | 0.21043 | 0.295239 | yes |
| P1 | 400 | (2,32) | double_path | 0.135962 | 0.135433 | 0.726179 | yes |
| P1 | 400 | (2,32) | double_path_base | 0.135962 | 0.135433 | 0.726179 | yes |
| P1 | 400 | (2,32) | path | 0.135624 | 0.135082 | 0.738889 | yes |
| P1 | 400 | (2,32) | path_base | 0.135624 | 0.135082 | 0.738889 | yes |
| P1 | 400 | (3,10) | double_path | 0.324057 | 0.313002 | 2.01158 | yes |
| P1 | 400 | (3,10) | double_path_base | 0.053089 | 0.0523493 | 1.40005 | yes |
| P1 | 400 | (3,10) | path | 16.7471 | 16.7714 | 0.321877 | yes |
| P1 | 400 | (3,10) | path_base | 0.807329 | 0.806016 | 0.300066 | yes |

### G-4

> **G-4 (speed, reported):** wall time per repetition (1,200 points) for `double_path_base`, P1, d=400, cells (3,16), (3,32), (4,6), (4,10), (4,16), on GPU and on CPU (numpy, same machine, report core count); peak GPU memory. CPU times that would exceed 2 hours may be extrapolated from one chunk; label them as extrapolated.

**Outcome: REPORTED.**

| cell | GPU seconds | CPU seconds | CPU extrapolated | speed-up | peak GPU GiB |
|---|---:|---:|---|---:|---:|
| (3,16) | 6.11839 | 267.441 | no | 43.7111 | 0.675204 |
| (3,32) | 51.8174 | 2182.25 | no | 42.1143 | 0.426577 |
| (4,6) | 6.16449 | 296.515 | no | 48.1005 | 0.832938 |
| (4,10) | 44.5209 | 2292.59 | no | 51.4946 | 0.741776 |
| (4,16) | 261.397 | 14062.3 | yes | 53.7967 | 0.900001 |

### G-5

> **G-5 (exploratory, only if G-2 and G-3 pass):** on GPU, n=4 with M in {6, 10, 16}, methods `path`, `double_path`, `double_path_base`; P1 and MR; d in {100, 400}; 10 repetitions. Report skill and all three costs next to the best n=2 and n=3 results from Study F. No verdict.

**Outcome: REPORTED.**

| PDE | d | M | method | skill | generator calls | Gaussian x d | seconds |
|---|---:|---:|---|---:|---:|---:|---:|
| MR | 100 | 6 | double_path | 7.19024 | 2.6928e+06 | 4.76568e+09 | 3.71723 |
| MR | 100 | 6 | double_path_base | 0.217151 | 2.982e+06 | 4.76568e+09 | 3.71397 |
| MR | 100 | 6 | path | 15557.4 | 1.2672e+06 | 1.4652e+09 | 1.14495 |
| MR | 100 | 10 | double_path | 1.30964 | 1.1784e+07 | 3.49812e+10 | 23.7837 |
| MR | 100 | 10 | double_path_base | 0.0812708 | 1.25532e+07 | 3.49812e+10 | 23.8574 |
| MR | 100 | 10 | path | 3193.05 | 5.424e+06 | 1.05972e+10 | 7.14743 |
| MR | 100 | 16 | double_path | 0.328723 | 4.67328e+07 | 2.22814e+11 | 152.586 |
| MR | 100 | 16 | double_path_base | 0.0446321 | 4.8654e+07 | 2.22814e+11 | 152.829 |
| MR | 100 | 16 | path | 746.205 | 2.12352e+07 | 6.69408e+10 | 45.8791 |
| MR | 400 | 6 | double_path | 32.322 | 2.6928e+06 | 1.90627e+10 | 11.0557 |
| MR | 400 | 6 | double_path_base | 0.347557 | 2.982e+06 | 1.90627e+10 | 11.0818 |
| MR | 400 | 6 | path | 520973 | 1.2672e+06 | 5.8608e+09 | 3.31399 |
| MR | 400 | 10 | double_path | 5.49816 | 1.1784e+07 | 1.39925e+11 | 82.7643 |
| MR | 400 | 10 | double_path_base | 0.0938191 | 1.25532e+07 | 1.39925e+11 | 82.9485 |
| MR | 400 | 10 | path | 106839 | 5.424e+06 | 4.23888e+10 | 24.9733 |
| MR | 400 | 16 | double_path | 1.03955 | 4.67328e+07 | 8.91256e+11 | 507.024 |
| MR | 400 | 16 | double_path_base | 0.04383 | 4.8654e+07 | 8.91256e+11 | 511.831 |
| MR | 400 | 16 | path | 25022.8 | 2.12352e+07 | 2.67763e+11 | 157.768 |
| P1 | 100 | 6 | double_path | 2.30493 | 2.6928e+06 | 4.76568e+09 | 2.13476 |
| P1 | 100 | 6 | double_path_base | 0.0837428 | 2.982e+06 | 4.76568e+09 | 2.13007 |
| P1 | 100 | 6 | path | 4954.99 | 1.2672e+06 | 1.4652e+09 | 0.653302 |
| P1 | 100 | 10 | double_path | 0.441907 | 1.1784e+07 | 3.49812e+10 | 11.5368 |
| P1 | 100 | 10 | double_path_base | 0.0406141 | 1.25532e+07 | 3.49812e+10 | 11.5755 |
| P1 | 100 | 10 | path | 1023.11 | 5.424e+06 | 1.05972e+10 | 3.508 |
| P1 | 100 | 16 | double_path | 0.139687 | 4.67328e+07 | 2.22814e+11 | 74.9127 |
| P1 | 100 | 16 | double_path_base | 0.023356 | 4.8654e+07 | 2.22814e+11 | 74.9457 |
| P1 | 100 | 16 | path | 263.818 | 2.12352e+07 | 6.69408e+10 | 22.6596 |
| P1 | 400 | 6 | double_path | 15.6655 | 2.6928e+06 | 1.90627e+10 | 5.89278 |
| P1 | 400 | 6 | double_path_base | 0.145 | 2.982e+06 | 1.90627e+10 | 5.87676 |
| P1 | 400 | 6 | path | 256887 | 1.2672e+06 | 5.8608e+09 | 1.79063 |
| P1 | 400 | 10 | double_path | 2.93801 | 1.1784e+07 | 1.39925e+11 | 43.6164 |
| P1 | 400 | 10 | double_path_base | 0.0566342 | 1.25532e+07 | 1.39925e+11 | 43.7292 |
| P1 | 400 | 10 | path | 55257.2 | 5.424e+06 | 4.23888e+10 | 13.1554 |
| P1 | 400 | 16 | double_path | 0.545541 | 4.67328e+07 | 8.91256e+11 | 256.759 |
| P1 | 400 | 16 | double_path_base | 0.0301918 | 4.8654e+07 | 8.91256e+11 | 258.23 |
| P1 | 400 | 16 | path | 13020.3 | 2.12352e+07 | 2.67763e+11 | 77.6754 |

## Best Study-F n=2/n=3 references

| PDE | d | n | cell/method | skill | generator calls | Gaussian x d | seconds |
|---|---:|---:|---|---:|---:|---:|---:|
| MR | 100 | 2 | M=512, path | 0.145218 | 614400 | 6.2976e+10 | 1930.7 |
| MR | 100 | 3 | M=32, double_path_base | 0.0354132 | 3.8412e+06 | 3.60077e+10 | 3457.46 |
| MR | 400 | 2 | M=512, path | 0.112974 | 614400 | 2.51904e+11 | 6053.43 |
| MR | 400 | 3 | M=32, double_path_base | 0.0291549 | 3.8412e+06 | 1.44031e+11 | 11678 |
| P1 | 100 | 2 | M=256, path | 0.122785 | 307200 | 1.57594e+10 | 992.36 |
| P1 | 100 | 3 | M=32, double_path_base | 0.0260918 | 3.8412e+06 | 3.60077e+10 | 1955.46 |
| P1 | 400 | 2 | M=512, path | 0.124583 | 614400 | 2.51904e+11 | 5402.28 |
| P1 | 400 | 3 | M=32, double_path_base | 0.0268048 | 3.8412e+06 | 1.44031e+11 | 8711.37 |

## Failures and deviations

- Mechanical failures: none.
- GPU and CPU G-3 rows use independent native random-number implementations and are not paired; agreement is assessed only by the preregistered Welch rule.
- G-3 repetitions are numbered 100–129 so all CPU draws are new relative to Study F.
- No post-hoc thresholds were introduced.

## Provenance

- Frozen preregistration commit: `aa5b9b1fd0d0a5ee0310b69895080c00c476cd35`.
- Analysis code commit: `cd5b9f22e94f101c306d2d3aebc8fea9c4738784`.
- Solver/run code commit: `89aeb537c23e52e6f42b28089a9b64eec3b59e52`.
- GPU: NVIDIA GeForce RTX 5070 Ti; driver 610.47; CUDA 12.8; PyTorch 2.11.0+cu128.
- FP64 throughput class: consumer GPU; reduced-rate FP64.
- CPU timing host: 16 logical cores; NumPy timings use the same machine as the GPU timings.
- G-3 rows: CPU 960/960, GPU 960/960.
- G-4 timing rows: 5/5; G-5 rows: 360/360 across 36/36 summary cells.
- Pairing: GPU methods within a cell share the preregistered device seed; GPU and CPU G-3 rows are intentionally unpaired. G-2 replayed every recorded primary and auxiliary NumPy draw.
- Arithmetic: float64 only; point sets are the frozen 1,200 points and verdicts use the frozen 960-point test split.

## Figures

- `results/gpu_mlp/figures/speedup.png`
- `results/gpu_mlp/figures/n4_skill.png`
