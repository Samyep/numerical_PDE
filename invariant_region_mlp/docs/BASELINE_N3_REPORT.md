# Correction-level baseline at n=3: pre-registered report

## Outcome table

| Criterion | Verdict |
|---|---|
| F-1 | PASS (12/12) |
| F-2 | PASS (12/12) |
| F-3 | PASS |
| F-4 | REPORTED (no verdict by design) |

## Pre-registered criteria (verbatim)

### F-1

> **F-1 (baseline repairs `path` at n=3):** `path_base` mean skill <= 0.5 x `path` mean skill in at least 10 of the 12 (PDE, d, M) n=3 combinations.

**Verdict: PASS (12/12; required 10).**

| PDE | d | M | ratio | criterion met |
|---|---:|---:|---:|---|
| P1 | 100 | 16 | 0.0442866 | yes |
| P1 | 100 | 24 | 0.0428394 | yes |
| P1 | 100 | 32 | 0.0429435 | yes |
| P1 | 400 | 16 | 0.0438045 | yes |
| P1 | 400 | 24 | 0.041608 | yes |
| P1 | 400 | 32 | 0.0401543 | yes |
| MR | 100 | 16 | 0.0873141 | yes |
| MR | 100 | 24 | 0.0831219 | yes |
| MR | 100 | 32 | 0.0810961 | yes |
| MR | 400 | 16 | 0.0722564 | yes |
| MR | 400 | 24 | 0.0700688 | yes |
| MR | 400 | 32 | 0.0682733 | yes |

### F-2

> **F-2 (baseline complements double):** `double_path_base` mean skill <= `double_path` mean skill in at least 10 of the 12 combinations.

**Verdict: PASS (12/12; required 10).**

| PDE | d | M | ratio | criterion met |
|---|---:|---:|---:|---|
| P1 | 100 | 16 | 0.376181 | yes |
| P1 | 100 | 24 | 0.506629 | yes |
| P1 | 100 | 32 | 0.598289 | yes |
| P1 | 400 | 16 | 0.221509 | yes |
| P1 | 400 | 24 | 0.307354 | yes |
| P1 | 400 | 32 | 0.42487 | yes |
| MR | 100 | 16 | 0.323269 | yes |
| MR | 100 | 24 | 0.446533 | yes |
| MR | 100 | 32 | 0.554971 | yes |
| MR | 400 | 16 | 0.215392 | yes |
| MR | 400 | 24 | 0.313421 | yes |
| MR | 400 | 32 | 0.392323 | yes |

### F-3

> **F-3 (n=3 below the n=2 floor at d=400):** for each PDE at d=400, the minimum over n=3 cells and methods {`double_path`, `double_path_base`} of mean skill is <= 0.9 x the minimum over all n=2 `path` cells (M = 8 ... 512, including reused rows). PASS if it holds for both PDEs. Report the same quantity at d=100 (no verdict).

**Verdict: PASS.**

| PDE | d | best n=2 path floor | best n=3 double family | ratio n3/n2 | verdict scope | met |
|---|---:|---|---|---:|---|---|
| P1 | 100 | M=256; 0.122785 | double_path_base M=32; 0.0260918 | 0.2125 | reported only | NA |
| P1 | 400 | M=512; 0.124583 | double_path_base M=32; 0.0268048 | 0.215155 | yes | yes |
| MR | 100 | M=512; 0.145218 | double_path_base M=32; 0.0354132 | 0.243863 | reported only | NA |
| MR | 400 | M=512; 0.112974 | double_path_base M=32; 0.0291549 | 0.258066 | yes | yes |

### F-4

> **F-4 (reported, no verdict):** frontiers (lower envelope of mean skill vs cost) of the n=2-only `path` family and of the combined family {`path`, `path_base`, `double_path`, `double_path_base`} over n in {2,3}, on each of the three cost axes, at 10 log-spaced cost levels where both exist; the attaining cell at each level; 1,000-draw bootstrap over reps.

**Verdict: REPORTED (reported, no pass/fail criterion).**

All 120 attaining-cell records and paired bootstrap intervals are in `results/baseline_n3/frontier_levels.csv`.

| axis | PDE | d | levels | combined better levels |
|---|---|---:|---:|---:|
| gaussian_coordinates | MR | 100 | 10 | 3 |
| gaussian_coordinates | MR | 400 | 10 | 3 |
| gaussian_coordinates | P1 | 100 | 10 | 4 |
| gaussian_coordinates | P1 | 400 | 10 | 3 |
| generator_calls | MR | 100 | 10 | 0 |
| generator_calls | MR | 400 | 10 | 0 |
| generator_calls | P1 | 100 | 10 | 0 |
| generator_calls | P1 | 400 | 10 | 0 |
| wall_time | MR | 100 | 10 | 2 |
| wall_time | MR | 400 | 10 | 2 |
| wall_time | P1 | 100 | 10 | 3 |
| wall_time | P1 | 400 | 10 | 2 |

## Predictions (reported, no verdict)

- `path_base_oracle <= path_base`: 0/12 cells; the oracle was instead 2.47549–5.55397x the baseline skill (lower is better).
- Root level-1 integrand variance baseline/path: median 0.087463, range 0.0625764–0.127348; reused M=16 path rows lack this new diagnostic and were not rerun.
- For path, M x skill / d ranges from 0.410679 to 0.688943; cell values are retained in `analysis_summary.json`.

## Failures and deviations

- Mechanical failures: none.
- No post-hoc thresholds were introduced.
- Existing rows were reused without rerunning. Consequently, the newly requested root-integrand variance is unavailable for reused M=16 path/double-path and n=2 rows; this does not enter F-1 through F-4.

## Provenance

- Frozen preregistration commit: `aa5b9b1fd0d0a5ee0310b69895080c00c476cd35`.
- Analysis code commit: `cd5b9f22e94f101c306d2d3aebc8fea9c4738784`.
- Solver/run code commit: `89aeb537c23e52e6f42b28089a9b64eec3b59e52`.
- Parent path-frontier commit: `a42e2c5cf3b3ad78ba702d42f3de6ab559150f0d`.
- Environment: Windows-11-10.0.26200-SP0; Python 3.12.11 | packaged by conda-forge | (main, Jun  4 2025, 14:29:09) [MSC v.1943 64 bit (AMD64)]; NumPy 2.3.4; 16 logical CPUs, 4 workers.
- Rows: 880 / 880; repetitions: 10 per cell.
- Reused rows: 280; new rows: 600.
- Primary-draw fingerprint audit: PASS; 360 paired chunks, 0 mismatches; records SHA-256 `9e1c655af0660c37d09a711ae103785e1e878931966f80c8f5626f95a25564e2`.
- Non-finite state values: 0; non-finite generator values: 0.
- Test points: frozen 1,200 points with the frozen 20% validation split; all criteria use the 960 test points.
- Arithmetic: float64; bootstrap: 1,000 paired repetition draws with seed 20261201.

## Figures

- `results/baseline_n3/figures/skill_n3.png`
- `results/baseline_n3/figures/frontier_generator_calls.png`
- `results/baseline_n3/figures/frontier_gaussian_coordinates.png`
- `results/baseline_n3/figures/frontier_wall_time.png`
