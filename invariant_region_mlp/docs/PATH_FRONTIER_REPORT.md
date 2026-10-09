# Pathwise-only equal-cost frontier: pre-registered report

## Outcome table

| Criterion | Verdict |
|---|---|
| E-1 | pathwise alone suffices at the frontier |
| E-2 | PASS |

## Pre-registered criteria (verbatim)

### E-1

> Let `best_double` be the lower envelope over `double` and `double_path`. "Double adds value at the frontier" if `best_double <= 0.9 x path` at >= 8 of the 10 cost levels, for at least 3 of the 4 (PDE, d) problems. Otherwise report "pathwise alone suffices at the frontier".

**Verdict: pathwise alone suffices at the frontier**

| Problem | qualifying levels | required | problem criterion |
|---|---:|---:|---|
| P1_d100 | 1 | 8 | no |
| P1_d400 | 0 | 8 | no |
| MR_d100 | 0 | 8 | no |
| MR_d400 | 0 | 8 | no |

### E-2

> the `path` frontier is <= 0.8 x the `raw` frontier at >= 8 of 10 cost levels, for each of the 4 problems.

**Verdict: PASS**

| Problem | qualifying levels | required | problem criterion |
|---|---:|---:|---|
| P1_d100 | 10 | 8 | yes |
| P1_d400 | 10 | 8 | yes |
| MR_d100 | 10 | 8 | yes |
| MR_d400 | 10 | 8 | yes |

## E-3: dimension robustness (reported, no verdict)

The median budget is the conventional numeric median of the ten frozen D-4 cost levels. The prediction column is reported mechanically but is not a criterion.

| PDE | method | median cost | d=100 skill (cell) | d=400 skill (cell) | ratio 400/100 | prediction | met |
|---|---|---:|---|---|---:|---|---|
| P1 | path | 114379 | 0.127271 ([2, 64]) | 0.128167 ([2, 64]) | 1.00704 | <= 1.3 | yes |
| P1 | best_double | 114379 | 0.127419 ([2, 64]) | 0.128314 ([2, 64]) | 1.00702 | <= 1.3 | yes |
| P1 | raw | 114379 | 0.531974 ([2, 64]) | 2.41553 ([2, 64]) | 4.5407 | >= 2 | yes |
| MR | path | 114379 | 0.15065 ([2, 64]) | 0.117099 ([2, 64]) | 0.777292 | <= 1.3 | yes |
| MR | best_double | 114379 | 0.152368 ([2, 64]) | 0.118223 ([2, 64]) | 0.775905 | <= 1.3 | yes |
| MR | raw | 114379 | 0.791221 ([2, 64]) | 2.70661 ([2, 64]) | 3.4208 | >= 2 | yes |

## Frontier-attaining cells (reported, no verdict)

Each entry is `method (n,M); mean skill`. Wall-time selections and bootstrap intervals are retained in `frontier_levels.csv`.

| PDE | d | level | calls | raw | path | double | double_path | best_double |
|---|---:|---:|---:|---|---|---|---|---|
| MR | 100 | 0 | 9600 | raw (2,8); 7.63145 | path (2,8); 0.201307 | double (2,8); 0.397213 | double_path (2,8); 0.210569 | double_path (2,8); 0.210569 |
| MR | 100 | 1 | 16515.5 | raw (2,8); 7.63145 | path (2,8); 0.201307 | double (2,8); 0.397213 | double_path (2,8); 0.210569 | double_path (2,8); 0.210569 |
| MR | 100 | 2 | 28412.5 | raw (2,16); 3.639 | path (2,16); 0.170604 | double (2,16); 0.220947 | double_path (2,16); 0.176337 | double_path (2,16); 0.176337 |
| MR | 100 | 3 | 48879.8 | raw (2,32); 1.72188 | path (2,32); 0.156803 | double (2,32); 0.169124 | double_path (2,32); 0.160043 | double_path (2,32); 0.160043 |
| MR | 100 | 4 | 84090.8 | raw (2,64); 0.791221 | path (2,64); 0.15065 | double (2,64); 0.154229 | double_path (2,64); 0.152368 | double_path (2,64); 0.152368 |
| MR | 100 | 5 | 144666 | raw (2,64); 0.791221 | path (2,64); 0.15065 | double (2,64); 0.154229 | double_path (2,64); 0.152368 | double_path (2,64); 0.152368 |
| MR | 100 | 6 | 248878 | raw (2,128); 0.330488 | path (2,128); 0.14804 | double (2,128); 0.149192 | double_path (2,128); 0.148968 | double_path (2,128); 0.148968 |
| MR | 100 | 7 | 428160 | raw (2,128); 0.330488 | path (2,128); 0.14804 | double (2,128); 0.149192 | double_path (2,128); 0.148968 | double_path (2,128); 0.148968 |
| MR | 100 | 8 | 736590 | raw (2,128); 0.330488 | path (2,128); 0.14804 | double (2,128); 0.149192 | double_path (2,128); 0.148968 | double_path (2,128); 0.148968 |
| MR | 100 | 9 | 1.2672e+06 | raw (2,128); 0.330488 | path (2,128); 0.14804 | double (2,128); 0.149192 | double_path (2,128); 0.148968 | double_path (2,128); 0.148968 |
| MR | 400 | 0 | 9600 | raw (2,8); 23.4181 | path (2,8); 0.166812 | double (2,8); 0.518894 | double_path (2,8); 0.171999 | double_path (2,8); 0.171999 |
| MR | 400 | 1 | 16515.5 | raw (2,8); 23.4181 | path (2,8); 0.166812 | double (2,8); 0.518894 | double_path (2,8); 0.171999 | double_path (2,8); 0.171999 |
| MR | 400 | 2 | 28412.5 | raw (2,16); 11.3103 | path (2,16); 0.13573 | double (2,16); 0.22441 | double_path (2,16); 0.139282 | double_path (2,16); 0.139282 |
| MR | 400 | 3 | 48879.8 | raw (2,32); 5.528 | path (2,32); 0.122593 | double (2,32); 0.142867 | double_path (2,32); 0.124586 | double_path (2,32); 0.124586 |
| MR | 400 | 4 | 84090.8 | raw (2,64); 2.70661 | path (2,64); 0.117099 | double (2,64); 0.121079 | double_path (2,64); 0.118223 | double_path (2,64); 0.118223 |
| MR | 400 | 5 | 144666 | raw (2,64); 2.70661 | path (2,64); 0.117099 | double (2,64); 0.121079 | double_path (2,64); 0.118223 | double_path (2,64); 0.118223 |
| MR | 400 | 6 | 248878 | raw (2,128); 1.29533 | path (2,128); 0.115049 | double (2,128); 0.116344 | double_path (2,128); 0.115589 | double_path (2,128); 0.115589 |
| MR | 400 | 7 | 428160 | raw (2,128); 1.29533 | path (2,128); 0.115049 | double (2,128); 0.116344 | double_path (2,128); 0.115589 | double_path (2,128); 0.115589 |
| MR | 400 | 8 | 736590 | raw (2,128); 1.29533 | path (2,128); 0.115049 | double (2,128); 0.116344 | double_path (2,128); 0.115589 | double_path (2,128); 0.115589 |
| MR | 400 | 9 | 1.2672e+06 | raw (2,128); 1.29533 | path (2,128); 0.115049 | double (2,128); 0.116344 | double_path (2,128); 0.115589 | double_path (2,128); 0.115589 |
| P1 | 100 | 0 | 9600 | raw (2,8); 5.28036 | path (2,8); 0.180082 | double (2,8); 0.309268 | double_path (2,8); 0.180755 | double_path (2,8); 0.180755 |
| P1 | 100 | 1 | 16515.5 | raw (2,8); 5.28036 | path (2,8); 0.180082 | double (2,8); 0.309268 | double_path (2,8); 0.180755 | double_path (2,8); 0.180755 |
| P1 | 100 | 2 | 28412.5 | raw (2,16); 2.49711 | path (2,16); 0.146369 | double (2,16); 0.176468 | double_path (2,16); 0.146808 | double_path (2,16); 0.146808 |
| P1 | 100 | 3 | 48879.8 | raw (2,32); 1.17274 | path (2,32); 0.133058 | double (2,32); 0.140313 | double_path (2,32); 0.133387 | double_path (2,32); 0.133387 |
| P1 | 100 | 4 | 84090.8 | raw (2,64); 0.531974 | path (2,64); 0.127271 | double (2,64); 0.128795 | double_path (2,64); 0.127419 | double_path (2,64); 0.127419 |
| P1 | 100 | 5 | 144666 | raw (2,64); 0.531974 | path (2,64); 0.127271 | double (2,64); 0.128795 | double_path (2,64); 0.127419 | double_path (2,64); 0.127419 |
| P1 | 100 | 6 | 248878 | raw (2,128); 0.214672 | path (2,128); 0.124132 | double (2,128); 0.124278 | double_path (2,128); 0.124246 | double_path (2,128); 0.124246 |
| P1 | 100 | 7 | 428160 | raw (2,128); 0.214672 | path (2,128); 0.124132 | double (2,128); 0.124278 | double_path (2,128); 0.124246 | double_path (2,128); 0.124246 |
| P1 | 100 | 8 | 736590 | raw (2,128); 0.214672 | path (2,128); 0.124132 | double (2,128); 0.124278 | double_path (2,128); 0.124246 | double_path (2,128); 0.124246 |
| P1 | 100 | 9 | 1.2672e+06 | raw (2,128); 0.214672 | path (2,128); 0.124132 | double (2,128); 0.124278 | double_path (3,16); 0.0938583 | double_path (3,16); 0.0938583 |
| P1 | 400 | 0 | 9600 | raw (2,8); 20.9844 | path (2,8); 0.180416 | double (2,8); 0.489339 | double_path (2,8); 0.181078 | double_path (2,8); 0.181078 |
| P1 | 400 | 1 | 16515.5 | raw (2,8); 20.9844 | path (2,8); 0.180416 | double (2,8); 0.489339 | double_path (2,8); 0.181078 | double_path (2,8); 0.181078 |
| P1 | 400 | 2 | 28412.5 | raw (2,16); 10.1476 | path (2,16); 0.148842 | double (2,16); 0.223541 | double_path (2,16); 0.149322 | double_path (2,16); 0.149322 |
| P1 | 400 | 3 | 48879.8 | raw (2,32); 4.95707 | path (2,32); 0.134847 | double (2,32); 0.149351 | double_path (2,32); 0.1351 | double_path (2,32); 0.1351 |
| P1 | 400 | 4 | 84090.8 | raw (2,64); 2.41553 | path (2,64); 0.128167 | double (2,64); 0.130179 | double_path (2,64); 0.128314 | double_path (2,64); 0.128314 |
| P1 | 400 | 5 | 144666 | raw (2,64); 2.41553 | path (2,64); 0.128167 | double (2,64); 0.130179 | double_path (2,64); 0.128314 | double_path (2,64); 0.128314 |
| P1 | 400 | 6 | 248878 | raw (2,128); 1.14386 | path (2,128); 0.127009 | double (2,128); 0.127348 | double_path (2,128); 0.127104 | double_path (2,128); 0.127104 |
| P1 | 400 | 7 | 428160 | raw (2,128); 1.14386 | path (2,128); 0.127009 | double (2,128); 0.127348 | double_path (2,128); 0.127104 | double_path (2,128); 0.127104 |
| P1 | 400 | 8 | 736590 | raw (2,128); 1.14386 | path (2,128); 0.127009 | double (2,128); 0.127348 | double_path (2,128); 0.127104 | double_path (2,128); 0.127104 |
| P1 | 400 | 9 | 1.2672e+06 | raw (2,128); 1.14386 | path (2,128); 0.127009 | double (2,128); 0.127348 | double_path (2,128); 0.127104 | double_path (2,128); 0.127104 |

## Best n=3 versus n=2 (reported, no verdict)

| PDE | d | method family | best n=2 | best n=3 | n3/n2 |
|---|---:|---|---|---|---:|
| P1 | 100 | path | path M=128; 0.124132 | path M=16; 2.64376 | 21.298 |
| P1 | 100 | best_double | double_path M=128; 0.124246 | double_path M=16; 0.0938583 | 0.755426 |
| P1 | 400 | path | path M=128; 0.127009 | path M=16; 10.3223 | 81.2722 |
| P1 | 400 | best_double | double_path M=128; 0.127104 | double_path M=16; 0.165191 | 1.29965 |
| MR | 100 | path | path M=128; 0.14804 | path M=16; 4.30589 | 29.086 |
| MR | 100 | best_double | double_path M=128; 0.148968 | double_path M=16; 0.15368 | 1.03163 |
| MR | 400 | path | path M=128; 0.115049 | path M=16; 12.1578 | 105.675 |
| MR | 400 | best_double | double_path M=128; 0.115589 | double_path M=16; 0.185809 | 1.60749 |

## Provenance

- Frozen preregistration commit: `8912ae69c0cd4c7bf53f3e2624ba93aa864771dd`.
- Analysis code commit: `e258255eb2ef7bcccfe2323102a524de8387f2cc`.
- Frozen double-estimator source/results commit: `50c3c7d8c83db707cd6c6878835132547a7d2fab`.
- New-row code commits: `e258255eb2ef7bcccfe2323102a524de8387f2cc`.
- Combined rows: 1920 / 1920; missing: 0.
- Path rows: 480 / 480; reused: 160; newly computed: 320.
- New-row aggregate worker wall time: 26064.9 seconds.
- Non-finite state values: 0; non-finite generator values: 0.
- Paired primary-draw audit: 480 pairs, 0 fingerprint mismatches, verdict PASS.
- The paired audit traces one complete registered chunk for every task identity. It does not rerun or overwrite any existing 1,200-point result row.
- Test points: the frozen 1,200 points per PDE/dimension and frozen 20% validation split; verdicts use the 960-point test subset.
- Arithmetic: float64 throughout; 10 repetitions per cell; bootstrap: 1,000 paired repetition draws with seed 20261201.
- Primary cost: mean generator calls. Mean wall time is carried as the secondary cost in both frontier output tables.

## Figure

- `results/path_frontier/figures/frontier.png`

## Verification

- The focused new-and-dependent test suite passed: 27 tests passed (`path_frontier`, `double_estimator`, `mechanism_suite`, and `expert_iteration`).
- Artifact audit: 1,920 unique combined rows, no duplicates, exactly 10 repetitions per registered cell, 320 new atomic JSON artifacts, and no missing tasks.
- Numerical audit: no non-finite state or generator values; every row uses float64, base seed 20261201, 1,200 fixed points, and the registered dimension-dependent chunk size.
- Pairing audit: every one of the 480 path/raw row identities agrees on point fingerprint, generator calls, base seed, chunk size, and primary seed scheme. The independent traced-chunk audit additionally found 0/480 primary-draw fingerprint mismatches and 0/480 draw-count mismatches.
- Cost audit: all 40 evaluation rows use the exact frozen D-4 cost levels; all requested bootstrap intervals are present.
- The frontier figure was visually inspected at original resolution for clipping, malformed log axes, unreadable labels, and missing methods; no rendering issue was found.
- No existing module, manuscript file, or existing result was modified.
