# Pathwise-only equal-cost frontier: pre-registered study

Repository `Samyep/numerical_PDE`, subproject `invariant_region_mlp`. Branch from `ir-mlp-double-estimator` (commit `50c3c7d8` or later) and create branch `ir-mlp-path-frontier`. You are authorized to run long experiments locally and use all available CPU/GPU resources reasonably, but do not delete existing results or rewrite the manuscript.

## Hard rules

1. Do not modify any existing module or result. New code goes in `invariant_region_mlp/experiments/path_frontier/` and must **import** the frozen double-estimator code (`experiments/double_estimator/`) unchanged; results go in `results/path_frontier/`, the report in `docs/PATH_FRONTIER_REPORT.md`.
2. Freeze this document first: commit it verbatim as `experiments/path_frontier/PREREGISTRATION.md` before any implementation or run, and record the commit hash in the report.
3. Use exactly the double-estimator protocol: base seed `20261201`, the same primary seed `SeedSequence([20261201, d, n, M, rep, chunk_index])`, the same chunk size, the same fixed 1,200 test points with the 20% validation split, float64, 10 repetitions. The new `path` rows must be paired with the existing `raw`, `double` and `double_path` rows (verify that the primary-draw fingerprints of `path` and `raw` match, as in the previous audit).
4. Verdicts are mechanical; report every failure; no post-hoc thresholds.

## Question

In the double-estimator study the equal-cost frontier of the repaired MLP is attained at n=2, where the pathwise terminal gradient and the double estimator repair the same noise. The `path` method was not run on the equal-cost grid. Does the double estimator add anything **at the frontier**, or does pathwise alone suffice?

## Runs

Method `path` only (pathwise terminal gradient, raw correction levels; the existing `path` implementation), on the equal-cost grid of the double-estimator study:

- PDEs: P1 and MR (k=10, scale=4, T=0.1), d in {100, 400};
- cells: `n=2: M in {8,16,32,64,128}`, `n=3: M in {4,6,10,16}`, `n=4: M in {3,4,6}`;
- 10 repetitions. Reuse every existing row (`raw`, `double`, `double_path`, and `path` where it already exists in the core grid); do not rerun them.

Cost axis: mean generator calls (primary), wall time (secondary), exactly as in D-4. Frontier: lower envelope of mean test skill against cost; evaluate at the same 10 cost levels as D-4; 1,000-draw bootstrap over repetitions.

## Pre-registered criteria

- **E-1 (does double add at the frontier?)** Let `best_double` be the lower envelope over `double` and `double_path`. "Double adds value at the frontier" if `best_double <= 0.9 x path` at >= 8 of the 10 cost levels, for at least 3 of the 4 (PDE, d) problems. Otherwise report "pathwise alone suffices at the frontier".
- **E-2 (pathwise vs raw):** the `path` frontier is <= 0.8 x the `raw` frontier at >= 8 of 10 cost levels, for each of the 4 problems.
- **E-3 (dimension robustness of the frontier):** for each of `path`, `best_double`, `raw`, report the ratio frontier(d=400) / frontier(d=100) at the median cost level. Prediction (no verdict): <= 1.3 for `path` and `best_double`, >= 2 for `raw`.
- **Reported, no verdict:** the cell at which each frontier is attained at each cost level, and the best `path` and `best_double` skill at n=3 versus n=2.

## Outputs
```bash
results/path_frontier/  rows.csv  frontier.csv  frontier_levels.csv  analysis_summary.json  figures/frontier.png
docs/PATH_FRONTIER_REPORT.md
```

One figure: frontiers of raw, path, double, double\_path (log cost axis), one panel per (PDE, d). The report starts with an outcome table of E-1 and E-2, then each criterion verbatim with its verdict, then the E-3 table, then provenance.
