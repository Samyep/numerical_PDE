# Joint nonseparable PDE extension report

> Prospective comparison frozen at preregistration commit `0f1485703`. The
> PDE family had prior Tangent-Actor-only evidence; the query banks and all
> FH-MLP/Deep-BSDE comparisons here are new. This report does not modify the manuscript.

## Outcome

| Criterion | Verdict | Mechanical result |
|---|---:|---|
| J-1 | PASS | 5/5 dimensions: [20, 50, 100, 200, 400] |
| J-2 | PASS | 5/5 dimensions: [20, 50, 100, 200, 400] |
| J-3 | PASS | 3/5 dimensions: [20, 50, 100] |
| J-4 | PASS | d400/d20 = 1.05195 |

## Primary RMSE

| d | Tangent-Actor | FH-MLP (n=3,M=16) | Deep BSDE | Best |
|---:|---:|---:|---:|---|
| 20 | 0.00118191 | 0.0113507 | 0.00258583 | Tangent-Actor |
| 50 | 0.00106482 | 0.0256061 | 0.0037551 | Tangent-Actor |
| 100 | 0.00103796 | 0.051534 | 0.0348293 | Tangent-Actor |
| 200 | 0.00121731 | 0.102422 | NUMERICAL FAILURE | Tangent-Actor |
| 400 | 0.00124331 | 0.213023 | NUMERICAL FAILURE | Tangent-Actor |

## Reference and provenance

- Reference gate: **VALIDATED**; 256² vs 128² value RMSE `1.00079e-05`; PDE residual RMSE `5.09455e-07`.
- Driver: `f(u,z)=u-u^3-0.5*||z+u v||^2`, with `nabla_z f_u=-v`.
- Methods use frozen budgets; truth enters only post-training evaluation.
- Artifact audit: **VALIDATED**; raw JSON files `636`, prediction files `631`.
- Code commit at aggregation: `b7dd884225936a0d67f770efce1207bc478d4f7c`.

## Full FH grid

All `30` dimension/configuration rows are in `results/single_point_joint_extension/fh_grid.csv`; no truth-based cell was selected.

## Deviations and failures

Deep BSDE produced nonfinite predictions and is recorded as a numerical failure at dimensions `[200, 400]`. These cells are excluded from the finite-RMSE count for J-3; no replacement configuration was run. After every preregistered seed in a cell had produced a nonfinite result, the remaining queries in that failed cell were not attempted. Details: d=200: first failures [{'seed': 0, 'point': 4}, {'seed': 1, 'point': 0}, {'seed': 2, 'point': 10}]; attempted 17/96 runs; d=400: first failures [{'seed': 0, 'point': 0}, {'seed': 1, 'point': 1}, {'seed': 2, 'point': 0}]; attempted 4/96 runs
