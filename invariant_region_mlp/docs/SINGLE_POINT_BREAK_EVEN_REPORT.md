# Matched-Accuracy Single-Point Amortization Break-Even

## Outcome

| Target | Method | Frozen configuration | Final status | Final pooled RMSE | Measured break-even |
|---|---|---|---|---:|---|
| epsilon_1 (2.5e-04) | Tangent-Actor | mcv_32 | FINAL_ACCURACY_PASSED | 0.000248 | reference curve |
| epsilon_1 (2.5e-04) | Full-History MLP | n3_M16 | FINAL_ACCURACY_PASSED | 0.0002262 | NOT_REACHED_WITHIN_TESTED_BUDGET |
| epsilon_1 (2.5e-04) | Deep BSDE | updates_4000 | FINAL_ACCURACY_PASSED | 0.0001879 | (8, 16] queries |
| epsilon_2 (1.0e-04) | Tangent-Actor | mcv_256 | FINAL_ACCURACY_PASSED | 7.401e-05 | NOT_EVALUATED_COMPUTE_CAP |
| epsilon_2 (1.0e-04) | Full-History MLP | n4_M12 | FINAL_ACCURACY_PASSED | 9.183e-05 | NOT_EVALUATED_COMPUTE_CAP |
| epsilon_2 (1.0e-04) | Deep BSDE | -- | TARGET_NOT_REACHED | -- | TARGET_NOT_REACHED |

The headline crossing is the first **actually executed** common query-count
prefix at which mean baseline time is no smaller than mean Tangent-Actor total
time.  No line fitted to one query and no interpolated crossing is used.
The strict-target accuracy result is retained even when its runtime comparison
is marked `NOT_EVALUATED_COMPUTE_CAP`; no reduced-replicate curve is substituted.
Tangent-Actor's primary final RMSE is a mechanical pass but a narrow one:
it lies only `0.806%` below the frozen threshold.  The
reference-uncertainty audit remains two orders of magnitude smaller than that
threshold; the narrow margin is not hidden or reclassified.

## Measured amortization result

- **epsilon_1, Full-History MLP:** `NOT_REACHED_WITHIN_TESTED_BUDGET`; at the largest common measured prefix (`262144` queries), TA used `3320.803` s and the baseline used `578.146` s (baseline/TA = `0.1741`).
- **epsilon_1, Deep BSDE:** observed between `8` and `16` queries; at the first crossing, TA used `801.253` s and the baseline used `1185.948` s.
- **epsilon_2, Full-History MLP:** `NOT_EVALUATED_COMPUTE_CAP`.
- **epsilon_2, Deep BSDE:** `TARGET_NOT_REACHED`.

These statements stop at the largest prefix actually executed by all three
replicates.  In particular, a missing crossing is not converted into an
asymptotic claim by fitting per-query slopes.

## Frozen calibration choices

- **epsilon_1:** `tangent_actor` = `mcv_32` (calibration RMSE 0.0002332); `fh_mlp_standard` = `n3_M16` (calibration RMSE 0.0002316); `deep_bsde` = `updates_4000` (calibration RMSE 0.0001988).
- **epsilon_2:** `tangent_actor` = `mcv_256` (calibration RMSE 9.925e-05); `fh_mlp_standard` = `n4_M12` (calibration RMSE 7.983e-05); `deep_bsde` = `None` (calibration RMSE 0.0001547).

Calibration used 32 new queries.  Final accuracy used a disjoint 64-query
bank and could not change any configuration.  Tangent-Actor is trained once;
Full-History MLP reruns its recursion at every point; native Deep BSDE trains
a new scalar-$Y_0$ model at every point.

## Reference and correctness gates

- Reference gate: **VALIDATED**.
- Reference uncertainty RMS: `1.955633e-07` (limit `5e-6`).
- Maximum reference uncertainty: `1.628558e-06` (limit `1e-5`).
- 200-state PDE/terminal audit: **PASS**.
- H100 numerical/source preflight: **VALIDATED** on `NVIDIA H100 80GB HBM3`.
- Final artifact audit: **VALIDATED**.

The frozen prose contained one unavoidable set-theoretic conflict: it required
both calibration and final banks (and the previous bank) to contain the exact
origin while also calling all banks disjoint.  The mandated origin is retained
and reported as the sole overlap; every stochastic row is disjoint.  Bank
audit status: `VALIDATED_WITH_MANDATED_ORIGIN_EXCEPTION`.

## Timing and accounting

Primary times are synchronized method-active times on NVIDIA H100 80GB GPUs.
Queue time, reference construction, truth evaluation, plotting, and artifact
upload are excluded from the curves.  The task records currently account for
`15.362` allocated H100 task-hours across `145` tasks,
including `17` recorded failed tasks.  Detailed rows are in
`results/single_point_break_even/resource_accounting.csv`.

For Tangent-Actor, every plotted total contains its measured Actor training,
reference-free selection, checkpoint loading, and the executed query prefix.
For Full-History MLP it contains all recursive inference.  For Deep BSDE it
contains all fresh per-query training.  The timing bank has 524,288 Gaussian
states and deliberately contains no truth array.

## Provenance

- Preregistration commit: `498589cc2`.
- Frozen implementation commit: `ae6ff1e1cfec4dd49db168481e5a31d7fc5671ef`.
- Analysis code commit: `64343a6d4e43d3796a1bdbfe0c08f60663e56736`.
- Frozen selection hash: `311804c418ece6276a1fb5f08c5afdb31c1f36d8455a6c5184f56660c53f2f4c`.
- Study revision: `single_point_break_even_v1`.
- Figure: `results/single_point_break_even/figures/matched_accuracy_break_even.pdf`.
- H100 archive (checkpoints, timing bank, raw predictions, and Slurm logs):
  `/storage/ice1/4/8/ywang4280/ir-mlp-single-point-break-even-worktree/invariant_region_mlp/results/single_point_break_even/`.

Large checkpoints and the truth-free 524,288-query timing bank remain in that
ICE archive rather than being placed in Git.  Their hashes and the hashes of
every reported prediction are retained in the machine-readable records.

The optional second PDE is not part of this report.  It may begin only after
this Allen--Cahn artifact audit is validated and under a separately frozen
addendum.
