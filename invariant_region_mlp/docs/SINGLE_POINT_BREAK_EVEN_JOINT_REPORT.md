# Joint Nonseparable Matched-Accuracy Amortization Break-Even

## Outcome

| Target | Method | Frozen configuration | Final status | Final pooled RMSE | Measured break-even |
|---|---|---|---|---:|---|
| epsilon_1 (6.0e-02) | Tangent-Actor | direct | FINAL_ACCURACY_PASSED | 0.0037644 | reference curve |
| epsilon_1 (6.0e-02) | Full-History MLP | n2_M16 | FINAL_ACCURACY_PASSED | 0.017657 | OBSERVED_BREAK_EVEN |
| epsilon_1 (6.0e-02) | Deep BSDE | -- | TARGET_NOT_REACHED | -- | TARGET_NOT_REACHED |
| epsilon_2 (2.0e-02) | Tangent-Actor | direct | FINAL_ACCURACY_PASSED | 0.0037644 | reference curve |
| epsilon_2 (2.0e-02) | Full-History MLP | n2_M32 | FINAL_ACCURACY_PASSED | 0.011861 | NOT_EVALUATED_COMPUTE_CAP |
| epsilon_2 (2.0e-02) | Deep BSDE | -- | TARGET_NOT_REACHED | -- | TARGET_NOT_REACHED |

The crossing is the first actually executed common query prefix at which the
mean baseline total is no smaller than the mean Tangent-Actor total.  No
single-query extrapolation, curve fit, interpolation, ensemble prediction, or
post-final configuration change is used.

## Measured amortization result

- **epsilon_1, Full-History MLP:** observed in `(524288, 1048576]` queries; at the first measured crossing TA used `1420.284` s and the baseline `1521.897` s.
- **epsilon_1, Deep BSDE:** `TARGET_NOT_REACHED`.
- **epsilon_2, Full-History MLP:** `NOT_EVALUATED_COMPUTE_CAP`.
- **epsilon_2, Deep BSDE:** `TARGET_NOT_REACHED`.

## Frozen calibration choices

- **epsilon_1:** `tangent_actor` = `direct` (calibration RMSE 0.0046623); `fh_mlp_standard` = `n2_M16` (calibration RMSE 0.020764); `deep_bsde` = `None` (calibration RMSE 0.016322).
- **epsilon_2:** `tangent_actor` = `direct` (calibration RMSE 0.0046623); `fh_mlp_standard` = `n2_M32` (calibration RMSE 0.013974); `deep_bsde` = `None` (calibration RMSE 0.016322).

Calibration used 32 new queries and final accuracy used a disjoint 64-query
bank, except for the mandated shared origin.  Tangent-Actor is trained once;
FH-MLP re-executes its recursion per point; native Deep BSDE trains a fresh
scalar-$Y_0$ model per point.  Secondary timing is deliberately omitted when
it would require a different selected configuration, preserving the frozen
18 H100-hour ceiling.

## Correctness and accounting

- Independent reference: **VALIDATED**; fresh-bank primary/audit RMSE `1.059356e-05`.
- H100/source preflight: **VALIDATED** on `NVIDIA H100 80GB HBM3`.
- Artifact audit: **VALIDATED**.
- Allocated H100 task-hours, including recorded failures: `5.8382` / `18`.

Method-active time includes Actor label construction, fitting, reference-free
selection and inference for TA; recursive inference for FH-MLP; and complete
fresh per-query training for Deep BSDE.  The 1,048,576-state timing bank has
no truth array.  Queueing, reference construction, truth evaluation, hashing,
plotting and upload are excluded from method curves.

## Provenance

- Frozen addendum commit: `8a18360be907da36b8b875c4fdf6ba0b0cc1007d`.
- Imported validated-source commit: `b7dd884225936a0d67f770efce1207bc478d4f7c`.
- Frozen implementation commit: `a39b7192a691bfff8329fc101f5db15027ee2d7f`.
- Frozen selection hash: `a14f2c78d1dae61912f1820eb60e884c3aa9f9856031b22ce7b433e9236b0e38`.
- Analysis commit: `7111bd0de808f40801ff677ee90774407f859b5b`.
- Figure: `results/single_point_break_even_joint/figures/joint_matched_accuracy_break_even.pdf`.

Large checkpoints and the truth-free timing bank remain in the ICE archive;
their hashes and every reported prediction hash are retained in compact
records.  This is the final preregistered extension; no third PDE is started.
