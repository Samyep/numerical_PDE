# High-Accuracy FH-MLP Amortization Break-Even

## Outcome

| Target | TA-direct final RMSE | TA-FK-MCV final RMSE | FH-MLP config | FH final RMSE | Measured break-even |
|---|---:|---:|---|---:|---|
| primary (0.01) | 0.00333844 | 0.00190042 | n2_M128 | 0.00851991 | (524,288, 1,048,576] |
| secondary (0.005) | 0.00333844 | 0.00190042 | -- | -- | PAIR_NOT_ACCURACY_QUALIFIED |

The previous audited RMSE <= 0.06 result remains unchanged: its observed bracket was (524,288, 1,048,576] queries. The new brackets above use only actually executed common prefixes; no interpolation or single-query extrapolation enters the verdict.

## Frozen calibration

| Target | Method | Calibration status | Selected configuration | Calibration RMSE |
|---|---|---|---|---:|
| primary | ta_direct | SELECTED | direct | 0.00320842 |
| primary | ta_fk_mcv | SELECTED | mcv_32 | 0.00204158 |
| primary | fh_mlp_standard | SELECTED | n2_M128 | 0.00973479 |
| secondary | ta_direct | SELECTED | direct | 0.00320842 |
| secondary | ta_fk_mcv | SELECTED | mcv_32 | 0.00204158 |
| secondary | fh_mlp_standard | TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET | -- | 0.00859731 |

## Measured cumulative timing

- **primary:** observed bracket `(524,288, 1,048,576]`; at the first crossing TA used `22.839 +/- 0.051` minutes and FH-MLP used `25.854 +/- 0.084` minutes.
  At `1,048,576` queries, TA used `11.7%` less measured cumulative computation.
- **secondary:** `PAIR_NOT_ACCURACY_QUALIFIED`.

## Correctness, accounting, and interpretation

- Artifact audit: **VALIDATED**.
- Actual allocation: **3.564 / 30 H100 GPU-hours**, including failed allocations and retries recorded by Slurm.
- Reference refinement RMSE: `1.017e-05`.
- FH batching audit: **NOT_ADOPTED**. A multi-query call changes the random-stream assignment and therefore is not an exact same-query, same-sample replay of the audited per-point recursion. No faster batched implementation was substituted.
- TA training and reference-free checkpoint selection are charged in every cumulative curve even though byte-identical audited checkpoints were reused; new experimental allocation is reported separately.
- The main timing curve uses `TA-direct`. The current manuscript describes FK-MCV as the final evaluation estimator, so FK-MCV accuracy and timing are explicitly separated rather than relabelled as direct inference.
- A missing crossing means `NOT_REACHED_WITHIN_TESTED_BUDGET`; it is not converted into an interpolated claim.

## Recorded failures and scheduler deviations

| Stage | Slurm allocation | State | Elapsed (s) | H100 GPU-hours |
|---|---|---|---:|---:|
| preflight | 6633981 | FAILED | 62 | 0.017222 |
| profile | 6633999_[0-27%8] | CANCELLED by 3335748 | 0 | 0.000000 |
| final_accuracy | 6634793_1 | FAILED | 6 | 0.001667 |
| final_accuracy | 6634793_2 | FAILED | 3 | 0.000833 |
| final_accuracy | 6634793_3 | FAILED | 3 | 0.000833 |
| final_accuracy_retry_2 | 6635275_3 | TIMEOUT | 129 | 0.035833 |

- The first preflight failed on an over-broad prior/prior bank-overlap assertion; the fresh-bank separation rule was corrected and revalidated before calibration.
- Three final-array allocations were rejected by the frozen commit guard after audit-only commits were synchronized to the run checkout. They produced no numerical records; the same manifest indices were rerun at the frozen numerical commit.
- The FH final bundle wrote all five complete records and its complete task-status artifact before a deliberately shortened two-minute scheduler limit expired during environment teardown. Slurm therefore records `TIMEOUT`; every prediction file, finite-value check, input hash and prediction hash passed the final artifact audit.

## Provenance

- Preregistration commit: `60c04015682a5f77ff03fd829f575d059620a1bd`.
- Report-generation commit: `bde292577e2e89aeef524d05b0ea87ed5feff6c4`.
- Frozen selection hash: `d8f701c82d613432545fb0f232666b48e01a38a32716f5e316b5084f8f30334f`.
- Raw timing: `results/high_accuracy_fhmlp_break_even/timing_replicates.csv`.
- Main figure: `results/high_accuracy_fhmlp_break_even/figures/high_accuracy_break_even_primary.pdf`.
