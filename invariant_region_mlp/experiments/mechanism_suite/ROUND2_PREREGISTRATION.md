<!--
Frozen before any round-2 implementation or computation.
Source attachment SHA-256: 9996F09EC64196CAA4BB7E10493B375358645C4B5283261AE9B5B99A2A376B62
-->

# Mechanism suite, round 2: confirmatory re-run with corrected pre-registration

Agent instructions. Branch from `ir-mlp-mechanism-suite` (HEAD `37666b23`), suggested name `ir-mlp-mechanism-suite-r2`.
Do **not** modify the manuscript, `FullHistoryMLP`, or any round-1 result. Extend the existing
`experiments/mechanism_suite/` code; write results to `results/mechanism_suite_r2/` and the report to
`docs/MECHANISM_SUITE_R2_REPORT.md`.

## 0. Why a round 2, and the rules that keep it honest

Round 1 applied its pre-registered criteria mechanically. Some failures were real; others came from errors in
the round-1 specification:

| Round-1 issue | Type | Round-2 fix |
|---|---|---|
| C2 inequality `box/oracle <= 0.5 raw/oracle` contradicts its own parenthetical "half of the excess error" | spec error | use `Gc = (raw - method)/(raw - oracle) >= 0.5` |
| C3 counted tuned shrink-to-centre as a "certificate prior"; it uses the data | spec error | shrink-to-centre becomes an alternative retraction (R5), not a C3 comparator |
| C3 compared box with `illegal_a0.9`, which is within noise of box | spec error | illegal tighter sets restricted to `a in {0.25, 0.5, 0.75}` |
| C3 was evaluated jointly over all comparators in a cell | ambiguous wording | evaluated per comparator class |
| P4 primary method was the box hull of a 1-D segment | spec error | P4 primary is a segment-type certificate (R3) |
| P4 localization prediction assumed homoscedastic noise | spec error | test localization under injected homoscedastic noise (R4) |
| C2 fails in shallow, well-sampled cells for P2/P3 | **real finding** | now a regime hypothesis tested on new configurations (R1) |
| P4 certificate `|psi_s| <= 1` too loose; stronger suppression wins | **real finding** | test a tighter, still PDE-derived certificate (R3) |

Rules:
1. **Fresh randomness.** Base seed `20261107` (round 1: 20261006). Regenerate test points with a new seed and a new 20% validation split. Problem definitions (equations, parameters, the direction `w` for P1/P4) are unchanged.
2. **Round-1 verdicts stand.** The round-2 report must print the round-1 outcome table unchanged next to the round-2 table and state explicitly that round-2 criteria were written after round 1 was seen.
3. **New configurations** are added so the regime hypothesis is tested out of sample, not only on fresh draws of the same cells.
4. Mechanical verdicts as before; report every failure; no post-hoc thresholds.
5. Value accuracy only. Gradient errors are reported but no gradient claim is tested.

## 1. Shared definitions

- Skill = `RMSE(u_hat - u*) / std(u*)` on the test split. Paired seeds `SeedSequence([20261107, d, n, M, rep, chunk_index])`, identical chunk size for all methods in a study.
- `Gc(method) = (raw - method)/(raw - oracle_state)` on cell means. Also report `method/oracle_state`, per-rep and repetition-averaged skill.
- **Primary method**, fixed in advance: P1 box; P2/P3 box; P4 `tight_segment` (section R3). Also report P1 segment and P4 `segment` and `box` in every table.
- Tuning: every tuned factor chosen by mean validation skill within its `(PDE, d, n, M)` cell, then frozen for the test split.
- Comparator classes for C3:
  - **A, certificate only (data-free):** `centre`.
  - **B, generic suppression (certificate-free):** best tuned of `shrink_c` (c in {0.1, 0.25, 0.5, 0.75}), `z_zero`, `f_zero`.
  - **C, illegal tighter set:** best tuned of `illegal_a`, a in {0.25, 0.5, 0.75}. For P4 the illegal set is the primary certificate shrunk about its own centre.
- Paired win: strictly lower test skill in the same repetition. 10 reps per cell unless stated.

## R1. C2 as a regime hypothesis (stabiliser claim)

Hypothesis: the certified retraction removes most of the noise damage where raw MLP is unstable, and little where raw MLP is already well sampled.

- PDEs: P1, P2_a4, P2_a8, P3_rho1, P3_rho2, P4.
- d: P1 {20, 50, 100}; P2 {20, 50}; P3 {20, 50}; P4 {20, 50}.
- Deep stratum (n >= 4): `(4,3), (4,4), (4,6), (5,2)`. Shallow stratum (n <= 3): `(3,6), (3,10), (3,16), (2,32)`. New relative to round 1: `(4,4), (5,2), (3,16), (2,32)`.
- Methods: raw, primary, oracle_state, oracle_z, plus P1 segment and P4 segment/box.
- Pre-registered criteria:
  - **C2-deep:** primary has `Gc >= 0.5` in >= 75% of deep cells, for every PDE.
  - **C2-shallow (boundary prediction):** for P2/P3, primary has `Gc < 0.5` in >= 50% of shallow cells. For P1 and P4, `Gc >= 0.5` in >= 75% of shallow cells.
  - **C2-graded:** pooled over all P2/P3 cells, Spearman correlation between `Gc` and `log(raw/oracle_state)` is >= 0.5.

## R2. C3 with corrected comparator classes (geometry, not suppression or prior)

- Cells: `(3,6), (3,10), (4,3), (4,6)` at the two smallest d per PDE. All six PDE variants.
- Methods: primary, all class A/B/C comparators, raw; P1 also segment; P4 also segment and box.
- Pre-registered criteria, evaluated **separately per class**:
  - **C3-A:** primary wins >= 7/10 vs centre in >= 75% of cells.
  - **C3-B:** primary wins >= 7/10 vs best class-B control in >= 75% of cells.
  - **C3-C:** primary wins >= 7/10 vs best class-C control in >= 75% of cells.
- Predictions: P1–P3 pass A, B, C. P4 with `tight_segment` passes A, B, C. P4 with the round-1 `segment` is predicted to fail C (round-1 replication).
- The "certificate-prior dominated" flag now refers to class A only.

## R3. P4 with a tighter, PDE-derived certificate (tightness hypothesis)

Round 1 showed that `|psi_s| <= 1` is loose on the test distribution and that stronger-than-certified suppression wins. Hypothesis: a tighter certificate that is still derived from the PDE removes that advantage.

Derivation (write it in the report and verify numerically). P4 reduces to `psi_t + psi_ss - lambda_f |psi_s| = 0`, `psi(T,s) = log cosh(beta s)/beta`. Let `v = psi_s`. Formally `v_t + v_ss + b v_s = 0` with `b = -lambda_f sign(v)`, `|b| <= lambda_f`, `v(T,s) = tanh(beta s)`. By Feynman–Kac, `v(t,s) = E[tanh(beta S_T)]` with `dS = b dt + sqrt(2) dW`. Because `tanh` is increasing, comparison with the constant-drift problems `b = +-lambda_f` gives

```
v_minus(t,s) <= psi_s(t,s) <= v_plus(t,s),
v_pm(t,s) = E tanh( beta (s +- lambda_f (T-t) + sqrt(2 (T-t)) xi) ),  xi ~ N(0,1)
```

Compute `v_pm` with 80-node Gauss–Hermite.

- **tight_segment:** `c = (w . z)/sqrt2`, clip c to `[v_minus(t,s), v_plus(t,s)]`, return `sqrt2 c w`.
- **tightness family**, all valid: interval `[(1-theta)(-1) + theta v_minus, (1-theta)(+1) + theta v_plus]`, theta in {0, 0.25, 0.5, 0.75, 1}. theta=0 is the round-1 segment; theta=1 is tight_segment.
- G5 containment test: 1e5 points from the test distribution plus 1e4 points with t near T and s near 0. Report the maximum violation of the reference `psi_s` outside `[v_minus, v_plus]`.
  The bound is **sharp**: where `psi_s` keeps one sign, the drift is exactly `-+lambda_f` and `psi_s` lies on an edge of the interval. Reference discretisation error can therefore appear as a small violation. Rule: compute the violation on at least three successive reference grids. R3 proceeds only if the violation decreases under refinement and the finest-grid violation is below 1e-5. Otherwise stop R3 and report.
  (Pre-check done before handing over: an independent first-order upwind solve gave violations 3.8e-4, 1.4e-4, 5.1e-5 at ds = 0.02, 0.01, 0.005 for T-t = 0.1, and 1.5e-4, 4.0e-5, 9.9e-6 for T-t = 0.25; interval at s=0 is [-0.13, 0.13], [-0.24, 0.24], [-0.35, 0.35] for T-t = 0.1, 0.25, 0.5.)
- Implementation: never shrink the interval by a safety margin; project exactly onto `[v_minus, v_plus]`.
- Runs: R1 and R2 cells for P4, plus the tightness family at `(3,6), (4,6)`, d in {20, 50}, 10 reps.
- Pre-registered criteria:
  - **C-tight-1:** test skill decreases monotonically in theta (mean over reps) in >= 75% of the tightness-family cells.
  - **C-tight-2:** tight_segment passes C3-A, C3-B and C3-C (R2 thresholds).

## R4. Rectification localization under homoscedastic noise (corrected mechanism test)

- P4, configurations `(4,3), (3,6)`, d in {20, 50}, 10 reps.
- Dose rows as in round-1 E2 (`s in {0.5, 1, 2, 4}`, noise added to the exact z, scale = certificate half-width).
- For every generator call record `s = w . x_child` and the dose-induced generator error `f(z* + eps) - f(z*)`. Bin s as in round 1.
- Additionally, for raw E1 rows, report the normalised bias `|bias| / RMS(z_hat - z*)` per s-bin.
- Pre-registered criteria:
  - **C-loc-1:** for every dose >= 1, the central bin `[-0.25, 0.25)` has the largest mean |dose-induced generator error|.
  - **C-loc-2:** in raw E1 rows, the normalised bias is largest in the central bin.

## R5. Hard projection versus soft certificate-informed contraction (descriptive with direction)

- P2/P3, all R1 cells. Compare primary (hard box) with tuned `shrink_centre_c`, c in {0.25, 0.5, 0.75}.
- Directional prediction: soft has lower mean skill than hard in >= 50% of shallow cells, and hard has lower mean skill than soft in >= 75% of deep cells.
- Report the mean difference with a paired bootstrap 95% interval per cell. No C-claim depends on R5.

## R6. Dimension stability with a pre-registered threshold

- P1, `(4,6)`, d in {20, 50, 100, 200, 400} (d=400 is new), 10 reps; methods raw, box, segment, oracle_state, centre.
- Criterion **C-dim:** `max_d(box/oracle) / min_d(box/oracle) <= 1.25`.

## R7. What is left between box and oracle (descriptive)

For all R1 cells, report per-rep skill, repetition-averaged skill, and their difference for box and oracle_state
(bias vs variance), plus mean signed generator error at child states for box and oracle_state.
Question to answer in text: is the residual box/oracle gap mainly bias introduced or left by the projection?

## Priority, compute and outputs

Run order: R3 (with its G5 test) -> R2 -> R1 -> R4 -> R6 -> R5 -> R7. R5 and R7 reuse R1 rows wherever possible.
Round 1 used about 42 worker-hours; round 2 should be of similar size. If compute is limited, drop d=50 from R1 for P2/P3 before dropping any new configuration.

```
results/mechanism_suite_r2/
  r1_regime.csv  r2_controls.csv  r2_tuning_choices.json  r3_tightness.csv  r3_containment.json
  r4_localization.csv  r5_soft_vs_hard.csv  r6_dimension.csv  r7_gap.csv
  analysis_summary.json  figures/
docs/MECHANISM_SUITE_R2_REPORT.md
```

The report starts with two outcome tables: round 1 (unchanged) and round 2. It then lists every pre-registered criterion above verbatim with its verdict, followed by all failures.
