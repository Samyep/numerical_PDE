# Main high-dimensional pilot on RTX 5070 Ti

## Outcome

| PDE | Decision | Evidence |
|---|---|---|
| B6 full-rank HJB | **GO_MAIN** | Mature reference and artifacts pass reuse audit; nonlinear skill 0.263 (100D) / 0.242 (400D), effective rank 29.0 / 77.1, stable learned iterations. |
| Dense Fourier joint | **GO_MAIN** | Exact joint $u,z$ reference; nonlinear skill 0.274, effective rank 23.7, and a2/a0 RMSE ratio 0.6896 under the frozen two-update pilot. |
| Published coupled LQG | **NEEDS_METHOD_FIX** | Reference is validated, but the problem is nearly heat-only (skill 0.0155) and learned policies run away; the reference-free normalized selector chooses the worst checkpoint. Even after stabilization it is a supporting, not headline, benchmark. |

The decisive result is mixed but useful: Dense Fourier supplies the missing joint value–gradient nonlinear case, while the recreated published LQG exposes an out-of-domain policy/selector instability rather than a success case.

## Main summary

| PDE | d | Effective rank | Nonlinear skill | a0 RMSE | a1 RMSE | a2 RMSE | a4 RMSE | MC contribution at final shown | Total GPU time |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| B6 full-rank HJB | 100 | 28.99 | 0.263 | 0.00165 | 0.001343 | 0.001337 | 0.001334 | 0.00135 | 1017.4 s |
| B6 full-rank HJB | 400 | 77.14 | 0.2418 | 0.003873 | 0.003386 | 0.003063 | 0.002881 | 0.002549 | 1004.8 s |
| Dense Fourier joint | 100 | 23.73 | 0.2744 | 0.002918 | 0.002257 | 0.002012 | — | 0.001767 | 294.8 s |
| Published coupled LQG | 100 | 15.2 | 0.01546 | 0.05961 | 20.18 | 1.334e+05 | — | 8480 | 166.5 s |

B6 uses `REUSED_MATURE_RESULT`: 8,192 states, 128 paths/label, 32 steps, 4×384 network, 5,000 optimizer steps per update and four updates. Its timing is therefore not an equal-budget comparison with the two-update pilots.

## Mathematical and artifact correctness

- Validation gate: **VALIDATED**. Dense terminal identity, analytic gradient, PDE residual, tangent derivatives/reconstruction, terminal-aware anchoring, and plain-FK/MCV consistency all pass.
- LQG: the tridiagonal matrix is SPD; the rescaled Gauss–Laguerre reference agrees across 96/144/192/256 nodes at roundoff, its analytic gradient agrees with finite differences, the PDE residual is at roundoff, and an independent 262,144-sample Monte Carlo check is within 1.66 standard errors.
- Prediction/checkpoint/B6 reuse audit: **PASS**. Every saved prediction hash was checked and each reported RMSE, relative $L^2$, bias and MC contribution was recomputed from the `.npz` arrays.
- No truth entered labels, fitting, early stopping, or checkpoint selection. Truth arrays were loaded only after the selector had completed.

## Dense Fourier joint PDE

This is an analytically integrable high-dimensional stress test, not a claim of a generic closed-form nonlinear PDE. It has genuine joint dependence through $f(u,z)=\tfrac12\tanh(u)\lVert z\rVert^2$ and uses the required anchored value and gradient heads.

| Actor | Mean RMSE | Pooled two-replicate RMSE (diagnostic) | MC contribution | Direct value-head RMSE | Direct gradient-vector RMSE |
|---|---:|---:|---:|---:|---:|
| a0 | 0.00291839 | 0.00258964 | 0.0018852 | 0.185361 | 0.126415 |
| a1 | 0.002257 | 0.00189626 | 0.00178871 | 0.0262264 | 0.061565 |
| a2 | 0.00201249 | 0.00159883 | 0.00176687 | 0.0253169 | 0.0579381 |

The frozen reference-free selector chose **a2**. a2 reduces mean solution RMSE by 31.0% versus a0; pooling the two already-saved independent predictions gives a larger 38.3% reduction. Because all three solution RMSEs are near their MC contributions, the pooled figures are diagnostics rather than a new budget or headline criterion.

## Recreated published coupled LQG

This is labelled `RECREATED_PUBLISHED_BENCHMARK`, not an official exact reproduction. Section 3.3 of the paper uses $[0,0.5]\times B^d$; the official repository defaults to $T=1$ and its coefficient arrays depend on batch shape. We use the paper interval, the unit-ball law, and one globally frozen coefficient draw.

| Actor | Mean RMSE | Pooled RMSE | MC contribution | Direct gradient-vector RMSE | Selector absolute component RMSE | Selector normalized score |
|---|---:|---:|---:|---:|---:|---:|
| a0 | 0.0596127 | 0.0590395 | 0.0116199 | 2.28807 | 0.233521 | 54.4355 |
| a1 | 20.1828 | 20.1814 | 0.377801 | 1.51652 | 1.65133 | 1.00446 |
| a2 | 133400 | 133206 | 8479.7 | 13.1837 | 11419.8 | 0.989766 |

The run is finite but unstable. Policy evaluation necessarily visits states outside the unit-ball training distribution; the evidence is consistent with uncontrolled off-domain learned drift, and the frozen linear evaluation becomes enormous. More importantly, self-normalization hides the scale runaway: absolute fixed-point discrepancy grows from 0.234 to 1.65 to 11,420, while normalized scores fall from 54.4 to 1.004 to 0.990, so the selector chooses a2. This is a method/selector failure discovered by the pilot; the recorded choice and failure are preserved unchanged.

The benchmark itself is weakly nonlinear on the paper distribution: heat-only RMSE is 0.0127, only 0.0155 solution standard deviations. Thus even a repaired version should be supporting evidence rather than the headline difficult PDE.

## B6 reuse audit

- Mature source commit: `281977883e0a5ce4f05895878680b582620461be`.
- Frozen B6 bank hashes: 100D `84b2132607bc6b17e44adbea838915f8ce657813f6a0987e34deca308660f1c1`; 400D `0df59abc155e90aa674532760578dc2563cbdb4aade97a849961c8a8edd9f2f9`.
- The selected calibration rows, raw evaluation records, seed-0 training records, selector choice a4, and available selected checkpoint hashes were cross-checked by the read-only adapter.

## Compute and deviations

- Dense new pilot: 294.8 RTX 5070 Ti seconds for both updates, selector, and all six 512-state/1,024-path evaluations.
- LQG new pilot: 166.5 RTX 5070 Ti seconds on the same accounting basis.
- Total new training and inference: 461.3 GPU seconds (0.128 GPU h), far below the 3 h cap.
- No CUDA OOM, chunk fallback, timeout, retry, nonfinite sample, hyperparameter search, extra seed, 200D/400D new training, FH-MLP grid, or formal strong-baseline training occurred.
- LQG a1/a2 instability and the selector scale pathology are failures, not budget truncations or omitted rows.

## Answers and recommended next stage (not run)

1. **Dense Fourier does exercise joint nonlinearity.** Its heat-only skill is 0.274 and driver RMS is 0.265; both value and gradient heads learn substantially.
2. **It is meaningfully high-effective-dimensional.** Effective rank is 23.7 and 50.8% of centred gradient variance lies outside the top ten directions. The two-update learned Actor improves solution RMSE, though the final estimate is MC-floor limited.
3. **Published LQG is close to heat-only on the paper distribution.** It is not a strong-nonlinearity headline benchmark, and this pilot additionally finds a serious policy/selector stability issue.
4. **B6 and Dense Fourier are complementary main cases.** B6 supplies stable full-rank quadratic HJB geometry; Dense supplies dense, joint $u,z$ nonlinearity. LQG should appear only as an external stress/failure case until fixed.
5. **Proposed formal suite:** B6 and Dense at 100D/200D/400D with three seeds; compare Tangent-Actor against DPI, Deep BSDE, PINN+SDGD, heat-only and a0 under measured same-GPU budgets. Planning estimate: 10–12 GPU h for Tangent-Actor training/audits and 24–40 GPU h for all strong baselines. Reserve at most 2 GPU h for a separately pre-registered LQG stabilization pilot (bounded/off-domain actor control plus a scale-aware selector) before deciding whether to include it.

## Provenance

- Frozen protocol commit: `50908e13560649025aefb874a625e55ecf76edd0`.
- Pilot runner commit: `9b22d5ad1d394bcd00143683e7044157881ab199`.
- Study revision: `main_highdim_pilot_5070_v1`.
- Base/highdim source commit: `281977883e0a5ce4f05895878680b582620461be`.
- Official SCaSML repository audited at commit `c2412b942b1720084accf9ac9651c6085e073919`.
- Source reconciliation: [ICLR 2026 paper](https://proceedings.iclr.cc/paper_files/paper/2026/file/cf8353a0efc3644454a5591a7cb55c17-Paper-Conference.pdf) and [official SCaSML repository](https://github.com/Francis-Fan-create/SCaSML).
- Hardware/software details, parameter hashes, bank hashes, seeds, and checkpoint/prediction hashes are stored in `test_bank_manifest.json`, `validation.json`, `artifact_audit.json`, and the raw JSON records.

## Figures

- `results/main_highdim_pilot/figures/actor_iteration_rmse.png`
- `results/main_highdim_pilot/figures/effective_dimension.png`
- `results/main_highdim_pilot/figures/nonlinearity_strength.png`
- `results/main_highdim_pilot/figures/pilot_accuracy_compute.png`
