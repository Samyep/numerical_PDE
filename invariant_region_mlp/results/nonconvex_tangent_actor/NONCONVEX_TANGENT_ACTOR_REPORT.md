# Nonconvex–Nonconcave Tangent-Actor Report

## Decision: **CONDITIONAL GO**

The experiment evaluates the unchanged undamped Terminal-Aware Tangent-Actor update with alpha=1 martingale control variates. Checkpoints are selected by an independently replicated, noise-corrected fixed-point consistency score before the numerical reference is loaded.

The current decision remains conditional because:

- accuracy advantage over all frozen baseline endpoints not established.
- d=400 pilot did not achieve a material 10% true-RMSE reduction from a0.

## Frozen selected results

| family | d | lambda | budget | training labels | seed | Actor | Policy Value RMSE | Relative L2 | systematic / MC RMSE | Value-Actor RMSE | Gradient-Actor vector RMSE | MCV variance reduction | ESS |
|---|---:|---:|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| coupled | 100 | 0.2 | pilot | mcv | 0 | a2 | 0.00198365 | 0.00422741 | 0.00177607 / 0.000883417 | 0.0275067 | 0.0630903 | 7.44× | 0.9504 |
| primary | 100 | 0.0 | pilot | mcv | 0 | a2 | 0.00278866 | 0.00556405 | 0.00253653 / 0.00115874 | 0.0283038 | 0.0525457 | 11.99× | 1.0000 |
| primary | 100 | 0.05 | pilot | mcv | 0 | a2 | 0.00141132 | 0.00283456 | 0.00100269 / 0.000993202 | 0.0218316 | 0.0519344 | 14.65× | 0.9984 |
| primary | 100 | 0.2 | main | mcv | 0 | a4 | 0.000684943 | 0.00137593 | 0.000559518 / 0.000395077 | 0.00421344 | 0.007139 | 51.75× | 0.9733 |
| primary | 100 | 0.2 | main | mcv | 1 | a4 | 0.000617239 | 0.00123992 | 0.000476799 / 0.000391979 | 0.00467583 | 0.00684501 | 51.79× | 0.9732 |
| primary | 100 | 0.2 | main | mcv | 2 | a4 | 0.000666913 | 0.00133971 | 0.000538077 / 0.000394013 | 0.0043375 | 0.00680163 | 50.99× | 0.9733 |
| primary | 100 | 0.2 | pilot | mcv | 0 | a2 | 0.00150872 | 0.00303075 | 0.00120383 / 0.000909409 | 0.0259386 | 0.0607912 | 9.76× | 0.9736 |
| primary | 100 | 0.2 | pilot | plain | 0 | a2 | 0.00151712 | 0.00304762 | 0.00117708 / 0.000957148 | 0.0274475 | 0.0710434 | 8.78× | 0.9740 |
| primary | 100 | 0.5 | pilot | mcv | 0 | a2 | 0.00954707 | 0.018995 | 0.00947334 / 0.00118423 | 0.0357353 | 0.0970053 | 2.22× | 0.8513 |
| primary | 400 | 0.2 | pilot | mcv | 0 | a2 | 0.00629103 | 0.0124982 | 0.0060322 / 0.00178595 | 0.079095 | 0.153572 | 2.82× | 0.9781 |

## Actor-iteration trajectories

Arrows list `a0 -> a1 -> ...` in training order. Value and Gradient errors are post-freeze audits; the fixed-point trajectory is the reference-free ranking audit used for selection.

| run | linear-policy mean (diagnostic) | Policy Value RMSE | Gradient-Actor vector RMSE | fixed-point inconsistency |
|---|---|---|---|---|
| coupled_d100_k2_lam0p2_eta1_pilot_mcv_s0 | 0.4657 → 0.4653 → 0.4644 | 0.006456 → 0.00573 → 0.001984 | 0.1792 → 0.1103 → 0.06309 | 2.822 → 1.051 → 0.2341 |
| primary_d100_k2_lam0_eta0_pilot_mcv_s0 | 0.4936 → 0.4931 → 0.4929 | 0.002322 → 0.002414 → 0.002789 | 0.123 → 0.07502 → 0.05255 | 0.7629 → 0.4335 → 0.2652 |
| primary_d100_k2_lam0p05_eta0_pilot_mcv_s0 | 0.4903 → 0.4898 → 0.4893 | 0.00274 → 0.003883 → 0.001411 | 0.1348 → 0.08068 → 0.05193 | 0.8613 → 0.3733 → 0.08492 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s0 | 0.4938 → 0.4918 → 0.4917 → 0.4917 → 0.4917 | 0.006031 → 0.0007646 → 0.0006983 → 0.0006922 → 0.0006849 | 0.1619 → 0.02285 → 0.01172 → 0.008576 → 0.007139 | 1.878 → 0.03702 → 0.008418 → 0.005203 → 0.003538 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s1 | 0.4937 → 0.4918 → 0.4917 → 0.4917 → 0.4917 | 0.005903 → 0.0007173 → 0.0006283 → 0.0006267 → 0.0006172 | 0.1619 → 0.02374 → 0.01138 → 0.008586 → 0.006845 | 2.048 → 0.04733 → 0.01456 → 0.006221 → 0.003001 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s2 | 0.4938 → 0.4918 → 0.4918 → 0.4917 → 0.4917 | 0.006031 → 0.0007314 → 0.0006864 → 0.0006605 → 0.0006669 | 0.1619 → 0.02212 → 0.011 → 0.008293 → 0.006802 | 2.227 → 0.03981 → 0.009683 → 0.004841 → 0.003108 |
| primary_d100_k2_lam0p2_eta0_pilot_mcv_s0 | 0.4937 → 0.4927 → 0.4921 | 0.005779 → 0.00505 → 0.001509 | 0.1619 → 0.102 → 0.06079 | 1.483 → 0.7351 → 0.1451 |
| primary_d100_k2_lam0p2_eta0_pilot_plain_s0 | 0.4938 → 0.4927 → 0.4921 | 0.006345 → 0.005956 → 0.001517 | 0.1619 → 0.1026 → 0.07104 | 2.442 → 1.998 → 0.7264 |
| primary_d100_k2_lam0p5_eta0_pilot_mcv_s0 | 0.5007 → 0.4997 → 0.5007 | 0.02424 → 0.02455 → 0.009547 | 0.2055 → 0.1332 → 0.09701 | 8.067 → 2.624 → 0.8083 |
| primary_d400_k2_lam0p2_eta0_pilot_mcv_s0 | 0.4986 → 0.4986 → 0.4985 | 0.006416 → 0.006354 → 0.006291 | 0.1622 → 0.1581 → 0.1536 | 2.55 → 2.482 → 2.41 |

## MCV label-noise control

For the identical frozen `a0`, state bank, path count, and Euler grid, independent label replications give:

| label | plain replication RMSE | MCV replication RMSE | RMSE reduction |
|---|---:|---:|---:|
| Value | 0.0131571 | 0.00728776 | 1.81× |
| Gradient vector | 0.027432 | 0.0150871 | 1.82× |

The same-budget frozen training endpoints separate label quality from final Value accuracy:

| training labels | selected Actor | Value RMSE | Gradient-vector RMSE | fixed-point inconsistency |
|---|---|---:|---:|---:|
| MCV | a2 | 0.00150872 | 0.0607912 | 0.145066 |
| plain FK | a2 | 0.00151712 | 0.0710434 | 0.726448 |

MCV lowers independent Value/Gradient label-replication RMSE by 1.81×/1.82×. At the selected pilot endpoint, its Gradient RMSE is lower by a factor of 1.17 and its fixed-point score is lower by 5.01×, while the final Value RMSE ratio is only 1.006×. Thus this budget establishes a training-quality benefit, but not a material Value-accuracy advantage.

## Accuracy versus curvature strength

The canonical seed-0 pilot sweep changes only lambda and retains the same full-ambient architecture and budget.

| lambda | selected Actor | a0 Value RMSE | selected Value RMSE | true reduction from a0 | Gradient-vector RMSE | MCV variance reduction | ESS |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | a2 | 0.00232208 | 0.00278866 | -20.1% | 0.0525457 | 11.99× | 1.0000 |
| 0.05 | a2 | 0.00273961 | 0.00141132 | 48.5% | 0.0519344 | 14.65× | 0.9984 |
| 0.2 | a2 | 0.00577928 | 0.00150872 | 73.9% | 0.0607912 | 9.76× | 0.9736 |
| 0.5 | a2 | 0.0242392 | 0.00954707 | 60.6% | 0.0970053 | 2.22× | 0.8513 |

The lambda=0 control is a retained negative signal: fixed-point consistency and Gradient error improve, yet the observed selected Value RMSE worsens. Its RMSE gap is smaller than the reported MC-error scale, so this finite-bank ordering is not claimed as statistically resolved. At lambda=0.05 and 0.5 the first update can also worsen Value before the second update improves it. These observations rule out a monotone-improvement interpretation.

## Reference-free selection audit

The post-freeze true-best column is an audit only; it was unavailable during checkpoint selection.

| run | selected | post-freeze true best | match | selected excess RMSE | within-run rank correlation |
|---|---|---|---:|---:|---:|
| coupled_d100_k2_lam0p2_eta1_pilot_mcv_s0 | a2 | a2 | yes | 0 | 1.000 |
| primary_d100_k2_lam0_eta0_pilot_mcv_s0 | a2 | a0 | no | 0.000466579 | -1.000 |
| primary_d100_k2_lam0p05_eta0_pilot_mcv_s0 | a2 | a2 | yes | 0 | 0.500 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s0 | a4 | a4 | yes | 0 | 1.000 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s1 | a4 | a4 | yes | 0 | 1.000 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s2 | a4 | a3 | no | 6.40634e-06 | 0.900 |
| primary_d100_k2_lam0p2_eta0_pilot_mcv_s0 | a2 | a2 | yes | 0 | 1.000 |
| primary_d100_k2_lam0p2_eta0_pilot_plain_s0 | a2 | a2 | yes | 0 | 1.000 |
| primary_d100_k2_lam0p5_eta0_pilot_mcv_s0 | a2 | a2 | yes | 0 | 0.500 |
| primary_d400_k2_lam0p2_eta0_pilot_mcv_s0 | a2 | a2 | yes | 0 | 1.000 |

Across canonical MCV runs, the frozen selector matches the post-freeze best candidate in 77.8% of runs. This supports the score as a useful diagnostic in the tested regime, not as an error bound.

## Euler-step sensitivity of frozen selections

Each entry is an independently seeded MCV evaluation with the same number of states and paths. The spread therefore includes both discretization and finite-MC variation.

| run | Actor | 32-step RMSE | 64-step RMSE | 128-step RMSE |
|---|---|---:|---:|---:|
| coupled_d100_k2_lam0p2_eta1_pilot_mcv_s0 | a2 | 0.00260464 | 0.00195547 | 0.00183591 |
| primary_d100_k2_lam0_eta0_pilot_mcv_s0 | a2 | 0.00246529 | 0.00250321 | 0.00309355 |
| primary_d100_k2_lam0p05_eta0_pilot_mcv_s0 | a2 | 0.00179946 | 0.00150893 | 0.0015301 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s0 | a4 | 0.00111405 | 0.000670375 | 0.000378794 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s1 | a4 | 0.00111911 | 0.000621932 | 0.000366243 |
| primary_d100_k2_lam0p2_eta0_main_mcv_s2 | a4 | 0.00118006 | 0.000623245 | 0.000360269 |
| primary_d100_k2_lam0p2_eta0_pilot_mcv_s0 | a2 | 0.00189639 | 0.00173641 | 0.00121496 |
| primary_d100_k2_lam0p2_eta0_pilot_plain_s0 | a2 | 0.00193218 | 0.00165894 | 0.00166688 |
| primary_d100_k2_lam0p5_eta0_pilot_mcv_s0 | a2 | 0.0102023 | 0.00863508 | 0.00848722 |
| primary_d400_k2_lam0p2_eta0_pilot_mcv_s0 | a2 | 0.00685241 | 0.0061489 | 0.00639588 |

## Frozen endpoint comparison

This table uses the reference-free selected main Actor, the predeclared highest Full-History endpoint `n=3,M=64`, and the PINN final step. No baseline row is selected by test RMSE.

| method | frozen configuration | finite / total seeds | Value RMSE mean ± sd | mean total seconds |
|---|---|---:|---:|---:|
| Tangent-Actor + MCV | reference-free FP-selected Actor | 3 / 3 | 0.000656365 ± 3.51e-05 | 2910.8 |

## Ambient-dimension scaling diagnosis

| d | params / train point | a1 selected step | a1 label-replication Gradient RMSE | selected Value RMSE | selected active / inactive Gradient RMSE |
|---:|---:|---:|---:|---:|---:|
| 100 | 254.7 | 250 | 0.015087 | 0.00150872 | 0.060729 / 0.0027515 |
| 400 | 367.3 | 1 | 0.017067 | 0.00629103 | 0.15318 / 0.010965 |

The d=400 Gradient-label replication noise is only 1.13× the d=100 value, while the selected error remains almost entirely in the two active directions and the inactive-output energy is negligible. Training after the first step nearly interpolates the training bank but worsens validation. The observed failure is therefore most consistent with full-ambient sample complexity/generalization under many irrelevant coordinates, rather than exponential-weight instability, MC-label noise alone, or spurious inactive Gradient output.

The predeclared expansion rule therefore stopped before a coupled d=400 run or a costly d=400 budget increase. This is a retained scaling limitation, not a missing positive run.

The 10% d=400 material-improvement threshold is used only for the final research GO/CONDITIONAL-GO judgment. It was introduced after freezing all training and checkpoint-selection rules and did not influence any Actor choice.

## Interpretation

Completed 10 training runs and 10 post-freeze evaluations. The rank correlation between the reference-free consistency score and post-freeze true Value error is 0.9552123552123553.

No monotone policy-improvement theorem or nonnegative Bregman-regret identity is used for these drivers. Linear-policy values are diagnostics only. The experiment supports claims only for smooth bounded-curvature drivers in the tested locally stable regime.

Every completed undamped run remained finite and above the ESS safety floor, so the predeclared trigger for a separate damping ablation was not met. No damped result is mixed into the primary method.

The Actor always receives the full ambient state. Active coordinates and reference gradients are used only in the post-freeze audit. These are ambient 100D/400D tests with a two-dimensional reference structure, not evidence for arbitrary effective dimension.

The independent comparison uses raw Full-History MLP and a plain soft-terminal PINN. Native scalar Deep BSDE does not estimate the required random-time solution function, the available Deep Picard adaptation embeds our terminal-aware Actor, and the available DPI implementation is tied to a control-form HJB; they are excluded rather than modified into favorable but method-mismatched baselines.

## Answers to the research questions

1. The undamped generic tangent update is accurate on the primary and coupled d=100 nonconvex drivers, including a genuinely nonseparable `(u,z)` driver. This is empirical local-regime evidence, not global convergence.
2. Iteration is numerically stable in every completed run, but true Value error is not monotone and additional updates matter differently across curvature strengths.
3. The independent fixed-point score usually ranks candidates well, but the lambda=0 control has an observed finite-bank selection mismatch whose gap is not resolved beyond MC error; the score cannot replace a theorem or calibrated validation bound.
4. MCV reduces label-replication noise and becomes a strong inference control variate once the Actor is accurate. With weak Actors it can increase inference variance.
5. Nonseparable coupling remains solvable at d=100, with lower variance-reduction benefit and higher error than the primary high-budget run.
6. The d=400 pilot is stable but gains only marginal accuracy; current evidence points to full-ambient statistical generalization as the scaling bottleneck.
7. The resulting paper claim should remain restricted to smooth bounded-curvature, locally stable semilinear drivers with careful reference-free selection.

## Artifacts

- [Actor iteration behavior](figures/actor_iteration_behavior.png)
- [Fixed-point consistency versus true error](figures/fixed_point_vs_true_error.png)
- [Curvature-region coverage](figures/curvature_region_coverage.png)
- [Accuracy versus curvature strength](figures/curvature_strength_accuracy.png)
- [MCV variance reduction](figures/mcv_variance_reduction.png)
- [Baseline accuracy and cost](figures/baseline_accuracy_cost.png)
- `ACTOR_CONSISTENCY_DIAGNOSTICS.csv`
- `ACTOR_ITERATION_RESULTS.csv`
- `SELECTION_RELIABILITY.csv`
- `EULER_SENSITIVITY.csv`
- `MCV_VARIANCE_RESULTS.csv`
- `BASELINE_COMPARISON.csv`
- `COMPUTE_ACCOUNTING.csv`
- `REPRODUCIBILITY_MANIFEST.json`
