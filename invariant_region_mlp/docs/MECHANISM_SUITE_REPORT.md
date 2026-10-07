# Mechanism benchmark suite: noise entering a nonlinear generator

## Outcome

| PDE | C1 | C2 | C3 | C3 joint-cell fraction | certificate-prior dominated | centre better fraction | shrink-centre better fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | PASS | PASS | FAIL | 0.5000 | no | 0.0000 | 0.1250 |
| P2_a4 | PASS | FAIL | FAIL | 0.5000 | yes | 0.0000 | 0.3750 |
| P2_a8 | PASS | FAIL | PASS | 0.7500 | no | 0.0000 | 0.2500 |
| P3_rho1 | PASS | FAIL | FAIL | 0.6250 | yes | 0.0000 | 0.3750 |
| P3_rho2 | PASS | FAIL | FAIL | 0.6250 | yes | 0.0000 | 0.3750 |
| P4 | PASS | PASS | FAIL | 0.0000 | yes | 0.0000 | 0.7500 |

Across the 6 admitted problem variants, C1 is supported for 6/6, C2 for 2/6, and C3 for 1/6. The suite therefore gives broad support to nonlinear-generator noise damage, but not to universal half-error removal by box retraction or to the claim that the box advantage generally requires the data rather than a tuned certificate prior.

The table applies the pre-registered thresholds mechanically; a failed criterion is not softened by qualitative trends. Lower skill is better, and a constant predictor has skill 1.

## Protocol and provenance

All computations used float64, 1,200 fixed points with a fixed 20% validation split, paired tree seeds `SeedSequence([base_seed,d,n,M,rep,chunk_index])`, and identical chunk sizes within each study. E0 uses 8 at d=20 and 4 otherwise; E1/E3/E5/E6 use 16 at d=20 and 4 otherwise; E2 uses 16/8/4 at d=20, d=50–100, and d=200 respectively. The unchanged `FullHistoryMLP` from `experiments/active_vb_high_budget/vb_mlp_methods.py` supplies the recursion; only the state immediately handed to f is transformed, and the root is never clipped. Base seed: 20261006.

The gate snapshot records code commit `514fe7a83527ed96e9090fe7a7f3d5a0a25ba939` on branch `ir-mlp-mechanism-suite`, Python 3.12.11, NumPy 2.3.4, SciPy 1.18.1, and 16 logical CPUs. The unchanged recursion source has SHA-256 `f4581babdd5ceed8354337ee7783112eb6b2590facc2a0115e7c02ccd50990b8`.

P2, P3, and N4 use the inherited VB test geometry with 1,000 interior and 200 face-boundary points; P1, P4, and N3 use t uniform on [0,T) and x uniform on [-1,1]^d. E0 containment independently samples 100,000 points from the corresponding distribution, including the VB boundary stratum.

### Problem matrix

| ID | parameters | d | certificate |
| --- | --- | --- | --- |
| P1 | ridge LSE; sigma=sqrt(2), T=0.25, kappa=1 | 20,50,100,200 | PDE-derived segment, box, ball, span |
| P2_a4 | VB-a; a=4, sigma=0.5, T=0.5 | 20,50,100 | u/sign PDE-derived; z cap solution-informed |
| P2_a8 | VB-a; a=8, sigma=0.5, T=0.5 | 20,50,100 | u/sign PDE-derived; z cap solution-informed |
| P3_rho1 | Burgers-Fisher; a=4, rho=1, sigma=0.5, T=0.5 | 20,50 | u/sign PDE-derived; z cap solution-informed |
| P3_rho2 | Burgers-Fisher; a=4, rho=2, sigma=0.5, T=0.5 | 20,50 | u/sign PDE-derived; z cap solution-informed |
| P4 | norm HJB; beta=2, lambda_f=1, T=0.5 | 20,50,100 | PDE-derived segment, box, ball, span |
| N3 | 5-direction LSE; strength=0.5, T=0.5 | 100 | PDE-derived convex hull, box, ball, span |
| N4 | published VB; a=1, sigma=sqrt(2), T=0.5 | 20,50,100 | u/sign PDE-derived; z cap solution-informed |

### Study matrix

| study | design | purpose |
| --- | --- | --- |
| E0 | all candidates; (3,6),(4,6); 3 reps | hard admission gates |
| E1 | 4 (n,M) cells; all admitted d; 10 reps | oracle/channel decomposition |
| E2 | (4,3),(3,6); six doses; 10 reps | one-shot causal dose response |
| E3 | 4 cells; two smallest d; 10 reps; validation-only tuning | geometry and suppression controls |
| E4 | derived from untuned/tuned E3 geometry rows | certificate component ablation |
| E5 | N3/N4; (3,6),(4,6); 10 reps | negative-control predictions |
| E6 | P1; (4,6); four d; 10 reps | box/oracle dimension stability |

The supplied attachment did not contain the referenced `reference_code/` directory. P1–P3 and N4 were implemented directly from the equations in the protocol; N3 is explicitly marked as a canonical mathematical reconstruction. This missing bundle is a reproducibility limitation, not hidden.
During E0 optimization, 11 preliminary P4 artifacts using FITPACK's direct pointwise derivative evaluator were preserved under `raw/superseded/`, and 1 paired comparison artifact under `raw/derivative_validation/`; all are excluded from every result table. The accepted artifacts use the audited cached derivative spline (maximum discrepancy from the direct spline derivative is recorded in G1 and the [paired runtime audit](../results/mechanism_suite/reference_cache/derivative_evaluator_audit.json)).

P4 used beta=2, lambda_f=1, T=0.5. Lambda_f=1 was the stronger pre-listed suggestion and was fixed before inspecting any MLP result; the lambda_f=0.5 alternatives were not exhaustively screened. The suggested T=0.25 setting was screened out because G2=0.1206<0.15; T=0.5 gave G2=0.2018. Its six-grid monotone upwind reference has common-node Richardson disagreement 8.189e-08, stored on 513×9601 points. The unextrapolated finest-pair difference is 2.376e-05; because the monotone scheme is first order, the formal reference and G1 comparison use successive common-node third-Richardson estimates, and both numbers are reported.
For P4 d=20, the spline-differential residual is 3.350e-07, the worst value in the three-step independent fourth-order FD audit is 6.009e-08, and the central gradient-check error is 8.884e-10.

P4 had no dimension list in the supplied protocol, so E0/E1/E2 use d={20,50,100}; the compute-limited E3 panel uses d={20,50}. This choice was fixed before MLP results were inspected.

## E0: hard gates

| PDE | G1 | G2 | G3 | G4 | G5 | G6 | verdict | failed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N3 | PASS | FAIL | PASS | PASS | PASS | PASS | exclude | G2 |
| N4 | PASS | PASS | FAIL | PASS | PASS | PASS | exclude | G3 |
| P1 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| P2_a4 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| P2_a8 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| P3_rho1 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| P3_rho2 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| P4 | PASS | PASS | PASS | PASS | PASS | PASS | admit | none |
| N1 | historical—not rerun | FAIL | historical—not rerun | historical—not rerun | historical—not rerun | historical—not rerun | historical negative control | G2 |
| N2 | historical—not rerun | FAIL | historical—not rerun | historical—not rerun | historical—not rerun | historical—not rerun | historical negative control | G2 |
| N5 | historical—not rerun | historical—not rerun | historical—not rerun | FAIL | historical—not rerun | historical—not rerun | historical negative control | G4 |

N3 was predicted in advance to fail only G2; N4 was predicted to fail only G3. Only primary PDEs passing all six gates enter C1–C3.
[Full gate diagnostics, containment counts, certificate labels, and reference audits](../results/mechanism_suite/gates.json)

## Pre-registered acceptance criteria (verbatim)

- **C1 supported** for a PDE if `D >= 0.5` in at least 3 of 4 headline configurations at each admitted d, and the E2 dose curve for raw input increases monotonically in s.
- **C2 supported** if `box/oracle <= 0.5 * raw/oracle` (equivalently box removes at least half of the excess error) in at least 3 of 4 configurations.
- **C3 supported** if box wins at least 7/10 paired reps against each of: best tuned ->0 control, best tuned illegal set, centre, best tuned shrink_centre, in at least 75% of E3 cells.
- Any PDE where centre or shrink_centre beats box in more than 25% of cells must be reported as "certificate-prior dominated" for that regime.
- Negative controls must behave as predicted; if not, report and explain.

For C2, the written ratio inequality is algebraically stricter than the parenthetical ‘half of excess error’ statement. Verdicts use the explicit pre-registered inequality `box/oracle <= 0.5 * raw/oracle`; `Gc` is reported separately so the alternative excess-gap reading remains visible.

## E1: oracle decomposition

Per-cell entries below are means over 10 paired repetitions. `D=(raw-oracle)/raw`; `Gc=(raw-box)/(raw-oracle)`. The CSV also contains every per-repetition skill, median skill, repetition-averaged skill, gradient relative L2, work counters, violations, activations, and non-finite counts.

| pde | d | n | M | raw | raw_median | raw_repavg | box | box_repavg | segment | oracle_state | oracle_state_repavg | oracle_z | oracle_u | raw_gradrel | box_gradrel | oracle_state_gradrel | damage_share_D | gap_closed_Gc | box_over_oracle |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | 20 | 3 | 6 | 57.7891 | 56.1704 | 50.2406 | 0.1775 | 0.1312 | 0.2239 | 0.0714 | 0.0221 | 0.0714 | 57.7891 | 362.8944 | 1.1439 | 0.9126 | 0.9988 | 0.9982 | 2.4864 |
| P1 | 20 | 3 | 10 | 13.0137 | 12.9776 | 12.4226 | 0.1401 | 0.1154 | 0.1495 | 0.0385 | 0.0118 | 0.0385 | 13.0137 | 59.4850 | 0.7659 | 0.5409 | 0.9970 | 0.9922 | 3.6357 |
| P1 | 20 | 4 | 3 | 5.398e+07 | 4.545e+07 | 2.588e+07 | 0.2389 | 0.1491 | 0.3868 | 0.0988 | 0.0300 | 0.0988 | 5.398e+07 | 8.510e+08 | 1.6335 | 1.0997 | 1.0000 | 1.0000 | 2.4184 |
| P1 | 20 | 4 | 6 | 1.038e+05 | 9.142e+04 | 7.190e+04 | 0.1429 | 0.1093 | 0.1948 | 0.0291 | 0.0092 | 0.0291 | 1.038e+05 | 1.162e+06 | 0.8326 | 0.3730 | 1.0000 | 1.0000 | 4.9106 |
| P1 | 50 | 3 | 6 | 544.4517 | 544.1792 | 495.5439 | 0.1787 | 0.1504 | 0.2205 | 0.0702 | 0.0233 | 0.0702 | 544.4517 | 5179.8731 | 1.5937 | 1.4699 | 0.9999 | 0.9998 | 2.5460 |
| P1 | 50 | 3 | 10 | 111.9600 | 111.5037 | 108.1736 | 0.1493 | 0.1354 | 0.1485 | 0.0379 | 0.0123 | 0.0379 | 111.9600 | 779.4540 | 1.0385 | 0.8687 | 0.9997 | 0.9990 | 3.9420 |
| P1 | 50 | 4 | 3 | 1.639e+10 | 1.451e+10 | 6.869e+09 | 0.2155 | 0.1599 | 0.3804 | 0.0984 | 0.0306 | 0.0984 | 1.639e+10 | 3.271e+11 | 2.1238 | 1.7618 | 1.0000 | 1.0000 | 2.1909 |
| P1 | 50 | 4 | 6 | 1.804e+07 | 1.653e+07 | 1.429e+07 | 0.1409 | 0.1225 | 0.1945 | 0.0287 | 0.0094 | 0.0287 | 1.804e+07 | 2.878e+08 | 1.0378 | 0.6003 | 1.0000 | 1.0000 | 4.9111 |
| P1 | 100 | 3 | 6 | 3476.5716 | 3479.8847 | 3195.9656 | 0.1817 | 0.1559 | 0.2241 | 0.0720 | 0.0213 | 0.0720 | 3476.5716 | 4.448e+04 | 2.1013 | 2.0107 | 1.0000 | 1.0000 | 2.5230 |
| P1 | 100 | 3 | 10 | 686.3247 | 686.1092 | 663.5225 | 0.1529 | 0.1416 | 0.1465 | 0.0392 | 0.0127 | 0.0392 | 686.3247 | 6356.6407 | 1.3486 | 1.1862 | 0.9999 | 0.9998 | 3.8979 |
| P1 | 100 | 4 | 3 | 9.118e+11 | 7.095e+11 | 4.761e+11 | 0.2093 | 0.1588 | 0.3880 | 0.1017 | 0.0313 | 0.1017 | 9.118e+11 | 2.807e+13 | 2.7619 | 2.4087 | 1.0000 | 1.0000 | 2.0587 |
| P1 | 100 | 4 | 6 | 1.212e+09 | 1.191e+09 | 1.053e+09 | 0.1445 | 0.1297 | 0.1938 | 0.0294 | 0.0095 | 0.0294 | 1.212e+09 | 2.528e+10 | 1.3095 | 0.8202 | 1.0000 | 1.0000 | 4.9162 |
| P1 | 200 | 3 | 6 | 2.430e+04 | 2.435e+04 | 2.266e+04 | 0.1822 | 0.1651 | 0.2122 | 0.0680 | 0.0205 | 0.0680 | 2.430e+04 | 4.600e+05 | 2.7858 | 2.8865 | 1.0000 | 1.0000 | 2.6782 |
| P1 | 200 | 3 | 10 | 4601.0372 | 4581.7137 | 4465.1033 | 0.1608 | 0.1536 | 0.1422 | 0.0378 | 0.0118 | 0.0378 | 4601.0372 | 6.460e+04 | 1.7468 | 1.7060 | 1.0000 | 1.0000 | 4.2562 |
| P1 | 200 | 4 | 3 | 6.652e+13 | 5.905e+13 | 3.675e+13 | 0.2038 | 0.1696 | 0.3768 | 0.0952 | 0.0304 | 0.0952 | 6.652e+13 | 3.097e+15 | 3.5414 | 3.4629 | 1.0000 | 1.0000 | 2.1396 |
| P1 | 200 | 4 | 6 | 1.127e+11 | 1.118e+11 | 9.967e+10 | 0.1463 | 0.1375 | 0.1869 | 0.0283 | 0.0095 | 0.0283 | 1.127e+11 | 3.793e+12 | 1.5684 | 1.1802 | 1.0000 | 1.0000 | 5.1670 |
| P2_a4 | 20 | 3 | 6 | 0.2472 | 0.2469 | 0.0826 | 0.1894 | 0.1682 | — | 0.0821 | 0.0255 | 0.0920 | 0.3040 | 2.3988 | 1.0350 | 1.0881 | 0.6677 | 0.3499 | 2.3061 |
| P2_a4 | 20 | 3 | 10 | 0.1444 | 0.1439 | 0.0579 | 0.1525 | 0.1423 | — | 0.0452 | 0.0136 | 0.0540 | 0.1718 | 1.4185 | 0.6848 | 0.6273 | 0.6867 | -0.0823 | 3.3724 |
| P2_a4 | 20 | 4 | 3 | 0.6280 | 0.6267 | 0.2002 | 0.2149 | 0.1699 | — | 0.1144 | 0.0359 | 0.1331 | 0.7539 | 6.6701 | 1.4335 | 1.3842 | 0.8179 | 0.8042 | 1.8794 |
| P2_a4 | 20 | 4 | 6 | 0.2270 | 0.2282 | 0.0748 | 0.1278 | 0.1115 | — | 0.0340 | 0.0110 | 0.0423 | 0.2663 | 2.3516 | 0.6974 | 0.4467 | 0.8501 | 0.5140 | 3.7562 |
| P2_a4 | 50 | 3 | 6 | 0.4767 | 0.4809 | 0.2127 | 0.3536 | 0.3293 | — | 0.1290 | 0.0407 | 0.1769 | 0.6240 | 6.6636 | 2.4413 | 2.4008 | 0.7293 | 0.3542 | 2.7402 |
| P2_a4 | 50 | 3 | 10 | 0.3134 | 0.3131 | 0.1918 | 0.3220 | 0.3115 | — | 0.0734 | 0.0240 | 0.1264 | 0.3905 | 4.0770 | 1.6818 | 1.4177 | 0.7658 | -0.0359 | 4.3879 |
| P2_a4 | 50 | 4 | 3 | 2.0038 | 1.9922 | 0.6555 | 0.3761 | 0.3269 | — | 0.1569 | 0.0486 | 0.2436 | 2.1564 | 34.4993 | 3.2017 | 2.9007 | 0.9217 | 0.8813 | 2.3976 |
| P2_a4 | 50 | 4 | 6 | 0.6455 | 0.6476 | 0.2668 | 0.3072 | 0.2924 | — | 0.0528 | 0.0149 | 0.0932 | 0.7529 | 10.0664 | 1.6264 | 0.9815 | 0.9182 | 0.5707 | 5.8219 |
| P2_a4 | 100 | 3 | 6 | 0.7861 | 0.7836 | 0.4845 | 0.5659 | 0.5358 | — | 0.1838 | 0.0594 | 0.3446 | 1.1263 | 13.0464 | 4.5912 | 4.1845 | 0.7662 | 0.3656 | 3.0797 |
| P2_a4 | 100 | 3 | 10 | 0.5703 | 0.5707 | 0.4293 | 0.5372 | 0.5231 | — | 0.1082 | 0.0328 | 0.2932 | 0.7661 | 8.1225 | 3.2672 | 2.4864 | 0.8102 | 0.0715 | 4.9640 |
| P2_a4 | 100 | 4 | 3 | 8.8679 | 8.1794 | 3.1155 | 0.6070 | 0.5416 | — | 0.2192 | 0.0716 | 0.4158 | 4.7182 | 199.5128 | 6.2940 | 4.9475 | 0.9753 | 0.9552 | 2.7690 |
| P2_a4 | 100 | 4 | 6 | 1.7323 | 1.7180 | 0.9036 | 0.5323 | 0.5155 | — | 0.0764 | 0.0241 | 0.1945 | 1.6574 | 37.5498 | 3.2378 | 1.7035 | 0.9559 | 0.7247 | 6.9701 |
| P2_a8 | 20 | 3 | 6 | 0.3007 | 0.3041 | 0.1024 | 0.2085 | 0.1851 | — | 0.0859 | 0.0263 | 0.0959 | 0.3739 | 3.3160 | 1.3853 | 1.3294 | 0.7144 | 0.4292 | 2.4278 |
| P2_a8 | 20 | 3 | 10 | 0.1757 | 0.1744 | 0.0691 | 0.1686 | 0.1567 | — | 0.0481 | 0.0146 | 0.0572 | 0.2121 | 1.9473 | 0.9058 | 0.7607 | 0.7261 | 0.0557 | 3.5031 |
| P2_a8 | 20 | 4 | 3 | 0.7778 | 0.7779 | 0.2458 | 0.2387 | 0.1909 | — | 0.1187 | 0.0367 | 0.1382 | 0.9507 | 9.5822 | 1.9280 | 1.7152 | 0.8474 | 0.8180 | 2.0105 |
| P2_a8 | 20 | 4 | 6 | 0.2754 | 0.2781 | 0.0895 | 0.1470 | 0.1285 | — | 0.0359 | 0.0113 | 0.0443 | 0.3281 | 3.3266 | 0.9321 | 0.5458 | 0.8697 | 0.5360 | 4.0967 |
| P2_a8 | 50 | 3 | 6 | 0.5721 | 0.5756 | 0.2565 | 0.4288 | 0.4023 | — | 0.1444 | 0.0449 | 0.1973 | 0.7561 | 8.2551 | 3.0222 | 2.7027 | 0.7476 | 0.3350 | 2.9700 |
| P2_a8 | 50 | 3 | 10 | 0.3705 | 0.3705 | 0.2235 | 0.3901 | 0.3785 | — | 0.0827 | 0.0265 | 0.1406 | 0.4691 | 5.0172 | 2.0639 | 1.5927 | 0.7769 | -0.0682 | 4.7190 |
| P2_a8 | 50 | 4 | 3 | 2.4594 | 2.4744 | 0.8061 | 0.4672 | 0.4153 | — | 0.1738 | 0.0535 | 0.2694 | 2.6467 | 43.8958 | 3.9404 | 3.2943 | 0.9293 | 0.8716 | 2.6879 |
| P2_a8 | 50 | 4 | 6 | 0.7726 | 0.7726 | 0.3153 | 0.3849 | 0.3691 | — | 0.0590 | 0.0167 | 0.1028 | 0.9149 | 12.4775 | 1.9765 | 1.1065 | 0.9237 | 0.5433 | 6.5287 |
| P2_a8 | 100 | 3 | 6 | 0.8951 | 0.8890 | 0.5446 | 0.6527 | 0.6197 | — | 0.2009 | 0.0643 | 0.3740 | 1.2922 | 15.9937 | 5.6797 | 4.7212 | 0.7756 | 0.3491 | 3.2497 |
| P2_a8 | 100 | 3 | 10 | 0.6513 | 0.6476 | 0.4894 | 0.6180 | 0.6026 | — | 0.1190 | 0.0367 | 0.3169 | 0.8811 | 9.9411 | 4.0199 | 2.7938 | 0.8172 | 0.0627 | 5.1913 |
| P2_a8 | 100 | 4 | 3 | 9.7218 | 9.1189 | 3.4604 | 0.7167 | 0.6484 | — | 0.2395 | 0.0791 | 0.4485 | 5.5308 | 228.8916 | 7.6374 | 5.6028 | 0.9754 | 0.9497 | 2.9926 |
| P2_a8 | 100 | 4 | 6 | 1.9415 | 1.9360 | 1.0103 | 0.6314 | 0.6138 | — | 0.0832 | 0.0261 | 0.2112 | 1.9081 | 45.0949 | 3.9133 | 1.9219 | 0.9571 | 0.7050 | 7.5888 |
| P3_rho1 | 20 | 3 | 6 | 0.2738 | 0.2731 | 0.0960 | 0.2010 | 0.1763 | — | 0.0974 | 0.0303 | 0.1102 | 0.3291 | 2.8292 | 1.2720 | 1.3945 | 0.6442 | 0.4128 | 2.0630 |
| P3_rho1 | 20 | 3 | 10 | 0.1652 | 0.1648 | 0.0728 | 0.1663 | 0.1548 | — | 0.0549 | 0.0164 | 0.0675 | 0.1882 | 1.7059 | 0.8182 | 0.8136 | 0.6675 | -0.0096 | 3.0268 |
| P3_rho1 | 20 | 4 | 3 | 0.6965 | 0.6926 | 0.2195 | 0.2244 | 0.1738 | — | 0.1292 | 0.0406 | 0.1501 | 0.8129 | 7.7757 | 1.6951 | 1.7324 | 0.8144 | 0.8322 | 1.7364 |
| P3_rho1 | 20 | 4 | 6 | 0.2630 | 0.2636 | 0.0883 | 0.1380 | 0.1216 | — | 0.0402 | 0.0131 | 0.0498 | 0.2981 | 2.8185 | 0.7657 | 0.5724 | 0.8473 | 0.5610 | 3.4348 |
| P3_rho1 | 50 | 3 | 6 | 0.5112 | 0.5136 | 0.2385 | 0.3558 | 0.3297 | — | 0.1403 | 0.0446 | 0.1958 | 0.6483 | 7.2439 | 2.6314 | 2.6892 | 0.7256 | 0.4189 | 2.5361 |
| P3_rho1 | 50 | 3 | 10 | 0.3423 | 0.3427 | 0.2171 | 0.3250 | 0.3137 | — | 0.0798 | 0.0260 | 0.1413 | 0.4043 | 4.4901 | 1.7878 | 1.5912 | 0.7668 | 0.0659 | 4.0712 |
| P3_rho1 | 50 | 4 | 3 | 2.0701 | 2.1031 | 0.6661 | 0.3723 | 0.3181 | — | 0.1690 | 0.0523 | 0.2551 | 2.2206 | 37.1055 | 3.4042 | 3.2333 | 0.9184 | 0.8931 | 2.2031 |
| P3_rho1 | 50 | 4 | 6 | 0.6927 | 0.6988 | 0.2834 | 0.3006 | 0.2854 | — | 0.0574 | 0.0161 | 0.1008 | 0.7841 | 11.0447 | 1.6587 | 1.0994 | 0.9171 | 0.6172 | 5.2364 |
| P3_rho2 | 20 | 3 | 6 | 0.3034 | 0.3024 | 0.1146 | 0.2163 | 0.1852 | — | 0.1139 | 0.0356 | 0.1370 | 0.3531 | 3.3151 | 1.5951 | 1.7076 | 0.6246 | 0.4597 | 1.8989 |
| P3_rho2 | 20 | 3 | 10 | 0.1903 | 0.1908 | 0.0932 | 0.1805 | 0.1656 | — | 0.0651 | 0.0193 | 0.0892 | 0.2037 | 2.0435 | 1.0308 | 1.0031 | 0.6577 | 0.0786 | 2.7707 |
| P3_rho2 | 20 | 4 | 3 | 0.7672 | 0.7656 | 0.2393 | 0.2384 | 0.1763 | — | 0.1454 | 0.0458 | 0.1797 | 0.8711 | 8.9385 | 2.0598 | 2.0913 | 0.8104 | 0.8504 | 1.6394 |
| P3_rho2 | 20 | 4 | 6 | 0.3009 | 0.3017 | 0.1022 | 0.1477 | 0.1296 | — | 0.0467 | 0.0152 | 0.0632 | 0.3290 | 3.3284 | 0.8902 | 0.7009 | 0.8450 | 0.6026 | 3.1660 |
| P3_rho2 | 50 | 3 | 6 | 0.5457 | 0.5476 | 0.2669 | 0.3631 | 0.3335 | — | 0.1515 | 0.0484 | 0.2231 | 0.6691 | 7.8281 | 2.9000 | 2.9783 | 0.7224 | 0.4632 | 2.3970 |
| P3_rho2 | 50 | 3 | 10 | 0.3732 | 0.3751 | 0.2453 | 0.3317 | 0.3188 | — | 0.0864 | 0.0280 | 0.1645 | 0.4158 | 4.9151 | 1.9639 | 1.7651 | 0.7685 | 0.1446 | 3.8397 |
| P3_rho2 | 50 | 4 | 3 | 2.1511 | 2.1507 | 0.6835 | 0.3741 | 0.3110 | — | 0.1815 | 0.0561 | 0.2788 | 2.2719 | 40.2214 | 3.7304 | 3.5672 | 0.9156 | 0.9022 | 2.0618 |
| P3_rho2 | 50 | 4 | 6 | 0.7328 | 0.7446 | 0.2960 | 0.2952 | 0.2785 | — | 0.0621 | 0.0174 | 0.1141 | 0.8096 | 12.0059 | 1.7661 | 1.2178 | 0.9153 | 0.6524 | 4.7562 |
| P4 | 20 | 3 | 6 | 2.3888 | 2.4022 | 2.2741 | 0.4795 | 0.4096 | 0.3666 | 0.2265 | 0.0703 | 0.2265 | 2.3888 | 3.8989 | 1.2028 | 0.8417 | 0.9052 | 0.8830 | 2.1174 |
| P4 | 20 | 3 | 10 | 1.3421 | 1.3445 | 1.2934 | 0.3545 | 0.3244 | 0.2181 | 0.1186 | 0.0385 | 0.1186 | 1.3421 | 1.7260 | 0.7238 | 0.4629 | 0.9116 | 0.8072 | 2.9890 |
| P4 | 20 | 4 | 3 | 5.7989 | 5.7916 | 5.3874 | 0.5614 | 0.4170 | 0.5963 | 0.3316 | 0.1058 | 0.3316 | 5.7989 | 12.7253 | 1.6598 | 1.1527 | 0.9428 | 0.9580 | 1.6932 |
| P4 | 20 | 4 | 6 | 1.5720 | 1.5607 | 1.4911 | 0.2888 | 0.2462 | 0.2325 | 0.0923 | 0.0290 | 0.0923 | 1.5720 | 2.8735 | 0.7647 | 0.3429 | 0.9413 | 0.8672 | 3.1277 |
| P4 | 50 | 3 | 6 | 6.0909 | 6.0907 | 5.9617 | 0.5837 | 0.5321 | 0.3710 | 0.2267 | 0.0693 | 0.2267 | 6.0909 | 15.4012 | 1.9752 | 1.3563 | 0.9628 | 0.9391 | 2.5749 |
| P4 | 50 | 3 | 10 | 3.4878 | 3.4901 | 3.4415 | 0.4625 | 0.4439 | 0.2184 | 0.1175 | 0.0373 | 0.1175 | 3.4878 | 6.2197 | 1.1576 | 0.7469 | 0.9663 | 0.8976 | 3.9364 |
| P4 | 50 | 4 | 3 | 22.1838 | 22.1965 | 21.1828 | 0.6281 | 0.5259 | 0.5911 | 0.3286 | 0.1036 | 0.3286 | 22.1838 | 87.2040 | 2.6504 | 1.8524 | 0.9852 | 0.9863 | 1.9117 |
| P4 | 50 | 4 | 6 | 6.2837 | 6.3107 | 6.1660 | 0.3639 | 0.3388 | 0.2321 | 0.0929 | 0.0284 | 0.0929 | 6.2837 | 12.5077 | 1.1959 | 0.5552 | 0.9852 | 0.9562 | 3.9176 |
| P4 | 100 | 3 | 6 | 11.5472 | 11.5121 | 11.3659 | 0.6372 | 0.5981 | 0.3468 | 0.2189 | 0.0698 | 0.2189 | 11.5472 | 43.5891 | 2.7191 | 1.8349 | 0.9810 | 0.9631 | 2.9110 |
| P4 | 100 | 3 | 10 | 6.7805 | 6.7870 | 6.7183 | 0.5361 | 0.5229 | 0.2076 | 0.1117 | 0.0361 | 0.1117 | 6.7805 | 18.0113 | 1.5671 | 1.0072 | 0.9835 | 0.9364 | 4.7981 |
| P4 | 100 | 4 | 3 | 61.9022 | 61.9585 | 59.5093 | 0.6774 | 0.5963 | 0.5662 | 0.3194 | 0.1026 | 0.3194 | 61.9022 | 375.7326 | 3.5422 | 2.5119 | 0.9948 | 0.9942 | 2.1210 |
| P4 | 100 | 4 | 6 | 18.3503 | 18.3605 | 18.1084 | 0.4281 | 0.4121 | 0.2214 | 0.0887 | 0.0288 | 0.0887 | 18.3503 | 62.6378 | 1.5357 | 0.7485 | 0.9952 | 0.9814 | 4.8276 |

[Full E1 rows](../results/mechanism_suite/e1_oracle_decomposition.csv) · [Oracle bars](../results/mechanism_suite/figures/e1_oracle_decomposition.png)

### Child u/z noise correlation

The correlation is reported rather than suppressing the known P2 `oracle_u` anomaly. The z scalar is the mean coordinate error, proportional to the sum-z error used by the product driver.

| pde | box | oracle_state | oracle_u | oracle_z | raw | segment |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | -0.0447 | -0.0459 | -0.0449 | -0.0459 | -0.0449 | -0.0451 |
| P2_a4 | -0.1711 | -0.1670 | -0.1342 | -0.1745 | -0.1580 | — |
| P2_a8 | -0.2369 | -0.2375 | -0.1977 | -0.2441 | -0.2238 | — |
| P3_rho1 | -0.2093 | -0.2051 | -0.1799 | -0.2104 | -0.1968 | — |
| P3_rho2 | -0.1625 | -0.1583 | -0.1358 | -0.1632 | -0.1505 | — |
| P4 | -4.508e-04 | -5.090e-04 | 5.214e-04 | -5.090e-04 | 5.214e-04 | -4.798e-04 |

Replacing only u can remove covariance-driven cancellation while leaving the noisy sum-z channel intact, but the correlation is dimension-dependent and need not stay negative. The following dimension-resolved table exposes the sign and magnitude together with skill and generator bias. Where the correlation weakens or reverses, covariance cancellation cannot by itself explain the adverse `oracle_u` result.

| pde | d | method | test_skill | child_u_z_noise_correlation | generator_bias | generator_mae |
| --- | --- | --- | --- | --- | --- | --- |
| P2_a4 | 20 | oracle_state | 0.0689 | -0.3548 | 0.0000 | 0.0000 |
| P2_a4 | 20 | oracle_u | 0.3740 | -0.3351 | 0.0050 | 0.2707 |
| P2_a4 | 20 | oracle_z | 0.0804 | -0.3582 | -0.0465 | 0.0557 |
| P2_a4 | 20 | raw | 0.3116 | -0.3475 | -0.0584 | 0.2418 |
| P2_a4 | 50 | oracle_state | 0.1030 | -0.1559 | 0.0000 | 0.0000 |
| P2_a4 | 50 | oracle_u | 0.9809 | -0.1202 | 0.0328 | 0.5739 |
| P2_a4 | 50 | oracle_z | 0.1600 | -0.1633 | -0.1754 | 0.1918 |
| P2_a4 | 50 | raw | 0.8598 | -0.1452 | -0.1619 | 0.5300 |
| P2_a4 | 100 | oracle_state | 0.1469 | 0.0097 | 0.0000 | 0.0000 |
| P2_a4 | 100 | oracle_u | 2.0670 | 0.0528 | 0.0120 | 0.9978 |
| P2_a4 | 100 | oracle_z | 0.3120 | -0.0020 | -0.4096 | 0.4318 |
| P2_a4 | 100 | raw | 2.9891 | 0.0187 | -0.3573 | 0.9410 |
| P2_a8 | 20 | oracle_state | 0.0721 | -0.4471 | 0.0000 | 0.0000 |
| P2_a8 | 20 | oracle_u | 0.4662 | -0.4206 | 0.0251 | 0.2758 |
| P2_a8 | 20 | oracle_z | 0.0839 | -0.4499 | -0.0403 | 0.0463 |
| P2_a8 | 20 | raw | 0.3824 | -0.4351 | -0.0384 | 0.2362 |
| P2_a8 | 50 | oracle_state | 0.1150 | -0.2182 | 0.0000 | 0.0000 |
| P2_a8 | 50 | oracle_u | 1.1967 | -0.1767 | 0.0800 | 0.6074 |
| P2_a8 | 50 | oracle_z | 0.1776 | -0.2246 | -0.1628 | 0.1752 |
| P2_a8 | 50 | raw | 1.0437 | -0.2035 | -0.1206 | 0.5364 |
| P2_a8 | 100 | oracle_state | 0.1607 | -0.0474 | 0.0000 | 0.0000 |
| P2_a8 | 100 | oracle_u | 2.4030 | 0.0040 | 0.1012 | 1.0786 |
| P2_a8 | 100 | oracle_z | 0.3377 | -0.0579 | -0.4007 | 0.4188 |
| P2_a8 | 100 | raw | 3.3024 | -0.0328 | -0.2965 | 0.9691 |

## E2: causal one-shot dose response

At each generator evaluation the child estimate is replaced by exact state plus the registered perturbation, so a dose-induced state error never enters a later nonlinear generator (the resulting estimator contributions can still aggregate linearly). The P3 u-channel uses scale 0.5, the half-width of its certified interval [0,1], because the protocol requested that channel but did not state a scale. Raw-dose monotonicity is evaluated separately for both registered configurations and every admitted dimension in `analysis_summary.json`.

[Dose rows](../results/mechanism_suite/e2_dose.csv) · [Method-comparison curves](../results/mechanism_suite/figures/e2_dose_response.png) · [All raw-dose dimensions](../results/mechanism_suite/figures/e2_raw_dose_all_dimensions.png)

## E3: geometry or suppression?

Every tuned factor was chosen solely by mean validation skill within its (PDE,d,n,M) cell, then frozen for the test split. The table gives strict paired box wins out of 10; ties are not wins.
C3 is evaluated jointly: a cell passes only when box records at least 7/10 wins against all four registered comparators in that same cell; at least 75% of a PDE's E3 cells must pass. Per-comparator marginal fractions are retained in `analysis_summary.json` as diagnostics.
The certificate-prior-dominated flag uses cell-level mean test skill: a prior control ‘beats’ box when its mean is lower.

| pde | d | n | M | box_wins_vs_best_shrink | comparator_best_shrink | box_mean_minus_best_shrink_mean | box_wins_vs_best_illegal | comparator_best_illegal | box_mean_minus_best_illegal_mean | box_wins_vs_centre | comparator_centre | box_mean_minus_centre_mean | box_wins_vs_best_shrink_centre | comparator_best_shrink_centre | box_mean_minus_best_shrink_centre_mean | box_wins_vs_raw | comparator_raw | box_mean_minus_raw_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | 20 | 3 | 6 | 10 | shrink_c0.25 | -0.3051 | 10 | illegal_a0.9 | -0.0037 | 10 | centre | -0.0796 | 5 | shrink_centre_c0.25 | -0.0017 | 10 | raw | -57.6116 |
| P1 | 20 | 3 | 10 | 10 | shrink_c0.5 | -0.1129 | 10 | illegal_a0.9 | -0.0084 | 10 | centre | -0.1122 | 10 | shrink_centre_c0.25 | -0.0425 | 10 | raw | -12.8736 |
| P1 | 20 | 4 | 3 | 10 | shrink_c0.1 | -0.2879 | 0 | illegal_a0.75 | 0.0091 | 10 | centre | -0.0260 | 10 | shrink_centre_c0.25 | -0.8703 | 10 | raw | -5.398e+07 |
| P1 | 20 | 4 | 6 | 10 | shrink_c0.25 | -0.3514 | 10 | illegal_a0.9 | -0.0042 | 10 | centre | -0.1086 | 10 | shrink_centre_c0.25 | -0.0498 | 10 | raw | -1.038e+05 |
| P1 | 50 | 3 | 6 | 10 | shrink_c0.25 | -0.1961 | 10 | illegal_a0.9 | -0.0047 | 10 | centre | -0.0802 | 10 | shrink_centre_c0.25 | -0.1734 | 10 | raw | -544.2730 |
| P1 | 50 | 3 | 10 | 10 | shrink_c0.25 | -0.3164 | 10 | illegal_a0.9 | -0.0089 | 10 | centre | -0.1050 | 0 | shrink_centre_c0.25 | 0.0238 | 10 | raw | -111.8107 |
| P1 | 50 | 4 | 3 | 10 | shrink_c0.1 | -0.3041 | 1 | illegal_a0.9 | 9.281e-04 | 10 | centre | -0.0513 | 10 | shrink_centre_c0.25 | -129.5461 | 10 | raw | -1.639e+10 |
| P1 | 50 | 4 | 6 | 10 | shrink_c0.25 | -0.2769 | 10 | illegal_a0.9 | -0.0071 | 10 | centre | -0.1125 | 10 | shrink_centre_c0.25 | -0.2455 | 10 | raw | -1.804e+07 |
| P2_a4 | 20 | 3 | 6 | 10 | shrink_c0.75 | -0.0436 | 10 | illegal_a0.9 | -0.0185 | 10 | centre | -0.0489 | 0 | shrink_centre_c0.5 | 0.0084 | 10 | raw | -0.0577 |
| P2_a4 | 20 | 3 | 10 | 10 | shrink_c0.75 | -0.0308 | 10 | illegal_a0.9 | -0.0190 | 10 | centre | -0.0768 | 0 | shrink_centre_c0.75 | 0.0246 | 0 | raw | 0.0082 |
| P2_a4 | 20 | 4 | 3 | 10 | shrink_c0.5 | -0.1475 | 10 | illegal_a0.9 | -0.0137 | 10 | centre | -0.0353 | 10 | shrink_centre_c0.25 | -0.0080 | 10 | raw | -0.4131 |
| P2_a4 | 20 | 4 | 6 | 10 | shrink_c0.75 | -0.0715 | 10 | illegal_a0.9 | -0.0146 | 10 | centre | -0.0982 | 10 | shrink_centre_c0.5 | -0.0191 | 10 | raw | -0.0992 |
| P2_a4 | 50 | 3 | 6 | 10 | shrink_c0.75 | -0.0358 | 10 | illegal_a0.9 | -0.0126 | 10 | centre | -0.0466 | 6 | shrink_centre_c0.5 | -0.0039 | 10 | raw | -0.1232 |
| P2_a4 | 50 | 3 | 10 | 2 | shrink_c0.75 | 0.0046 | 10 | illegal_a0.9 | -0.0143 | 10 | centre | -0.0664 | 0 | shrink_centre_c0.75 | 0.0276 | 3 | raw | 0.0086 |
| P2_a4 | 50 | 4 | 3 | 10 | shrink_c0.5 | -0.2108 | 10 | illegal_a0.9 | -0.0040 | 10 | centre | -0.0313 | 10 | shrink_centre_c0.25 | -0.0423 | 10 | raw | -1.6277 |
| P2_a4 | 50 | 4 | 6 | 10 | shrink_c0.75 | -0.0646 | 10 | illegal_a0.9 | -0.0036 | 10 | centre | -0.0707 | 10 | shrink_centre_c0.5 | -0.0236 | 10 | raw | -0.3383 |
| P2_a8 | 20 | 3 | 6 | 10 | shrink_c0.75 | -0.0553 | 10 | illegal_a0.9 | -0.0158 | 10 | centre | -0.0732 | 10 | shrink_centre_c0.5 | -0.0070 | 10 | raw | -0.0922 |
| P2_a8 | 20 | 3 | 10 | 10 | shrink_c0.75 | -0.0295 | 10 | illegal_a0.9 | -0.0168 | 10 | centre | -0.1042 | 0 | shrink_centre_c0.75 | 0.0149 | 9 | raw | -0.0071 |
| P2_a8 | 20 | 4 | 3 | 10 | shrink_c0.5 | -0.1608 | 10 | illegal_a0.9 | -0.0095 | 10 | centre | -0.0531 | 10 | shrink_centre_c0.25 | -0.0217 | 10 | raw | -0.5391 |
| P2_a8 | 20 | 4 | 6 | 10 | shrink_c0.75 | -0.0755 | 10 | illegal_a0.9 | -0.0116 | 10 | centre | -0.1222 | 10 | shrink_centre_c0.5 | -0.0323 | 10 | raw | -0.1284 |
| P2_a8 | 50 | 3 | 6 | 10 | shrink_c0.75 | -0.0160 | 9 | illegal_a0.9 | -0.0014 | 10 | centre | -0.0919 | 10 | shrink_centre_c0.75 | -0.0361 | 10 | raw | -0.1433 |
| P2_a8 | 50 | 3 | 10 | 0 | shrink_c0.75 | 0.0389 | 10 | illegal_a0.9 | -0.0044 | 10 | centre | -0.1195 | 0 | shrink_centre_c0.75 | 0.0352 | 0 | raw | 0.0196 |
| P2_a8 | 50 | 4 | 3 | 10 | shrink_c0.5 | -0.2020 | 8 | illegal_a0.75 | -0.0044 | 10 | centre | -0.0621 | 10 | shrink_centre_c0.25 | -0.0622 | 10 | raw | -1.9923 |
| P2_a8 | 50 | 4 | 6 | 10 | shrink_c0.75 | -0.0466 | 10 | illegal_a0.75 | -0.0099 | 10 | centre | -0.1157 | 10 | shrink_centre_c0.5 | -0.0321 | 10 | raw | -0.3878 |
| P3_rho1 | 20 | 3 | 6 | 10 | shrink_c0.75 | -0.0410 | 10 | illegal_a0.9 | -0.0178 | 10 | centre | -0.0346 | 0 | shrink_centre_c0.5 | 0.0100 | 10 | raw | -0.0728 |
| P3_rho1 | 20 | 3 | 10 | 10 | shrink_c0.75 | -0.0157 | 10 | illegal_a0.9 | -0.0194 | 10 | centre | -0.0564 | 0 | shrink_centre_c0.75 | 0.0270 | 5 | raw | 0.0011 |
| P3_rho1 | 20 | 4 | 3 | 10 | shrink_c0.5 | -0.1442 | 10 | illegal_a0.9 | -0.0136 | 10 | centre | -0.0224 | 10 | shrink_centre_c0.25 | -0.0081 | 10 | raw | -0.4721 |
| P3_rho1 | 20 | 4 | 6 | 10 | shrink_c0.75 | -0.0721 | 10 | illegal_a0.9 | -0.0161 | 10 | centre | -0.0777 | 10 | shrink_centre_c0.5 | -0.0144 | 10 | raw | -0.1250 |
| P3_rho1 | 50 | 3 | 6 | 10 | shrink_c0.75 | -0.0462 | 10 | illegal_a0.9 | -0.0123 | 10 | centre | -0.0356 | 8 | shrink_centre_c0.5 | -0.0092 | 10 | raw | -0.1554 |
| P3_rho1 | 50 | 3 | 10 | 3 | shrink_c0.75 | 0.0038 | 10 | illegal_a0.9 | -0.0143 | 10 | centre | -0.0518 | 0 | shrink_centre_c0.75 | 0.0143 | 10 | raw | -0.0173 |
| P3_rho1 | 50 | 4 | 3 | 10 | shrink_c0.5 | -0.2229 | 10 | illegal_a0.9 | -0.0055 | 10 | centre | -0.0226 | 10 | shrink_centre_c0.25 | -0.0457 | 10 | raw | -1.6978 |
| P3_rho1 | 50 | 4 | 6 | 10 | shrink_c0.75 | -0.0890 | 10 | illegal_a0.9 | -0.0067 | 10 | centre | -0.0584 | 10 | shrink_centre_c0.5 | -0.0264 | 10 | raw | -0.3921 |
| P3_rho2 | 20 | 3 | 6 | 10 | shrink_c0.75 | -0.0424 | 10 | illegal_a0.9 | -0.0166 | 10 | centre | -0.0241 | 1 | shrink_centre_c0.5 | 0.0083 | 10 | raw | -0.0871 |
| P3_rho2 | 20 | 3 | 10 | 10 | shrink_c0.75 | -0.0092 | 10 | illegal_a0.9 | -0.0187 | 10 | centre | -0.0412 | 0 | shrink_centre_c0.5 | 0.0236 | 10 | raw | -0.0098 |
| P3_rho2 | 20 | 4 | 3 | 10 | shrink_c0.5 | -0.1465 | 10 | illegal_a0.9 | -0.0123 | 10 | centre | -0.0123 | 10 | shrink_centre_c0.25 | -0.0121 | 10 | raw | -0.5288 |
| P3_rho2 | 20 | 4 | 6 | 10 | shrink_c0.75 | -0.0791 | 10 | illegal_a0.9 | -0.0165 | 10 | centre | -0.0587 | 10 | shrink_centre_c0.5 | -0.0151 | 10 | raw | -0.1532 |
| P3_rho2 | 50 | 3 | 6 | 10 | shrink_c0.75 | -0.0586 | 10 | illegal_a0.9 | -0.0119 | 10 | centre | -0.0265 | 7 | shrink_centre_c0.25 | -0.0027 | 10 | raw | -0.1826 |
| P3_rho2 | 50 | 3 | 10 | 5 | shrink_c0.75 | -0.0015 | 10 | illegal_a0.9 | -0.0139 | 10 | centre | -0.0393 | 0 | shrink_centre_c0.5 | 0.0078 | 10 | raw | -0.0415 |
| P3_rho2 | 50 | 4 | 3 | 10 | shrink_c0.25 | -0.2572 | 10 | illegal_a0.9 | -0.0062 | 10 | centre | -0.0148 | 10 | shrink_centre_c0.25 | -0.0509 | 10 | raw | -1.7770 |
| P3_rho2 | 50 | 4 | 6 | 10 | shrink_c0.75 | -0.1150 | 10 | illegal_a0.9 | -0.0092 | 10 | centre | -0.0473 | 10 | shrink_centre_c0.25 | -0.0157 | 10 | raw | -0.4376 |
| P4 | 20 | 3 | 6 | 0 | shrink_c0.5 | 0.0382 | 0 | illegal_a0.5 | 0.1924 | 10 | centre | -0.3633 | 0 | shrink_centre_c0.5 | 0.0382 | 10 | raw | -1.9093 |
| P4 | 20 | 3 | 10 | 0 | shrink_c0.5 | 0.0993 | 0 | illegal_a0.5 | 0.1288 | 10 | centre | -0.4743 | 0 | shrink_centre_c0.5 | 0.0993 | 10 | raw | -0.9876 |
| P4 | 20 | 4 | 3 | 5 | shrink_c0.25 | -4.862e-04 | 0 | illegal_a0.5 | 0.1839 | 10 | centre | -0.3162 | 5 | shrink_centre_c0.25 | -4.862e-04 | 10 | raw | -5.2375 |
| P4 | 20 | 4 | 6 | 0 | shrink_c0.5 | 0.0206 | 0 | illegal_a0.75 | 0.0983 | 10 | centre | -0.5371 | 0 | shrink_centre_c0.5 | 0.0206 | 10 | raw | -1.2832 |
| P4 | 50 | 3 | 6 | 0 | shrink_c0.25 | 0.1895 | 0 | illegal_a0.5 | 0.2993 | 10 | centre | -0.2573 | 0 | shrink_centre_c0.25 | 0.1895 | 10 | raw | -5.5072 |
| P4 | 50 | 3 | 10 | 0 | shrink_c0.25 | 0.0214 | 0 | illegal_a0.5 | 0.2391 | 10 | centre | -0.3647 | 0 | shrink_centre_c0.25 | 0.0214 | 10 | raw | -3.0254 |
| P4 | 50 | 4 | 3 | 0 | shrink_c0.25 | 0.0965 | 0 | illegal_a0.5 | 0.2633 | 10 | centre | -0.2443 | 0 | shrink_centre_c0.25 | 0.0965 | 10 | raw | -21.5557 |
| P4 | 50 | 4 | 6 | 10 | shrink_c0.25 | -0.0974 | 0 | illegal_a0.5 | 0.1494 | 10 | centre | -0.4601 | 10 | shrink_centre_c0.25 | -0.0974 | 10 | raw | -5.9197 |

[Control rows](../results/mechanism_suite/e3_controls.csv) · [Frozen choices](../results/mechanism_suite/e3_tuning_choices.json) · [Win heatmap](../results/mechanism_suite/figures/e3_win_count_heatmap.png)

## E4: certificate ablation

| pde | ablation_role | mean | median |
| --- | --- | --- | --- |
| P1 | ball | 0.5110 | 0.5173 |
| P1 | box | 0.1730 | 0.1616 |
| P1 | segment | 0.2374 | 0.2090 |
| P1 | span_only | 573.1421 | 2.3496 |
| P2_a4 | best_validation_looser | 0.2646 | 0.2727 |
| P2_a4 | box | 0.2554 | 0.2614 |
| P2_a4 | sign_only | 6.7031 | 2.3369 |
| P2_a8 | best_validation_looser | 0.3362 | 0.3414 |
| P2_a8 | box | 0.3042 | 0.3114 |
| P2_a8 | sign_only | 8.4160 | 3.0027 |
| P3_rho1 | best_validation_looser | 0.2641 | 0.2762 |
| P3_rho1 | box | 0.2604 | 0.2631 |
| P3_rho1 | sign_only | 6.8812 | 2.5164 |
| P3_rho2 | best_validation_looser | 0.2681 | 0.2759 |
| P3_rho2 | box | 0.2684 | 0.2687 |
| P3_rho2 | sign_only | 6.9924 | 2.6749 |
| P4 | ball | 0.8172 | 0.8737 |
| P4 | box | 0.4653 | 0.4690 |
| P4 | segment | 0.3533 | 0.2976 |
| P4 | span_only | 0.4884 | 0.3943 |

[Ablation rows](../results/mechanism_suite/e4_ablation.csv) · [Ablation bars](../results/mechanism_suite/figures/e4_certificate_ablation.png)

## E5: negative controls

| pde | d | n | M | box | centre | oracle_state | raw |
| --- | --- | --- | --- | --- | --- | --- | --- |
| N3 | 100 | 3 | 6 | 0.4624 | 0.0875 | 0.0868 | 28.8970 |
| N3 | 100 | 4 | 6 | 0.2424 | 0.0373 | 0.0355 | 4615.2420 |
| N4 | 20 | 3 | 6 | 1.1499 | 1.0718 | 0.3623 | 1.3467 |
| N4 | 20 | 4 | 6 | 0.9364 | 0.9005 | 0.1482 | 5.1142 |
| N4 | 50 | 3 | 6 | 1.4935 | 1.4424 | 0.4708 | 1.7971 |
| N4 | 50 | 4 | 6 | 1.3803 | 1.2913 | 0.1937 | 11.8517 |
| N4 | 100 | 3 | 6 | 1.4883 | 1.4700 | 0.5372 | 1.8812 |
| N4 | 100 | 4 | 6 | 1.4397 | 1.4013 | 0.2167 | 8.6007 |

N3's statement “centre >= box” is interpreted as predictive performance, hence centre skill <= box skill. N4 passes its prediction only when oracle-state skill remains above 0.10 in every tested cell.

N3 prediction: PASS. N4 prediction: PASS.

Historical no-new-run controls: N1 is documented in [HJB_LIFE_OR_DEATH_ABLATION.md](HJB_LIFE_OR_DEATH_ABLATION.md) and [hjb_life_or_death_summary.json](../results/hjb_life_or_death_summary.json); N2 in [FUNDING_LIFE_OR_DEATH_ABLATION.md](FUNDING_LIFE_OR_DEATH_ABLATION.md) and [funding_life_or_death_summary.json](../results/funding_life_or_death_summary.json); N5 in [batchir_negative_controls.md](batchir_negative_controls.md). Their supplied gate diagnoses (G2 for N1/N2, G4 for N5) were cited, not rerun.

[Negative-control rows](../results/mechanism_suite/e5_negative.csv)

## E6: P1 dimension sweep

| d | box | centre | oracle_state | raw | segment | box_over_oracle |
| --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.1429 | 0.2515 | 0.0291 | 1.038e+05 | 0.1948 | 4.9106 |
| 50 | 0.1409 | 0.2535 | 0.0287 | 1.804e+07 | 0.1945 | 4.9111 |
| 100 | 0.1445 | 0.2469 | 0.0294 | 1.212e+09 | 0.1938 | 4.9162 |
| 200 | 0.1463 | 0.2392 | 0.0283 | 1.127e+11 | 0.1869 | 5.1670 |

No numerical threshold for dimension-independence was pre-registered. The observed mean-skill `box/oracle_state` ratio ranges from 4.9106 to 5.1670; relative range=0.0515, coefficient of variation=0.0221, and max/min=1.0522.

[Dimension rows](../results/mechanism_suite/e6_dimension.csv) · [Dimension plot](../results/mechanism_suite/figures/e6_dimension.png)

## P4 rectification mechanism

The table pools raw E1 generator calls and weights each s-bin by its number of calls; the central bin [-0.25,0.25) is the pre-specified neighborhood of the rectification point psi_s=0.

| s_low | s_high | count | generator_bias | generator_mae |
| --- | --- | --- | --- | --- |
| −∞ | -2.0000 | 339022 | -2.1192 | 2.1295 |
| -2.0000 | -1.0000 | 4224271 | -2.0806 | 2.0881 |
| -1.0000 | -0.2500 | 14689181 | -1.7558 | 1.7593 |
| -0.2500 | 0.2500 | 14821990 | -1.4712 | 1.4726 |
| 0.2500 | 1.0000 | 15449336 | -1.7505 | 1.7540 |
| 1.0000 | 2.0000 | 4775593 | -2.0950 | 2.1028 |
| 2.0000 | ∞ | 384607 | -2.1159 | 2.1263 |

The bias is negative in every populated bin, matching the direction expected when noisy gradients inflate the norm inside the negative norm driver.
The central-bin raw generator bias/MAE is -1.4712/1.4726; the largest outer-bin bias magnitude is 2.1192 (MAE 2.1295) on [−∞,-2.0000). Thus the pre-stated localization prediction that the Jensen damage is largest near psi_s=0 is not supported by the raw recursive E1 diagnostic.
The raw recursive diagnostic combines Jensen rectification with state-dependent recursive variance, so a failed localization pattern does not contradict the controlled E2 dose response; it does reject the stronger claim that raw E1 damage is maximized at the rectification point.

[P4 s-bin plot](../results/mechanism_suite/figures/p4_generator_bias_by_s_bin.png)

## Failures and adverse findings

- N3 failed gate(s) G2 and was excluded from C1–C3.
- N4 failed gate(s) G3 and was excluded from C1–C3.
- N1 failed gate(s) G2 and was excluded from C1–C3.
- N2 failed gate(s) G2 and was excluded from C1–C3.
- N5 failed gate(s) G4 and was excluded from C1–C3.
- P1: pre-registered C3 criterion failed (joint passing-cell fraction=0.5000; required >=0.75).
- P2_a4: pre-registered C2 criterion failed (d=20: 1/4 cells; d=50: 2/4 cells; d=100: 2/4 cells).
- P2_a4: pre-registered C3 criterion failed (joint passing-cell fraction=0.5000; required >=0.75).
- P2_a4: certificate-prior dominated in the E3 regime (centre fraction=0.0000, shrink-centre fraction=0.3750).
- P2_a8: pre-registered C2 criterion failed (d=20: 1/4 cells; d=50: 2/4 cells; d=100: 2/4 cells).
- P3_rho1: pre-registered C2 criterion failed (d=20: 1/4 cells; d=50: 2/4 cells).
- P3_rho1: pre-registered C3 criterion failed (joint passing-cell fraction=0.6250; required >=0.75).
- P3_rho1: certificate-prior dominated in the E3 regime (centre fraction=0.0000, shrink-centre fraction=0.3750).
- P3_rho2: pre-registered C2 criterion failed (d=20: 2/4 cells; d=50: 2/4 cells).
- P3_rho2: pre-registered C3 criterion failed (joint passing-cell fraction=0.6250; required >=0.75).
- P3_rho2: certificate-prior dominated in the E3 regime (centre fraction=0.0000, shrink-centre fraction=0.3750).
- P4: pre-registered C3 criterion failed (joint passing-cell fraction=0.0000; required >=0.75).
- P4: certificate-prior dominated in the E3 regime (centre fraction=0.0000, shrink-centre fraction=0.7500).
- P4: the pre-stated rectification-localization prediction was not observed (central-bin |generator bias|=1.4712; largest outer-bin |generator bias|=2.1192).

## Work accounting

Across E1/E2/E3/E5/E6: 22,240 method-repetition rows (2,818 exact artifact reuses), 6,934,488,000 non-duplicated generator calls, 6,957,888,000 recursively evaluated states, 80,428,224,000 MLP stochastic samples, 45,279,648,000 additional registered dose-normal variates, and 41.77 non-duplicated summed worker-hours. E0-sourced E5/E6 artifacts count once here because E0 is outside this table. Rows with any non-finite state/generator count: 0.

## Reproduction

```powershell
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage all --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.analyze_suite
```

The runner is resumable and validates existing per-repetition artifacts before skipping them.
