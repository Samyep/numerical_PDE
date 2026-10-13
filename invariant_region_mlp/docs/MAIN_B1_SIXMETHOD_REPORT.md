# Main B1 six-method benchmark report

> **2026-10-13 follow-up:** This frozen table uses the original small actor. After increasing actor label quality and training to cross-budget convergence, Actor--MLP is best for value prediction in 12/15 cells and all 12 cells with `k >= 5`. Gradient prediction does not match the best neural gradient baseline. See [the matched follow-up benchmark](MAIN_B1_CONVERGED_ACTOR_VS_BASELINES.md). The original table below is retained for provenance.

**Numerical answer:** No. Complete Actor--MLP is best in 0/15 cells; it beats Vanilla in 15/15 and Double in 13/15. The best of the other five methods wins 15/15 cells. Total measured experiment use is 3.496 GH200 h.

All entries below are actual held-out results at `t=0`, `x~N(0,I_d)`. MLP-family methods use `n=3, M=16` and 10 evaluation seeds. Neural baselines use five independent training seeds. Values are RMSE with a 95% normal-approximation interval across seed-level RMSEs.
Ours trains the one frozen actor per `(d,k)` required by the protocol, then uses 10 independent recursive evaluation seeds. Its interval therefore measures MLP evaluation variation conditional on that actor; it does not measure actor-retraining variation.

## Main Table 1 — accuracy

| d | k | MLP | Double | DPI | Deep BSDE | PINN+SDGD | Ours |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1 | 714.3 ± 3.1 | 0.6381 ± 0.021 | 0.1356 ± 0.026 | 0.188 ± 0.0092 | **0.06497 ± 0.0067** | 0.1807 ± 0.0022 |
| 100 | 5 | 244.5 ± 2.4 | 0.2596 ± 0.0053 | **0.1179 ± 0.0076** | 0.5914 ± 0.027 | 0.1232 ± 0.026 | 0.1719 ± 0.0019 |
| 100 | 10 | 220.6 ± 1.2 | 0.2408 ± 0.0063 | 0.161 ± 0.014 | 0.8006 ± 0.044 | **0.1022 ± 0.0084** | 0.1747 ± 0.0017 |
| 100 | 20 | 124 ± 0.99 | 0.1574 ± 0.0042 | 0.1616 ± 0.014 | 0.7894 ± 0.041 | **0.1125 ± 0.0064** | 0.1317 ± 0.00073 |
| 100 | 40 | 83.21 ± 0.74 | 0.1162 ± 0.0022 | 0.1893 ± 0.0079 | 0.668 ± 0.004 | **0.07866 ± 0.0042** | 0.1268 ± 0.00052 |
| 200 | 1 | 4548 ± 22 | 0.9688 ± 0.031 | 0.1226 ± 0.01 | 0.2314 ± 0.01 | **0.07322 ± 0.0074** | 0.2216 ± 0.0012 |
| 200 | 5 | 1676 ± 10 | 0.3876 ± 0.0092 | 0.14 ± 0.012 | 0.8349 ± 0.023 | **0.1318 ± 0.021** | 0.2728 ± 0.0045 |
| 200 | 10 | 1531 ± 8.8 | 0.3542 ± 0.009 | **0.1449 ± 0.013** | 0.9493 ± 0.022 | 0.1663 ± 0.03 | 0.3237 ± 0.0071 |
| 200 | 20 | 862.4 ± 5.7 | 0.2211 ± 0.0078 | 0.1983 ± 0.021 | 0.8393 ± 0.01 | **0.1237 ± 0.015** | 0.1836 ± 0.00061 |
| 200 | 40 | 541.4 ± 4.7 | 0.1539 ± 0.0021 | 0.2124 ± 0.012 | 0.654 ± 0.0068 | **0.08023 ± 0.0055** | 0.1864 ± 0.00038 |
| 400 | 1 | 3.291e+04 ± 1.9e+02 | 1.794 ± 0.072 | 0.1557 ± 0.029 | 0.2914 ± 0.022 | **0.08632 ± 0.013** | 0.3302 ± 0.0021 |
| 400 | 5 | 1.269e+04 ± 85 | 0.7047 ± 0.019 | 0.1597 ± 0.02 | 1.091 ± 0.026 | **0.1318 ± 0.027** | 0.3426 ± 0.0019 |
| 400 | 10 | 1.199e+04 ± 76 | 0.6426 ± 0.022 | 0.1901 ± 0.0095 | 1.103 ± 0.024 | **0.1788 ± 0.018** | 0.3544 ± 0.0015 |
| 400 | 20 | 5955 ± 37 | 0.3539 ± 0.0082 | 0.1994 ± 0.011 | 0.7856 ± 0.0096 | **0.1253 ± 0.016** | 0.2778 ± 0.00056 |
| 400 | 40 | 4145 ± 30 | 0.258 ± 0.0081 | 0.2459 ± 0.014 | 0.645 ± 0.0066 | **0.087 ± 0.0067** | 0.2472 ± 0.00019 |

## Comparative result

Ours/Vanilla geometric-mean RMSE ratio: 0.0001536; Ours/Double: 0.6151; Ours/DPI: 1.347; Ours/Deep BSDE: 0.3575; Ours/PINN+SDGD: 2.094.

## Cost and protocol

| Method | Training GPU h | Full-test GPU h | Total GPU h | Parameters | Peak GiB |
|---|---:|---:|---:|---:|---:|
| MLP | 0.0000 | 0.0063 | 0.0063 | 0 | 52.26 |
| Double | 0.0000 | 0.0138 | 0.0138 | 0 | 53.73 |
| Ours | 0.0029 | 0.0098 | 0.0128 | 42,468–119,568 | 56.42 |
| DPI | 0.5933 | 0.0000 | 0.5933 | 840,705–994,305 | 1.70 |
| PINN+SDGD | 1.8857 | 0.0001 | 1.8858 | 223,745–300,545 | 0.38 |
| Deep BSDE | 0.6269 | 0.0000 | 0.6269 | 105,061–220,561 | 0.12 |

Summed measured method-cell work is 3.139 GPU h. Scheduler-accounted total experiment use is 3.496 GH200 h. `main_cost_table.csv` separates training and full-test inference. Actor label generation and training are included under Ours; independent recursive branches are included for Double.

Deep BSDE is a global/amortized protocol trained over the same initial-state distribution; its `Z(t,x)` is learned independently. DPI uses control-variate BEL gradient labels and warm-started Picard rounds. PINN+SDGD uses two independent coordinate subsets, giving an unbiased gradient for the squared full residual, and hard-enforces the terminal condition.

The final DPI pass retains the authors' HJB value-network architecture (four width-512 ELU hidden layers) and at most 5,120 optimizer steps. Under the fixed allocation it uses 8 rather than 40 Picard rounds and 32/128 rather than 4,096 training/validation Monte Carlo samples per label. This profile was frozen from the public configuration and Stage A timing, without using held-out Stage B error. The preliminary 4x128 tanh pass is archived under `superseded/` and is excluded from all final tables.

The compute plot compares the measured cost and accuracy of the frozen primary configurations. The methods do not consume equal time, and no additional budget sweep was run, so this is not presented as an equal-spend frontier.

## Correctness and calibration

Stage 0 passed the terminal identity, analytic and autodiff PDE residual, gradient, active Hessian rank, P1 equivalence, diffusion scaling, actor identity, SDGD unbiasedness, BSDE sign, and CPU/GPU identical-draw recursion checks. Reference-only calibration froze the coefficient scale at 2.0 before solver evaluation; no cell was flagged as weakly nonlinear.

Analytic values and gradients were used only for Stage 0 and held-out diagnostics/evaluation. They were not training targets, regularizers, or checkpoint-selection signals.

## Supplementary diagnostics

All six methods completed all 15 cells; total non-finite value predictions: 0. The median held-out actor gradient RMSE is 0.9583; the median noisy training-label RMSE against the post-training analytic diagnostic is 2.342. Across 15 cells, actor label generation plus training took 10.588 GPU seconds.

`stage_b_results.csv` records normalized value RMSE, gradient RMSE, relative gradient error, value bias with intervals, value variation across runs, and non-finite counts for every method-cell.

## Artifacts

- `stage_b_results.csv`, `main_accuracy_table.csv`, `main_cost_table.csv`, and `main_cost_summary.csv`
- `actor_training.csv` and `pde_calibration.csv`
- `figures/accuracy_vs_effective_dimension.*`
- `figures/accuracy_vs_total_compute.*`
- `figures/gradient_accuracy.*`

## Method references

- [Deep Picard Iteration paper](https://arxiv.org/abs/2409.08526) and [authors' implementation](https://github.com/DeepOptimalControl/DeepPicardIteration).
- [Deep BSDE paper](https://doi.org/10.1073/pnas.1718942115) and [authors' implementation](https://github.com/frankhan91/DeepBSDE).
- [SDGD paper](https://arxiv.org/abs/2307.12306) and [authors' implementation](https://github.com/zheyuanhu01/SDGD_PINN).

No Stage C, component ablation, `(M,n)` sweep, or SCaSML run was performed.
