# Baseline paper-fidelity review

**Review status:** ready for independent review. All corrected seed-0 Stage-B
runs and the separate SDGD 250D fidelity anchor are complete. The dimension
sweep remains blocked by `REVIEW_GATE_STATUS.json` until an independent
scientific review explicitly approves continuation.

## Scope and decision rule

This audit asks whether each external baseline implements the published
method and reproduces a native paper result credibly. Source-code fidelity by
itself is insufficient. A custom-PDE run is reported separately as a
`faithful_transfer` or `adapted` result and cannot receive a paper-reproduction
PASS merely because it executes.

All new runs use the one fixed training seed `0`. Historical seed-1/2 runs and
superseded profiles remain in `HISTORICAL_AND_NONPRIMARY_RESULTS.csv`; they do
not enter primary rankings, means, checkpoint selection, or completion counts.
The original per-run training budgets are unchanged.

## Original papers and pinned sources

| Method | Paper evidence used | Pinned official source | Native configuration |
|---|---|---|---|
| DPI and PINN-HTE | Han, Hu, Long, and Zhao, *Deep Picard Iteration for High-Dimensional Nonlinear PDEs*, arXiv:2409.08526v2, Table 1, Section 5.3, Figures 4 and 9 | `DeepOptimalControl/DeepPicardIteration` at `e8616c2f531c9f95299b246a13c1180b8f02822a` | `scripts/hjb/base_100d_T1.0_w0.1_0.yaml` and `scripts/hjb/pinn_100d_T1.0_v16.yaml` |
| SDGD-PINN | Hu, Shukla, Karniadakis, and Kawaguchi, *Tackling the Curse of Dimensionality with Physics-Informed Neural Networks*, arXiv:2307.12306v6, Sections 5.2--5.3, Tables 4--5 | `zheyuanhu01/SDGD_PINN` at `c597cc4fe9caac0b3ca118ea84453405f25ff833` | `Section_5.3/HJB_Log.py`; `Section_5.2/FBSNNs_PINN.py` plus its PDE subclasses |
| Deep BSDE | Han, Jentzen, and E, *Solving High-Dimensional Partial Differential Equations Using Deep Learning*, PNAS 2018 / arXiv:1707.02568v3, Figures 2--3 | `frankhan91/DeepBSDE` at `e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` | `configs/hjb_lq_d100.json` and `configs/allencahn_d100.json` |

The checked-out source revisions exactly match these hashes. The frozen
matrix in `experiments/tangent_actor_icml_final_suite/BASELINE_FIDELITY_MATRIX.md`
records the transfer classification for every suite cell.

## Native numerical reproduction

The DPI paper reports Figure-4 bar plots rather than a numerical table for the
100D, `T=1` OU/GMM case. Values labeled “figure read” below are approximate
visual readings from that plot and are not represented as exact paper data.
Its bars summarize three paper seeds; our primary result is the mandated
single seed 0.

| Method | Original paper setting | Original published performance | Our reproduced native seed-0 performance | Fidelity verdict |
|---|---|---|---|---|
| DPI | OU/GMM HJB, 100D, `T=1`; distributed value/gradient task | rMAE about `0.032`, g-rMAE about `0.34` (Figure 4, figure read) | rMAE `0.0394264`, g-rMAE `0.312475`; relative L2 `0.0569723` | **PASS**: value is close and gradient is slightly better than the plotted paper bar |
| PINN-HTE | same OU/GMM HJB task, HTE-16 | rMAE about `0.063`, g-rMAE about `0.50` (Figure 4, figure read) | rMAE `0.0874229`, g-rMAE `0.582427`; relative L2 `0.124724` | **CONDITIONAL**: source configuration is faithful, but both errors are worse than the paper plot and the single-seed gap remains unresolved |
| Deep BSDE | HJB-Log, 100D, scalar `u(0,0)` | relative error `0.17%`, five-run summary (Figure 2) | absolute error `0.0103037`, relative error `0.224473%` | **PASS** jointly with the Allen row: the fixed seed is close to the paper result |
| Deep BSDE | Allen--Cahn, 100D, scalar `u(0,0)` | relative error `0.30%`, five-run summary (Figure 3) | absolute error `0.000372480`, relative error `0.705704%` | **PASS** with historical context: the mandated seed 0 is high, while the preserved three-seed reproduction mean is `0.32039%` and no seed was selected by accuracy |
| SDGD-PINN | HJB-Log, 250D, distributed relative L2 | `6.147e-3`, `210` minutes, five seeds (Table 5) | relative L2 `9.07227e-3`; RMSE `4.97144e-2`; `204.39` training minutes | **CONDITIONAL**: the source algorithm and runtime reproduce closely, but the fixed seed-0 error is `47.6%` above the five-seed paper value |

DPI/PINN native metrics were recomputed on the frozen 960-state subset of the
authors' evaluation law; rMAE and g-rMAE use the exact paper definitions.
Deep BSDE remains a single-initial-state solver and appears only in the native
point table. Its full training cost is charged for each solved initial state.

## Architecture, budget, terminal treatment, and sampling audit

| Method/profile | Paper/source budget | Actual corrected runner | Match assessment |
|---|---|---|---|
| DPI OU/GMM | PISGradNet; 40 Picard rounds; `M=N=4096`; 16 epochs/round; batch 512; double; Adam `1e-3`; gradient supervision | Identical source objects and YAML; four 512-unit ELU blocks plus the source scalar time network; hard terminal construction | Exact algorithm/configuration; only trusted-checkpoint loading and noninteractive logging compatibility changes |
| DPI generic semilinear transfer | Table-1 `K=80`, `M=N=4096`, 16 epochs/round; 4x128 ELU; batch 512; Adam `1e-3`; source gradient loss | Same, with the paper's `K=80` overriding the repository YAML's inconsistent `K=100` | Faithful method transfer; not a native paper reproduction |
| PINN-HTE OU/GMM | source PISGradNet; HTE-16; 20,000 epochs; hard terminal; 4x512 ELU; Adam `1e-3` | Identical YAML and source network | Exact algorithm/configuration |
| PINN-HTE generic semilinear transfer | source Burgers PINN; HTE-16; 30,000 epochs; soft terminal weight 10; 4x128 ELU | Identical source profile with only equation/horizon/dimension tensors changed | Faithful method transfer |
| SDGD HJB | 3x1024 tanh hidden layers; hard terminal; 10,000 Adam epochs; linear `1e-3 -> 0`; 100 residual points; PGD-20; 100 training coordinates and 10 adversarial coordinates | Same; model/optimizer/scheduler/collocation/RNG checkpointing only | Exact source algorithm/configuration; numerical reproduction is CONDITIONAL because of the endpoint accuracy gap |
| SDGD generic semilinear | 3x1024 tanh; 10,000 Adam epochs; exponential gamma `0.9995`; 100 trajectories x 20 time steps; two independent 10-coordinate trace estimators at every time slice; terminal value/gradient weights 20/1 | Same, including fresh coordinate subsets at each time slice | Faithful transfer; no HJB adversarial profile or hard terminal is imported |
| Deep BSDE HJB / Allen | scalar trainable `Y0`, free `Z0`, 19 nonshared time networks, two width-110 hidden layers, 20 intervals, float64, clipped terminal loss; 2,000 / 4,000 optimizer steps | Audited PyTorch compatibility port with the same architecture, schedules, batch 64, validation 256, and update count | Framework compatibility change only; numerical native reproduction is credible |

The common function test law is independent of training and frozen before the
new runs. Native-paper laws and metrics are also retained separately. No
external baseline receives Tangent-Actor's terminal residual, tangent update,
martingale control variate, rank projection, oracle coordinates, or analytic
nonterminal labels.

## Corrected fidelity defects and preserved negative results

The review found two material mismatches before the corrected metrics were
admitted to the primary table:

1. The earlier DPI/PINN transfer applied the OU/GMM-specific PISGradNet profile
   to every PDE. Corrected runs use that profile only on OU/GMM and use the
   paper's generic semilinear profile elsewhere.
2. The earlier SDGD transcription used ten trace coordinates for both HJB
   training and PGD. The paper and pinned HJB source prescribe 100 coordinates
   for training and ten for PGD. The generic Section-5.2 solver is a distinct
   soft-terminal SDE-trajectory protocol, not the HJB adversarial solver.

All earlier files and checkpoints are preserved. Their rows are labeled
`diagnostic_legacy_hjb_profile` or
`diagnostic_invalid_coordinate_budget` in
`HISTORICAL_AND_NONPRIMARY_RESULTS.csv`. They cannot support a Tangent-Actor
ranking. The pinned HJB source uses two independent coordinate sets for its
`method=2` estimator even though the final paper calls the method Algorithm 3;
we preserve the pinned source semantics and disclose the paper/source naming
difference.

## Compatibility changes

- PyTorch 2.6 changed `torch.load` to `weights_only=True`; DPI restores the old
  trusted-checkpoint behavior for checkpoints generated locally by the pinned
  source. No tensor or numerical update changes.
- Rich terminal progress and online W&B finalization are disabled under Slurm.
  Evaluation callbacks, optimization, data, and checkpoints are unchanged.
- Deep BSDE moves from unavailable TensorFlow/Keras to the previously audited
  PyTorch port. The BSDE recursion, time-specific networks, batch
  normalization, loss, dtype, and schedules are unchanged.
- SDGD adds exact resume state because its monolithic official loop exceeds
  practical queue windows. Model, optimizer, scheduler, collocation points,
  and all RNG states are restored at epoch boundaries.
- The nonsmooth MLP-failure transfers use PyTorch's minimum-norm subgradient at
  zero and are classified `adapted`; they are not smooth-driver theorem tests.

## Convergence and undertraining evidence

`BASELINE_CONVERGENCE.csv` and
`fidelity_figures/baseline_native_convergence.pdf` collect the complete native
seed-0 histories. DPI decreases through its Picard rounds. PINN-HTE finishes
above the paper's plotted accuracy and remains CONDITIONAL rather than being
used as evidence of method superiority. Deep BSDE's HJB and Allen terminal
validation losses reach stable low plateaus at the fixed source budgets.

The corrected SDGD 250D run completed all 10,000 epochs. Its online reporting
relative error moved from `8.555e-2` at epoch 0 to `9.192e-3` at epoch 9,900;
the independently recomputed final-bank relative L2 is `9.072e-3`. The curve
is noisy and attains lower intermediate values, but the protocol uses the
published fixed endpoint and never selects a checkpoint by test accuracy.
The `204.39` minute measured training time agrees with the paper's `210`
minutes, so under-budget training does not explain the remaining accuracy gap.

## Completed 100D whole-function results

Every entry below is the fixed seed-0 endpoint. Cells show `Value RMSE /
relative L2 / training seconds`. DPI and PINN-HTE use their specialized native
OU/GMM profile only on OU/GMM and the paper's generic semilinear profile
elsewhere. SDGD uses the HJB adversarial profile on HJB-Log and OU/GMM and the
Section-5.2 SDE profile on the other smooth PDEs. These transfers receive no
paper-reproduction verdict.

| PDE | DPI | PINN-HTE | SDGD-PINN |
|---|---:|---:|---:|
| HJB-Log | `0.133935 / 0.029215 / 447.5` | `0.053621 / 0.011696 / 303.5` | `0.039884 / 0.008700 / 12069.3` |
| OU/GMM | `17.6638 / 0.056972 / 627.4` **native** | `38.6698 / 0.124724 / 767.2` **native** | `62.4838 / 0.201533 / 14701.3` |
| Joint concave | `0.002142 / 0.003606 / 1244.3` | `0.003226 / 0.005433 / 507.1` | `0.023332 / 0.039288 / 5396.1` |
| Joint convex | `0.001857 / 0.004238 / 1234.0` | `0.003267 / 0.007458 / 500.0` | `0.012780 / 0.029173 / 5485.1` |
| Nonconvex coupled | `0.001699 / 0.003621 / 4011.6` | `0.003068 / 0.006539 / 1240.8` | `0.008877 / 0.018918 / 6280.9` |
| MLP failure | `0.010021 / 0.009241 / 261.4` **adapted** | `0.039939 / 0.036830 / 317.3` **adapted** | `0.066871 / 0.061665 / 5625.0` **adapted** |
| Allen-Cahn | `0.001067 / 0.023284 / 414.5` | `0.000948 / 0.020696 / 280.8` | `0.002692 / 0.058777 / 5798.6` |

The OU/GMM SDGD transfer is especially poor on the distributed test law even
though its native-point error is small. The HJB source profile samples its
source-prescribed collocation law, which differs substantially from the broad
OU/GMM evaluation law; this is reported as a transfer limitation rather than
evidence that the published SDGD method fails on its own benchmark.

Deep BSDE remains separate because it solves only one initial state. Its 100D
absolute errors are `0.010304` on HJB-Log and `0.000372` on Allen-Cahn
(native), then `71.0813`, `0.000381`, `0.001252`, `0.002180`, and `0.066461`
on OU/GMM, joint concave, joint convex, nonconvex coupled, and MLP failure.
The final MLP-failure entry is `adapted`; the other five nonnative entries are
faithful single-point method transfers.

## Unsupported or incomplete comparisons

- Deep BSDE has no whole-function row because the original method learns one
  scalar initial value per training run.
- A source-native SDGD HJB result at 100D does not exist in Table 5; suite 100D
  is a faithful dimension transfer. The separate 250D run is the native
  fidelity anchor.
- Generic DPI/PINN/SDGD transfers have no published target-PDE number and are
  judged for method fidelity and convergence, not paper-accuracy reproduction.
- The MLP-failure norm driver is nonsmooth at zero. All neural baseline results
  there are explicitly adapted and do not establish smooth-PDE convergence.

## Review gate

**Current state: PENDING INDEPENDENT REVIEW.** The corrected 100D tables,
convergence evidence, and compute estimate are complete. The post-gate sweep
contains 175 method cells and is estimated at `139.53` GH200 GPU-hours, with a
wide extrapolation interval of `81.67--284.82` GPU-hours. No dimension launcher
may run unless `REVIEW_GATE_STATUS.json` is explicitly changed to `approved`.
