# Baseline source and fidelity audit

| Method | Pinned official source | Native configuration | Our intervention |
|---|---|---|---|
| Deep BSDE | `frankhan91/DeepBSDE` at `e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` | `hjb_lq_d100.json`, `allencahn_d100.json` | Audited PyTorch 2.10 compatibility port because TensorFlow is absent on the GH200 image |
| Deep Picard Iteration | `DeepOptimalControl/DeepPicardIteration` at `e8616c2f531c9f95299b246a13c1180b8f02822a` | `scripts/hjb/base_100d_T1.0_w0.1_0.yaml` | Fixed seeds, wall/memory logging, final state export |
| PINN-HTE | Same DPI commit | `pinn_100d_T1.0_v16.yaml` | Fixed seeds, wall/memory logging, final state export |
| D-DBSDE | Same DPI commit | `diffusion_100d_T1.0.yaml` | Fixed seeds, wall/memory logging, final state export |
| SDGD-PINN | `zheyuanhu01/SDGD_PINN` at `c597cc4fe9caac0b3ca118ea84453405f25ff833` | `Section_5.3/HJB_Log.py` defaults | Supply the omitted test file from an independent exact reference; execute source unmodified with `runpy`; export final state afterward |

## Deep BSDE port

The port retains the native scalar free `Y0`, free `Z0`, 19 nonshared
time-slice networks, input/hidden/output batch normalization, two 110-unit ReLU
hidden layers, division of learned time-slice `Z` by dimension, 20-step Euler
BSDE, clipped quadratic terminal loss, float64 arithmetic, Adam epsilon,
published learning rates, batch sizes, validation sizes, and 2,000/4,000-step
schedules.  It does not amortize `Y0(x)`, add `g` or `grad g`, constrain values,
or use MCV.  The framework change and all measured outputs are labelled as a
compatibility port rather than an untouched TensorFlow run.

## DPI fixed instance recovery

The official HJB directory tracks
`mean_100d_ms=1.0_vs=2.0_5.pt` and
`pi_100d_ms=1.0_vs=2.0_5.pt`, but not the variance tensor that its constructor
loads in the same `try` block.  Running it as checked out would enter the
fallback and randomly replace the tracked mean and weights.  We preserve the
tracked tensors byte-for-byte and create five `2I` covariance matrices, exactly
the `var_scale=2` fallback specified by the source.  SHA-256 hashes and the
reason are in `prepared_inputs_manifest.json`.

The DPI family keeps its native PISGradNet terminal construction.  That is an
official method feature and is not removed merely because our Actor also has a
terminal construction.  No MCV or Tangent-Actor feature is inserted.

## SDGD test data recovery

The official script requires `HJB_Log_100_grad.pkl`, but the file is absent
from the repository.  It enters only periodic L2 reporting, never the training
loss.  We generate 20,000 states from the script's stated distribution and
label them with independently converged Hopf-Cole quadrature.  The first 1,200
states are the frozen common bank, which permits an exact shared 960-state
test subset.  The official script still performs its native per-epoch test
reporting, and that cost is charged.

