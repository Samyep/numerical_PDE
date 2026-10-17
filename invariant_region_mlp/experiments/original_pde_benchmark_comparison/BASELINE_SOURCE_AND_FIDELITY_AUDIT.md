# Baseline source and fidelity audit

| Method | Pinned official source | Native configuration | Our intervention |
|---|---|---|---|
| Deep BSDE | `frankhan91/DeepBSDE` at `e76fed80b1995daf4a0dd6d3e2cf64a0931dabda` | `hjb_lq_d100.json`, `allencahn_d100.json` | Audited PyTorch 2.10 compatibility port because TensorFlow is absent on the GH200 image |
| Deep Picard Iteration | `DeepOptimalControl/DeepPicardIteration` at `e8616c2f531c9f95299b246a13c1180b8f02822a` | `scripts/hjb/base_100d_T1.0_w0.1_0.yaml` | Fixed seeds, wall/memory logging, final state export |
| PINN-HTE | Same DPI commit | `pinn_100d_T1.0_v16.yaml` | Fixed seeds, wall/memory logging, final state export |
| D-DBSDE | Same DPI commit | `diffusion_100d_T1.0.yaml` | Fixed seeds, wall/memory logging, final state export; execute the standalone solver once as documented below |
| SDGD-PINN | `zheyuanhu01/SDGD_PINN` at `c597cc4fe9caac0b3ca118ea84453405f25ff833` | `Section_5.3/HJB_Log.py` defaults | Supply the omitted test file from an independent exact reference; split its otherwise uncheckpointed loop at exact epoch boundaries; export final state afterward |

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

The D-DBSDE YAML inherits `PICARD.N=40` from the DPI base.  In the upstream
runner, the `Diffusion` branch trains a fresh standalone D-DBSDE network and
returns without updating Picard state; the outer loop therefore repeats the
entire 15,000-epoch run 40 times, discards the first 39 networks, and exports
only the last.  The repository README describes the YAML as one D-DBSDE run.
We set only `PICARD.N=1` for this launcher defect.  The PISGradNet, terminal
construction, loss, batch size, 50 SDE intervals, step size, optimizer, and
15,000 epochs are unchanged.  Both the inherited value and this correction
are recorded in every run artifact; no result is selected across repetitions.

The upstream custom Rich progress callback raises `IndexError` before the
first batch under this cluster's noninteractive Slurm stream.  The launcher
removes that display callback and sets Lightning's `enable_progress_bar=False`.
Evaluation callbacks, model summaries, TensorBoard/W&B logging, data, labels,
losses, precision, optimizers, and iteration counts remain unchanged.

The pinned repository also predates PyTorch 2.6.  Its Picard checkpoints store
the equation object and call `torch.load` without a `weights_only` argument.
PyTorch 2.6 changed that default to `True`, causing iteration 2 to reject the
locally generated iteration-1 checkpoint.  The launcher restores the former
`weights_only=False` default only when upstream omits the keyword.  This is a
checkpoint-format compatibility change and does not alter any numerical step.

## SDGD test data recovery

The official script requires `HJB_Log_100_grad.pkl`, but the file is absent
from the repository.  It enters only periodic L2 reporting, never the training
loss.  We generate 20,000 states from the script's stated distribution and
label them with independently converged Hopf-Cole quadrature.  The first 1,200
states are the frozen common bank, which permits an exact shared 960-state
test subset.  The official script still performs its native per-epoch test
reporting, and that cost is charged.

The monolithic 10,000-epoch script exceeds the queue's two-hour wall limit and
has no resume hook.  The production launcher executes the official class
definitions verbatim and transcribes only the body of `train_adam` into
checkpointable epoch chunks.  It persists the model, Adam and scheduler
states, current collocation points, recorded losses, and NumPy, Torch CPU, and
all CUDA RNG states.  Before the full run, an independent uninterrupted
100-epoch official execution and the checkpointable 100-epoch execution
produced bitwise-identical values for every parameter tensor (maximum absolute
difference exactly zero).  This is recorded in
`sdgd_resumable_equivalence_audit.json`; the split changes neither an update nor
the official terminal ansatz.

The script writes its trajectory as `.xlsx` after training.  The cluster base
environment lacked that optional writer, so `openpyxl==3.1.5` is installed in
an isolated dependency directory.  This changes only post-training file
serialization; the official source, model, optimizer, losses, samples, and
reporting frequency remain unchanged.

## Fixed-configuration reproduction status

The official configurations were run once for each amortized DPI-family
method and with three independent seeds for each native Deep BSDE problem.
These are fixed reproductions rather than baseline hyperparameter searches.

| Method and native task | Seeds | Accuracy | Mean training time |
|---|---:|---:|---:|
| Deep BSDE, HJB-Log at `u(0,0)` | 3 | absolute error `0.0096785 ± 0.0027611` | 65.15 s |
| Deep BSDE, Allen–Cahn at `u(0,0)` | 3 | absolute error `0.00016911 ± 0.00017982` | 127.69 s |
| DPI, OU/GMM common 960-state bank | 1 | Value RMSE `17.6638` | 627.42 s |
| PINN-HTE, OU/GMM common 960-state bank | 1 | Value RMSE `38.6698` | 767.21 s |
| D-DBSDE, OU/GMM common 960-state bank | 1 | Value RMSE `30.7247` | 453.05 s |

The uncertainty after `±` is the sample standard deviation across seeds, not
a confidence interval.  DPI, PINN-HTE, and D-DBSDE have one official run each;
their table entries must not be interpreted as seed-stability estimates.  The
formal comparison report also keeps DPI's distributed function task separate
from diagnostic values at `x=0`, which is not that benchmark's published
native task.

