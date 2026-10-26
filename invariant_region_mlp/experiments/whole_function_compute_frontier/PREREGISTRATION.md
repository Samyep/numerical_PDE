# Two-PDE 100D Accuracy--Compute Frontier (frozen design)

Frozen before any new frontier test error was computed. The authoritative
machine-readable version is `EXPERIMENT_PROTOCOL.json`, whose SHA-256 is
recorded in `configuration_manifest.json`.

## Scope

- PDEs: B6 Full-Rank HJB and Dense Fourier, dimension 100 only.
- Methods: Tangent-Actor (TA), DPI, and PINN-HTE.
- Training seed: 0.
- Evaluation: the existing frozen 960-state whole-function bank for each PDE.
- Hardware: NVIDIA H100 80GB for all new numerical training and evaluation.
- No new SDGD-PINN job may be submitted.
- Total allocation cap: 18 H100 GPU-hours, including smoke tests, failures,
  retries, evaluation, and engineering overhead.

The primary comparison is the observed 5--20 minute region. Genuine training
may continue to 30 active minutes. Times are measured, not nominal. Missing
times are not interpolated, and an old checkpoint is not repeated at a later
time to manufacture coverage.

## Frozen configurations

| Method | C0 | C1 | C2 |
|---|---|---|---|
| TA | 4 original actor updates, 128 paths/label, 5000 fit steps; at most 8 unchanged updates | 256 paths/label, 5000 fit steps, 4 updates | 128 paths/label, 10000 fit steps, 4 updates |
| DPI | original 16 epochs/Picard and 4096 integral samples; iteration 80 retained, at most 320 genuine rounds | 32 epochs/Picard, 4096 integral samples, 80 rounds | 16 epochs/Picard, 8192 integral samples, 80 rounds |
| PINN-HTE | HTE-16, batch 512, original epoch 30000 retained | HTE-32, batch 512 | HTE-16, batch 1024 |

All other original network, objective, optimizer, sampling, terminal-loss, and
PDE-specific mathematical choices are unchanged. C2's upstream keys are
`DATA.kwargs.n_estimate_integral` for DPI and `TRAIN.BATCH_SIZE` for PINN.
The proposed PINN batch must pass a truth-free H100 smoke test before formal
training; failure would require a documented preregistration amendment before
any test result is opened.

## Checkpoints and stopping

Observation targets are 2, 5, 8, 10, 15, 20, 25, and 30 active minutes.
Checkpoints occur only at algorithmically valid boundaries: completed TA actor
updates, completed DPI Picard rounds, or completed PINN optimizer steps. The
original four-update/80-round/30000-epoch checkpoints are always retained.

TA and PINN may stop after at least 15 active minutes when their last three
eligible independent validation measurements improve by less than 1% relative
to the first of those three. TA uses the frozen noise-corrected fixed-point
inconsistency. PINN uses a fixed 256-state truth-free bank and

`sqrt(mean(PDE residual^2) + 10 mean(terminal error^2))`

with 16 fixed Hutchinson probes. DPI has no source-compatible independent
validation stop and therefore uses its iteration limit or the 30-minute cap.
No stop or extension decision may use frozen test truth.

TA-C0's audited checkpoints and phase times may be reused because continuation
at actor boundaries is exact. DPI-C0 and PINN-C0 must be rerun as single
instrumented trajectories because their historical artifacts do not preserve
all state needed for exact continuation. Previously completed work is charged,
not treated as free.

## Truth gate and metrics

All 18 training trajectories must reach a terminal training status before a
global marker opens access to nonterminal frozen truth. Until then, evaluation
functions abort. Truth is forbidden from training, configuration choice,
selection, stopping, and extension.

Primary metric:

`RMSE = sqrt(mean((prediction - truth)^2))`.

Secondary metrics are `sum(abs(error))/sum(abs(truth))` (RelMAE) and relative
L2. TA evaluation is fixed at 1024 paths, 64 Euler steps, its existing
martingale control variate, and two independent inference replicates.

The main x-axis is synchronized cumulative active training plus causal
selection. Training includes TA label generation and fitting, DPI data
generation plus fitting, and PINN trace/residual/optimizer work. Inference,
job-wall, and allocated-GPU time are separate fields. C0 is primary; C1/C2 are
sensitivity analyses and are never collapsed into a test-selected envelope.

## Execution order and audit

First run a non-result smoke test. Then run all six C0 trajectories before the
twelve C1/C2 trajectories. Only after all 18 stop, open the truth gate and
evaluate every prespecified checkpoint, including poor ones.

The final audit checks monotone time, completion-before-timestamp semantics,
checkpoint and prediction hashes, inclusion of TA label cost, truth isolation,
original-budget C0 consistency, explicit numerical failures, scheduler
reconciliation, shared protocols, and absence of SDGD jobs. Historical C0
consistency is preregistered as

`abs(new RMSE - old RMSE) <= max(0.05, 0.50 * old RMSE)`.

Any failed audit remains a failure in the report.
