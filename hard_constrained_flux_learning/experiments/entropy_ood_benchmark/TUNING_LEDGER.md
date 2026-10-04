# Tuning and implementation ledger

This file records successful, failed, interrupted, and corrective runs.  Test
results were not used to choose a checkpoint or hyperparameter.

## Euler stress screen

- Existing HLLC, MUSCL-HLLC, residual FNO, and three HCFL checkpoints were
  evaluated without retraining in the initial exploratory screen.  Four of the
  five named cases had appeared in earlier local work, so this is not labelled
  a blind or preregistered test.
- The 64-cell HCFL completed all five cases and beat both 64-cell finite-volume
  baselines in mean error.  The recorded condition for increasing the
  training grid was therefore not met; no higher-grid rescue run was made.

## RoeNet adaptation

- A first optimization attempt was interrupted at update 3500 after detecting
  that checkpoint selection applied a relative-error threshold even when the
  validation completion count changed.  That logic could prefer an error
  improvement over physical completion and was rejected before stress-test
  evaluation.
- An independent audit found that the first claimed lexicographic
  implementation still applied a relative-improvement threshold to a scalar
  score near 100.  This could skip a better error at identical completion.
- The corrected run uses an exact tuple ordering: validation completion first,
  completed-only NRMSE second.  The relative threshold affects plateau
  scheduling only, never checkpoint retention.
- The corrected run selected update 2000 (NRMSE `0.0490284`) and stopped at
  update 12000 after four learning-rate reductions.  Validation completion was
  134/136.  The two failures remain in the denominator.
- The adaptation uses the official 64-wave architecture idea but matched HCFL
  data, normalization, and a regularized inverse.  It is never labelled an
  exact paper reproduction.

## Official data-free PINN

- No PINN hyperparameter was tuned locally.  The authors' Sod and Lax L-NN2
  state dictionaries are loaded directly.
- Continuous predictions are integrated into cell averages with fixed
  order-eight Gauss-Legendre quadrature.  Center sampling was not used.
- These are native per-instance checkpoints; no unsupported zero-shot OOD
  result is reported.

## 2-D SWE data and smoke checks

- The split sizes, seeds, OOD intervals, horizons, task geometry, and metric
  order were written to `PROTOCOL.md`/`protocol.json` before model training.
- Pre-result smoke checks covered constant-state preservation, x/y flux shapes,
  Roe reconstruction, hard projection, reference positivity, and official
  operator input/output shapes.
- Short 20-update and 100-update HCFL smoke runs were used only to catch
  runtime errors and estimate cost.  Their checkpoints were overwritten and
  are absent from all result tables.
- A two-epoch FNO/clawFNO smoke run likewise checked the official source
  interface; both checkpoints were overwritten by the formal run.

## Corrective implementation findings

1. During the first full four-cell HCFL run, a validation check exposed a
   `1.44e-6` fully-discrete entropy excess after float64 beta bisection was cast
   back to float32.  The run was stopped.  A rounded-state postcheck was added:
   any remaining violation sets that batch member to the HLL endpoint.  The
   optimizer path and data were unchanged, and the formal run restarted from
   update zero.
2. The first locked-test evaluation stopped before writing results because the
   black-box conservation audit summed `(X,channel)` instead of `(Y,X)`.  The
   diagnostic axes were fixed; checkpoints and predictions were not changed.
3. Saved-frame trapezoidal boundary diagnostics are approximate and can show a
   positive coarse-frame entropy balance even when every HCFL internal update
   passes.  The report therefore separates exact flux-level safety metrics
   from black-box saved-frame diagnostics.
4. Independent review found that entropy had been evaluated after an Euler or
   SWE trajectory was already nonphysical.  Such values are thermodynamically
   undefined.  Failed trajectories now keep their failure/minimum-state and
   conservation diagnostics, while entropy and TV/error quantities are `NA`.
5. Operator selection was made explicitly completion-first.  The FNO winner
   completes 20/20 validation trajectories.  clawFNO completes 0/20 at every
   epoch; epoch 258 is retained only as an inadmissible diagnostic checkpoint,
   not reported as a physical validation winner.

## Formal 2-D HCFL selection

Both arms use width 72, feasibility weight `1e-3`, Adam learning rate `5e-4`,
the same trajectories, and the same deployment wrapper.  HLL is absent from
the training loss and appears only as the safe deployment endpoint.

| Arm | Parameters | Best update | Stop update | Validation completion | Best validation NRMSE |
|---|---:|---:|---:|---:|---:|
| 4-cell | 6,411 | 2,000 | 5,000 | 20/20 | 0.189533 |
| 6-cell | 6,843 | 2,000 | 5,000 | 20/20 | **0.188898** |

Both stopped at a validation plateau after reaching the minimum learning rate;
neither hit the 8,000-update run cap or the global 50,000-update protocol cap.
The six-cell arm was selected using validation before evaluating the test splits.  Its
advantage over four cells is only about 0.3%, so this is not presented as a
statistically significant stencil ranking.

## Formal neural-operator selection

- Vanilla FNO: best epoch 477, completed all 500 cosine-schedule epochs,
  validation primitive relative L2 `0.0153913`.
- clawFNO: no physically admissible validation checkpoint (0/20 completed at
  every epoch).  The diagnostic epoch 258 stopped at epoch 358 by fixed
  patience; its all-prediction number `162.911` is not a physical error metric.
- Both use 2,464,683 parameters, the official source classes, and the published
  radial-dam modes/width/lr/weight decay.  The same local train/validation data
  are used, so both remain labelled official-architecture adaptations.

## Higher-grid rule

No higher-grid HCFL was launched.  Euler-64 met both accuracy and safety goals.
The 2-D operator benchmark is defined at the official 32-grid deployment
resolution and HCFL halves HLL error on ID and doubled-horizon tests.  Raising
only HCFL's deployment grid after seeing the 32-grid results would break
the matched comparison.  A future 64-grid 2-D study should be committed as a
new locked design with all baselines rerun at the same resolution.
