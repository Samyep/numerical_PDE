# 1D Euler reference-precision ablation

> Historical fixed-budget screen: the 1,100-update models were not shown to
> converge. Their metrics remain as an audit trail, but their weights are not
> retained; this comparison must be rerun with validation convergence before
> drawing a model conclusion.

## Pre-run hypothesis

> The remaining 1D Euler accuracy gap is partly caused by numerical diffusion
> and discretization error in the 512-cell Rusanov + SSP-RK2 trajectory
> teacher. Training the same direct-vector HLLC-HCFL model on 2048-cell
> references, conservatively restricted to the same 64-cell grid, should
> improve moderate-OOD and several canonical rollout errors when both arms are
> evaluated against 2048-cell references, without harming ordinary-ID accuracy
> or the existing hard safety guarantees.

This hypothesis and the decision rule below were recorded before running the
experiment.

## Single changed scientific factor

Only the fine reference grid used to generate training targets changes:

- `reference_512`: 512-cell Rusanov + SSP-RK2 teacher;
- `reference_2048`: 2048-cell Rusanov + SSP-RK2 teacher.

Both are conservatively averaged to the same 64-cell learning grid. To make
the comparison paired, every initial condition is first generated on the
2048-cell grid. The 512-cell initial condition is its conservative four-cell
average. Consequently, both arms have the same 64-cell initial state exactly;
only the fine-grid evolution used to form later targets differs.

## Held fixed

- 580 independent training trajectories: 220 ordinary, 260 broad random, and
  100 random extreme;
- 16 snapshots and `dt_snapshot = 4e-4`;
- direct-vector HLLC correction architecture and width;
- model initialization and minibatch-index sequence;
- baseline-derived normalization for both arms;
- Adam settings, batch size, gradient clipping, and 1,100 update steps;
- hard Tadmor projection, local conservative admissibility limiter, and
  trajectory-wise fully-discrete entropy safeguard;
- 64-cell learning grid and float32 stored tensors;
- evaluation cases, horizons, error normalization, and safety metrics.

The common evaluation targets use the 2048-cell reference so that the
low-resolution arm is not evaluated against its own teacher.

## Evaluation suite

- 90 ordinary ID trajectories;
- 90 broad-random in-support trajectories;
- 90 moderate-OOD smooth trajectories with frequencies 4--6;
- Sod, Lax, collision, strong-pressure jump, and near-vacuum expansion.

The random suites retain the existing 16-snapshot horizon (`t_final=0.006`).
Canonical cases retain the existing 64-snapshot horizon (`t_final=0.0252`).

## Seed-0 screening rule

The 2048-cell arm passes the screen only if:

1. ordinary-ID and moderate-OOD NRMSE each degrade by no more than 5%;
2. at least three of five canonical cases improve;
3. mean canonical NRMSE improves by at least 5%;
4. no canonical case degrades by more than 10%; and
5. every rollout remains admissible with zero measured Tadmor violations.

Only a passing seed-0 result should receive multi-seed confirmation.

## Run

```powershell
python run_precision_ablation.py --seed 0
```
