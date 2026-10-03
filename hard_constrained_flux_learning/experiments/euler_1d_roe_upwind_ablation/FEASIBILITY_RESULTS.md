# Central + nonnegative Roe: proposal-feasibility result

Seed: 0.  This experiment changes only the training objective of the existing
Central + nonnegative Roe model.  All three arms have 6,627 parameters and use
the identical hard-projected flux for every PDE update:

```text
F_update = hard_entropy_projection(F_raw)
loss = trajectory_loss + lambda_feas * mean(relu(r_raw)^2)
r_raw = (v_R - v_L)^T F_raw - (psi_R - psi_L)
```

No HLLC+Roe model was retrained or modified.

## Convergence and in-distribution feasibility

| lambda_feas | Best / stop update | Validation rollout NRMSE | Raw violation rate | Mean positive residual squared | Hard-projection rate |
|---:|---:|---:|---:|---:|---:|
| 0 | 14,800 / 16,300 | **0.0290998** | 4.0858% | 0.116126 | 3.8553% |
| 1e-4 | 14,800 / 17,800 | **0.0290703** | 3.8912% | 0.096673 | 3.6592% |
| 1e-3 | 9,800 / 13,800 | 0.0320910 | **3.5156%** | **0.087125** | **3.2844%** |

Both feasibility arms met the independent-validation plateau rule before the
50,000-update cap.  On teacher-forced 64-cell validation states, the auxiliary
loss behaves as intended.  Relative to the control, `1e-4` reduces the raw
violation rate by 4.76% and mean squared violation by 16.75%, while changing
validation NRMSE by only -0.10%.  The stronger `1e-3` penalty reduces those
feasibility metrics further but increases validation NRMSE by 10.28%.

## Held-out 64-cell canonical cases

| lambda_feas | Five-case mean rollout NRMSE |
|---:|---:|
| 0 | **0.0470929** |
| 1e-4 | 0.0508625 |
| 1e-3 | 0.0522568 |

The light penalty is already 8.00% worse than the control on the canonical
mean, despite matching the aggregate validation metric.

## Zero-shot 512-cell deployment

All checkpoints were trained on 64-cell trajectories and applied unchanged to
512 cells.  Metrics use strict native HLLC-2048 conservatively restricted to
512 only for scoring.

| lambda_feas | Completed | Mean rollout NRMSE | Mean final NRMSE | Projection rate | Mean TV ratio | Excess extrema |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 5/5 | **0.0497970** | **0.0619062** | **10.5247%** | **1.0430** | 74 |
| 1e-4 | 5/5 | 0.0599907 | 0.0735079 | 10.5611% | 1.0578 | 80 |
| 1e-3 | 5/5 | 0.0540494 | 0.0651243 | 11.9882% | 1.0482 | 73 |

Relative to the control, mean 512-cell rollout error rises by 20.47% for
`1e-4` and 8.54% for `1e-3`.  Neither penalty reduces the projection rate on
the autoregressive 512-cell trajectories.  The light penalty has the largest
near-vacuum degradation and more excess extrema; the stronger penalty is
smoother than the light arm but remains less accurate than the control.

## Conclusion

The proposed loss restores a raw entropy-normal training signal and measurably
makes teacher-forced 64-cell proposals more feasible.  That mechanism is real.
However, the gain does not transfer to autoregressive 512-cell states: rollout
projection reliance is unchanged or higher, and both tested weights reduce
cross-resolution accuracy.  The likely failure mode is state-distribution and
resolution shift, not loss nonconvergence.

As a strict loss ablation against the zero-penalty nonnegative control, the
result is negative for aggregate accuracy: neither penalty wins at 512 cells.
That is not the same as saying that the complete feasibility-trained method is
unusable.  The `lambda_feas=1e-3` model completes all five cases, has a mean
TV ratio of `1.0482` and 73 excess extrema, and visibly avoids the severe
oscillation of Central + signed Roe (mean TV ratio `1.6577`, 790 excess
extrema).  Its mean rollout/final NRMSE (`0.0540494` / `0.0651243`) is also
3.94% / 5.60% lower than HLLC + Roe correction in this seed.

We therefore retain **Central + nonnegative Roe + proposal-feasibility loss
(`lambda_feas=1e-3`)** as one of the two effective methods, alongside
**HLLC + Roe correction**.  The lighter `1e-4` arm is not retained.  This
classification records the good complete method without overstating causal
evidence: versus the zero-penalty nonnegative control, `1e-3` changes total
excess extrema only from 74 to 73 and slightly worsens aggregate TV/range
metrics.  More autoregressive and multi-resolution feasibility training would
be needed to show that the penalty itself supplies an additional general
improvement.

## Artifacts

- `results/upwind_feasibility512_seed0.png`
- `results/upwind_feasibility512_seed0.json`
- `results/upwind_feasibility512_seed0.csv`
- `results/scientific_integrity_audit_seed0.json`
- Two converged feasibility-loss checkpoints and their training curves
