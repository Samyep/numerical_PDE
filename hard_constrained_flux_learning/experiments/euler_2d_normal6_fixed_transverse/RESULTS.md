# Results: six-cell normal HCFL with fixed transverse transport

## Verdict

The proposed dimension-consistent construction is numerically viable and
passes the complete safety audit.  Fixed transverse Roe transport reduces
excess density total variation, especially for oblique Riemann and contact /
shear data, but it does **not** materially improve aggregate NMAE over the same
trained checkpoint with transverse transport disabled.  It also does not beat
the native PyClaw Roe-64 baseline.

This is therefore a useful method and ablation result, not evidence that the
current 2-D HCFL is the final accuracy winner.

## Protocol

- Learned proposal: central physical flux plus nonnegative entropy-fixed Roe
  dissipation, using six normal cells and one network shared by x and y faces.
- Fixed component: Clawpack-style transverse increment-wave transport
  (`transverse_waves=1`) with zero trainable parameters.
- Constraint: proposal feasibility loss during training and hard interface
  entropy projection before every PDE update.
- Data: 512 by 512 reference trajectories conservatively restricted to 64 by
  64; 192 training, 36 validation, and 36 disjoint test trajectories across
  six initial-condition families.
- Model size: 7,348 trainable parameters, versus 10,804 for the earlier
  18-cell learned-transverse-context model.
- Seeds: 0, 1, and 2.  The cap was 50,000 optimizer updates; all runs stopped at
  update 12,000 because validation had plateaued at the minimum learning rate.
  The selected checkpoints occurred at updates 1,500, 1,500, and 1,000.

## Held-out accuracy

All values are means over the held-out trajectories; neural-model values also
average the three seeds.

| Method | NMAE | Density MAE | x-velocity MAE | y-velocity MAE | Pressure MAE | Signed density-TV error | Completion |
|---|---:|---:|---:|---:|---:|---:|---:|
| PyClaw Roe-64 | 0.011282 | 0.007118 | 0.004035 | 0.003759 | 0.004259 | -3.64% | 100% |
| HLLC-64 | 0.023791 | 0.013758 | 0.008761 | 0.008039 | 0.009598 | -11.43% | 100% |
| HCFL-18 learned transverse context | 0.014732 | 0.008614 | 0.005342 | 0.004981 | 0.005984 | +12.41% | 100% |
| HCFL-6, same checkpoint with fixed transverse disabled | 0.014507 | 0.008509 | 0.005247 | 0.004897 | 0.005901 | +13.39% | 100% |
| **HCFL-6 + fixed transverse Roe** | **0.014508** | **0.008597** | **0.005192** | **0.004884** | **0.005929** | **+11.07%** | **100%** |

Relative to the earlier 18-cell model, the new complete method has 1.52% lower
mean NMAE and wins 28 of 36 paired cases.  The paired bootstrap 95% interval
for the NMAE difference is `[-4.88e-4, 1.74e-5]`, however, so this small gain
is not conclusive.  Relative to PyClaw Roe-64, its NMAE is 28.59% higher and it
wins only 8 of 36 paired cases.  It is substantially better than HLLC-64,
winning 35 of 36 paired cases.

## Causal transverse ablation

The on/off comparison below uses the **same trained weights**; only the fixed
transverse term is disabled at inference.  It therefore isolates the effect of
the fixed transport more cleanly than the comparison to the 18-cell model.

| Family | NMAE without transverse | NMAE with fixed transverse | Relative NMAE change | TV error without | TV error with |
|---|---:|---:|---:|---:|---:|
| All | 0.0145071 | 0.0145077 | +0.004% | +13.39% | +11.07% |
| Oblique Riemann | 0.0160575 | 0.0159921 | -0.41% | +30.54% | +20.76% |
| Contact / shear | 0.0057072 | 0.0052768 | -7.54% | +25.29% | +24.05% |
| Oblique quadrant | 0.0334193 | 0.0337865 | +1.10% | +4.80% | +4.24% |
| Radial interface | 0.0154955 | 0.0156702 | +1.13% | +1.77% | +1.34% |
| Colliding waves | 0.0156690 | 0.0156722 | +0.02% | +18.01% | +16.16% |
| Smooth packet | 0.0006944 | 0.0006487 | -6.59% | -0.09% | -0.14% |

The aggregate paired NMAE difference is `5.82e-7`, with bootstrap 95%
interval `[-1.51e-4, 1.39e-4]`.  Thus the fixed term's defensible benefit in
this run is lower roughness, not an aggregate accuracy gain.

## Constraint and safety audit

All automated checks passed:

- 36/36 test trajectories completed for every new-model seed.
- Minimum observed density and pressure were 0.1553 and 0.0731.
- The positivity fallback was never activated.
- The fully discrete entropy fallback was activated in 32 of 8,844 batch
  substeps (0.362%); the maximum recorded entropy balance remained negative,
  at `-8.58e-7`.
- Maximum hard-projection interface residual: `4.29e-8`.
- Maximum relative conservation closure: `2.30e-9`.
- The transverse solver has zero trainable parameters.
- Learned Roe multipliers stayed in `[0.00119, 1.97869]` and the learned flux
  was not the standard Roe flux.

The hard interface projection is a structural guarantee of the implemented
face proposal.  The fully discrete entropy and positivity statements above
are empirical results on this held-out suite plus the deployment safety
wrapper; they are not presented as an unconditional theorem for arbitrary
states and time steps.

## Interpretation

The experiment supports the design choice requested here: a 2-D extension can
retain the 1-D learned object exactly (six normal cells, shared x/y weights)
and add genuinely 2-D transport as a fixed, parameter-free numerical module.
The exact 1-D reduction and x/y rotation symmetry are covered by tests.

It also identifies the remaining issue.  The fixed transverse term is only
about `3.5e-4` of the normal-flux L1 magnitude, while the learned mean Roe
multipliers are 0.158--0.206 instead of the standard Roe value 1.  Together
with the remaining +11.07% density-TV excess, this is consistent with the
normal learned dissipation being too small in some wave families.  This is a
diagnosis to test next, not proof from this experiment alone.

## Artifacts

- `results/AUDIT.json`: machine-readable pass/fail audit and paired statistics.
- `results/test_case_metrics.csv`: all per-case metrics.
- `results/test_summary.json`: complete aggregate accuracy and safety data.
- `results/final_time_density_heatmaps.png`: deterministic first held-out case
  from each family (not hand-selected).
- `results/final_time_density_error_heatmaps.png`: corresponding absolute
  density errors.
- `results/test_nmae_by_family.png` and
  `results/test_density_tv_by_family.png`: family-wise accuracy and roughness.
- `results/validation_convergence.png`: three-seed validation curves.
