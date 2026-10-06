# Training-data-only enrichment: results

## Outcome

Doubling the training set with six stronger two-dimensional families did **not**
improve the retained HCFL method on the locked test set.  The enriched model is
safe and converged, but its held-out rollout NMAE is 3.20% higher than the
original-data model.  This enrichment should therefore not replace the
original 192-trajectory training set.

## Controlled comparison

The architecture, 7,348 trainable parameters, shared x/y six-cell stencil,
central physical flux, nonnegative Roe dissipation, proposal-feasibility loss,
hard interface projection, fixed Roe transverse transport, optimizer, loss
weights, update cap, validation rule, and 64 by 64 deployment grid were held
fixed.  Input normalization, conserved-variable loss scales, and NMAE scales
were also frozen to the original training archive.

The only changed factor was the sampled training archive:

- original: 192 trajectories, 32 from each of six families;
- enriched: all original trajectories plus 32 crossed-Riemann, 32 corrugated
  contact, 32 shock/contact-triple, 32 shock-vortex, 32 double-blast, and 32
  four-way-collision trajectories, for 384 total;
- all added targets: PyClaw 5.10 Roe with `transverse_waves=2` at 512 by 512,
  conservatively restricted to 64 by 64;
- validation and test archives: unchanged byte-for-byte.

This experiment tests the combined enrichment.  It does not identify data
quantity and the six new families as separate causal factors.

## Validation convergence

All three runs stopped by the preregistered validation-plateau rule at the
minimum learning rate; none reached the 50,000-update cap.

| Seed | Best update | Best validation NMAE | Stop update |
|---:|---:|---:|---:|
| 0 | 6,000 | 0.018327 | 15,000 |
| 1 | 6,500 | 0.018130 | 15,500 |
| 2 | 1,500 | 0.018426 | 12,000 |

## Locked held-out test

Metrics below average 36 test cases and three validation-selected seeds (108
HCFL rollouts per training archive).  Lower is better except that signed TV
error is ideally near zero.

| Metric | Original 192 | Enriched 384 | Relative change |
|---|---:|---:|---:|
| NMAE | 0.014508 | 0.014972 | +3.20% |
| Density MAE | 0.008597 | 0.008754 | +1.83% |
| x-velocity MAE | 0.005192 | 0.005429 | +4.56% |
| y-velocity MAE | 0.004884 | 0.005090 | +4.22% |
| Pressure MAE | 0.005929 | 0.006035 | +1.79% |
| Signed final density-TV error | +11.07% | +13.25% | farther from zero |

The per-case, seed-averaged NMAE difference (enriched minus original) is
+0.000465.  A 20,000-resample paired bootstrap gives a 95% interval of
[+0.000040, +0.000923], and the enriched model wins only 12 of 36 cases.
Thus the overall degradation is not explained by a single unlucky seed.

By family, the enriched model improves colliding-waves NMAE from 0.015672 to
0.015161 (3.26% lower), but worsens the other five families.  The largest
relative loss occurs on the smooth packet (0.000649 to 0.001701).  Contact/
shear also becomes more oscillatory: its signed density-TV error grows from
+24.05% to +38.21%.  The added interactions therefore help one targeted regime
while diluting or conflicting with accuracy on the original distribution.

For context, unchanged PyClaw Roe-64 obtains NMAE 0.011282 on the same test
archive, so enrichment does not close the remaining gap to that comparator.

## Constraint and numerical audit

All 108 enriched-model test rollouts complete.  Across the three seeds:

- minimum density is 0.113832 and minimum pressure is 0.074297;
- the positivity fallback is never activated;
- maximum fully discrete entropy balance is 4.804e-7, within the implemented
  5e-7 numerical tolerance;
- maximum hard-projection interface residual is 3.501e-7;
- maximum conservation closure is 2.897e-9.

The data audit verifies finite values, positive density and pressure, the exact
original-data prefix, fixed saved times, and unchanged validation/test hashes.
All 19 numerical-invariant regression tests pass.

## Decision

The narrow hypothesis "the current error is mainly caused by too few or too
simple training trajectories" is not supported by this equal-weight doubling.
The original checkpoint remains the better main result.  A future data study
would need to isolate quantity from family mix, for example with a fixed-size
replacement ablation or a validation-prespecified stratified sampler; doing so
would be a new experiment, not a reinterpretation of this result.

## Artifacts

- `results/AUDIT.json`: hashes, checks, paired bootstrap, MAE, NMAE, and TV.
- `results/validation_data_ablation.png`: all validation curves and selected
  checkpoints.
- `results/test_nmae_data_ablation.png`: overall and family NMAE.
- `results/test_mae_data_ablation.png`: dimensional primitive-variable MAE.
- `results/test_tv_data_ablation.png`: signed density-TV error.
- `results/final_time_density_heatmaps.png`: representative final density.
- `results/final_time_density_error_heatmaps.png`: matched density-error maps.
- `results/enriched_training_examples.png`: the six added training families.
