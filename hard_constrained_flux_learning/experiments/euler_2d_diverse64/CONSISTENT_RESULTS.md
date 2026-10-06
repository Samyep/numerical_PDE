# Method-consistent fixed-64 2-D Euler result

This is the controlled test of the direct 2-D extension of the retained 1-D
HCFL flux:

```text
central physical flux
+ nonnegative entropy-fixed Roe dissipation multipliers
+ proposal feasibility loss
+ hard interface entropy projection.
```

The comparator is the earlier `flat18` HLLC plus signed-Roe-correction
ablation.  The two neural methods use the same 18-cell input, 10,804
parameters, high-resolution data, train/validation/test splits, minibatch
stream, optimizer and stopping rules, trajectory substeps, feasibility
weight, hard projection, and deployment safety wrapper.  The controlled
factor is the raw flux parameterization.

## Protocol

- Reference: PyClaw 5.10 four-wave 2-D Euler Roe, `transverse_waves=2`, on a
  `512 x 512` grid and conservatively restricted to `64 x 64`.
- Data: 192 train, 36 validation, and 36 held-out test trajectories; six
  balanced initial-condition families and disjoint seeds.
- Horizon: `t=0.1`, saved every `0.005`.
- Training: six differentiable PDE substeps per saved interval, three model
  seeds, 50,000-update cap, and validation-plateau learning-rate reductions.
- Selection: physical completion first and validation rollout NMAE second.
  The test split was first evaluated only after all checkpoints were frozen.
- Data hashes: train
  `a99ffa2bae301902a14792406dcbbf0b599c375a4989580e19e7ca05fa559f4f`,
  validation
  `7e449077c17a113b8930abed3bb1a249ea147a525ab9023a618b810236f3225c`,
  and test
  `ed877ba706a48d10266790da9919cf2d5b102bb55dd9a5cc696871aaad8bdea0`.

All three runs stopped at update 12,000 with
`validation_plateau_at_minimum_learning_rate`:

| seed | best update | best validation NMAE |
|---:|---:|---:|
| 0 | 500 | 0.0180050 |
| 1 | 1,000 | 0.0177766 |
| 2 | 1,000 | 0.0179150 |

The full curves are in
[`results/consistent_validation_convergence.png`](results/consistent_validation_convergence.png).

## Held-out accuracy

All values below use every saved rollout time and every one of the 36 unseen
test trajectories.  Neural values average the three independently trained
seeds.

| method | rollout NMAE |
|---|---:|
| PyClaw Roe-64 | **0.0112824** |
| Central + nonnegative Roe-18 | 0.0147320 |
| HLLC + signed Roe-18 | 0.0151796 |
| HLLC-64 | 0.0237910 |

The method-consistent model improves NMAE by **2.95%** relative to the old
signed-correction model, but remains **30.6% worse** than PyClaw Roe-64.  It
therefore does not support a claim that HCFL beats the strong same-grid Roe
solver.

The paired 36-case comparison against the signed model gives a mean NMAE
difference of `-4.476e-4`, a fixed-seed bootstrap 95% interval of
`[-8.570e-4, -6.305e-5]`, and 21/36 case wins.  Against Roe-64 the difference
is `+3.450e-3`, the interval is `[+2.283e-3, +4.717e-3]`, and HCFL wins only
8/36 cases.  Against HLLC-64 it wins 35/36 cases.

Component MAE shows that the aggregate improvement is not uniform:

| method | density MAE | x-velocity MAE | y-velocity MAE | pressure MAE |
|---|---:|---:|---:|---:|
| PyClaw Roe-64 | **0.007118** | **0.004035** | **0.003759** | **0.004259** |
| Central + nonnegative Roe-18 | 0.008614 | 0.005342 | 0.004981 | 0.005984 |
| HLLC + signed Roe-18 | 0.008592 | 0.005577 | 0.005165 | 0.006268 |
| HLLC-64 | 0.013758 | 0.008761 | 0.008039 | 0.009598 |

Relative to the signed model, the new parameterization improves the two
velocity MAEs and pressure MAE by about 4.2%, 3.6%, and 4.5%, respectively,
while density MAE is about 0.25% worse.  Family-level NMAE is lower for five
of six families and marginally higher for colliding waves.  See
[`results/consistent_test_nmae_by_family.png`](results/consistent_test_nmae_by_family.png).

## Entropy, positivity, and conservation audit

- Every test trajectory completed for every seed; the minimum observed
  density and pressure were `0.1357` and `0.0740`.
- Positivity fallback activations: **0**.
- Fully discrete entropy fallback activations: **2 / 8,852** evaluated batch
  substeps (`0.0226%`).
- Maximum interface projection residual: `1.424e-7`.
- Maximum relative conservation closure: `3.386e-9`.
- Maximum fully discrete entropy balance: `-1.502e-5` (admissible).
- The unprojected proposal violates the audited interface entropy inequality
  on 5.3--5.9% of interfaces, so the hard projection is materially active;
  the proposal alone is not entropy admissible.

This supports the narrow claim that the deployed update is conservative,
positive on this test set, and entropy admissible within the audit tolerance.
It does not establish positivity without the deployment safety wrapper for
arbitrary states.

## What the network learned

The learned Roe multipliers obey their design range: global observed minimum
`0.00324`, maximum `1.97623`.  However, their per-seed means are only
`0.182--0.209`, and almost every multiplier differs from one.  The network is
therefore removing roughly 80% of standard Roe dissipation, not learning a
small perturbation of Roe.

The relative response to the transverse rows is only `0.84--1.17%`.  Despite
receiving 18 cells, this network has mostly learned to ignore transverse
information rather than reproduce PyClaw's transverse-wave coupling.

## Oscillation finding and decision

The signed final density-TV error is **+12.41%** for the new method versus
**+11.89%** for the old signed model.  The fixed, preselected heatmaps show
parallel ringing around oblique Riemann, quadrant, and colliding-wave
features for both neural models:

- [`results/consistent_final_time_density_heatmaps.png`](results/consistent_final_time_density_heatmaps.png)
- [`results/consistent_final_time_density_error_heatmaps.png`](results/consistent_final_time_density_error_heatmaps.png)
- [`results/consistent_test_density_tv_by_family.png`](results/consistent_test_density_tv_by_family.png)

Therefore nonnegative Roe multipliers plus feasibility loss are a modest and
statistically supported accuracy improvement over the signed correction, but
they do **not** solve the 2-D oscillation problem and do **not** beat Roe-64.
Nonnegativity prevents a wave from becoming explicitly antidiffusive, yet a
multiplier near zero can still remove nearly all stabilizing dissipation.
Hard entropy projection enforces entropy admissibility but is not a TVD or
monotonicity constraint.

The next controlled design should target the two diagnosed failures rather
than merely enlarge the same MLP: impose a meaningful dissipation floor or a
shock-sensor-dependent lower bound, and introduce an explicit transverse
fluctuation/corner-coupling path with its own supervised or trajectory signal.
The present model and plots should be retained as the method-consistent
ablation, not promoted as the final 2-D result.

Machine-readable details are in
[`results/consistent_AUDIT.json`](results/consistent_AUDIT.json) and
[`results/consistent_test_summary.json`](results/consistent_test_summary.json).
