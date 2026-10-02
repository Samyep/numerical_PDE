# Seed-0 result: validation-converged 1D Euler audit

## Outcome

All six arms met the pre-registered validation-plateau criterion before the
20,000-update cap. The old 1,100-update budget was not convergence for any arm:
the selected checkpoints occurred between 6,100 and 12,500 updates, and
stopping occurred between 10,000 and 16,500 updates.

The main conclusion changes. The fixed-budget evidence does not support claims
that Conv1d or wave-coverage training is intrinsically worse. The learned-
dissipation and characteristic parameterizations improve substantially when
given enough optimization and now lead the seed-0 comparison.

## Convergence summary

Checkpoint selection used only the independent common validation set. The five
canonical cases below remained final tests and were never used for scheduling,
selection, or stopping.

| arm | best / stop update | validation at 1,100 | best validation | validation change | canonical mean at 1,100 | canonical mean converged | canonical change |
|---|---:|---:|---:|---:|---:|---:|---:|
| direct / broad | 8,300 / 10,000 | 0.050717 | 0.049430 | -2.54% | 0.078174 | 0.088057 | +12.64% |
| invariant / broad | 8,300 / 11,500 | 0.049484 | 0.047696 | -3.61% | 0.080315 | 0.089859 | +11.88% |
| characteristic / broad | 10,400 / 13,900 | 0.043154 | 0.037819 | -12.36% | 0.084122 | 0.075569 | -10.17% |
| dissipation / broad | 12,500 / 16,500 | 0.043570 | **0.031557** | **-27.57%** | 0.080006 | **0.062575** | **-21.79%** |
| CNN / broad | 6,100 / 11,100 | 0.050852 | 0.049599 | -2.46% | 0.077533 | 0.079270 | +2.24% |
| direct / wave | 12,200 / 15,200 | 0.050816 | 0.049490 | -2.61% | 0.077122 | 0.088924 | +15.30% |

Negative change is an improvement. The canonical mean averages Sod, Lax,
collision, strong-pressure, and near-vacuum expansion.

## Architecture result at convergence

Relative to the converged direct/broad arm:

| arm | validation difference | canonical-mean difference |
|---|---:|---:|
| invariant / broad | -3.51% | +2.05% |
| characteristic / broad | -23.49% | -14.18% |
| dissipation / broad | **-36.16%** | **-28.94%** |
| CNN / broad | +0.34% | -9.98% |

The dissipation model is best on ordinary ID, broad-random ID, Lax, collision,
strong-pressure, near-vacuum, and the five-case canonical mean. The
characteristic model is narrowly best on moderate high-frequency OOD, while
the invariant model is narrowly best on Sod.

Therefore:

- the prior fixed-budget architecture ranking is retired;
- learned dissipation is the strongest seed-0 candidate and characteristic
  correction is second;
- CNN is effectively tied with direct on the selection metric (0.34% worse)
  and is 9.98% better on the canonical mean, so “CNN is worse” is unsupported;
- this remains a single-seed screening result and needs multi-seed confirmation.

## Wave-coverage result at convergence

With the direct architecture fixed, direct/wave differs from direct/broad by:

- +0.12% on validation rollout NRMSE;
- -0.52% on ordinary ID;
- -0.42% on broad-random ID;
- +0.58% on moderate high-frequency OOD;
- +0.98% on the canonical mean.

These are small, mixed differences. The correct conclusion is not that
wave-coverage data is bad; it is that this *replacement* design provides no
material overall advantage over broad-random training after both are allowed
to converge. Additive or interface-targeted wave coverage remains open.

## Validation convergence is not severe-OOD optimization

Every selected checkpoint improves the common validation rollout metric, but
the canonical mean worsens for direct, invariant, CNN, and wave arms. It
improves strongly for the characteristic and dissipation arms. This is not a
contradiction: the canonical problems are deliberately held-out distribution
shifts, whereas the validation set samples the declared training families.

Consequently, “trained to convergence” should mean convergence of the declared
validation objective, not automatic optimality on every severe OOD problem.
Future model selection should retain a clean ordinary validation set and add a
separate, independently parameterized stress-validation family; the named
canonical tests must remain final tests.

## Safety and reproducibility

- Every converged rollout remained admissible.
- Measured Tadmor violation rate was zero for every arm and split.
- The minimum observed density was `0.002176`; minimum pressure stayed at or
  above the declared `1e-5` floor.
- The exact 1,100-update direct/broad and direct/wave state dictionaries are
  bit-for-bit identical to the prior wave-coverage experiment (maximum tensor
  difference `0.0`).
- Boundary treatment stayed periodic. No learned or physical-boundary closure
  was introduced into this optimization audit.

The convergence trace, per-split metrics, pairwise comparisons, metadata,
checkpoints, and figure are stored in [`results/`](results/).

## Best checkpoint versus high-resolution FVM

The validation-selected winner, `dissipation_broad`, was rolled out on the five
canonical initial conditions without retraining or selecting on them. Its
64-cell predictions were compared with strict periodic 2048-cell Rusanov +
SSP-RK2 FVM trajectories. The profile figure is intentionally restricted to
two series: HCFL-64 and 64 point samples taken from the original 2048-cell
solution. Each group of 32 consecutive fine cells contributes its near-center
cell (`16 + 32*i`, zero-based), with no averaging, and the resulting values are
connected by ordinary lines. Quantitative errors elsewhere still use the
conservatively averaged FVM-2048 series so the reported coarse-cell metric
retains its finite-volume meaning.

| case | HCFL-64 NRMSE | FVM-64 NRMSE | HCFL reduction |
|---|---:|---:|---:|
| Sod | **0.01898** | 0.02746 | 30.88% |
| Lax | **0.08221** | 0.12329 | 33.33% |
| collision | **0.10936** | 0.25047 | 56.34% |
| strong pressure | **0.13796** | 0.21944 | 37.13% |
| near-vacuum expansion | **0.10755** | 0.23829 | 54.87% |
| five-case mean | **0.09121** | 0.17179 | **46.91%** |

The simplified profile directly compares the learned 64-cell values with an
unaveraged, stride-sampled view of the fine solution. It is a visual diagnostic,
not the conservative finite-volume error definition used by the tables and
error map.

An additional resolution-transfer profile compares native FVM-512 with the
same validation-selected dissipation checkpoint deployed on 512 cells. The
checkpoint is not retrained: it remains the model learned on 64-cell states,
and eight CFL-matched HCFL substeps are used per saved interval. Native
FVM-2048 is retained as a light background reference. This is a visual transfer
diagnostic, not a trained-HCFL-512 result; see
`results/fvm512_vs_hcfl512_with_fvm2048_seed0.png`.

The same checkpoint scored a canonical mean of `0.06257` against the 512-cell
FVM reference used by the convergence audit, but `0.09121` against FVM-2048, a
45.76% increase. Thus the absolute error assessment is still sensitive to
reference resolution.

Every converged arm was therefore evaluated on the same FVM-2048 trajectories:

| converged arm | FVM-2048 canonical mean | change from direct/broad |
|---|---:|---:|
| dissipation / broad | **0.09121** | **-20.23%** |
| characteristic / broad | 0.10083 | -11.82% |
| CNN / broad | 0.10569 | -7.57% |
| direct / broad | 0.11435 | reference |
| direct / wave | 0.11536 | +0.89% |
| invariant / broad | 0.11597 | +1.42% |

Dissipation remains the best arm on the FVM-2048 five-case mean and on four of
five individual cases; invariant is narrowly best on Sod. This strengthens the
seed-0 architecture result, but it is a held-out post-hoc evaluation rather
than a new checkpoint-selection rule.

See `results/best_vs_fvm_profiles_seed0.png`,
`results/best_vs_fvm_error_maps_seed0.png`, and
`results/best_vs_fvm_summary_seed0.json`. The all-arm comparison is in
`results/all_arms_vs_fvm_seed0.csv`.
