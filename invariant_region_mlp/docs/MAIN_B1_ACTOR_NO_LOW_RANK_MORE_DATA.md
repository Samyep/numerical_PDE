# No-projection actor with twice as many labeled states

## Result

Doubling the raw-label data and allowing a larger optimization budget does not
improve the no-projection actor at `(d,k)=(400,40)`.  Its Stage-B value RMSE is
**0.12104 ± 0.00482**, statistically indistinguishable at this resolution from
the first no-projection result, **0.12016 ± 0.00335**.  Their paired mean
difference is `+0.00088 ± 0.00263` over the same ten MLP repetitions.

The enlarged run remains **49.7% worse** than the projected Actor--MLP and
**39.1% worse** than PINN+SDGD, the best previous value baseline in this cell.
All ten paired seed-level errors are higher than those of the projected actor.

| Method | Training points | Value RMSE | Ratio to projected |
|---|---:|---:|---:|
| Actor--MLP, no projection, enlarged data | 65,536 | 0.12104 ± 0.00482 | 1.497 |
| Actor--MLP, no projection, first trial | 32,768 | 0.12016 ± 0.00335 | 1.486 |
| Plain actor at the original level-5 checkpoint | 32,768 | 0.10743 ± 0.00196 | 1.328 |
| Actor--MLP, learned rank-40 projection | 32,768 | **0.08087 ± 0.00080** | 1.000 |
| PINN+SDGD | — | 0.08700 ± 0.00667 | 1.076 |

Intervals are 95% normal approximations across seed-level RMSEs.  Comparisons
among the three Actor--MLP rows use identical test points and MLP random-number
seeds.

## Frozen enlarged protocol

The run keeps the plain ambient-input, ambient-output actor with four
width-384 SiLU hidden layers and 751,888 trainable parameters.  It contains no
label projection, output projection, or `label_basis`.  It restarts from the
original unprojected level-5 checkpoint so that it does not inherit the first
trial's downstream degradation.

The original 32,768 training and 8,192 validation states are combined with
independent sets of the same sizes, yielding 65,536 training and 16,384
validation states.  Every added raw label uses `n=3`, `M=20`, and ten Monte
Carlo replicas.  The base and added training-point hashes are respectively
`9c2713ec86580fe8c7cd95d5a169f50fe4c3b40489a7e263c523026da0fd6843` and
`a919e8400be1fe03530c16224c3f0be02270ff829774c5044053579a7103581e`.

Three training stages use learning rates `1e-4`, `2.5e-5`, and `5e-6`.  Their
predeclared maximum budgets total 26,000 steps, versus 17,000 in the first
confirmation trial.  All stages reached the frozen validation-and-parameter
plateau rule after 4,700, 3,000, and 2,100 steps, so forcing the remaining
maximum steps was unnecessary under the selection protocol.

The final checkpoint is evaluated once on the same 960 held-out points, ten
Actor--MLP repetitions, `n=3`, `M=16`, and `B1|actor_mlp` seed namespace.  The
test-point hash remains
`c38372e2258615471a85deb55de71c289673fcb3d673be682d5a6d36aca184d8`.
Only value accuracy is evaluated; the run never constructs `exact_z` or uses a
gradient metric for training, selection, or reporting.

## Optimization behavior and interpretation

The enlarged validation-label RMSE decreases from 1.76618 to 1.76603 and then
1.76582.  This is slightly below the first no-projection trial's final 1.77505,
but it does not transfer to Stage-B value accuracy.  The enlarged result is
also 12.7% worse than the original level-5 checkpoint, with higher RMSE in all
ten paired seeds.

The direct conclusion is limited but clear: at this cell, twice as many
independently labeled states and a larger allowed optimizer budget do not
replace the learned output projection.  A plausible explanation is that more
states leave the per-label ambient orthogonal noise unchanged, whereas the
rank-40 projection removes that noise before optimization.  This explanation
is an inference from the ablation rather than a separate noise-scaling test.

## Execution and artifacts

Slurm job `3373792` completed on GH200 node `gh034` in 7 minutes 23 seconds
with exit code `0:0`.  New label construction took 399.85 seconds, actor
optimization took 17.67 seconds, and the complete Python runner took 426.60
seconds.

- Scalar result:
  `../results/main_b1_sixmethods/actor_convergence/no_low_rank_more_data/result_d400_k40.json`
- Local label caches, checkpoints, and value predictions:
  `../results/main_b1_sixmethods/actor_convergence/no_low_rank_more_data/`
- Independent run record and scheduler logs:
  `../../../actor_no_low_rank_ablation_20261013/`
