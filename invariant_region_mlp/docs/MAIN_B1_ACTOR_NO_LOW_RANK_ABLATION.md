# Actor output-projection ablation at `(d,k)=(400,40)`

## Result

Removing the actor's learned rank-`k` output projection eliminates its value
advantage in this focused setting.  The plain ambient-output Actor--MLP reaches
a mean value RMSE of **0.12016 ± 0.00335**, compared with
**0.08087 ± 0.00080** for the projected actor and **0.08700 ± 0.00667** for
PINN+SDGD, the best previous value baseline in this cell.  Thus the unprojected
actor is 48.6% worse than the projected actor and 38.1% worse than the best
previous baseline.

The comparison with the projected actor is paired over the same ten MLP
repetitions.  All 10/10 seed-level RMSEs are higher without the projection; the
mean paired increase is **0.03929 ± 0.00311**.  The uncertainty intervals are
95% normal approximations across seed-level RMSEs.

| Method | Value RMSE | Ratio to projected | Ratio to best previous |
|---|---:|---:|---:|
| Actor--MLP, no output projection | 0.12016 ± 0.00335 | 1.486 | 1.381 |
| Actor--MLP, learned rank-40 projection | **0.08087 ± 0.00080** | 1.000 | 0.930 |
| PINN+SDGD | 0.08700 ± 0.00667 | 1.076 | 1.000 |

## Matched protocol

This is a one-cell structural ablation using the most demanding existing
unprojected pilot, `(d,k)=(400,40)`.  The actor remains an ambient-input
four-hidden-layer, width-384 SiLU network with 751,888 trainable parameters,
but its 400-dimensional output is used directly.  No label projection, output
projection, or learned subspace appears in the ablated model.  Here "no low
rank" refers to the actor output constraint; the frozen benchmark PDE retains
its specified effective dimension `k=40`.

The run starts from the protocol-v2 plain actor's level-5 checkpoint.  It
reuses the cached raw `n=3`, `M=20`, ten-replica labels on 32,768 training and
8,192 validation points, then applies the same level 6--8 confirmation
learning rates (`1e-4`, `2.5e-5`, and `5e-6`) selected by the protocol-v3
projected actor.  Final evaluation uses the same 960 held-out points, ten
Actor--MLP repetitions, `n=3`, `M=16`, and `B1|actor_mlp` random-number
namespace.  The held-out point hash is
`c38372e2258615471a85deb55de71c289673fcb3d673be682d5a6d36aca184d8`.

Only value accuracy is evaluated.  The runner never constructs `exact_z` and
does not report direct-actor or recursive gradient accuracy.

## Optimization behavior

The unprojected actor's raw-label validation RMSE changes from 1.78018 at the
starting level-5 checkpoint to 1.77543, 1.77512, and 1.77505 after confirmation
levels 6, 7, and 8.  Despite that small improvement against noisy ambient
labels, Stage-B value RMSE increases from **0.10743 ± 0.00196** at the old
level-5 checkpoint to **0.12016 ± 0.00335** after level 8.  Lower learning-rate
confirmation therefore does not recover the projected actor's performance in
this setting.

This result is intentionally limited to the one previously selected setting;
it does not estimate a full 15-cell no-projection performance surface.  At the
high-dimensional cell where the structural ablation is most stringent, the
learned output projection is necessary for the reported value win.

## Execution and artifacts

Slurm job `3373686` completed on GH200 node `gh042` in 41 seconds with exit
code `0:0`.  Neural confirmation training took 16.31 seconds, and the complete
runner took 24.95 seconds after process startup.

- Scalar result:
  `../results/main_b1_sixmethods/actor_convergence/no_low_rank_ablation/result_d400_k40.json`
- Local checkpoints and value predictions:
  `../results/main_b1_sixmethods/actor_convergence/no_low_rank_ablation/`
- Independent run record and scheduler logs:
  `../../../actor_no_low_rank_ablation_20261013/`
