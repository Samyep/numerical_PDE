# Euler direct-512 validation, seed 0

## Setup

- Method: symmetric four-cell `Central + nonnegative Roe + feasibility`.
- Parameters: 6,411.
- Training targets: periodic HLLC-2048 + SSP-RK2 trajectories,
  conservatively restricted to 512 finite-volume cell averages.
- Split: 580 training trajectories and 136 independent validation
  trajectories. Named deployment cases are excluded from both splits.
- Deployment comparison: native HLLC-2048 reference, zero-network Roe-512,
  and HCFL-512 under the same interface projection and outer safety stack.
- Seed count: one. These numbers establish the seed-0 result, not uncertainty
  across random seeds.

The float32 reference-generation audit against a float64 CPU trajectory gave
a normalized RMS difference of `1.99e-7` and maximum absolute state difference
of `1.91e-6`.

## Convergence and network activity

Training stopped at update 39,500 after the validation metric plateaued at the
minimum learning rate, `3e-6`. The selected checkpoint is update 30,500:

| quantity | value |
|---|---:|
| initial validation rollout NRMSE | 0.0428443 |
| selected validation rollout NRMSE | 0.0206336 |
| relative reduction | 51.84% |
| minimum validation density | 2.317e-4 |
| minimum validation pressure | 8.203e-3 |

The learned multiplier did not collapse to the Roe identity or to a zero
correction. Its mean absolute deviation from one is `0.6385`, its mean is
`0.3661`, and its observed range is `[2.74e-4, 1.9983]`.

## Held-out deployment results

Rollout NRMSE against HLLC-2048 restricted to 512 cell averages:

| case | Roe-512 | HCFL-512 | HCFL reduction |
|---|---:|---:|---:|
| Sod | 0.01030 | 0.00441 | 57.2% |
| Lax | 0.04812 | 0.01576 | 67.2% |
| Collision | 0.05209 | 0.01350 | 74.1% |
| Strong pressure | 0.16581 | 0.05357 | 67.7% |
| Near-vacuum | 0.02222 | 0.00584 | 73.7% |
| Contact exits left | 0.00944 | 0.01008 | -6.8% |
| Contact exits right | 0.00944 | 0.00796 | 15.7% |
| Pressure wave exits right | 0.13013 | 0.07574 | 41.8% |

Across the first four cases, mean rollout NRMSE falls by 68.43% and final-time
NRMSE by 70.17%. Across the four boundary-interaction cases, the corresponding
reductions are 41.83% and 60.67%.

This accuracy improvement has an oscillation caveat. HCFL has more detected
excess significant extrema: 30 versus 3 in the first group and 20 versus 2 in
the boundary-interaction group. Mean normalized final TV excess is worse in
the first group (`0.0514` versus `0.0236`) but better in the second (`0.0367`
versus `0.0460`). Thus the result supports a large accuracy improvement, not a
claim that every smoothness diagnostic improves.

## Removing `F_low`

The learned proposal and feasibility loss do not use `F_low`; it appears only
in the outer local-positivity and fully-discrete entropy fallback. The same
selected checkpoint was therefore redeployed with both fallback blends
removed while retaining the interface entropy projection.

- Seven of eight named cases completed without a fully-discrete entropy
  increase. On several easy cases this version was numerically identical, or
  nearly identical, to the complete solver. This supports the observation that
  `F_low` is usually inactive and is not the source of the accuracy gain.
- The nonperiodic near-vacuum case failed at snapshot 29, internal substep 226:
  minimum density became `-2.437e-6`. Before failure, 37 of 226 fully-discrete
  checks violated the entropy balance (16.37%), with maximum increase
  `5.94e-2`.
- With the complete safety stack, the same near-vacuum HCFL rollout completed
  with minimum density `1.202e-2`, minimum pressure `9.257e-4`, and rollout
  NRMSE `0.00584`. Its fully-discrete limiter was active on 7.54% of checks and
  reached `beta = 0`.

Conclusion: removing `F_low` is a useful simplified empirical variant for
benign cases, but it is not equivalent to the hard-safe method and cannot
replace `F_low` when positivity and fully-discrete entropy are claimed as
guarantees. A clean paper design can keep `F_low` outside the neural network as
a rarely active safety fallback, rather than presenting it as part of the
learned flux architecture.

## Artifacts

- `results/euler_roe512_vs_hcfl512_with_hllc2048_seed0.png`
- `results/validation_convergence_train512_seed0.png`
- `results/deployment_roe_vs_hcfl_train512_seed0.csv`
- `results/no_f_low_ablation_seed0.csv`
- `results/summary_roe_vs_hcfl_train512_seed0.json`
- `results/central_nonnegative_feas_s4_train512_converged_best_seed0.pt`
