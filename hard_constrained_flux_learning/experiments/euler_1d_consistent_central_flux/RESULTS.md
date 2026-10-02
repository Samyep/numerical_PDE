# Results: exactly consistent central flux + learned correction

Seed: 0.  The model was selected only by the independent validation rollout
NRMSE.  Canonical Riemann cases and the 512-cell transfer test were not used
for training, early stopping, or checkpoint selection.

As in the controlled comparator runs, the training targets are strict
512-cell Rusanov + SSP-RK2 trajectories conservatively restricted to 64 cell
averages.  The network is trained and validated on those 64-cell trajectories.

## Model tested

At interface `i + 1/2`, the raw proposal is

```text
F_hat = 0.5 * (F(U_i) + F(U_{i+1}))
        + 0.18 * ||(P_{i+1} - P_i) / sigma_P||_2
          * scale * tanh(MLP(P_{i-2:i+2}))
```

The jump norm contains no epsilon.  Therefore its learned term is exactly
zero whenever `U_i == U_{i+1}`, even when the other three stencil states are
different.  The architecture has 6,627 trainable parameters, matching the
controlled neural comparators.  It uses neither HLLC nor a Roe decomposition
inside the raw proposal.  The common Tadmor projection and rollout safety
stack remain enabled.

## Convergence and mechanism checks

- Training met the preregistered validation-plateau rule at update 15,700.
- The selected checkpoint is update 15,200, with validation rollout NRMSE
  `0.0457603363`.
- The 50,000-update value was a cap, not the selected training length.
- Maximum raw and projected error in the varied-outer-stencil test of
  `F_hat(U,U) = F(U)` is exactly `0.0`.
- The learned correction did not collapse to zero: its RMS is 2.517% of the
  central-flux RMS, and its flux-divergence RMS is 22.344% of the central
  flux-divergence RMS.
- The hard Tadmor projection changes 23.646% of validation interfaces; its
  RMS change is 1.651% of the raw-flux RMS.

For comparison, the selected validation rollout NRMSE values are:

| Model | Validation rollout NRMSE |
|---|---:|
| Direct complete flux | 0.0280469 |
| HLLC + Roe-dissipation correction | 0.0315568 |
| Central flux + exactly consistent correction | 0.0457603 |
| HLLC + direct vector correction | 0.0494303 |

Thus restoring an analytic central flux and exact consistency improves on the
HLLC direct-vector-correction validation score, but does not match the Roe
dissipation parameterization or the direct complete-flux fit.

## 64-cell held-out rollouts

Mean NRMSE over the three random-distribution splits:

| Model | Mean NRMSE |
|---|---:|
| Direct complete flux | 0.0113292 |
| HLLC + Roe-dissipation correction | 0.0129886 |
| Central flux + exactly consistent correction | 0.0144642 |
| HLLC + direct vector correction | 0.0170641 |

Mean NRMSE over the five canonical Riemann cases:

| Model | Mean NRMSE |
|---|---:|
| HLLC + Roe-dissipation correction | 0.0625750 |
| HLLC + direct vector correction | 0.0880568 |
| Direct complete flux | 0.0975063 |
| Central flux + exactly consistent correction | 0.110243 |

The consistent-central model helps relative to direct complete flux on Sod,
collision, and near-vacuum at 64 cells, but it is not uniformly better and is
substantially worse on Lax and strong-pressure cases.  Roe-dissipation is
better on all five canonical cases.

## Unchanged deployment on 512 cells

The 64-cell-trained checkpoint was applied convolutionally to 512 cells
without retraining or parameter changes.  NRMSE is computed against a strict
native HLLC-2048 trajectory restricted conservatively to 512 cell averages;
the plotted HLLC-2048 curve itself remains on its native 2048-point grid.

| Method | Completed cases | Five-case mean rollout NRMSE |
|---|---:|---:|
| HLLC-512 | 5/5 | 0.0742369 |
| Roe-dissipation HCFL-512 | 5/5 | 0.0562689 |
| Consistent-central HCFL-512 | 5/5 | 0.0938677 |

The new model completes every case, unlike the earlier direct complete-flux
model's near-vacuum cost failure.  However, its mean error is 26.44% higher
than native HLLC-512 and 66.82% higher than Roe-dissipation HCFL-512.  It is
worse than native HLLC-512 on each of the five cases.  The plotted profiles
show shock-localized oscillations, especially in Lax, collision, strong
pressure, and near-vacuum.  Its pressure approaches the `1e-5` safety floor in
the latter three difficult cases, and hard projection acts on roughly 42--44%
of interfaces for collision and near-vacuum.

## Interpretation

The experiment separates consistency from useful dissipation:

1. Exact `F_hat(U,U) = F(U)` fixes a real structural defect of unconstrained
   complete-flux prediction.  This architecture also avoids the previous
   catastrophic 512-cell failure, but the central base and jump gate were
   introduced together, so that stability change cannot be attributed to
   consistency alone.
2. Consistency constrains only the diagonal `U_L = U_R`; it does not determine
   the correct upwind response when `U_L != U_R`.
3. The arithmetic central flux has no built-in shock dissipation.  The neural
   correction must learn the entire dissipative mechanism from 64-cell data,
   and that learned mechanism does not transfer cleanly to 512-cell shock
   widths.
4. Entropy projection and positivity/global-entropy fallback ensure the stated
   inequalities and admissibility, but do not guarantee non-oscillatory or
   high-accuracy profiles.

For this seed and protocol, the result supports retaining the Roe-wave
dissipation structure.  It does not support replacing it with an unrestricted
vector correction around the arithmetic central flux.  Multi-seed replication
would still be required before making a statistical claim.

## Artifacts

- `results/central_consistent_64cell_comparison_seed0.png`
- `results/central_consistent512_vs_hllc_and_roe_with_hllc2048_seed0.png`
- `results/report_seed0.json`
- `results/central_consistent512_vs_hllc_and_roe_with_hllc2048_seed0.json`
- `results/central_consistent_broad_converged_best_seed0.pt`
