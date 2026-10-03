# Results: Roe coordinates and automatic upwinding

Seed: 0.  All checkpoints below satisfied the same independent-validation
plateau rule.  Canonical Riemann problems and 512-cell deployment were never
used for checkpoint selection.

## Controlled setup

- Training data: 580 strict 512-cell Rusanov + SSP-RK2 trajectories,
  conservatively restricted to 64 cell averages.
- Validation data: 136 disjoint trajectories generated from different seeds.
- Network: five primitive cells, two width-72 hidden layers, 6,627 trainable
  parameters.
- Maximum training cap: 50,000 updates; every retained run stopped earlier by
  validation plateau at the minimum learning rate.
- Shared hard safety stack: Tadmor half-space projection, local admissibility
  limiter, and fully-discrete total-entropy limiter.
- Boundary condition: periodic.

The two central-Roe models start from the same entropy-fixed classical Roe
flux.  The only output-map difference is:

```text
signed:       d = 1 + 2*tanh(z),  -1 < d < 3
auto-upwind:  d = 1 +   tanh(z),   0 < d < 2
```

Both use

```text
F_hat = 0.5*(F_L + F_R) - 0.5*R@(d*|lambda|_fix*alpha)
alpha = solve(R, U_R - U_L)
```

and therefore satisfy `F_hat(U,U)=F(U)` exactly.  The nonnegative multiplier
fixes each characteristic wave's dissipative direction, but the common safety
stack is still required for positivity and entropy guarantees.

## Convergence

| Model | Best update | Stop update | Validation rollout NRMSE |
|---|---:|---:|---:|
| Central + signed Roe multipliers | 8,800 | 13,800 | **0.0245126** |
| Direct complete flux, physical coordinates | 20,000 | 23,000 | 0.0280469 |
| Central + nonnegative auto-upwind | 14,800 | 16,300 | **0.0290998** |
| Direct complete flux, Roe coordinates | 8,000 | 13,000 | 0.0305313 |
| HLLC + signed Roe correction | 12,500 | 16,500 | 0.0315568 |
| Central + unrestricted vector correction | 15,200 | 15,700 | 0.0457603 |

All three new checkpoint metrics were independently recomputed after loading
from disk and matched their recorded minima exactly.

## Learned mechanism

| Diagnostic on validation states | Signed Roe | Auto-upwind Roe |
|---|---:|---:|
| Minimum multiplier | -0.8859 | 0.0000662 |
| Maximum multiplier | 2.5624 | 1.9971 |
| Mean multiplier | 0.2341 | 0.5086 |
| Negative multiplier fraction | **25.326%** | **0%** |
| Mean left/contact/right multiplier | 0.162 / 0.462 / 0.079 | 0.386 / 0.681 / 0.458 |
| Hard-projection intervention rate | 9.511% | 3.855% |

The signed model obtains part of its accuracy by learning anti-diffusive wave
coefficients.  The nonnegative model also learns substantially less
dissipation than classical Roe (`d=1`), but never reverses its direction.

The direct Roe-coordinate complete-flux arm does not acquire consistency merely
from changing coordinates: its maximum raw equal-interface error is `4.7478`.
Roe decomposition alone is therefore an output basis, not an upwind or
consistency constraint.

## Held-out 64-cell results

| Model | Three random splits mean | Five canonical cases mean |
|---|---:|---:|
| Central + signed Roe multipliers | **0.0118338** | **0.0417283** |
| Central + nonnegative auto-upwind | 0.0129260 | **0.0470929** |
| HLLC + signed Roe correction | 0.0129886 | 0.0625750 |
| Direct complete flux, physical coordinates | 0.0113292 | 0.0975063 |
| Direct complete flux, Roe coordinates | 0.0133854 | 0.1087450 |
| Central + unrestricted vector correction | 0.0144642 | 0.110243 |

Both structured central-Roe models improve markedly on the previous best
canonical mean.  The coordinate-only complete-flux model does not.

## Zero-shot deployment on 512 cells

The same 64-cell-trained checkpoint is applied convolutionally to 512 cells
without retraining or parameter changes.  Metrics use strict native HLLC-2048
conservatively restricted to 512 cell averages; the plotted reference curve
itself remains on all 2048 native cells.

### Per-case rollout NRMSE

| Case | Old HLLC + Roe | Signed central-Roe | Auto-upwind central-Roe | Roe-complete |
|---|---:|---:|---:|---:|
| Sod | 0.00682 | 0.00984 | **0.00570** | 0.01305 |
| Lax | **0.03076** | 0.03199 | 0.03356 | 0.13628 |
| Collision | 0.05115 | 0.04339 | **0.04038** | 0.27538 |
| Strong pressure | 0.12973 | **0.11596** | 0.13005 | 0.41742 |
| Near-vacuum | 0.06288 | 0.04552 | **0.03930** | 0.21492 |

### Aggregate result

| Method | Completed | Mean rollout NRMSE | Mean final-snapshot NRMSE |
|---|---:|---:|---:|
| Central + signed Roe multipliers | 5/5 | **0.0493395** | 0.0641927 |
| Central + nonnegative auto-upwind | 5/5 | **0.0497970** | **0.0619062** |
| Old HLLC + signed Roe correction | 5/5 | 0.0562689 | 0.0689903 |
| Native HLLC-512 | 5/5 | 0.0742369 | 0.0913162 |
| Central + unrestricted vector correction | 5/5 | 0.0938677 | 0.114931 |
| Direct complete flux in Roe coordinates | 5/5 | 0.211411 | 0.347707 |

Relative to the old HLLC+Roe model, mean rollout error falls by 12.31% for the
signed model and 11.50% for automatic upwinding.  The signed mean is only 0.92%
lower than auto-upwind, while auto-upwind has a 3.56% lower final-snapshot
mean, wins Sod/collision/near-vacuum, and never uses a negative wave multiplier.
The earlier direct physical-complete model is omitted from this table because
its near-vacuum 512-cell run hit the preregistered safe-step cost limit, making
an all-five-case mean invalid.

The safety diagnostics also distinguish them:

| 512-cell diagnostic averaged over five cases | Signed Roe | Auto-upwind Roe |
|---|---:|---:|
| Hard-projection intervention rate | 18.036% | 10.525% |
| Local admissibility intervention rate | 0.0113% | 0.0575% |
| Fully-discrete entropy intervention rate | 0% | 0% |
| Minimum pressure over all cases | `1.11e-5` | `1.92e-3` |

The signed model approaches the `1e-5` pressure floor in near-vacuum, whereas
automatic upwinding keeps substantially more pressure margin.  The direct
Roe-complete model is oscillatory, needs hard projection on 47.70% of 512-cell
interfaces on average, and also approaches the pressure floor.

## Oscillation-aware interpretation

Rollout NRMSE alone hides the main weakness of the signed central-Roe arm.  An
additional diagnostic over the final density, velocity, and pressure profiles
of all five 512-cell cases gives:

| Method | Mean TV / reference TV | Excess local extrema | Mean normalized range violation | Interpretation |
|---|---:|---:|---:|---|
| HLLC + Roe correction | **1.0351** | **55** | 1.017% | Retain: reliable HLLC-anchored method |
| Central + signed Roe | 1.6577 | 790 | 3.655% | Reject as a main method: severe oscillation |
| Central + nonnegative Roe, no feasibility loss | 1.0430 | 74 | **0.786%** | Architectural control |
| Central + nonnegative Roe + feasibility (`1e-3`) | 1.0482 | 73 | 1.385% | Retain: feasibility-trained alternative |

The successful second configuration is the **nonnegative Roe model trained
with proposal-feasibility loss**, not the signed model.  Relative to signed
Roe, its mean TV ratio falls from `1.6577` to `1.0482`, excess extrema fall
from 790 to 73 (90.8%), and mean normalized range violation falls from 3.655%
to 1.385% (62.1%).  The visibly severe signed-Roe oscillation is therefore no
longer present in the retained feasibility-trained configuration.

The controlled comparison against the zero-penalty nonnegative model needs a
more precise interpretation: `lambda_feas=1e-3` changes aggregate excess
extrema only from 74 to 73 and does not improve aggregate TV or range
violation.  It does reduce selected visible oscillations (for example, Lax
density excess extrema fall from 13 to 7), but the present ablation cannot
attribute all of the improvement over signed Roe to the feasibility term
alone.  The nonnegative output map and feasibility training act together in
the retained method.

Therefore there are **two successful methods to carry forward**:

1. **HLLC + Roe correction**: the conservative reference architecture, with a
   physics-based Riemann-solver anchor and the fewest excess extrema here.
2. **Central + nonnegative Roe + proposal-feasibility loss
   (`lambda_feas=1e-3`)**: a distinct, effective feasibility-trained
   architecture whose final profiles avoid the severe signed-Roe oscillation.

These methods are complementary candidates.  The present one-seed experiment
does not justify declaring either one universally superior.

## Conclusion

The experiments answer the two proposed directions differently:

1. **Direct complete flux + Roe decomposition is not enough.**  An invertible
   basis change does not provide consistency, upwinding, or useful
   cross-resolution shock behavior.
2. **An analytic central flux plus Roe-wave assembly is useful.**  Both
   structured variants beat the previous learned method on the five-case
   512-cell mean.
3. **Central + nonnegative Roe + proposal-feasibility loss is an effective
   second method alongside the existing HLLC + Roe correction method.**  It
   retains exact consistency, nonnegative characteristic dissipation, and a
   hard-projected PDE update, while sharply suppressing the signed-arm
   oscillations.
4. **HLLC + Roe correction remains a main method**, rather than being replaced:
   it supplies the stronger physics-based baseline and has slightly lower
   oscillation counts in these tests.
5. The signed arm should remain only as an accuracy/anti-diffusion ablation;
   its low aggregate NRMSE is not sufficient to accept its visibly oscillatory
   solutions.

This is still a single-training-seed result over five named problems.  It is
strong evidence for the architecture, not yet a multi-seed statistical claim.

A focused follow-up added the raw-proposal feasibility loss while keeping this
nonnegative architecture and the hard-projected PDE update unchanged.  It
improved teacher-forced 64-cell proposal feasibility but worsened zero-shot
512-cell accuracy relative to the zero-penalty control.  Nevertheless, the
`1e-3` configuration remains a useful oscillation-controlled method and is the
second retained approach here; this is different from claiming that the
penalty wins its own ablation.  See `FEASIBILITY_RESULTS.md`.

## Artifacts

- `results/roe_upwind_summary_seed0.png`
- `results/roe_upwind512_with_hllc2048_seed0.png`
- `results/summary_seed0.json`
- `results/roe_upwind512_with_hllc2048_seed0.json`
- `results/scientific_integrity_audit_seed0.json`
- Three converged `*_converged_best_seed0.pt` checkpoints
