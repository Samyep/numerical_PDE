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

## Conclusion

The experiments answer the two proposed directions differently:

1. **Direct complete flux + Roe decomposition is not enough.**  An invertible
   basis change does not provide consistency, upwinding, or useful
   cross-resolution shock behavior.
2. **An analytic central flux plus Roe-wave assembly is useful.**  Both
   structured variants beat the previous learned method on the five-case
   512-cell mean.
3. **Nonnegative automatic upwinding is the preferred main method.**  Its mean
   error differs from signed by only 0.92% in this one seed, which is not enough
   evidence to favor signed over the stronger structural guarantee.  Its
   mechanism is physically interpretable, its final-time mean is lower, its
   projection burden is smaller, and its near-vacuum pressure margin is much
   better.
4. The signed arm should remain as an accuracy/anti-diffusion ablation, not be
   presented as the safer default.

This is still a single-training-seed result over five named problems.  It is
strong evidence for the architecture, not yet a multi-seed statistical claim.

## Artifacts

- `results/roe_upwind_summary_seed0.png`
- `results/roe_upwind512_with_hllc2048_seed0.png`
- `results/summary_seed0.json`
- `results/roe_upwind512_with_hllc2048_seed0.json`
- `results/scientific_integrity_audit_seed0.json`
- Three converged `*_converged_best_seed0.pt` checkpoints
