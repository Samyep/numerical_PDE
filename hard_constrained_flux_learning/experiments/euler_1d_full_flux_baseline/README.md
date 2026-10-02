# Direct full-flux baseline

## Pre-run question

Can the same five-cell MLP learn the complete Euler interface flux directly,
without receiving HLLC as an analytic base flux and without being restricted to
a correction parameterization?

## Single changed factor

The candidate maps the normalized five-cell primitive stencil directly to
three numbers interpreted as the complete shared interface flux:

```text
(rho, u, p) on five cells -> MLP(15, 72, 72, 3) -> F_NN
```

There is no `F_HLLC + correction`, jump gate, Roe basis, or analytic output
scale.  Everything after the proposal is held fixed: the hard Tadmor
half-space projection, local density/pressure limiter, fully-discrete entropy
line search, and conservative finite-volume divergence.

The controlled comparators are the converged `direct_broad` correction model
and the converged `dissipation_broad` model.  All three have 6,627 trainable
parameters and use the same broad training tensor, independent validation
tensor, initialization seed, minibatch schedule, optimizer, and stopping rule.

## Important identifiability check

Trajectory loss observes only

```text
F_(i+1/2) - F_(i-1/2)
```

so it cannot determine a spatially constant additive flux.  The experiment
therefore reports both rollout error and diagnostics that are easy to hide in
rollout error alone:

- raw and projected flux magnitude;
- spatial flux-divergence magnitude (to rule out collapse to a constant/zero
  update);
- hard-projection intervention rate and size;
- positivity and fully-discrete entropy limiter intervention rates;
- consistency error on constant states, `F_NN(U,U) - F_physical(U)`.

The last diagnostic is not used for training or checkpoint selection.

## Run

```powershell
python run_full_flux_baseline.py --self-test
python run_full_flux_baseline.py --seed 0
python plot_results.py --seed 0
python plot_full_flux512.py --seed 0
```

Only the validation-selected checkpoint is retained.  If the validation
plateau criterion is not reached before the 50,000-update cap, the checkpoint
is deleted and the run exits as unsuccessful.  This arm uses a larger cap than
the earlier audit because its validation error was still improving at update
20,000; the stopping rule itself is unchanged.
