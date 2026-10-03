# Euler 1D symmetric stencil ablation

This experiment compares 2-, 4-, and 6-cell neural input stencils for the
two retained HCFL flux proposals:

1. HLLC + learned Roe-coordinate correction.
2. Central + nonnegative Roe dissipation with proposal-feasibility loss
   `lambda_feas = 1e-3`.

For a flux at interface `i+1/2`, all new stencils are symmetric about the
interface:

| cells | primitive-state input |
|---:|---|
| 2 | `(P_i | P_{i+1})` |
| 4 | `(P_{i-1}, P_i | P_{i+1}, P_{i+2})` |
| 6 | `(P_{i-2}, P_{i-1}, P_i | P_{i+1}, P_{i+2}, P_{i+3})` |

The legacy five-cell stencil is intentionally not silently changed; its old
checkpoints remain loadable. It was `(P_{i-2}, P_{i-1}, P_i | P_{i+1},
P_{i+2})`, hence asymmetric about the interface.

All six new arms use the same seed, width, optimizer, 580 training
trajectories, 136 independent validation trajectories, and validation-plateau
checkpoint rule. Evaluation includes periodic 64-cell held-out sets,
periodic 512-cell canonical problems, and zero-shot transmissive nonperiodic
512-cell problems. In nonperiodic deployment the NN is evaluated only on
interior interfaces; both physical boundary fluxes remain analytic.

Run the invariant tests:

```powershell
python run_stencil_ablation.py --self-test
```

Run or resume the full seed-0 experiment:

```powershell
python run_stencil_ablation.py --phase all --seed 0 --resume
```

Recompute the independent checkpoint, indexing, boundary, activity, safety,
and entropy audit after evaluation:

```powershell
python audit_results.py --seed 0
```

The measured results and interpretation are recorded in `RESULTS.md`.

This is a single-seed screening ablation. A final architecture claim requires
confirmation over additional seeds.
