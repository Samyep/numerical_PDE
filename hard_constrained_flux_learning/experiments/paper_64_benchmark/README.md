# Paper-facing 64-cell benchmark

This directory turns the earlier single-seed architecture screens into a
paper-facing, matched 64-cell experiment for 1D Euler and homogeneous 1D SWE.

The common HCFL method is the symmetric four-cell
`central + nonnegative Roe dissipation + proposal-feasibility loss` model.
Its forward solver retains the hard Tadmor projection, the conservative
positivity limiter, and the fully-discrete entropy limiter.  Seed 0 is loaded
from the original controlled stencil experiments; seeds 1 and 2 are trained
here with the same validation-plateau stopping rule and a 50,000-update cap.
Nonconverged checkpoints are deleted rather than reported.

The Euler training set has 580 trajectories (220 ordinary, 260 broad, and
100 extreme) and the disjoint validation set has 136 trajectories (44
ordinary, 52 broad, 20 extreme, and 20 randomized wave-coverage cases).  The
SWE split is 580 training trajectories (220 ordinary, 160 broad, and 200
Froude-coverage) and 136 validation trajectories (44, 32, and 60,
respectively).  Euler training targets use the existing strict 512-cell
Rusanov reference generator; SWE targets use its 1024-cell HLL reference generator.
Both are conservatively restricted to 64 finite-volume averages.  The final
paper-facing test metric is stricter and is recomputed against independent
2048-cell HLLC (Euler) or HLL (SWE) trajectories.

The completed numerical results, failure accounting, and claim boundaries are
summarized in [`RESULTS.md`](RESULTS.md).  The comparison with representative
PINN, neural-operator, divergence-free, learned-discretization, and neural-FV
work is in [`LITERATURE_COMPARISON.md`](LITERATURE_COMPARISON.md).

Run the replicate training:

```powershell
python train_hcfl_replicates.py --seeds 1 2
python train_operator_baselines.py --systems euler swe --seed 0
python train_learned_flux_baselines.py --systems euler swe --seed 0
python run_64_benchmark.py --seed 0 --random-trajectories 24 --reference-cells 2048 --hcfl-seeds 0 1 2 --fno-seed 0 --learned-flux-seed 0
python audit_benchmark.py
```

The later benchmark phase compares the converged model on 64 cells against
native first-order and MUSCL finite-volume solvers and a matched periodic FNO
baseline.  A same-data four-cell learned numerical flux supplies a direct
conservative neural-FV negative control: it preserves conservation through a
flux difference, but has no positivity or entropy constraint.  Literature
comparisons are kept separate from same-task numerical
comparisons because PINNs, neural operators, divergence-free surrogates, and
learned-flux methods solve materially different learning problems.

The benchmark also contains the requested low-order-anchor ablation.  It uses
the identical converged network and hard Tadmor proposal, but never computes or
blends an `F_low`: a proposed forward-Euler update is accepted only when it is
positive and its total mathematical entropy does not increase; otherwise its
time step is halved and retried.  Accuracy, substep count, rejected steps, and
both entropy diagnostics are reported separately from the retained full safety
stack.
