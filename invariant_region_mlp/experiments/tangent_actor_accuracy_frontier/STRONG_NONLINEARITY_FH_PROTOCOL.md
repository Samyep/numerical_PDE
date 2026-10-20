# Paired high-curvature Full-History MLP pilot (no training)

## Goal

Measure raw Full-History MLP accuracy and runtime against the **already frozen**
Tangent-Actor + MCV checkpoint on a demonstrably more strongly nonlinear PDE.
The original primary `lambda=0.2, eta=0` has a strong FH baseline that is
cheaper near RMSE 5e-4; moving to `lambda=0.5, eta=0` increases the formal
maximum `|f_uu|` from `0.2*(2*pi)^2` to `0.5*(2*pi)^2` (2.5x). Compare on
the same input and independent truth arrays, and retain a negative result.

Here `f(u,z) = lambda*sin(2*pi*u) - ||z+eta*u*e||^2/2`; the default stronger
*value* nonlinearity is `lambda=0.5, eta=0`, `d=100`, `k=2`, with the
previously frozen **pilot** `a2` (64 Euler steps, 512 paths, 512 held-out states).
The known Actor pilot RMSE is `0.00954707` at ~520.55 s total including
training+selection+inference. This is a pilot, not the mature four-update Actor.

A separate, already-completed `lambda=0.2, eta=1` coupled run tests
nonseparable `u,z` interactions and has existing FH `n=4,M=96` evidence.
Do not confuse these distinct drivers or compare their RMSEs as paired runs.

## Preregistered FH controls

The supplementary sweep is frozen to:
`(n,M) = (3,32), (3,64), (3,96), (3,128), (4,16), (4,32), (4,64), (4,96)`.
Begin with the cheaper cases, and run the expensive `(4,96)` only if the
remaining allocated-GPU budget permits. Use seed 0 as primary; optional
additional seeds must be predeclared and all included, not selected by test RMSE.
Do not use test accuracy to choose the next budget.

Full-History is the original recursive estimator plus terminal-memory tiling;
it has **no** Actor, MCV, oracle projection, or terminal-aware learned correction.
The model and reference are unchanged. The test bank is loaded from the
immutable original `post_freeze_primary_d100_k2_lam0p5_eta0_pilot_mcv_s0.npz`
and requires exactly 512 test states. The checker reloads the independent
2D numerical reference and verifies truth array identity to 1e-9.
It does not retrain/reselect the Actor or change the final-suite 140-cell grid.

## Instructions on the GH200 HPC workspace

```bash
# Use a scratch checkout of branch ir-mlp-strong-nonlinearity-fh-audit.
# Keep the original HPC results/nonconvex_tangent_actor/{predictions,references}
# available at the same project-relative paths. These large files are not
# stored in GitHub.

python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.strong_nonlinearity_fh_pair --dry-run

# Example single cell; the 1800-second shell cap is in addition to the
# script's own per-cell timeout. Always charge actual Slurm GPU allocation.
timeout -k 60s 1800s python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.strong_nonlinearity_fh_pair \
  --family primary --lam 0.5 --eta 0 --budget pilot --d 100 \
  --actor-seed 0 --baseline-seed 0 --n 3 --M 32 \
  --device cuda --max-seconds 1750

# Repeat for the remaining frozen configurations, respecting available
# Slurm GPU allocations. Independent configs may run concurrently ONLY on
# different allocated GPUs and if account concurrency rules permit it.
```

Results are written only to
`invariant_region_mlp/results/strong_nonlinearity_fh_paired/`, never into the
original prediction/checkpoint directories. Existing output paths are
refused rather than overwritten. Each run records the original test-bank hash,
truth hash, Actor source JSON SHA-256, FH `(n,M)`, measured 512-point time,
and 960-point **linear extrapolation** explicitly labeled as *not measured*.
An out-of-time experiment must not receive a fabricated RMSE.

## Required analysis after execution

1. Plot measured RMSE versus measured total GPU seconds for all 512 points.
2. Do **not** multiply the already-batched FH runtime by 960. For a 960-point
   extrapolation, multiply the measured 512-point inference time by 960/512.
3. Actor 960 extrapolation: keep frozen training+reference-free selection cost
   fixed and scale only inference seconds.
4. Report whether FH reaches the frozen Actor error `0.00954707`. If yes,
   compare minimum measured total times (with an explicit disclaimer that
   configurations were chosen before RMSE evaluation); if not, label
   `NOT_REACHED` within the tested budgets, not theoretical impossibility.
5. Do not mix seed-0 point estimates with three-seed means without labeling;
   retain all negative/failed/OOT rows and reference uncertainty.
6. Never interpret a change in `lambda` alone as higher *effective dimension*;
   this remains an intrinsic two-direction example.
7. The current machine running the ChatGPT review has no CUDA GPU or live
   Slurm access. The runner is **prepared and GitHub-pushed but not run on
   GH200**; do not claim new experimental performance until actual HPC artifacts
   are returned.

## Existing control (already measured, lambda=0.2, eta=1)

For the 100D coupled benchmark at 512 query states:
- FH `n=4, M=96`: mean RMSE ~1.63057e-3, ~886.5 GPU s;
- frozen mature Actor: 64-step RMSE ~8.09e-4; 512-step RMSE ~1.36e-4;
- mature Actor training+selection ~3292.2 GPU s per seed;
- frozen Actor 64-step inference ~25.1 s and 512-step inference ~193.9 s per 512 states.

At 960 query points the linear time extrapolation gives FH ~27.7 min, Actor
~55.7 min (64-step) or ~60.9 min (512-step), without changing measured 512-point
RMSE. A matched total-time crossover is around 1950–2450 queries **only
conditional on a linear scaling assumption**, and the configurations have
different accuracy. This is not a measured amortization curve.

Scientific goal: establish whether stronger nonlinearity shifts the
accuracy–compute crossover in favor of learned tangent linearization, not
manufacture a FH failure by choosing an underpowered MLP budget.
