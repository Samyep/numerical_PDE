# SWE: direct 512-cell training validation

This experiment tests the resolution-transfer explanation for the ringing
seen when a 64-cell-trained closure is deployed on 512 cells.  It trains one
model from scratch at 512 cells:

- `Central + nonnegative Roe + proposal-feasibility loss`
- symmetric 4-cell stencil `(i-1, i | i+1, i+2)`
- HLL-2048 reference trajectories conservatively restricted to 512 cells
- the same 580 training and 136 independently seeded validation trajectories
  as the 64-cell experiment
- hard Tadmor projection in every learned forward substep
- validation-rollout checkpoint selection with a 50,000-update cap

The original 64-cell training interval is `dt = 5e-4`, for which
`dt / dx = 0.032`.  A single 512-cell update over the same interval would
instead have `dt / dx = 0.256` and would not be a fair or generally stable
experiment.  The 512-cell forward map therefore uses eight differentiable
substeps, each with `dt / dx = 0.032`, while retaining the same physical
supervision interval and final training horizon.

Run:

```powershell
C:\Users\Yeping\miniforge3\envs\bsae\python.exe run_train512.py --phase all --resume
```

The script uses CUDA when available.  Cached reference trajectories are
local-only; the converged checkpoint, metrics, figures, and scientific
summary are retained.

Seed-0 findings are recorded in [`RESULTS.md`](RESULTS.md).
