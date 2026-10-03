# Paper-consistent 1D SWE HCFL

This experiment applies the same retained HCFL method used for 1D Euler to
the homogeneous 1D shallow-water equations.

- Conservative state: `U = (h, hu)`.
- Primitive input: `P = (h, u)` on symmetric 4-cell and 6-cell stencils.
- Roe waves: `u_tilde - c_tilde` and `u_tilde + c_tilde`.
- Arm 1: HLL + learned Roe-coordinate correction.
- Arm 2: central flux + nonnegative Roe dissipation, with proposal-feasibility
  loss `lambda_feas = 1e-3`.
- Shared forward safety: hard Tadmor projection, conservative water-depth
  limiter, and fully-discrete total-entropy limiter.
- Selection: independent validation rollout plateau, with a 50,000-update
  cap rather than a fixed training budget.

Each retained method is trained and evaluated twice. The stencils are

- 4-cell: `(i-1, i | i+1, i+2)`
- 6-cell: `(i-2, i-1, i | i+1, i+2, i+3)`

No 5-cell model is part of this experiment. All four arms use exactly the
same trajectories, validation split, optimizer policy, and deployment
safeguards.

The training distribution deliberately includes subcritical, transcritical,
and both left- and right-going supercritical regimes. This is necessary to
exercise the two possible upwind directions; the legacy SWE data were almost
entirely subcritical.

Run invariant tests:

```powershell
python run_swe_consistent.py --self-test
```

Run or resume seed 0:

```powershell
python run_swe_consistent.py --phase all --seed 0 --resume
```

After all four validation-converged checkpoints exist, run the held-out
512-cell deployment audit (periodic and transmissive nonperiodic, both against
an HLL-2048 reference):

```powershell
python evaluate_deployment.py --seed 0
```

The nonperiodic audit uses the neural model only at interior interfaces.
Constant-extrapolation boundary states and physical SWE boundary fluxes are
fixed by the solver; no boundary flux is learned.

Seed-0 convergence, accuracy, hard-safety, fallback, and oscillation findings
are recorded in [`RESULTS.md`](RESULTS.md). In particular, lower NRMSE does
not eliminate shock ringing; the TV/extrema audit must be read alongside the
accuracy plot.

This first experiment is flat-bed SWE. Bathymetry requires a separate
well-balanced hydrostatic reconstruction and is not claimed here.
