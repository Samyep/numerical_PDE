# Roe-coordinate and automatic-upwind ablation

## Preregistered question

Three remaining designs are tested under the same seed-0 data, parameter
count, optimizer, independent validation selection, hard Tadmor projection,
positivity limiter, and fully-discrete entropy limiter used by the existing
Euler comparison.

### 1. Direct complete flux in Roe coordinates

```text
c = MLP(five-point primitive stencil)
F_hat = normalized_Roe_basis(U_L, U_R) @ (local_flux_scale * c)
```

This arm predicts the entire flux and has no analytic central/HLLC anchor,
jump gate, or exact consistency.  It asks whether a Roe coordinate system alone
repairs direct complete-flux learning.

### 2. Central flux with signed Roe multipliers

```text
d = 1 + 2*tanh(MLP(stencil))             # -1 < d < 3
F_hat = 0.5*(F_L + F_R) - 0.5*R@(d*|lambda|_fix*alpha)
alpha = solve(R, U_R - U_L)
```

This arm starts from the classical entropy-fixed Roe flux but permits a learned
wave multiplier to become negative and reverse the dissipative direction.

### 3. Central flux with nonnegative automatic upwinding

```text
d = 1 + tanh(MLP(stencil))               # 0 < d < 2
F_hat = 0.5*(F_L + F_R) - 0.5*R@(d*|lambda|_fix*alpha)
```

This arm has the same initial classical Roe flux and parameter count, but the
learned total wave multipliers cannot become negative.  Nonnegativity fixes
the characteristic upwind direction; it does not by itself prove positivity
or entropy stability, so the shared hard safety stack remains active.

The acoustic eigenvalues use a Harten-style sonic entropy fix.  The contact
speed remains `abs(u_roe)` so a stationary contact is not artificially given
an acoustic dissipation floor.

Both central-Roe arms satisfy `F_hat(U,U)=F(U)` exactly because `alpha=0` for
an equal interface.  The direct-complete Roe-coordinate arm does not.

### 4. Nonnegative automatic upwinding with proposal feasibility learning

This arm keeps design 3 unchanged and continues to use the hard-projected flux
for every finite-volume update.  It adds a loss only on the raw proposal:

```text
r_raw = (v_R - v_L)^T F_raw - (psi_R - psi_L)
loss = trajectory_loss + lambda_feas * mean(relu(r_raw)^2)
F_update = hard_entropy_projection(F_raw)
```

Thus the forward entropy guarantee is identical to design 3, while an active
raw violation receives a gradient in the missing entropy-normal direction.
The existing `central_roe_upwind_broad` checkpoint is the sole architectural
control.  A second weight of `1e-4` is retained as a scale-sensitivity check;
both arms use the identical Central + nonnegative Roe model.  No HLLC+Roe arm
is retrained or modified in this experiment.

## Selection protocol

- 580 training and 136 disjoint validation trajectories.
- Fine targets: strict 512-cell Rusanov + SSP-RK2, conservatively restricted
  to 64 cell averages.
- 6,627 trainable parameters per model at width 72.
- Maximum 50,000 optimizer updates.
- Stop only after the existing validation-rollout plateau rule reaches its
  minimum learning rate.
- Retain only the minimum-validation-NRMSE checkpoint of a converged run.
- Canonical Riemann cases and 512-cell transfer are evaluation only.

## Run

```powershell
python run_roe_upwind_ablation.py --self-test
python run_roe_upwind_ablation.py --arm roe_complete_broad --seed 0
python run_roe_upwind_ablation.py --arm central_roe_signed_broad --seed 0
python run_roe_upwind_ablation.py --arm central_roe_upwind_broad --seed 0
python run_roe_upwind_ablation.py --arm central_roe_upwind_feas_broad --seed 0
python run_roe_upwind_ablation.py --arm central_roe_upwind_feas_light_broad --seed 0
python plot_roe_upwind512.py --seed 0
python plot_upwind_feasibility512.py --seed 0
python summarize_results.py --seed 0
python audit_results.py --seed 0
```

The focused feasibility-loss result and its cross-resolution accuracy tradeoff
are reported in `FEASIBILITY_RESULTS.md`.  The two retained effective methods
are Central + nonnegative Roe + proposal-feasibility loss (`lambda_feas=1e-3`)
and the existing HLLC + Roe correction.  The zero-penalty nonnegative model is
kept as the architectural control, while the signed central-Roe arm remains an
oscillatory ablation rather than a recommended method.  See `RESULTS.md` for
the TV/extrema evidence and the causal qualification.
