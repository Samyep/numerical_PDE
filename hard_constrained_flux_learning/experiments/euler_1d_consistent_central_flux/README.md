# Exactly consistent central-flux baseline

## Pre-run question

Does the failure of direct complete-flux learning come mainly from asking the
network to relearn the physical/central part of the Euler flux?  This arm puts
that part back analytically while retaining a direct vector correction and no
Roe decomposition:

```text
F_hat = 0.5 * (F(U_L) + F(U_R))
        + 0.18 * ||P_R - P_L||_normalized * scale * tanh(MLP(stencil))
```

The jump norm contains no epsilon.  Therefore the learned correction is
exactly zero whenever `U_L == U_R`, independently of the outer stencil, and

```text
F_hat(U, U) = F(U)
```

holds by construction before and after Tadmor projection.

## Controlled factors

Held fixed: five primitive cells as input, width 72, 6,627 trainable
parameters, broad seed-0 data, independent validation selection, optimizer,
learning-rate/plateau rule, hard Tadmor projection, local positivity limiter,
fully-discrete entropy limiter, and periodic finite-volume update.

Changed relative to `direct_broad`: HLLC is replaced by the arithmetic average
of the left/right physical fluxes.  Changed relative to `full_flux_broad`: the
physical central flux and exact consistency are restored analytically.

The 50,000-update cap is only an upper bound; training stops as soon as the
unchanged validation plateau rule is satisfied.  Only a converged,
validation-selected checkpoint is retained.

## Run

```powershell
python run_consistent_central_flux.py --self-test
python run_consistent_central_flux.py --seed 0
python plot_results.py --seed 0
python plot_central_consistent512.py --seed 0
```
