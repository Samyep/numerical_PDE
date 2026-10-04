# Hard projection + SCaSML: controlled defect-correction experiment

## What was tested

The recursive variable in SCaSML is a defect `(u_breve,z_breve)`. Hard projection is therefore applied to the **total** state

`u_total = u_hat + u_breve`, `z_total = z_hat + z_breve`,

and the solver returns `project(total) - surrogate`. We use corrected EBL normalization throughout.

Because the public SCaSML repository does not ship trained PINN checkpoints, this round isolates the interaction between defect correction and invariant-region projection using a fixed surrogate shared by both methods. For LCD/VB/DR/HJB the surrogate is the terminal function frozen in time. For the two finance benchmarks, where a frozen terminal payoff is too poor, the surrogate is the linearized Feynman--Kac value `E[g(X_T)]` computed by fixed antithetic Monte Carlo. Thus **the only difference between SCaSML and Hard-SCaSML is the recursive projection**.

## Main result

| PDE | dimensions | SCaSML -> Hard-SCaSML | conclusion |
|---|---|---|---|
| LCD | 10/20/30/60 | both essentially exact | neutral control |
| VB | 20/40/60/80 | hard is 10.5%, 61.9%, 65.8%, 81.2% worse | tight z projection conflicts with an already accurate defect correction |
| DR | 100/120/140/160 | change <0.05% | essentially neutral |
| HJB | 100/120/140/160 | error reduced by 49.2%, 44.4%, 42.0%, 38.6% | strong positive complementarity |
| Funding | 100D | MAE 0.1332 -> 0.1266 | small +5.0% preliminary gain; not significant with 10 seeds |
| Credit risk | 100D | MAE 0.6598 -> 0.6598 | exact negative control; projector never activates |

### HJB details

| d | surrogate rel-L2 | SCaSML | Hard-SCaSML | reduction |
|---:|---:|---:|---:|---:|
|100|0.942|1.811|**0.920**|49.2%|
|120|0.951|1.682|**0.936**|44.4%|
|140|0.956|1.607|**0.932**|42.0%|
|160|0.962|1.531|**0.940**|38.6%|

Here ordinary defect correction overshoots badly because the quadratic generator reuses noisy gradient defects. Projection prevents the total gradient from leaving the certified HJB region. It does not dramatically beat the surrogate itself in this controlled setup; its role is primarily stabilization.

### VB diagnostic

The bad VB result is specifically a gradient-projection effect. At 80D:

| method | rel-L2 |
|---|---:|
| SCaSML | 0.01014 |
| value-only | **0.00982** |
| 1x gradient ball | 0.03528 |
| joint 1x | 0.03527 |
| relaxed 2x gradient ball | 0.01102 |

So the correct conclusion is **not** that certified projection always helps. A tight total-gradient ball can increase finite-sample bias after defect correction even when the true state lies in the ball. Relaxing the ball nearly removes the damage.

### HJB diagnostic

At 140D:

| method | rel-L2 |
|---|---:|
| SCaSML | 1.596 |
| value-only | 1.043 |
| gradient-ball only | 0.957 |
| joint | **0.950** |

Both constraints help, with the gradient constraint contributing more.

## Scientific verdict

1. Hard projection and SCaSML are **strongly complementary on HJB-type quadratic gradient feedback**.
2. They are not universally complementary: VB is a clear counterexample at the certified 1x radius.
3. DR and LCD are useful neutral controls.
4. The previous credit-risk negative control remains negative under defect correction: no violations, no effect.
5. Funding shows only a small preliminary gain once the surrogate is already good; ten paired seeds are insufficient for a positive claim.
6. The paper should therefore narrow the combined-method claim to a **stability mechanism**, not a universal accuracy improvement.

## Important limitation

This is a controlled SCaSML defect-correction experiment, not a bit-for-bit reproduction of the paper's trained DeepXDE PINNs. The public repository contains training code but not the trained checkpoints. The experiment deliberately fixes the same surrogate for baseline and projected defect correction so that the incremental effect of projection is identifiable.
