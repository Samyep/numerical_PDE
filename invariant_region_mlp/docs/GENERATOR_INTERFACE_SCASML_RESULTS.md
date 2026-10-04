# Generator-interface projection with faithful SCaSML defect clipping

## Correction to the previous controlled prototype

The public `ScaSML_full_history.py` clips every returned defect state

`(u_breve,z_breve) <- clip((u_breve,z_breve), -equation.uncertainty, equation.uncertainty)`.

The first controlled Hard-SCaSML prototype accidentally replaced this SCaSML-specific defect clipping by the looser MLP-style clipping used in the standalone IR-MLP experiments. That made the recursive gradient defect much larger than in the actual SCaSML mechanism and exaggerated both the benefit on HJB and the damage on VB.

This round keeps the original SCaSML defect clipping in **every** method and changes only where an additional certified invariant-region projection is inserted.

## Compared methods

1. `SCaSML`: official-style defect clip only.
2. `Final-only`: same SCaSML recursion; project the total state only after inference.
3. `Generator-only`: keep the clipped SCaSML defect unchanged; immediately before evaluating the nonlinear generator, form total state `surrogate + defect`, project that total state, and evaluate `f` there.
4. `Recursive total`: add certified total-state projection at recursive return while retaining the official defect clip.

The same fixed surrogate is shared within each comparison. Main track uses corrected EBL normalization.

## Results

### Gradient-dependent nonlinear / VB

| d | SCaSML | generator-only | change | generator projection activation |
|---:|---:|---:|---:|---:|
|20|0.078613885|0.078613885|~0|0.486%|
|40|0.070402290|0.070402291|~0|0.328%|
|60|0.061174813|0.061174814|~0|0.263%|
|80|0.060820637|0.060820637|~0|0.226%|

The catastrophic VB degradation seen in the earlier prototype disappears completely once the official `uncertainty=1e-2` defect clipping is retained.

### HJB

Quick diagnostic with the public `uncertainty=0.1` defect clip:

| d | surrogate | SCaSML | generator-only | activation |
|---:|---:|---:|---:|---:|
|100|0.94133|0.93286|0.93286|0%|
|160|0.96255|0.95717|0.95717|0%|

The official defect clip is already tight enough in this controlled surrogate experiment that the certified HJB total-state projector never activates. Therefore the previous ~40--50% HJB Hard-SCaSML gain was not a faithful comparison to the public SCaSML implementation.

### Diffusion-reaction

| d | SCaSML | generator-only |
|---:|---:|---:|
|100|0.0115602|0.0115307|
|160|0.0190288|0.0189281|

Only small sub-percent changes; no strong claim.

### Linear convection-diffusion

LCD remains an exact neutral control: SCaSML, final-only, generator-only and recursive projection all give machine-precision error.

## Finance controlled experiments

The public SCaSML paper does not contain these two finance benchmarks, so there is no official `uncertainty` setting. Using the linearized Feynman--Kac surrogate and the same controlled defect solver:

- 100D funding, 6 paired seeds: SCaSML MAE 0.1299; generator-only MAE 0.1208 (~7% preliminary improvement).
- 100D credit risk: projector never activates and all methods are identical.

These remain exploratory until a principled SCaSML-style defect clipping scale is specified.

## Scientific takeaway

The user's concern was correct to test, but the faithful result is reassuring:

- **Do not remove or replace SCaSML's own defect clipping.**
- With that mechanism preserved, an additional generator-interface certified projector is mostly a no-op on the tested controlled surrogates and does not damage VB.
- The next decisive experiment is still the trained-PINN reproduction. A trained PINN may create total states that violate the certified region even though the defect itself is clipped to `uncertainty`. Only that experiment can show whether the certified guardrail adds value beyond the original SCaSML heuristic clip.

The combined-method claim should therefore be postponed until the official-style PINNs are retrained and a single checkpoint is shared between SCaSML and generator-projected SCaSML.
