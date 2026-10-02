# Scientific-integrity note: why HCFL-512 can beat HLLC-512

## Short answer

There is no teacher-performance paradox in the recorded experiment.  The
network was **not** given native 512-cell states and trained to reproduce a
native HLLC-512 step.  A 512-cell strict Rusanov simulation first generated a
fine trajectory, and every snapshot was then conservatively restricted to 64
finite-volume cell averages.  The network saw only those 64-cell states.  It
learned a local correction to an HLLC flux that approximates the evolution of
the restricted fine solution, i.e. a coarse-grid closure.

The later 512-cell comparison is zero-shot resolution transfer of that local
flux law.  Its comparator is a different operator, native HLLC-512, and its
common reference is HLLC-2048 conservatively restricted to 512 cells.  Nothing
in supervised learning requires the transferred closure to inherit the error
of native HLLC-512.

The current seed-0 result is good controlled evidence that the learned
correction is useful.  It is not yet a proof of general superiority.

## What was actually learned

Let `R_512->64` denote conservative restriction and let `S_Rus,512` denote one
saved-interval update of the strict 512-cell Rusanov solver.  The training
pairs have the form

```text
input  = R_512->64 U_n^512
target = R_512->64 S_Rus,512(U_n^512).
```

They are not pairs of native HLLC-512 inputs and outputs.  In particular,

```text
R_512->64 S_Rus,512(U) != S_HLLC,64(R_512->64 U)
```

in general.  Their difference is precisely a discretization/filtering closure
that contains information about fine-grid fluxes inside each coarse cell.
This is why a model operating on 64 degrees of freedom can improve on a native
64-cell solver without reconstructing all 512 fine values.

At deployment the selected dissipation model uses

```text
F_HCFL = F_HLLC - 1/2 R diag(d_theta |lambda|) R^{-1} Delta U,
```

where the network predicts bounded characteristic-wise multipliers `d_theta`
from a local five-cell stencil.  This form has three useful transfer
properties:

1. On a constant state, `Delta U = 0`, so the correction is exactly zero.
2. In a smooth region, `Delta U = dx U_x + O(dx^2)`, so the correction
   naturally shrinks when the grid is refined even though `dx` is not an
   explicit network input.
3. At a shock, `Delta U` remains order one and the network can change the
   characteristic dissipation, shock width, and phase error relative to HLLC.

The flux is translation equivariant and shared at every interface.  Deployment
on 512 cells uses eight updates per saved interval, preserving the training
ratio `dt/dx`; it does not interpolate a 64-vector into a 512-vector.

## Controls that rule out the obvious artificial explanations

The matched comparison freezes the selected checkpoint and evaluates every
method against the same conservative HLLC-2048-to-512 reference:

| method | five-case mean rollout NRMSE |
|---|---:|
| classical HLLC-512 + SSP-RK2 | 0.07424 |
| matched hard-safe HLLC-512, zero neural correction | 0.07191 |
| learned HCFL-512 | **0.05627** |

The learned model improves on the exactly matched zero-correction control by
21.75% on the mean and by 8.34%--36.33% on every one of the five canonical
cases.  The zero-correction control keeps the same HLLC proposal, grid,
time-update schedule, Tadmor projection, positivity path, and fully discrete
entropy path.  Therefore the remaining difference is the learned correction,
not a switch of base Riemann solver or safety implementation.

The positivity and fully discrete entropy fallback limiters did not intervene
on any of these learned 512-cell rollouts.  Maximum conservation drift was
below `1e-6`.  Hence the improvement is not produced by repeatedly falling
back to a hidden low-order solver.  The hard Tadmor projection remains part of
HCFL and is reported separately; its total relative flux RMS change was at
most 0.128% in a case.

## Leakage and selection audit

- Training generators use seeds `6000`, `7000`, and `8000` for the selected
  broad-data arm.  The alternative wave arm additionally uses `9000`.
- Validation uses independent seeds `14000`, `15000`, `16000`, and `17000`.
- The five named canonical states are constructed only in the final evaluation
  suite.  They are not used for gradients, learning-rate changes, early
  stopping, or checkpoint selection.
- Every retained model weight file is the minimum of its recorded independent
  validation-rollout curve, and every arm reached the declared plateau rule.
- Normalization statistics come from the training tensor only.
- The HLLC-2048 comparison is post-hoc: repository history records the frozen
  checkpoint before the matched comparison code and result artifacts.
- The model rollout receives only its current numerical state.  No reference
  trajectory or future value is passed to it.
- All methods in the matched table use the same initial state, final time,
  reference tensor, and error normalization.
- Exact SHA-256 initial-state checks find zero overlap between the regenerated
  training, validation, and canonical sets.

Run the executable audit with:

```powershell
python audit_scientific_integrity.py --seed 0
```

It writes `results/scientific_integrity_audit_seed0.json`, including split
fingerprints, validation-selection checks, checkpoint/result hashes, nonzero
correction diagnostics, matched-control results, limiter use, and conservation
checks.  Fixed 1,100-update states remain in the CSV audit trail because they
show why the old fixed-budget conclusion was invalid, but their nonconverged
weight files are intentionally not retained.

## What this result does and does not establish

The defensible claim is:

> For seed 0 and the five held-out canonical periodic tests, the frozen local
> dissipation correction transfers from 64-cell training states to 512 cells
> and reduces error by 21.75% relative to an otherwise identical zero-NN HLLC
> control, while retaining admissibility and conservation.

It would be too strong to claim that HCFL is generally better than high-order
classical solvers.  The evidence is still one training seed, five named tests,
and a finite HLLC-2048 reference.  A publication-level claim should require
multiple converged training seeds, a locked or private test generator, a
reference-resolution study (for example HLLC-4096 or converged WENO), and
matched MUSCL-HLLC/WENO accuracy-cost baselines.

If the goal changes to learning a native 512-cell closure, the appropriate
training target is a finer trajectory (for example 2048 or 4096 cells)
conservatively restricted to 512 cells.  That is a different experiment from
the current resolution-transfer test and should be labeled separately.
