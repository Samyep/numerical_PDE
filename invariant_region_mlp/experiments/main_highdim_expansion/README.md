# B6 and Dense-Fourier ICML expansion

This package runs the frozen 40-cell whole-function benchmark requested for
the B6 full-rank HJB and Dense Fourier joint PDE at
`d = 20, 50, 100, 200, 400`.  The four methods are Tangent-Actor + MCV, DPI,
PINN-HTE, and SDGD-PINN, all with training seed 0 and a shared 960-state bank
within each PDE/dimension pair.

The execution order is fixed before formal training:

1. `stage_a.slurm` calibrates the new B6 dimensions, freezes every parameter
   and test bank, runs mathematical and stochastic validation, runs non-result
   baseline smoke tests, characterizes the benchmarks, and freezes a runtime
   profile.
2. The 100D anchor runs all four methods for both PDEs.
3. Tangent-Actor fills the remaining four dimensions for both PDEs.
4. The three fixed source-faithful baseline transfers fill the remaining
   dimensions.

`slurm/auto_expansion_controller.sh` waits until the original seven-PDE suite
has no controller or Slurm job left, executes `sbatch --test-only` before every
actual submission, prevents duplicate suite jobs, and enforces both the
six-GH200-hour cell cap and separate 120-GH200-hour expansion cap.  Accuracy
never enters task ordering or model selection.

The formal Tangent-Actor row is the reference-free selected Actor evaluated at
64 Euler steps with inference replicate 0.  Replicate 1 and the 128/256-step
audits remain separate records.  Exact solutions are unavailable to training,
labels, stopping, and checkpoint selection.  They are opened after selection
for evaluation, and are used before training only for the frozen benchmark
calibration/characterization requested by the protocol.

DPI, PINN-HTE, and SDGD-PINN are labelled `faithful_transfer`.  Their source
budgets and architectures are fixed; no Tangent-Actor terminal or MCV feature
is added.  These rows do not claim to reproduce native experiments from the
source papers on these new PDEs.
