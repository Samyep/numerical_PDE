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

`slurm/auto_expansion_controller.sh` verifies that no original-suite job is
active and, by default, treats this checkout as an independent machine.  Set
`CHECK_OLD_SUITE=1` only on a machine that shares the original-suite queue.  It
executes `sbatch --test-only` before every actual submission, gives every cell
its own single-GPU Slurm job, prevents duplicate writers, and enforces both the
six-GPU-hour cell cap and separate 120-GPU-hour expansion cap.  Accuracy never
enters task ordering or model selection.

All machine paths are derived from the checkout.  `SOURCE_ROOT`, `DPI_ROOT`,
`EXPANSION_VENV`, and the module names can be overridden through the
environment.  `slurm/setup_ice.sh` installs the pinned DPI checkout and its
dependencies once; `slurm/common_env.sh` is then reused by Stage A and every
formal cell.  The ICE controller defaults to eight concurrent H100 cells, with
the account, partition, QOS, GPU type, and concurrency exposed as environment
variables.

Restarting the controller after a partial run validates and reuses an existing
completed Stage A rather than regenerating the frozen inputs. Completed formal
cells are likewise omitted from subsequent queues.

ICE allocations run a real CUDA arithmetic-and-synchronization preflight before
formal work. A narrowly classified `cudaErrorDevicesUnavailable`/driver-startup
failure is recorded as an immutable infrastructure event, charged to resource
accounting, and leaves the cell `CHECKPOINTED` so it can be retried. The
controller quarantines the affected Slurm node for the remainder of the run.
Application exceptions, numerical failures, CUDA OOM, and protocol-gate failures
remain terminal and are never converted into retries. Set
`AUTO_QUARANTINE_NODES=0` to disable automatic quarantine or provide a
comma-separated static list through `EXPANSION_EXCLUDE_NODES`.

These ICE controls do not change the frozen scientific protocol: parameter and
test banks, seeds, architectures, optimizer budgets, cell/global GPU-hour caps,
selection rules, and evaluation metrics are unchanged. Failed and retry
allocations—including jobs named `icml_retry_*`—remain in the 120-GPU-hour
accounting ledger.

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
