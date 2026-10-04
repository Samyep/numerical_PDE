# Entropy and OOD stress benchmark

This experiment is deliberately separate from the earlier matched random
64-cell benchmark.  Its exploratory design and implementation record is in
[`PROTOCOL.md`](PROTOCOL.md); it is not an immutable preregistration.

Run the Euler stress screen from the repository root with the CUDA-enabled
environment used by this project:

```powershell
python `
  numerical_PDE\hard_constrained_flux_learning\experiments\entropy_ood_benchmark\run_euler_stress.py
```

The first phase consumes only pre-existing checkpoints and writes a reference
cache, per-case metrics, an aggregate JSON file, and profile figures under
`results/frozen64`.

Train and evaluate the 2-D SWE comparison:

```powershell
python `
  numerical_PDE\hard_constrained_flux_learning\experiments\entropy_ood_benchmark\run_swe_radial.py `
  --phase all
```

The phases `data`, `hcfl`, `operators`, and `evaluate` can also be run
independently.  Formal output is written under `results/swe_radial`.  The
reference-grid refinement check is:

```powershell
python `
  numerical_PDE\hard_constrained_flux_learning\experiments\entropy_ood_benchmark\audit_swe_reference_convergence.py
```

The consolidated scientific interpretation, including negative results and
scope limitations, is in [`RESULTS.md`](RESULTS.md).  All tuning and corrective
runs are retained in [`TUNING_LEDGER.md`](TUNING_LEDGER.md).

Run the machine-readable integrity audit after reproducing the results:

```powershell
python `
  numerical_PDE\hard_constrained_flux_learning\experiments\entropy_ood_benchmark\audit_results.py
```

The audit checks the Euler and SWE case counts, Roe reconstruction, hard
projection residual, exact HCFL entropy tolerance, reference refinement, and
the pinned revisions of the external implementations.  Its output is written
to `results/AUDIT.json`.
