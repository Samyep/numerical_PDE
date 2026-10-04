# 2-D Euler Clawpack quadrant audit

This directory contains the prospective 64-grid HCFL experiment based on the
official PyClaw quadrant gallery problem.  Read `PROTOCOL.md` before running it
and `RESULTS.md` for the audited outcome.  The result improves the HLLC base
but does not outperform native PyClaw Roe-64; this limitation is retained.

The Clawpack dependency is isolated in a Linux container:

```text
docker build -f Dockerfile.clawpack -t hcfl-clawpack:5.10.0 ../../..
```

The large generated datasets and prediction arrays are intentionally not
versioned.  The small selected checkpoints, scripts, metadata, scalar results,
plots, and audit reports are versioned so a clone can reproduce the evaluation
without retraining.

The main stages are:

```text
python run_experiment.py audit-data
python run_experiment.py train --stencil 4 --seed 0
python run_experiment.py train --stencil 6 --seed 0
python run_experiment.py select
python run_experiment.py train --stencil 6 --seed 1
python run_experiment.py train --stencil 6 --seed 2
python run_experiment.py evaluate --stencil 6 --seeds 0 1 2
python audit_results.py
```
