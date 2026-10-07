# Mechanism benchmark suite

This directory implements the pre-registered benchmark for noise entering a
nonlinear MLP generator.  It reuses the active-VB `FullHistoryMLP` class
unchanged and overrides only the equation-generic state transformation and
generator hook.

Run from the `numerical_PDE` repository root:

```powershell
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e0 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e1 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e2 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e3 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e5 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_suite --stage e6 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.analyze_suite
```

`--stage all` runs the computational stages in order.  E4 is a certificate-
ablation view derived from the E3 artifacts, so it has no separate runner
stage.  Each repetition is
written atomically to `results/mechanism_suite/raw/`; valid existing artifacts
are skipped, so interrupted runs are resumable.  E1--E3 (and therefore the
derived E4 view) include only primary PDEs admitted by `gates.json`.

The supplied task attachment referred to a `reference_code/` directory, but
that directory was not present.  The implementation therefore follows the
written equations directly; N3 is explicitly recorded as a canonical
reconstruction of the specified multi-direction Hopf--Cole problem.

## Confirmatory round 2

The post-round-1 correction is frozen verbatim in
`ROUND2_PREREGISTRATION.md`.  It uses base seed 20261107, new point sets and a
new validation split, and writes only to `results/mechanism_suite_r2/`:

```powershell
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r3 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r2 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r1 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r4 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r6 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r5 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.run_round2 --stage r7 --workers 8
python -m invariant_region_mlp.experiments.mechanism_suite.analyze_round2
```

`--stage all` follows that exact order. R3 first constructs the 80-node
Gauss--Hermite P4 bound cache and applies the four-grid containment gate. If
the violation does not decrease under refinement or remains at least `1e-5`
on the finest fixed grid, the runner stops before any confirmatory P4 MLP
work. R5 and R7 reuse R1 data exactly where specified.
