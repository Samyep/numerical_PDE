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
