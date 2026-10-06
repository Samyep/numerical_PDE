# Allen--Cahn truncated-MLP recovery

This experiment establishes that Beck-style scalar truncated MLP is the
`[-r,r] x R^d` special case of Samplewise IR-MLP and tests two independent
implementations on the published Allen--Cahn setup of Becker et al. (2020).

The recovery track uses the published numerical choices `T=1`, diffusion
`sqrt(2)`, `f(u)=u-u^3`, `g(x)=1/(2+(2/5)||x||^2)`, and fixed truncation radius
`r=4`.  A separate `[0,1]` solution-side invariant interval is reported as a
diagnostic and is never conflated with the published numerical truncation.

From the repository root:

```powershell
python -m invariant_region_mlp.experiments.allen_cahn_recovery.run_recovery equivalence
python -m invariant_region_mlp.experiments.allen_cahn_recovery.run_recovery source
python -m invariant_region_mlp.experiments.allen_cahn_recovery.run_recovery pilot
python -m invariant_region_mlp.experiments.allen_cahn_recovery.run_recovery final
python -m invariant_region_mlp.experiments.allen_cahn_recovery.analyze_recovery
python -m unittest -v invariant_region_mlp.experiments.allen_cahn_recovery.test_allen_cahn_recovery
```

The scripts are resumable: each task writes an immutable stage/configuration
JSON artifact before the aggregate CSVs and figures are produced.
