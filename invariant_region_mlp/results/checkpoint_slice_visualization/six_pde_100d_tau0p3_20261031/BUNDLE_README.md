# Frozen 100D slice data at `t/T = 0.3`

This portable bundle contains every numerical array used to render the six
long frozen-checkpoint comparison figures below:

- `hjb_log`: 100D Radial Log-Terminal HJB
- `ou_gmm`: 100D OU-Gaussian Mixture HJB
- `joint_concave`: 100D Concave Cubic-Quadratic PDE
- `joint_convex`: 100D Convex Cubic-Quadratic PDE
- `nonconvex_coupled`: 100D Oscillatory Value-Gradient PDE
- `allen_cahn`: 100D Allen-Cahn

All PDEs use dimension 100, seed 0, and the same frozen 81-by-81 spatial slice
within each PDE. The four methods are `tangent_actor`, `dpi`, `sdgd_pinn`, and
`pinn_hte`. No checkpoint, training cache, retraining, or model selection is
included. Full-Rank Log-Sum-Exp and Fourier-Tanh are generated independently
on ICE and are intentionally not mixed into this DeltaAI bundle.

## Layout

```text
references/   six frozen physical reference slices
predictions/  24 method NPZ files with value/gradient/error arrays
metrics/      per-method RMSE metadata and completion records
audit/        checkpoint provenance and canonical-bank reproduction audits
figures/      PDF, 300-dpi PNG, and SVG previews
scripts/      portable renderer and JSON/CSV consolidator
```

`DATA_MANIFEST.json` records file sizes, SHA-256 hashes, NPZ shapes, dtypes,
finite-value checks, and direct error-array consistency checks. `SHA256SUMS`
covers every payload file except itself.

## Reproduce

```bash
python -m pip install numpy matplotlib
python scripts/render_long_pde_figures.py
python scripts/consolidate_results.py
```

The renderer performs no smoothing, denoising, interpolation, or checkpoint
selection. Solution panels use reference-derived limits; discrepancy panels
share logarithmic scales across the four methods within each PDE.
