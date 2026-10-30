# Frozen Slice Data: Three PDEs at `t/T = 0.3`

This portable bundle contains every numerical array used to render the three
long comparison figures for:

- `joint_concave`: 400D Jointly Concave;
- `nonconvex_coupled`: 100D Nonconvex--Nonconcave;
- `ou_gmm`: 100D OU/GMM HJB.

The four frozen methods are `tangent_actor`, `dpi`, `sdgd_pinn`, and
`pinn_hte`.  No checkpoint or training cache is included or required.  The
arrays are frozen post-inference outputs; rerendering performs no model
selection, smoothing, denoising, or interpolation.

## Directory layout

```text
references/   one frozen reference slice per PDE
predictions/  four frozen method predictions per PDE
metrics/      seed/budget metadata, RMSE values, and completion records
audit/        checkpoint provenance and canonical-bank reproduction audits
figures/      PDF, 300-dpi PNG, and SVG previews
scripts/      portable Matplotlib renderer and JSON/CSV consolidator
```

`DATA_MANIFEST.json` records file sizes, SHA-256 hashes, NPZ keys, shapes,
dtypes, finite-value checks, and RMSE values recomputed directly from the
arrays. `SHA256SUMS` covers every payload file except itself.

## NPZ schema

Each `references/<pde>_tau0p3.npz` contains the physical grid and exact or
independently validated reference. Important fields are:

```text
axis                 (81,)
s1, s2               (81, 81)
coordinates          (6561, 2)
t                     (6561,)
x                     (6561, d)
basis                 (d, 2)
x_perp                (d,)
truth                 (6561,) or (81, 81)
truth_gradient        (6561, 2) or (81, 81, 2)
time_fraction         scalar (= 0.3)
horizon               scalar
```

Each `predictions/<pde>_tau0p3_<method>.npz` contains:

```text
prediction            (81, 81)
gradient              (81, 81, 2)
value_error           (81, 81)
gradient_error        (81, 81, 2)
value_discrepancy     (81, 81) = abs(value_error)
gradient_discrepancy  (81, 81) = ||gradient_error||_2
```

The gradients are derivatives with respect to the displayed reduced
coordinates `(s1, s2)`. Gradient RMSE labels use the vector convention
`sqrt(mean(||gradient_error||_2^2))`.

## Reproduce the three figures

From the extracted bundle root:

```bash
python -m pip install numpy matplotlib
python scripts/render_long_pde_figures.py
```

The renderer writes the three PDF/PNG/SVG triplets back to `figures/`. To
render only selected PDEs, repeat `--pde`, for example:

```bash
python scripts/render_long_pde_figures.py --pde joint_concave --pde ou_gmm
```

To regenerate the consolidated CSV files:

```bash
python scripts/consolidate_results.py
```

## Scientific conventions

- All panels use the same 81-by-81 physical slice within a PDE.
- Solution colors use limits derived only from the reference surface.
- Value-discrepancy maps share one logarithmic scale across all four methods
  within a PDE; gradient-discrepancy maps share another.
- RMSE annotations are recomputed from the saved error arrays.
- `t/T = 0.3`; the actual time is `0.15` for the first two PDEs (`T=0.5`)
  and `0.30` for OU/GMM (`T=1`).
