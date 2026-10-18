# Reference accuracy audit

All references in this suite are independent of the learned solvers and are
opened only after a method's training and model selection are frozen. The
machine-readable test-bank and mathematical checks are in
`PREFLIGHT_AUDIT.json`; historical artifacts are protected by the hashes in
`HISTORICAL_REUSE_AUDIT.json`.

| PDE | Reference | Refinement or identity check | Current uncertainty proxy | Status |
|---|---|---|---:|---|
| HJB-Log | Hopf–Cole reduction with scaled Gauss–Laguerre quadrature | order 64 versus 128: max difference `3.55e-15`; order 128 versus 256: `3.55e-15`; order 128 versus adaptive quadrature on eight states: `5.33e-15` | `5.33e-15` on the audited queries | passed |
| OU/GMM | closed-form evolving diagonal Gaussian mixture | exact PDE residual RMS `8.28e-13`, maximum `1.19e-12`; NumPy/Torch value agreement `1.14e-13` | `1.19e-12` residual check | passed |
| Joint concave | periodic Fourier Strang solver, 2D | `n=64, dt=5e-4` versus `n=128, dt=2.5e-4`; value RMSE `8.45e-5`, gradient-vector RMSE `8.46e-5`; fine-grid PDE residual RMS `3.18e-7` | value RMSE `8.45e-5` | passed; comparisons below this scale are reference-limited |
| Joint convex | independent periodic solver applied to the original PDE; sign identity separately checked | same two grids; value RMSE `7.78e-5`, gradient-vector RMSE `7.77e-5`; fine-grid PDE residual RMS `4.00e-7` | value RMSE `7.78e-5` | passed; comparisons below this scale are reference-limited |
| Nonconvex coupled | independently refined 2D Fourier reference | old `n=128, dt=1.25e-4`, intermediate `n=192, dt=6.25e-5`, final `n=256, dt=3.125e-5`; intermediate-to-final random-query value RMSE `8.08e-6` | value RMSE `8.08e-6` | passed |
| MLP failure | exact Gaussian absolute-moment formula and its analytic first-coordinate derivative | exact terminal identity and minimum-norm inactive subgradient checks; tangent remainder at nonzero points `1.78e-15` | analytic floating-point evaluation | passed |
| Allen–Cahn, 100D | squared-radius PDE with monotone CN/Strang solver and Richardson extrapolation | successive native differences `3.87e-5`, `1.93e-5`, `9.64e-6`; final coarse `(ds=0.025,dt=1.25e-5)` and fine `(ds=0.0125,dt=6.25e-6)` | `9.64e-6` at `u(0,0)` | passed; invariant interval verified |

The joint-concave and joint-convex working solutions remain inside their
certified tangent domains: `[0.2, 0.8332]` and `[-0.8, -0.2]`, respectively.
The nonconvex coupled grid remains in `[0.2, 0.8]`. Dense lifted bases are
orthonormal to at worst `1.34e-15`, and the convex sign-transform identity is
exact to the tested floating-point precision.

Dimension-specific Allen–Cahn radial references are generated lazily on each
frozen common bank with two levels `(ds=0.1,dt=1e-4)` and
`(ds=0.05,dt=5e-5)`. Each generated artifact records the coarse/fine maximum
and RMS discrepancy and is included only after that audit succeeds. The older,
finer 100D reference above remains the native reproduction anchor. Dimensions
whose reference artifact has not yet been generated remain pending in the
execution manifest; they are not assigned an interpolated or neural truth.

Source evidence:

- `experiments/original_pde_benchmark_comparison/reference_accuracy_audit.json`
- `results/tangent_actor_icml_final_suite/references/joint_concave_k2_audit.json`
- `results/tangent_actor_icml_final_suite/references/joint_convex_k2_audit.json`
- `results/tangent_actor_accuracy_frontier/references/coupled_refinement_audit.json`
- `results/tangent_actor_icml_final_suite/PREFLIGHT_AUDIT.json`
