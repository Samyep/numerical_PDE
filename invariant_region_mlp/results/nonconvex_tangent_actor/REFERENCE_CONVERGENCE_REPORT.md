# Reference Convergence Report

References use a periodic Fourier pseudospectral discretization in the two active coordinates. The production solution uses Strang diffusion splitting and an explicit midpoint nonlinear step. A classical method-of-lines RK4 solve supplies an integration-scheme cross-check.

| Driver | production (n,dt) | fine→production Value RMSE | Gradient-vector RMSE | PDE residual RMSE | aliasing sensitivity | RK4 scheme difference | range |
|---|---|---:|---:|---:|---:|---:|---|
| coupled | (128,0.000125) | 3.152e-05 | 4.557e-05 | 2.877e-06 | 1.091e-15 | 2.485e-08 | [0.200000, 0.800000] |
| primary | (128,0.000125) | 3.358e-05 | 4.852e-05 | 1.794e-06 | 3.627e-15 | 9.470e-09 | [0.200000, 0.800000] |

Reference arrays are retained as reproducible local artifacts and are never imported by Actor training or checkpoint selection. Only post-freeze evaluation loads them.
