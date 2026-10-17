# Reference Convergence Report

References use a periodic Fourier pseudospectral discretization in the two active coordinates. The production solution uses Strang diffusion splitting and an explicit midpoint nonlinear step. A classical method-of-lines RK4 solve supplies an independent integration-scheme cross-check.

Every row uses four refinement levels `(n,dt)=(48,1e-3),(64,5e-4),(96,2.5e-4),(128,1.25e-4)`. The three reported refinement columns compare coarse, medium, and fine solutions with production.

| Driver | lambda | eta | coarse / medium / fine Value RMSE vs production | fine Gradient-vector RMSE | PDE residual RMSE | aliasing sensitivity | RK4 scheme difference | terminal Value / Gradient RMSE | range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| coupled | 0.2 | 1 | 2.060e-04 / 1.020e-04 / 3.152e-05 | 4.557e-05 | 2.877e-06 | 1.091e-15 | 2.485e-08 | 6.274e-05 / 8.965e-05 | [0.200000, 0.800000] |
| primary | 0 | 0 | 2.547e-04 / 1.299e-04 / 3.976e-05 | 5.690e-05 | 8.002e-07 | 4.456e-15 | 4.653e-09 | 6.309e-05 / 8.935e-05 | [0.200000, 0.800000] |
| primary | 0.05 | 0 | 2.457e-04 / 1.236e-04 / 3.840e-05 | 5.462e-05 | 9.370e-07 | 2.026e-15 | 5.476e-09 | 6.298e-05 / 8.914e-05 | [0.200000, 0.800000] |
| primary | 0.2 | 0 | 2.192e-04 / 1.086e-04 / 3.358e-05 | 4.852e-05 | 1.794e-06 | 3.627e-15 | 9.470e-09 | 6.274e-05 / 8.965e-05 | [0.200000, 0.800000] |
| primary | 0.5 | 0 | 1.811e-04 / 8.932e-05 / 2.804e-05 | 3.993e-05 | 5.005e-06 | 3.107e-15 | 2.166e-08 | 6.335e-05 / 8.929e-05 | [0.200000, 0.800000] |

The fine-to-production Value discrepancy is the conservative reference-error proxy used in post-freeze reports. It is below one tenth of the selected method RMSE in every reported run, including the highest-accuracy main run.

Reference arrays are retained as reproducible local artifacts and are never imported by Actor training or checkpoint selection. Only post-freeze evaluation loads them.
