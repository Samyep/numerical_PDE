# Nonseparable High-Budget Report

The unchanged undamped generic tangent update was trained on `f(u,z)=0.2 sin(2πu)-0.5||z+u e||²` with full ambient inputs. Three seeds use the same mature budget and reference-free selection rule as the primary experiment.

The independent periodic reference was refined to `n=256`, `dt=3.125e-05`. Its conservative previous-to-current Value RMSE is `8.080e-06`, and the computed solution remains in `[0,1]`.

| iteration | seeds | 64-step Value RMSE mean ± sd | Gradient-vector RMSE | Value-Actor RMSE | fixed-point score | MCV VR | ESS |
|---|---:|---:|---:|---:|---:|---:|---:|
| a0 | 3 | 0.0063739 ± 0.0001 | 0.179185 | 0.0986236 | 3.2562 | 1.94× | 0.9526 |
| a1 | 3 | 0.00103878 ± 4.5e-05 | 0.0253496 | 0.0150248 | 0.069572 | 27.64× | 0.9488 |
| a2 | 3 | 0.000897387 ± 3.4e-06 | 0.0121168 | 0.00753993 | 0.017514 | 43.04× | 0.9488 |
| a3 | 3 | 0.000868999 ± 4.3e-06 | 0.00849375 | 0.00544626 | 0.010776 | 47.93× | 0.9488 |
| a4 | 3 | 0.00086301 ± 7.6e-06 | 0.00694584 | 0.00470067 | 0.0061911 | 50.25× | 0.9488 |

Reference-free selection produced 3 frozen seeds. Their mean observed 64-step improvement from `a0` is `86.5%`. This benchmark supplies local empirical evidence for genuinely mixed `(u,z)` curvature, not a universal nonconvex convergence theorem.

## Selected-Actor Euler audit

| steps | RMSE mean ± sd | MC contribution |
|---:|---:|---:|
| 64 | 0.00080926 ± 5.3e-06 | 0.000168166 |
| 128 | 0.000422355 ± 3.3e-06 | 0.000124274 |
| 256 | 0.000229824 ± 5e-06 | 9.51093e-05 |
| 512 | 0.000136337 ± 4.7e-06 | 7.64973e-05 |

## Measured cost

Mean production training plus reference-free selection cost is `3292.2` GPU-seconds per seed; the post-freeze Actor/policy audit costs `629.3` GPU-seconds per seed. Queue time is excluded.

## Strong baseline conclusion

At the finest completed Actor evaluation (512 Euler steps), the three-seed mean RMSE is `0.000136337`. The strongest replicated raw Full-History endpoint is `n=4, M=96` at `0.00163057` (95% t interval `[0.00145455, 0.00180658]`). Tangent-Actor has the lower mean error: **True**.
