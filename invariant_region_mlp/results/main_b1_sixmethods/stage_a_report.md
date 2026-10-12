# B1 six-method Stage A report

**Verdict:** PASS. All six methods ran on all three pilot PDEs, all errors were finite, and checkpoint/resume artifacts were written.

Completed 18/18 method-cells.
Measured cumulative method time: 0.0446 GPU h.

| d | k | Vanilla | Double | DPI | Deep BSDE | PINN+SDGD | Ours |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 1 | 5135 | 2.802 | 0.3178 | 0.337 | 0.0673 | 0.2171 |
| 100 | 10 | 2153 | 0.9812 | 0.6085 | 1.076 | 0.1747 | 0.2475 |
| 100 | 20 | 1051 | 0.5586 | 0.5741 | 0.7625 | 0.1329 | 0.1399 |

The conservative serial Stage B estimate is 6.056 GH200 h with a 10 h wall limit.

The table reports value RMSE. Vanilla is numerically unstable at this shallow budget. PINN+SDGD is best at k=1 and k=10; Ours and PINN+SDGD are close at k=20. No configuration was removed or changed in response.

Missing cells are recorded as `N/A` and are never converted to numerical scores.
