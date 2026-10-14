# Experiment Log

- 2026-10-14: Created branch `ir-mlp-concave-semilinear-tangent-actor` from commit `96c324b1`.
- 2026-10-14: Implemented generic concave drivers, reusable convex sign transform, weighted Feynman–Kac rollout, full state-dependent VJP, bounded Value tangent head and normalized Value/Gradient supervision.
- 2026-10-14: CPU preflight passed 15 tests; independent periodic reference refinement and residual validation saved.
- 2026-10-14: Ran the fixed a0→a1→a2 pilot on one GH200 allocation. All a2 updates passed the reference-free paired rule, so a targeted second allocation added a3 without using test RMSE.
- 2026-10-14: Generated tables/figures from raw outputs and independently recomputed Value metrics from saved prediction arrays.
- Total allocation time was 31m09s (0.519 single-GH200 GPU-hours); there were no failed or retried jobs.

Scheduler accounting:

```text
3382289|COMPLETED|00:20:58|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3382289.batch|COMPLETED|00:20:58|0:0|540416K|cpu=8,gres/gpu=1,mem=64G,node=1
3382289.extern|COMPLETED|00:20:58|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3382289.0|COMPLETED|00:20:46|0:0|1690560K|cpu=8,gres/gpu=1,mem=64G,node=1
3382500|COMPLETED|00:10:11|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3382500.batch|COMPLETED|00:10:11|0:0|38.50M|cpu=8,gres/gpu=1,mem=64G,node=1
3382500.extern|COMPLETED|00:10:11|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3382500.0|COMPLETED|00:10:01|0:0|1517184K|cpu=8,gres/gpu=1,mem=64G,node=1
```
