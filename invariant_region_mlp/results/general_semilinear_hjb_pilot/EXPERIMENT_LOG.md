# Experiment Log

- 2026-10-14: Created branch `ir-mlp-general-semilinear-hjb-pilot` from the completed compute-frontier branch.
- 2026-10-14: Implemented the generic controlled-HJB interface, full-trajectory VJP labels, constrained Actor map, Riccati reference, and independent scalar finite-difference references.
- 2026-10-14: CPU validation completed: `11 passed`; scalar reference refinement artifacts saved.
- 2026-10-14: Slurm script passed local validation and scheduler `--test-only`; submitted job `3381749` to one GH200.
- 2026-10-14: Job `3381749` completed in 6m20s. Repeated validation/final evaluation as job `3381795` (3m48s) after correcting antithetic uncertainty to use pair means; policy Value estimates were unchanged.
- 2026-10-14: The reference-free a3 rule selected LQG and Quartic but not Bounded. Targeted job `3381820` completed in 2m00s.
- GPU job accounting at report generation:

```text
3381749|COMPLETED|00:06:20|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381749.batch|COMPLETED|00:06:20|0:0|363328K|cpu=8,gres/gpu=1,mem=64G,node=1
3381749.extern|COMPLETED|00:06:20|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381749.0|COMPLETED|00:06:10|0:0|2017600K|cpu=8,gres/gpu=1,mem=64G,node=1
3381795|COMPLETED|00:03:48|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381795.batch|COMPLETED|00:03:48|0:0|694.50M|cpu=8,gres/gpu=1,mem=64G,node=1
3381795.extern|COMPLETED|00:03:48|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381795.0|COMPLETED|00:03:36|0:0|1227456K|cpu=8,gres/gpu=1,mem=64G,node=1
3381820|COMPLETED|00:02:00|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381820.batch|COMPLETED|00:02:00|0:0|44480K|cpu=8,gres/gpu=1,mem=64G,node=1
3381820.extern|COMPLETED|00:02:00|0:0||billing=1000,cpu=8,gres/gpu=1,mem=64G,node=1
3381820.0|COMPLETED|00:01:50|0:0|1408896K|cpu=8,gres/gpu=1,mem=64G,node=1
```

- Generated all tables and plots from raw JSON/NPZ outputs, then independently recomputed Value metrics from the saved arrays.
