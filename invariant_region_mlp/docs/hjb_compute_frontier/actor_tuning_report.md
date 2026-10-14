# Actor tuning report

Status: GPU profiling pending.

The ambient-output Actor has no rank projection on either PDE. The declared study compares terminal-pathwise `q0`, one policy improvement `q1`, a matched fresh-label continuation `q1_refresh`, a true `q2`, and a true `q3`. Each controlled-policy label differentiates the full Euler trajectory with the source Actor frozen while preserving its state derivative.

The pilot uses 2,048 terminal-label points for `q0`, then 4,096 points, 256 paths per label, 16 Euler steps, and 4,000 optimizer steps for each subsequent phase. The final preset raises this to 4,096/256/6,000 for `q0` and 8,192/512/32/6,000 for policy-evaluation phases. A high-quality branch declares width 512, 16,384 points, 1,024 paths, 64 steps, and 12,000 optimizer steps for `q1` and its matched continuation. It will be launched only if measured profiling fits the 18-hour ceiling.

Per the user's explicit override, Actor candidates are ranked by test Relative L2/RMSE. The report will show both model-path cost and total exploration cost so this oracle selection is not presented as a deployable reference-free procedure. External baselines do not receive the same selection privilege.
