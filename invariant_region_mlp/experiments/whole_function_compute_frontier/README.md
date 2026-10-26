# Whole-function compute frontier

This isolated module implements the preregistered B6/Dense-Fourier 100D
accuracy--compute study. It imports the mature benchmark and solver code rather
than editing it.

The formal sequence is:

1. `run_frontier freeze`
2. truth-free `smoke`
3. train all six C0 trajectories
4. train all twelve C1/C2 trajectories
5. `run_frontier mark-ready` (the global truth gate)
6. evaluate all trajectories
7. collect Slurm accounting, aggregate, audit, plot, and report

`submit_all.sh` encodes this dependency chain. Outputs are written only under
`invariant_region_mlp/results/whole_function_compute_frontier/`.
