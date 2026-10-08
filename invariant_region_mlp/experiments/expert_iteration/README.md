# Expert-iteration preregistered study

`PREREGISTRATION.md` is the frozen protocol. The current authorized run is its
Study-A decision subset at `d = 20, 100, 160`; Study B is intentionally not
started.

From `invariant_region_mlp/experiments/expert_iteration/`:

```powershell
python run_study_a.py reference
python run_study_a.py train
python run_study_a.py run --phase decision
python run_study_a.py collect
python analyze.py
```

All stages are resumable. A completed reference, network checkpoint, or row is
left untouched. The `decision` phase includes the final-checkpoint primary cell,
the two earlier surrogate checkpoints needed by A-5, and the preregistered
exact-Laplacian control at `d=100`, seed 0. Secondary MLP cells are available via
`--phase secondary` but are outside the currently authorized decision run.

Results are written only under `invariant_region_mlp/results/expert_iteration/`;
the report is `invariant_region_mlp/docs/EXPERT_ITERATION_REPORT.md`.

