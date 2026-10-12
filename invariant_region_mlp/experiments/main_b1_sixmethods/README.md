# Main B1 six-method benchmark

This directory implements the isolated three-stage benchmark specified in the
task attachment. It does not run Stage C, an `(M,n)` sweep, component
ablations, or SCaSML.

The frozen PDE grid is `d in {100,200,400}` and `k in {1,5,10,20,40}` with
`T=0.25`, `sigma=sqrt(2)`, `K=2k+1`, and coefficient scale 2.0 for `k>1`.
The scale was selected once from exact-reference diagnostics before any solver
results. For `k=1`, the direction and coefficients exactly reproduce P1.

The six method identifiers are:

- `vanilla_mlp`: original centered-BEL full-history MLP.
- `double_mlp`: independent-replicate quadratic-driver MLP.
- `dpi`: Deep Picard Iteration with control-variate BEL gradient labels.
- `deep_bsde`: global/amortized Deep BSDE over `x0 ~ N(0,I_d)`.
- `pinn_sdgd`: hard-terminal PINN with two-sample SDGD.
- `actor_mlp`: pathwise terminal gradient, level-matched baseline, and a
  frozen ambient-input actor inside full-history MLP.

Run from the repository root with the cluster Python environment:

```bash
python -m pytest invariant_region_mlp/experiments/main_b1_sixmethods/test_correctness.py -q
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_stage0
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_stage_a --all
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_stage_b --all
```

The final DPI table uses the frozen resource-adapted profile in
`dpi_final_config.json`. It retains the authors' four-hidden-layer width-512
ELU HJB architecture and a maximum of 5,120 optimizer steps, with 8 Picard
rounds and 32/128 training/validation Monte Carlo samples per label to fit the
fixed allocation. Run its archival and resumable replacement pass with:

```bash
export B1_DPI_PROFILE_PATH="$PWD/invariant_region_mlp/experiments/main_b1_sixmethods/dpi_final_config.json"
python -m invariant_region_mlp.experiments.main_b1_sixmethods.run_final_dpi
```

The command copies the preliminary 4x128 tanh artifacts under
`results/main_b1_sixmethods/superseded/dpi_initial_stage_b/` before replacing
them. Checkpoints are reused only when their complete training signature
matches the final profile.

Every method-cell writes an atomic JSON record and a compressed prediction
archive, and neural methods checkpoint each training seed. Repeating a command
resumes by skipping completed records and loading existing checkpoints.
`--force` deliberately replaces the requested method-cell.

Exact values and gradients are called only by Stage 0 calibration and by
held-out evaluation after a neural checkpoint has been selected. They are not
passed to actor, DPI, Deep BSDE, or PINN training.
