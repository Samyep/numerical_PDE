# OU/GMM FH-MLP pre-run audit

**Status: BLOCKED_MISSING_ARTIFACTS.**

## Immutable artifacts

| Artifact | Present | SHA-256 check | Path |
|---|---|---|---|
| main_results | yes | PASS | `C:\Users\Yeping\workspace\HCFL\ir-mlp-ou-gmm-fh-time-accuracy-worktree\invariant_region_mlp\results\tangent_actor_icml_final_suite\MAIN_RESULTS.csv` |
| parameter_npz | no | MISSING | `C:\Users\Yeping\workspace\HCFL\ir-mlp-ou-gmm-fh-time-accuracy-worktree\invariant_region_mlp\results\tangent_actor_icml_final_suite\inputs\ou_gmm_d100.npz` |
| actor_prediction_npz | no | MISSING | `C:\Users\Yeping\workspace\HCFL\ir-mlp-ou-gmm-fh-time-accuracy-worktree\invariant_region_mlp\results\tangent_actor_icml_final_suite\paired_predictions\main_ou_gmm_d100_tangent_actor_s0_n64p0_p1024_mcv.npz` |

The deterministic generator reproduced the frozen 960-point array hash in memory, but this reconstruction was not saved or substituted for a missing original artifact.
The formal GPU sweep remains fail-closed until the original parameter and paired-prediction NPZ files are present.

## Mathematical checks

- `terminal_value_max_abs`: `0.0`
- `terminal_gradient_max_abs`: `0.0`
- `pde_residual_relative_max_200_states`: `4.631436141997569e-15`
- `weighted_ou_bel_relative_max`: `2.026980648162402e-16`
- `ou_variance_formula_max_abs`: `0.0`
- `solver_has_no_nonterminal_exact_oracle`: `True`
- `constant_minus_d_integrated_deterministically`: `True`
- `passed`: `True`

## Provenance

- Unit tests accepted for this FH audit: `True`.
- Unit-test result: `2 failed, 48 passed in 3.79s`.
- The two failures are confined to final-suite DPI tests requiring the absent optional `picard` package; the OU/GMM FH tests themselves pass.
- Git commit: `fffe0fca43fdb6a9f7adaf5eebf318780d7ac141`
- Python: `3.12.11`
- PyTorch: `2.11.0+cu128`
- CUDA device: `NVIDIA GeForce RTX 5070 Ti`
- No FH timing or RMSE point was produced while the artifact gate was closed.
