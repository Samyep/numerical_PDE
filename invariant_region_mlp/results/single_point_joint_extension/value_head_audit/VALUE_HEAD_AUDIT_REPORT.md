# Joint Nonseparable direct value-head audit

## Scope

This audit extracts `actor_r` from the four formal saved prediction NPZ files.
It does not retrain, reselect a checkpoint, load a model, or run FK/MCV.
All metrics use the original 32-state frozen bank. The result is one selected
seed-0 model, not a multi-seed mean and not a 960-state whole-function result.

- Worktree HEAD audited: `b7dd884225936a0d67f770efce1207bc478d4f7c`
- Formal origin tip: `24440691ac7beb9e4f2b8d09db6517ca1fee7248`
- Expected formal tip matched: `True`

## Results

| d | Actor | Direct-head RMSE | Direct-head RelMAE | Full TA FK-MCV RMSE | FH-MLP RMSE | Deep BSDE RMSE |
|---:|:---:|---:|---:|---:|---:|---:|
| 20 | a3 | 2.46e-3 | 3.61e-3 | 1.18e-3 | 1.14e-2 | 2.59e-3 |
| 50 | a4 | 2.97e-3 | 3.55e-3 | 1.06e-3 | 2.56e-2 | 3.76e-3 |
| 100 | a4 | 3.52e-3 | 3.83e-3 | 1.04e-3 | 5.15e-2 | 3.48e-2 |
| 200 | a4 | 6.69e-3 | 8.28e-3 | 1.22e-3 | 1.02e-1 | NF |
| 400 | a4 | 1.28e-2 | 1.57e-2 | 1.24e-3 | 2.13e-1 | NF |

Full TA (FK-MCV) is retained only as an auxiliary reference and is not the
direct value-head metric.

## Baseline-cell outcome

The direct value head wins **8/8** legal baseline cells
against FH-MLP and nonfailed Deep BSDE runs. Deep BSDE numerical failures are
reported as NF and are excluded from the denominator.

## Audit checks

- PASS: `all_five_value_heads_from_formal_selected_checkpoints`
- PASS: `no_test_rmse_used_for_actor_selection`
- PASS: `extraction_only_no_fk_mcv_or_rollout_called`
- PASS: `four_inference_replicas_have_exactly_identical_direct_heads`
- PASS: `all_predictions_match_original_frozen_banks`
- PASS: `all_banks_have_32_states`
- PASS: `all_direct_metrics_finite`
- PASS: `single_training_seed_zero_not_described_as_multiseed`
- PASS: `result_not_described_as_960_state_whole_function_accuracy`
- PASS: `forbidden_legacy_nonseparable_checkpoint_not_used`
- PASS: `formal_origin_tip_matches_expected`

## Replacement assessment

Suitable for replacing the Oscillatory Value-Gradient diagnostic: the direct head wins 8/8 legal FH-MLP/Deep-BSDE cells, while the comparison remains explicitly a selected seed-0, 32-state frozen-bank diagnostic.
