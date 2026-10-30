# Final Component Ablation Study for Tangent-Actor

Generated: 2026-10-29T22:30:28.462844+00:00

## Experimental definition and provenance

The preregistered 4 PDE × 4 method study used the mature 100D budget, four tangent updates, seeds 0/1/2, a frozen 960-state bank, 64 Euler steps, and 1,024 paths per state. Replicate 0 is primary and replicate 1 is an independent inference repeat. No test reference entered fitting, noisy-label checkpointing, or reference-free selection. The convex working solution and gradient were sign-converted back to the original physical variable before metrics.

Exact Full-TA actor-training reuse: 8 cells; newly trained Full-TA: 4 cells. A1 and A2 each contain 12 newly trained cells; A3 has no training.

The historical oscillatory bank had 512 states, so its old RMSE was not reused. The preregistered final-suite generation rule supplied one new frozen 960-state bank shared by all variants. Full actors were re-evaluated on it.

## Primary RMSE table

| Method | Radial Log-Terminal HJB | Concave Cubic–Quadratic | Convex Cubic–Quadratic | Oscillatory Value–Gradient |
|---|---:|---:|---:|---:|
| Tangent-Actor (Full) | 8.677e-04 ± 5.630e-05 | 3.415e-04 ± 9.784e-06 | 2.232e-04 ± 5.019e-06 | 8.579e-04 ± 4.671e-05 |
| w/o MCV | 3.038e-03 ± 5.583e-05 | 1.848e-03 ± 1.507e-05 | 1.563e-03 ± 7.414e-05 | 1.931e-03 ± 4.856e-05 |
| w/o Terminal Anchoring | 1.012e-03 ± 9.352e-05 | 3.555e-04 ± 1.102e-05 | 2.342e-04 ± 7.932e-06 | 9.005e-04 ± 7.225e-05 |
| w/o Tangent Updates | 2.753e-03 ± 3.352e-05 | 3.350e-03 ± 2.382e-05 | 5.140e-03 ± 1.636e-05 | 6.445e-03 ± 1.562e-06 |

For A0–A2, ± is the sample standard deviation over the three primary training-seed RMSEs. A3 is a deterministic terminal-extension actor: its center and ± summarize two independent inference-replicate RMSEs, not fabricated training seeds. The replicate predictions remain separate and were never averaged pathwise.

## Component interpretation

### MCV

- Radial Log-Terminal HJB: paired ablated-minus-full RMSE +2.170e-03 on average; mean ratio 3.51×.
- Concave Cubic--Quadratic: paired ablated-minus-full RMSE +1.506e-03 on average; mean ratio 5.41×.
- Convex Cubic--Quadratic: paired ablated-minus-full RMSE +1.339e-03 on average; mean ratio 7.00×.
- Oscillatory Value--Gradient: paired ablated-minus-full RMSE +1.073e-03 on average; mean ratio 2.25×.

### Terminal anchoring

- Radial Log-Terminal HJB: paired ablated-minus-full RMSE +1.446e-04 on average; mean ratio 1.17×.
- Concave Cubic--Quadratic: paired ablated-minus-full RMSE +1.403e-05 on average; mean ratio 1.04×.
- Convex Cubic--Quadratic: paired ablated-minus-full RMSE +1.106e-05 on average; mean ratio 1.05×.
- Oscillatory Value--Gradient: paired ablated-minus-full RMSE +4.265e-05 on average; mean ratio 1.05×.

### Iterative tangent learning

- Radial Log-Terminal HJB: fixed-A0 / Full mean-RMSE ratio 3.17×. A3 remains a fixed linearization and is not FH-MLP.
- Concave Cubic--Quadratic: fixed-A0 / Full mean-RMSE ratio 9.81×. A3 remains a fixed linearization and is not FH-MLP.
- Convex Cubic--Quadratic: fixed-A0 / Full mean-RMSE ratio 23.03×. A3 remains a fixed linearization and is not FH-MLP.
- Oscillatory Value--Gradient: fixed-A0 / Full mean-RMSE ratio 7.51×. A3 remains a fixed linearization and is not FH-MLP.

The formal curve below evaluates every stored Full-TA Actor iterate with the same 960-state bank, 64 Euler steps, 1,024 paths, and paired Brownian seed. The iteration diagnostic was executed on cuda; all five iterates for a seed use the same device-specific random stream. Each learned update receives the mature 5,000-step optimization budget, so A1 is the controlled one-update comparison rather than a prematurely stopped checkpoint.

| PDE | A0 | A1 | A2 | A3 | A4 | Selected Full vs A1 |
|---|---:|---:|---:|---:|---:|---:|
| Radial Log-Terminal HJB | 2.717e-03 ± 1.5e-05 | 1.539e-03 ± 1.3e-04 | 1.093e-03 ± 4.4e-05 | 9.822e-04 ± 9.4e-05 | 8.677e-04 ± 5.6e-05 | +43.6% |
| Concave Cubic--Quadratic | 3.340e-03 ± 2.3e-05 | 3.633e-04 ± 1.2e-05 | 3.463e-04 ± 1.0e-05 | 3.423e-04 ± 8.8e-06 | 3.417e-04 ± 9.9e-06 | +6.0% |
| Convex Cubic--Quadratic | 5.147e-03 ± 1.5e-05 | 2.387e-04 ± 4.0e-06 | 2.262e-04 ± 3.2e-06 | 2.240e-04 ± 3.5e-06 | 2.232e-04 ± 5.0e-06 | +6.5% |
| Oscillatory Value--Gradient | 6.424e-03 ± 3.1e-05 | 1.642e-03 ± 1.0e-04 | 9.878e-04 ± 1.6e-04 | 8.848e-04 ± 4.9e-05 | 8.579e-04 ± 4.7e-05 | +47.8% |

- Radial Log-Terminal HJB: the reference-free selected multi-update Actor reduces mean RMSE by 43.6% relative to A1; selected iteration(s) across seeds: A4.
- Concave Cubic--Quadratic: the reference-free selected multi-update Actor reduces mean RMSE by 6.0% relative to A1; selected iteration(s) across seeds: A3, A4.
- Convex Cubic--Quadratic: the reference-free selected multi-update Actor reduces mean RMSE by 6.5% relative to A1; selected iteration(s) across seeds: A4.
- Oscillatory Value--Gradient: the reference-free selected multi-update Actor reduces mean RMSE by 47.8% relative to A1; selected iteration(s) across seeds: A4.

These curves resolve whether the existing mature protocol benefits from later tangent updates. They do not establish that no alternative, longer single-update optimization schedule could attain the same accuracy; that stronger claim would require a separate compute-matched continuation control.

## Crossed MCV diagnostic

The existing `w/o MCV` Actors were trained and selected entirely with plain labels. Each frozen Actor was also evaluated with FK–MCV on exactly the same Brownian paths already stored by the formal evaluator. This yields the missing crossed cell—`w/o MCV training` with MCV inference—without retraining or test-based reselection. The comparison therefore separates MCV in Actor learning from MCV in final inference.

| PDE | Full (train + infer MCV) | w/o infer | w/o train | w/o all MCV | Infer ratio | Train ratio |
|---|---:|---:|---:|---:|---:|---:|
| Radial Log-Terminal HJB | 8.677e-04 ± 5.630e-05 | 3.032e-03 ± 4.755e-05 | 8.854e-04 ± 5.201e-05 | 3.038e-03 ± 5.583e-05 | 3.50× | 1.02× |
| Concave Cubic--Quadratic | 3.415e-04 ± 9.784e-06 | 1.845e-03 ± 1.480e-05 | 3.556e-04 ± 9.588e-06 | 1.848e-03 ± 1.507e-05 | 5.41× | 1.04× |
| Convex Cubic--Quadratic | 2.232e-04 ± 5.019e-06 | 1.562e-03 ± 7.526e-05 | 2.297e-04 ± 5.236e-06 | 1.563e-03 ± 7.414e-05 | 7.00× | 1.03× |
| Oscillatory Value--Gradient | 8.579e-04 ± 4.671e-05 | 1.913e-03 ± 5.211e-05 | 8.962e-04 ± 3.704e-05 | 1.931e-03 ± 4.856e-05 | 2.23× | 1.05× |

Seed- and replicate-level paired values are in `MCV_INFERENCE_ONLY_SEED_METRICS.csv`; the separate summary is in `MCV_INFERENCE_ONLY_SUMMARY.csv`. The filenames are retained for artifact compatibility, but both files now contain the full 2×2 training/inference crossing.

### MCV training-label effect

The table below supplements the direct crossed value comparison with Actor diagnostics. Gradient RMSE is post-selection only and never entered fitting or checkpoint selection. `MCV-evaluated ratio` compares the no-MCV-trained Actor against Full TA while holding the final MCV estimator fixed.

| PDE | Full Actor gradient RMSE | w/o-MCV Actor gradient RMSE | MCV-evaluated ratio |
|---|---:|---:|---:|
| Radial Log-Terminal HJB | 4.716e-02 | 4.761e-02 | 1.02× |
| Concave Cubic--Quadratic | 2.371e-03 | 5.250e-03 | 1.04× |
| Convex Cubic--Quadratic | 2.457e-03 | 5.414e-03 | 1.03× |
| Oscillatory Value--Gradient | 8.943e-03 | 1.240e-02 | 1.05× |

MCV-trained Actors have lower gradient RMSE on every PDE (about 1% for Log HJB, 55% for both joint cubic--quadratic PDEs, and 28% for the oscillatory PDE). With the final estimator held to FK–MCV, removing training MCV raises mean value RMSE by only about 2%–5% across the four PDEs. This is directionally consistent but small relative to three-seed variability in some columns; it must not be described as a universally resolved gain. In contrast, removing inference MCV raises RMSE by about 2.2×–7.0× and is the dominant measured MCV contribution under this budget.

## Costs, failures, and numerical resolution

Actual Slurm allocated GPU-hours across the recorded jobs: 22.064. The preregistered hard ceiling was 36. No variant-specific stabilization or post-test tuning was permitted.

Numerical failures recorded in the 16 summary cells: 0. Minimum exponential-weight ESS fraction: 0.9495.

| PDE | Reference uncertainty (value RMSE) | Terminal-anchoring mean delta |
|---|---:|---:|
| Radial Log-Terminal HJB | 1.514e-15 | +1.446e-04 |
| Concave Cubic--Quadratic | 8.454e-05 | +1.403e-05 |
| Convex Cubic--Quadratic | 7.784e-05 | +1.106e-05 |
| Oscillatory Value--Gradient | 8.080e-06 | +4.265e-05 |

The anchoring ablation is consistently worse in mean RMSE, but the concave and convex deltas are below their independent reference-resolution estimates and are therefore unresolved. The oscillatory delta exceeds its reference floor but is comparable to three-seed variability. Log HJB has negligible quadrature-reference uncertainty and shows the clearest anchoring effect, although it remains modest relative to the MCV and tangent-update effects. The MCV and fixed-A0 degradations are much larger than all listed reference floors.

The exact reference metadata for every seed is retained in `ABLATION_SEED_METRICS.csv`; the preflight reference audits are in `PREFLIGHT_AUDIT.json`. Small apparent rank changes comparable to reference resolution or independent inference-replicate fluctuation are treated as unresolved.

## Limitations

Three training seeds are sufficient for an ablation signal but not a universal variance claim. The formal A0–A4 curves use the same frozen 960-state banks, 64 Euler steps, 1,024 paths, and paired evaluation seeds as the primary study for all four PDEs. They compare stored per-update Actors under a common inference protocol; they do not prove that an arbitrarily longer-trained single update could never match later updates. These ablations neither establish universally lower variance than FH-MLP nor identify fixed A0 with FH-MLP.

## Artifacts

- `results/final_component_ablation/ABLATION_SEED_METRICS.csv`
- `results/final_component_ablation/ABLATION_SUMMARY.csv`
- `results/final_component_ablation/MCV_INFERENCE_ONLY_SUMMARY.csv`
- `results/final_component_ablation/ACTOR_ITERATION_SUMMARY.csv`
- `results/final_component_ablation/tables/final_component_ablation.tex`
- `results/final_component_ablation/tables/component_ablation_appendix.tex`
- `results/final_component_ablation/tables/mcv_inference_only.tex`
- `results/final_component_ablation/figures/component_ablation.pdf`
- `results/final_component_ablation/figures/actor_iteration_ablation.pdf`
