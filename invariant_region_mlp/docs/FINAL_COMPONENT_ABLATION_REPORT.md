# Final Component Ablation Study for Tangent-Actor

Generated: 2026-10-29T06:08:08.224309+00:00

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

## Inference-only MCV diagnostic

`Full TA w/o MCV Inference` is a supplemental diagnostic, not a fifth primary method. It freezes the Full TA Actor trained and selected with MCV, then replaces only the final FK–MCV estimator by plain FK on identical Brownian paths. In contrast, primary `w/o MCV` removes MCV from labels, reference-free selection, and inference.

| PDE | Full FK–MCV RMSE | Full w/o MCV inference RMSE | Plain / MCV |
|---|---:|---:|---:|
| Radial Log-Terminal HJB | 8.677e-04 ± 5.630e-05 | 3.032e-03 ± 4.755e-05 | 3.50× |
| Concave Cubic--Quadratic | 3.415e-04 ± 9.784e-06 | 1.845e-03 ± 1.480e-05 | 5.41× |
| Convex Cubic--Quadratic | 2.232e-04 ± 5.019e-06 | 1.562e-03 ± 7.526e-05 | 7.00× |
| Oscillatory Value--Gradient | 8.579e-04 ± 4.671e-05 | 1.913e-03 ± 5.211e-05 | 2.23× |

Seed- and replicate-level paired values are in `MCV_INFERENCE_ONLY_SEED_METRICS.csv`; the separate summary is in `MCV_INFERENCE_ONLY_SUMMARY.csv`. This comparison isolates inference variance reduction without conflating it with Actor-learning changes.

### MCV training-label effect

The table below separates Actor quality from final-estimator noise. Gradient RMSE is a post-selection diagnostic and never entered fitting or checkpoint selection. `Plain-value delta` is the primary w/o-MCV RMSE minus the same Full Actor's plain-FK RMSE, averaged at the PDE level.

| PDE | Full Actor gradient RMSE | w/o-MCV Actor gradient RMSE | Plain-value delta |
|---|---:|---:|---:|
| Radial Log-Terminal HJB | 4.716e-02 | 4.761e-02 | +5.700e-06 |
| Concave Cubic--Quadratic | 2.371e-03 | 5.250e-03 | +2.581e-06 |
| Convex Cubic--Quadratic | 2.457e-03 | 5.414e-03 | +2.885e-07 |
| Oscillatory Value--Gradient | 8.943e-03 | 1.240e-02 | +1.771e-05 |

MCV-trained Actors have lower gradient RMSE on every PDE (about 1% for Log HJB, 55% for both joint cubic--quadratic PDEs, and 28% for the oscillatory PDE). However, at 1,024 paths the plain-FK value RMSEs of the Full and w/o-MCV Actors differ by less than 1%; plain-estimator noise masks that Actor-quality improvement. Thus the decisive observed value-RMSE gain from MCV is inference variance reduction, while the Actor diagnostics provide separate evidence of cleaner supervision rather than a resolved plain-value improvement.

## Costs, failures, and numerical resolution

Actual Slurm allocated GPU-hours across the recorded jobs: 20.051. The preregistered hard ceiling was 36. No variant-specific stabilization or post-test tuning was permitted.

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

Three training seeds are sufficient for an ablation signal but not a universal variance claim. The oscillatory a0–a4 supporting curve uses the audited historical 512-state bank, while its primary 4×4 cells use the newly frozen 960-state bank; no values are mixed within a curve. Joint iteration curves have only the audited seed-0 evidence where additional seeds were unavailable. These ablations neither establish universally lower variance than FH-MLP nor identify fixed A0 with FH-MLP.

## Artifacts

- `results/final_component_ablation/ABLATION_SEED_METRICS.csv`
- `results/final_component_ablation/ABLATION_SUMMARY.csv`
- `results/final_component_ablation/MCV_INFERENCE_ONLY_SUMMARY.csv`
- `results/final_component_ablation/ACTOR_ITERATION_SUMMARY.csv`
- `results/final_component_ablation/tables/final_component_ablation.tex`
- `results/final_component_ablation/tables/mcv_inference_only.tex`
- `results/final_component_ablation/figures/component_ablation.pdf`
- `results/final_component_ablation/figures/actor_iteration_ablation.pdf`
