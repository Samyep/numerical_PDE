# Direct value-head full-dimensional audit

This is a frozen-artifact analysis: no model was retrained, no checkpoint was reselected, and no test value entered training or selection. The scope is the five PDEs whose audited Tangent-Actor architecture has a value output, at dimensions 20, 50, 100, 200, and 400. Gradient-only HJB Actors are excluded because they do not have a value head.

## Result

The direct head beats at least one independently trained baseline in 24/25 cells (96%), a majority of the three baselines in 17/25 cells (68%), and all three baselines in 6/25 cells (24%). Across all matched pairwise comparisons it wins 47/75 (62.7%): 8/25 against DPI, 15/25 against PINN-HTE, and 24/25 against SDGD-PINN.

The strongest defensible wording is therefore that the auxiliary value head is already a competitive whole-function solver without FK inference. The stronger claim that it outperforms almost all methods is not supported, primarily because DPI is better in 17/25 cells. The formal Tangent-Actor FK--MCV estimator remains better than the direct head in all 25 cells.

## Scope and comparability

- Every comparison within a PDE/dimension cell uses the exact same frozen   state bank and numerical reference.
- The analysis is seed 0, matching the official single-seed benchmark.
- The 100D Oscillatory value-head metric is recovered from the selected-a4   post-freeze record; its test time, states, and truth arrays are byte-identical   to the canonical 512-state benchmark artifact.
- Fourier--Tanh rows come from the audited ICE expansion result at Git   commit `d43b45bdcf068f8fb1e09cf2d9a7e1d5bfd5c038`.

## Files

- `VALUE_HEAD_FULL_DIMENSION.csv`: raw 25-cell comparison
- `VALUE_HEAD_SUMMARY.json`: aggregate win counts
- `SOURCE_AUDIT.json`: hashes and source provenance
- `tables/value_head_full_dimension_icml.tex`: compact ICML single-column table
- `tables/value_head_full_dimension.tex`: publication table
- `tables/value_head_full_dimension_preview.pdf`: rendered preview
