# Batch contraction experiments

This experiment asks whether IR-MLP's gain is specific to samplewise Euclidean projection or can be reproduced by batch-level variance contraction.

Methods:
- raw: no correction;
- samplewise_ir: project each child state independently;
- mpbc: mean-preserving batch contraction z_i' = zbar + alpha(z_i-zbar), using the largest common alpha that makes the batch feasible; if the empirical mean is itself infeasible, no all-feasible mean-preserving correction exists and the batch is left unchanged;
- uniform_shrink: z_i' = alpha z_i with one common alpha set by the most extreme violation in the batch.

All corrections are applied to recursive child states immediately before the nonlinear generator F. The root output is not clipped. Random paths are paired across methods.

Run the module or use the committed summary in results/batch_contraction/summary.json. See docs/batch_contraction_report.md for the interpretation.
