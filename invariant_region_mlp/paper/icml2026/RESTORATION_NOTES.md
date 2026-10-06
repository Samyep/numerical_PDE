# Restoration and evidence audit — 2026-10-05

## Deliverable

This edition starts from `Constraint_Aware_MLP_ICML2026_batchwise_extended.zip`.
The body remains seven pages; the complete PDF has 27 pages. The original 19
panels remain, 15 new distinct plots are added, and the Neufeld–Wu budget
panel appears twice (main and appendix): 34 unique panels / 35 placements.
Convergence theory and the related-work positioning are not expanded.

## What was actually run

Only the 100D, n=2, 36-setting overshoot diagnostic was newly executed.
It used ten paired repeats per M, six budgets, six valid radii, float64,
1,200 fixed evaluation points, and scaled 64-node Hopf–Cole reference
quadrature cross-checked at 128 nodes. The agreement is 1.776e-15 on those
points. Raw and corrected root values are saved, along with the matrix,
evaluation points, reference arrays, energies, and activity by repetition.
The original lost PNG was not recovered; the new reproduction is identified
as such in the main caption, appendix, metadata, and figure manifest.

The JSON summaries are independently recomputed from the NPZ with maximum
discrepancy zero. Descriptive Spearman correlations are 0.8862290862 and
0.9945945946. They are not independent-sample causal tests or new convergence
claims. Uniform random time concerns the n=2 value-channel experiment, not a
replacement for the gradient-weighted MLP theorem's sampling assumptions.

## Recovered, not rerun

1. Original two-budget Neufeld–Wu results (18 method records, 30 repeats)
   were recovered from the archival repository's parabolic archive. The
   existing plot and six-row table are unchanged; the plot is prominent in
   the main paper again.
2. Funding complete budget records and higher-depth geometry controls were
   recovered from the older source/results archive. They include n=4 and n=5
   and preserve the coordinate-box comparison. No Batch-IR series was added
   to those older draws. The original PNG work envelope is unchanged.
3. The controlled logistic surrogate-quality JSON has all 45 method rows
   across M=1,2,4 and five q settings, with 40 recorded repeats. All values
   and reported SDs are retained; the M=4 cases where IR loses to raw defect
   MLP are explicit. This is not trained-SCaSML performance.
4. Batchwise negative controls were fetched from the main repository's JSON.
   All 27 method records' displayed metrics are extracted without inventing
   per-replica arrays. Credit risk and Allen–Cahn are inactive; the linear
   control has an active gradient correction yet unchanged values.
5. The older batch report supplies the missing mean-feasibility and n=M=3
   Neufeld–Wu plots. These report-rounded values are not pooled with the new
   n=M=2 suite.
6. An elliptic neural-BSDE diagnostic is preserved outside the MLP claim.
   The generator diagnostic is mean absolute error, not MSE. Seed-0 panels
   are distinguished from the three-seed summary. The upstream /d heuristic
   remains better on terminal loss, and the joint ball does not improve the
   single-point value error. It is not sold as a new positive benchmark.

## Source integrity and limits

Main repository snapshot: Samyep/numerical_PDE @
771f26594aa7ee364d426ad29744b7716d76cdc3.
The old parabolic archive's Git blob SHA was verified after decoding:
60168f9e831fb30ba19ee233dcafd9980a6c8f07.
Per-source identifiers and extraction precision are in
`data/restored_evidence/PROVENANCE.json`.

The archive has two textual defects: a literal newline in an unused
Neufeld diagnostic key, and an invalid numeric token in the old,
superseded unscaled-reference violation JSON. The raw bytes are retained;
only the verified standard Neufeld metric fields are extracted with a
permissive JSON reader. The obsolete violation file is never used in a
current figure. Old reference-invalid HJB graphics stay in
`archive_superseded/`, explicitly not current evidence.

All 44 input graphical/data assets in the preservation manifest are
byte-for-byte unchanged. Figure names are not silently reused for new
stochastic results. The main paper keeps the individual-increment-variance
cases unfavorable to Batch-IR as well as all old ablations.

## Publication status

These are local revised manuscript artifacts. No GitHub write or new
repository commit was made during this restoration.
