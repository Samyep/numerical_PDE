# High-accuracy FH-MLP amortization benchmark

This isolated study compares the frozen Tangent-Actor and audited Full-History
MLP recursion on the 100D Joint Nonseparable PDE at RMSE thresholds `1e-2` and
`5e-3`.  Deep BSDE and all other baselines are out of scope.

The execution gates are:

1. `prepare` creates fresh calibration/final banks and a truth-free timing stream.
2. `preflight` validates the equation, sources, H100, and frozen TA artifacts.
3. `dispatch make-profile` runs the truth-free FH configuration pilot.
4. `profile` freezes the 30 H100-hour cost gate.
5. `dispatch make-calibration`, followed by `freeze`, selects configurations.
6. `dispatch make-final`, followed by `final_accuracy`, opens the disjoint final bank once.
7. `timing_grid` freezes deterministic timing refinements before formal timing.
8. `dispatch make-timing`, `analyze`, `accounting collect`, and `audit` finalize the study.

Every timing marker is an actually executed query prefix.  TA training and
reference-free selection are charged even though byte-identified checkpoints
are reused from the audited predecessor study.
