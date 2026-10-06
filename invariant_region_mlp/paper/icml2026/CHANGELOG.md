# Complete evidence restoration

See `RESTORATION_NOTES.md` and `FIGURE_MANIFEST.md`. All 19 incoming panels and 44 graphical/data assets are preserved; 15 new distinct figures are included. The old overshoot diagnostic is independently reproduced with raw arrays, not digitized. Main: 7 pages; full: 27 pages. No GitHub write was performed.

---

# Extended Batch-IR revision

- Imported results from commit 771f26594aa7ee364d426ad29744b7716d76cdc3.
- Added four new main-text panels and Appendix E protocols / seven tables.
- Kept all 15 earlier panels and all source assets.
- Distinguished uncorrected vs heuristic baselines, pooled increment variance vs averaged-correction variance, and value-only vs full-state checks.
- Related work and headline theory remain intact.

---

# Related-work and positioning revision — 2026-10-05

- Added independent Related Work immediately after the introduction.
- Replaced the old combined Related Work/Limitations section with a separate
  Limitations section, avoiding duplicate survey paragraphs.
- Rewrote novelty positioning and three contributions around certified
  value–gradient geometry, structural rescue, and mechanism-resolving evidence.
- Credited truncated MLP for scalar PDE-based truncation and the same
  nonexpansive modified-driver inheritance argument.
- Credited SCaSML's clipping, rather than treating clipping as a new ingredient.
- Distinguished hard-constrained network outputs and physical state-constrained
  reflected control from numerical solution–gradient constraints.
- Added the full-history gradient MLP reference and the reflected-control
  reference; corrected SCaSML and HardNet bibliographic metadata.
- Kept convergence in Appendix A and strengthened its norm/dimension caveat.
- Explicitly stated that the exact full-state MSE uses projection of the final
  gradient, as well as the recursive child returns. No new theorem was added.
- Kept all 15 visual panels, figure captions, graphic bytes, and data unchanged.
  The missing historic overshoot-energy original is still pending.
- No new numerical experiment; no GitHub push in this revision.

---

# Changes from the seven-page tight ICML draft

- Restored the recursive-interface diagram as a real vector diagram, not a
  boxed line of text; queued it early enough to appear at the method entry.
- Restored the nine-dimension counterexample state and value-only plots, with
  all correction variants and legends.
- Restored high-dimensional HJB and Neufeld–Wu outcome plots alongside the
  internal generator-error and recursive-variance diagnostics.
- Restored both batchwise tradeoff plots, including the funding mean-preserving
  control rather than presenting only selected table entries.
- Restored the depth sweep, independent MSE validation, primary RMSE check,
  funding-budget curve, violation-rate diagnostic, and Allen–Cahn control in
  the relevant appendix subsections. Main text points to them explicitly.
- Retained the numeric tables, moving redundant main-text tables into the
  experimental appendix instead of removing their information.
- Corrected a caption conflating 512-replica primary-sweep RMSE with 65,536-
  replica independent MSE. Replaced the erroneous independent-validation
  wording with two separate captions and a short protocol explanation.
- Corrected the parentheses in the HJB terminal condition to agree with its
  archived code and existing derivative formula. No experiment was rerun.
- Kept ordinary convergence inheritance in the appendix and batch-coupled
  statistical complexity as future work; no new headline theorem was added.
- Preserved the finance certificate's formal-justification caveat and the
  distinction between weighted and Euclidean projection metrics.
- Missing original overshoot-energy graphic is documented explicitly in
  `FIGURE_MANIFEST.md`; summary statistics were not used to invent its points.
