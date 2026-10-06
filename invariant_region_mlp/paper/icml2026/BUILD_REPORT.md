# Build and verification

- Input: `Constraint_Aware_MLP_ICML2026_batchwise_extended.zip`.
- Engine: pdfLaTeX, ICML 2026 review mode, unchanged font sizes and margins.
- Main body ends on page 7; references page 8; appendix pages 9–27.
- Final working-tree builds: no LaTeX warnings, undefined references/citations,
  missing graphics, overfull boxes, or underfull boxes.
- Figures: 34 distinct panels / 35 placements, including one inline TikZ diagram.
- 44 incoming graphical/data assets verified byte-for-byte unchanged.
- All 15 new distinct plot PDFs are referenced in the LaTeX source.
- 36-setting HJB summaries independently recomputed from raw arrays: discrepancy 0.
- PDF was rasterized and page layouts visually inspected, including all new panels.
- Source ZIP is also checked after independent extraction (see release check).

No external font files are shipped; normal PDF font embedding is used by LaTeX.
