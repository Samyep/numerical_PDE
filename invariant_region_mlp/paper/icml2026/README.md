# Constraint-Aware MLP — complete evidence edition

Upload this entire ZIP to Overleaf. Select `main.tex` and pdfLaTeX.
The project is self-contained for LaTeX compilation and uses the existing
ICML 2026 anonymous-review style. Do not remove the `figures`, `tables`, or
`sections` folders. No external figures or shell escape are required.

## Contents

- `main.tex`: seven-page main body, references, and original appendices.
- `sections/92_batchwise_extended.tex`: the prior Batch-IR extension.
- `sections/93_restored_evidence.tex`: overshoot reproduction, recovered
  funding and surrogate diagnostics, batch controls, and the bounded-scope
  elliptic diagnostic.
- `FIGURE_MANIFEST.md`: exact figure-to-page map, including old Figure 7.
- `RESTORATION_NOTES.md`: what was newly executed versus recovered.
- `data/restored_evidence/`: results, raw reproduction arrays, source bytes,
  provenance, and machine-checkable preservation/audit records.
- `archive_superseded/`: old invalid-reference evidence, retained but not
  shown as a current result.

Current PDF: 27 pages total; body ends on page 7, references on page 8.
There are 34 distinct panels (35 placements because one original panel
also appears in the appendix). All 19 incoming panels remain.

## Local build

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Check data and preservation (Python + numpy/scipy)

```sh
python scripts/audit_restored_evidence.py
```

## Regenerate only the new figures (also needs matplotlib)

```sh
python scripts/plot_restored_evidence.py
python scripts/make_restored_tables.py
```

## Independently rerun the new overshoot experiment

```sh
python scripts/reproduce_overshoot.py --out new_overshoot_run
```

The rerun command does not overwrite the archived input arrays. All other
added PDE results were recovered rather than rerun. The controlled surrogate
is not a trained SCaSML checkpoint, and the elliptic experiment is not MLP.

This edition has not been pushed to GitHub. The source experiments were read
at the commit recorded in `data/restored_evidence/PROVENANCE.json`.
