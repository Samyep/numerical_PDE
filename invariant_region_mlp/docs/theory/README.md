# One-step theory: gradient noise, generator bias, certified retraction

- `onestep_gradient_noise.tex` / `.pdf`: Lemma 1 (orthogonal gradient noise has trace (d-1)c_h/M),
  Prop. 2 (quadratic generator: exact Jensen bias), Prop. 3 (norm generator: bias of order sqrt(d/M)),
  Prop. 5 (certified retraction caps the generator error; law independent of d when the certificate lies in span(w)),
  Cor. 7 (one averaged nonlinear correction), Prop. 8 (ridge-linear generators never see the orthogonal noise).
- `onestep_check_table.tex`: generated from `results/theory_onestep/onestep_theory_check.json`.
- Reproduce: `cd experiments/theory_onestep && python verify_onestep_theory.py ../../results/theory_onestep` (about 4 min, 1 CPU),
  then `latexmk -pdf onestep_gradient_noise.tex` in this folder.

Scope: one generator evaluation / one Picard level with the terminal-block estimator under ridge structure.
The recursive amplification is documented empirically in the mechanism-suite reports, not proved here.
The manuscript is not modified.
