# Related-work and contribution audit

Revision: 2026-10-05. This is an editorial/source note, not an additional
contribution claim or a compiled appendix. The manuscript is `main.tex`.

## Positioning used in the paper

The contribution is the relation between PDE-certified value–gradient geometry
and nonlinear recursive feedback: a common correction interface, an exact
repair of a documented recurrence, and paired comparisons of separable and
batch-coupled corrections. Projection, clipping, and the modified-driver
inheritance argument are not presented as new primitives.

| Nearest work | Existing ingredient credited | Distinction retained here |
|---|---|---|
| Beck et al. (2020), truncated MLP | Scalar interval truncation inside the driver, with PDE bounds and a nonexpansive modified-driver argument | Structured value–gradient geometry and gradient-feedback rescue; batch-coupled versus separable corrections |
| Gradient-dependent MLP literature | Joint value/gradient stochastic solvers and conditional approximation/complexity results | Uses, rather than replaces, that numerical engine |
| Fan et al. (ICLR 2026), SCaSML | Learned-surrogate defect correction; experimental clipping | Studies feasible child-state geometry; does not claim a full trained-SCaSML comparison |
| Min and Azizan, HardNet v4 | Hard-constrained learned outputs | Recursive Monte Carlo states, rather than network capacity or training |
| Zhong (2025) | Constraints on a reflected controlled spatial process | Constraints on the numerical solution–gradient estimate |
| Hutzenthaler and Nguyen (2025) | The dimensionality counterexample and its low-dimensional target | Sufficient generator-compatible repair condition and exact projected MSE, with reduced MC retained as reference |

## Primary sources consulted

- Beck et al.: https://arxiv.org/html/1907.06729v1
  (Section 3, Eq. (76), Lemma 3.2); journal DOI
  https://doi.org/10.1515/jnma-2019-0074
- E et al. (2021): https://doi.org/10.1007/s42985-021-00089-5
- Hutzenthaler and Kruse (2020): https://doi.org/10.1137/17M1157015
- Hutzenthaler, Jentzen, and Kruse (2022):
  https://doi.org/10.1007/s10208-021-09514-y
- Neufeld and Wu: https://arxiv.org/abs/2310.12545
  and https://doi.org/10.1515/jnma-2024-0074
- Neufeld, Nguyen, and Wu (2025):
  https://doi.org/10.1016/j.jco.2025.101946
- SCaSML conference record:
  https://iclr.cc/virtual/2026/poster/10008446
  Experimental clipping was checked in the authors' manuscript:
  https://arxiv.org/html/2504.16172v3
- HardNet (version 4): https://arxiv.org/abs/2410.10807v4
- State-constrained drift control: https://arxiv.org/abs/2510.21607
- Dimensionality counterexample: https://arxiv.org/abs/2506.23969

## Bibliography corrections

The SCaSML entry now uses the ICLR conference title, *Physics-Informed
Inference Time Scaling for Solving High-Dimensional Partial Differential
Equations*. The HardNet entry explicitly cites version 4 (2025) and its
current two-author list, Youngjae Min and Navid Azizan, rather than mixing
metadata from an earlier version. The full-history gradient-MLP paper in
Foundations of Computational Mathematics and the reflected-control MLP paper
are added. The bibliography contains ten entries and is alphabetically ordered.

## Boundaries that must survive later revisions

- No “first projection/clipping in MLP” claim.
- No “first PDE-derived truncation” claim.
- No general final-error dominance follows from one-state projection geometry.
- Population target preservation alone is not statistical consistency.
- The exact MSE belongs to the specified orthogonal subspace projection,
  including its final gradient; arbitrary compatible retractions guarantee
  nonlinear collapse, not that same exact MSE.
- A dimension-uniform complexity claim requires the hypotheses of the selected
  MLP theorem in its own norm and coordinates.
- The batchwise general consistency/cost-to-accuracy problem remains future work.
- The structural certificate is known; there is no active-subspace discovery claim.
- A complete justification of the particular finance envelope remains pending.
- This is a targeted comparison of the closest retrieved primary works, not a
  proof that no other paper has any overlapping component.

## Asset preservation

Every existing figure environment, caption, label, external plot asset,
original graphic, archived graphic, and figure-data file is unchanged from
`Constraint_Aware_MLP_ICML2026_with_figures.zip`. There are still 15 panels in
ten figure groups. The earlier missing overshoot-energy original remains
pending; see `FIGURE_MANIFEST.md`. No new experiments were run.

This local revision has not been pushed to GitHub.
