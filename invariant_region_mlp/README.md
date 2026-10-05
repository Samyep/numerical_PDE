# Constraint-Aware Multilevel Picard Methods

**Canonical development location.** This directory is the maintained project. The older standalone \`Samyep/invariant-region-mlp\` repository is archival.

The project now studies a common principle rather than one particular projector:

> **Correct noisy recursive MLP states with a PDE-certified constraint retraction before they are reused by the nonlinear generator.**

For a sibling Monte Carlo batch \(Y=(Y_1,\ldots,Y_B)\) and certified feasible set \(\mathcal C\), we write
\[
\mathcal R_B:\;Y\mapsto \mathcal R_B(Y)\in\mathcal C^B,
\]
with \(\mathcal R_B(Y)=Y\) whenever the batch is already feasible. The correction is inserted inside the stochastic Picard recursion, immediately before evaluation of \(F\).

## Two algorithms

### 1. Samplewise IR-MLP

Each child state is corrected independently:
\[
[\mathcal R_{\rm sep}(Y)]_i=\Pi_{\mathcal C}(Y_i).
\]

This is the separable Euclidean nearest-point retraction. Its advantages are locality, no cross-sample coupling, Fejer/nonexpansive error control for convex sets, and a clean structural rescue theorem.

### 2. Batchwise Gauge IR (Batch-IR)

For a star-shaped/gauge feasible set
\[
\mathcal C=\{y:\gamma_{\mathcal C}(y)\le1\},
\]
define
\[
\Gamma(Y)=\max_i\gamma_{\mathcal C}(Y_i),\qquad
\alpha_B(Y)=\min\{1,\Gamma(Y)^{-1}\},
\]
and
\[
\mathcal R_{\rm batch}(Y)=\alpha_B(Y)Y.
\]

This is a coupled radial retraction: the entire sibling batch uses the largest feasible common scale. It preserves relative directions and ratios within the batch and has shown stronger finite-sample performance than samplewise projection on the tested HJB and nonlinear-funding settings, at the cost of cross-sample coupling and a more direct finite-sample bias.

A mean-preserving batch contraction is retained as a diagnostic baseline rather than a main method.

## Structural counterexample rescue

For the Hutzenthaler--Nguyen high-dimensional HJB counterexample
\[
u_t+\tfrac12\Delta u+
\left(\sum_{j=2}^d|\partial_{x_j}u|^2\right)^{1/2}=0,
\qquad
u(1,x)=|x_1|,
\]
the PDE/control structure certifies
\[
\nabla u_d(t,x)\in S_d:=\operatorname{span}(e_1).
\]
On this subspace,
\[
f_d(v)=0\qquad(v\in S_d).
\]

The exact rescue is therefore not tied to Euclidean projection. Any possibly batch-coupled recursive correction whose returned child gradients lie in \(S_d\) before \(f_d\) is evaluated makes every nonlinear correction vanish pathwise. Samplewise \(Q_d\) projection is one realization. A structural batchwise realization
\[
\widetilde v_i=\alpha_B Q_d v_i,\qquad 0\le\alpha_B\le1,
\]
has the same nonlinear collapse. Pure common scaling \(\alpha_Bv_i\) alone does **not** exact-rescue this counterexample because it does not remove inactive coordinates.

For samplewise \(Q_d\), at \((t,x)=(0,0)\) with \(N=m^n\),
\[
\mathbb E\|\widehat Y^{IR}-Y^*\|_2^2
=\frac{4-2/\pi}{N},
\]
independent of ambient dimension.

## Current evidence

- Counterexample: exact structural rescue, simulations through \(d=1000\), and final-output-only ablation.
- High-dimensional quadratic HJB: both recursive retractions strongly suppress nonlinear amplification; Batch-IR is stronger in the paired low-budget comparison.
- 100D nonlinear funding: Samplewise IR and Batch-IR both improve the baseline; Batch-IR is strongest in the paired comparison.
- Neufeld--Wu 100--300D gradient-dependent problems: Samplewise IR and Batch-IR are numerically identical in value error because both suppress the thresholded spurious generator activation.
- Linear convection-diffusion, credit-risk, and Allen--Cahn controls show little or no effect when the certified constraint is inactive.

See \`docs/constraint_aware_framework.md\` for the unified formulation, \`docs/batch_contraction_report.md\` for the paired batch study, \`docs/counterexample_rescue_theory.md\` for the rescue proof notes, and \`paper/\` for manuscript inserts.
