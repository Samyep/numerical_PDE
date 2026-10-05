# Unified constraint-aware MLP framework

## 1. Recursive constraint retractions

Let
\[
Y=(Y_1,\ldots,Y_B)
\]
denote a sibling Monte Carlo batch of recursive MLP child-state estimates and let \(\mathcal C(t,x)\) be a PDE-certified feasible state set. A **recursive constraint retraction** is a map
\[
\mathcal R_B:\mathcal Y^B\to\mathcal C^B
\]
that is the identity on feasible batches. It is applied to the noisy child states **before** the nonlinear generator \(F\).

This separates the structural requirement---return to a certified set before nonlinear reuse---from the statistical geometry used to perform the correction.

### Samplewise metric retraction

For closed convex \(\mathcal C\),
\[
[\mathcal R_{\rm sep}(Y)]_i=\Pi_{\mathcal C}(Y_i)
\]
solves
\[
\arg\min_{Z\in\mathcal C^B}\sum_i\|Z_i-Y_i\|_2^2.
\]
It is local and nonexpansive. If \(Y_i^*\in\mathcal C\),
\[
\|\Pi_{\mathcal C}Y_i-Y_i^*\|^2
\le
\|Y_i-Y_i^*\|^2-\operatorname{dist}(Y_i,\mathcal C)^2.
\]

### Batchwise gauge retraction

If
\[
\mathcal C=\{y:\gamma_{\mathcal C}(y)\le1\}
\]
is star-shaped about the chosen center, define
\[
\Gamma(Y)=\max_i\gamma_{\mathcal C}(Y_i),\qquad
\alpha_B=\min\{1,\Gamma(Y)^{-1}\},
\]
then
\[
[\mathcal R_{\rm batch}(Y)]_i=\alpha_BY_i.
\]
This uses the largest common feasible scale. It couples sibling samples and generally changes the empirical batch mean, but it preserves their relative directions and ratios.

The two algorithms therefore share the same feasible product set \(\mathcal C^B\), the same insertion point in MLP, and the same fixed-point target; they differ only in the retraction geometry.

## 2. Structural rescue principle

Let \(S\) be a PDE-certified set satisfying a generator-compatibility condition. In the strongest null-set case,
\[
F(y)=F(y^*)\qquad\forall y\in S.
\]
If a recursive correction, possibly batch-coupled, satisfies
\[
[\mathcal R_B(Y)]_i\in S
\]
for every child before \(F\) is evaluated, then the spurious part of the nonlinear recursive feedback is removed pathwise.

For the dimensionality counterexample,
\[
S_d=\operatorname{span}(e_1),\qquad
f_d(z)=\|z_{2:d}\|_2,
\]
so \(f_d|_{S_d}=0\). Hence every structurally compatible recursive retraction into \(S_d\) makes the nonlinear full-history corrections vanish.

Two realizations are
\[
z_i\mapsto Q_dz_i
\]
(samplewise) and
\[
z_i\mapsto\alpha_BQ_dz_i
\]
(batch-coupled structural contraction). Pure radial scaling \(z_i\mapsto\alpha_Bz_i\) is not sufficient because it does not force \(z_{i,2:d}=0\).

For the samplewise realization at the origin,
\[
\mathbb E\|\widehat Y-Y^*\|_2^2
=\frac{4-2/\pi}{N}.
\]
For a structural batchwise contraction with \(0\le\alpha_B\le1\), the same pathwise nonlinear collapse holds.

## 3. Empirical tradeoff

The paired batch study gives a complementary statistical message. For ordinary ball/ellipsoid constraints away from the exact counterexample, common-factor Batch-IR can outperform nearest-point Samplewise IR.

HJB, \(n=2,M=6\):
- \(d=20\): Samplewise IR 0.2412 vs Batch-IR 0.0966 relative \(L_2\).
- \(d=40\): Samplewise IR 0.1993 vs Batch-IR 0.0833.

100D nonlinear funding:
- \(M=10\): 0.5102 vs 0.4807 MAE.
- \(M=20\): 0.2965 vs 0.2539.
- \(M=40\): 0.1993 vs 0.1486.

Thus the manuscript does not claim Euclidean projection is universally optimal. The contribution is the broader principle of **constraint-aware recursive correction**, with two concrete retraction geometries exposing a structural-guarantee / finite-sample bias--variance tradeoff.
