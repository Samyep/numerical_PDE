# Audit notes: counterexample revision, 2026-10-04

This update adds new mathematics and newly executed experiments. It does not label historical tables as newly replicated.

## Corrections to earlier discussion

1. **Quantifiers.** Hutzenthaler--Nguyen Theorem 1.1 is `for every p>=0, every integer n>2p, liminf d^-p RMSE(U_(n,n)) = infinity`. A fixed n is not asserted to beat all polynomial exponents.
2. **Algorithm comparison.** The lower bound concerns the paper's unprojected full-history recursion on the diagonal n=m. Our scheme supplies certified symmetry and changes the recursion. It does not refute the lower bound or prove a lower bound for every possible use of MLP.
3. **Exact MSE requires output projection.** `(4-2/pi)/N` applies when the returned gradient is projected. Generator-only projection leaves `(d-1)/N` extra gradient MSE. Final-output-only projection cannot change the computed scalar value.
4. **Computational cost.** The O(d epsilon^-2) statement refers to a pruned implementation after proving every correction zero, with arbitrary N. The experiments traverse the full tree for matched-budget recursive comparisons. Diagonal n=m rounding may add a logarithmic overhead.
5. **Reduced MC is not hidden.** It equals our projected estimator pathwise on this intrinsically one-dimensional problem. We make no efficiency or novelty claim over exploiting known dimension reduction itself.
6. **Full-history MLP versus a noisy Picard cartoon.** The abstract conditionally unbiased-noise inequality is not automatically a convergence theorem for the nested MLP estimator. The rescue proof is a separate exact finite-tree argument.
7. **Norms and moving sets.** Euclidean Lipschitz continuity after projection does not automatically meet dimension-uniform maximum-norm assumptions. State-dependent projectors require spatial regularity and metric bookkeeping. We removed any implied automatic general complexity transfer.
8. **Finance metric.** The prior script rescales delta and maps back to z; this is a weighted projection, not generally an ordinary Euclidean projection onto the z ellipsoid. The paper now states the relevant metric.
9. **Historical correlations.** The earlier shared-data radius/budget sweep is descriptive. Its p-values should not be treated as independent-sample causal evidence.
10. **Bibliographic correction.** The SCaSML author order is Zexi Fan, Yan Sun, Shihao Yang, Yiping Lu, as in arXiv:2504.16172v3.

## Fresh empirical basis

- Original source: arXiv:2506.23969v1, pp. 2--3, Eq. (6).
- 36 configurations, dimensions 2--1000, diagonal n=m=1--4, 512 replicas each.
- Seven reported variants, all compared with paired input randomness; five recursive variants use the same full tree.
- 12 main-script checks and 17 independent validation checks pass.
- Independent moment validation: 65,536 replicas at each N=1,4,27,256.
- All main estimates finite; all recursive-subspace nonlinear corrections exactly zero.
- Main code and result files include seeds, precision, environment, scalar work counts, and source SHA-256.

## Historical evidence retained

HJB, funding, Neufeld--Wu, and negative-control tables were copied from the supplied `IR_MLP_Overleaf_updated.zip` and previous reports. They were NOT rerun in this update. The stable HJB reference numbers already in that draft remain the historical version, and older unscaled reference tables remain superseded. No full official SCaSML pipeline or new broad active-nonlinearity benchmark was run in this revision.