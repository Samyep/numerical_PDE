# Counterexample rescue: audited experiment report

## What was checked
The source is Hutzenthaler and Nguyen, arXiv:2506.23969v1 (2025), Theorem 1.1 and Eq. (6), pp. 2--3. Its precise quantifiers are: for every p>=0 and every integer n>2p, d^(-p) times the RMSE of U^d_(n,n)(0,0) tends to infinity (liminf). It does NOT say one fixed n outgrows every polynomial.

We implement that recursion with independent full-history branches, r=Uniform(0,1)^2, the centered terminal term, and the correct inverse-time density. The l=0 term is identically zero and skipped for every method. No clipping, shared-branch shortcut, antithetic sampling, or time floor is used for the raw baseline. Every compared recursive variant uses the same random tree. Inactive spatial positions need not be stored, but all d Brownian weights and noisy gradient coordinates are sampled and propagated.

## Protocol
36 configurations: d in {2,5,10,20,50,100,200,500,1000}, n=m in {1,2,3,4}, 512 independent estimator replicas each. The time and first-coordinate streams are also shared across dimensions. Thus an exactly flat empirical IR curve across d is expected under this coupling, not nine independent confirmations. Other coordinates use an independent random stream.

The modes are raw, recursive subspace, recursive unit ball, recursive coordinate box, generator-only subspace, final-only subspace, and dimension-reduced Monte Carlo. Both the ball and box contain the true gradient. The reduced MC control is intentionally included: on this family it equals full IR exactly.

## Main results: full state RMSE, n=m=3

| d | Raw Eq.(6) | Final-only subspace | Generator-only | Unit ball | Recursive subspace | Reduced MC |
|--:|--:|--:|--:|--:|--:|--:|
| 2 | 1.195661 | 0.939134 | 0.411586 | 0.776278 | 0.356391 | 0.356391 |
| 10 | 12.767196 | 5.411684 | 0.674330 | 1.358303 | 0.356391 | 0.356391 |
| 100 | 504.759377 | 66.088470 | 1.958497 | 1.426688 | 0.356391 | 0.356391 |
| 1000 | 16995.682733 | 697.972247 | 6.102105 | 1.433745 | 0.356391 | 0.356391 |

The exact expected IR state MSE is (4-2/pi)/m^n. It requires projection of the returned gradient, not only its use in f. At n=m=3 the theoretical RMSE is 0.3529442421. The empirical 512-replica RMSE is 0.356391 (bootstrap 95% interval below). Value-only RMSE also improves, so the result is not merely zeroing error coordinates in the final metric.

## Theory agreement and exact ablations

| n=m | N=m^n | Exact IR RMSE | Empirical IR RMSE | Bootstrap 95% interval |
|--:|--:|--:|--:|:--|
| 1 | 1 | 1.833952 | 1.797422 | [1.608331, 1.990312] |
| 2 | 4 | 0.916976 | 0.945043 | [0.874965, 1.017019] |
| 3 | 27 | 0.352944 | 0.356391 | [0.337274, 0.375606] |
| 4 | 256 | 0.114622 | 0.117464 | [0.111243, 0.123584] |

For generator-only projection, the value estimator is identical to full IR, but the returned inactive gradient coordinates remain noisy. Its exact full-state MSE is (d+3-2/pi)/N. For final-only projection, the value is identical to raw MLP: feedback corruption has already happened.

## Additional mathematical result
An orthogonal projection Q=BB^T, B^TB=I_r, commutes pathwise with the projected recursion whenever g_d(x)=g_r(B^Tx) and f_d(y,Bv)=f_r(y,v). The full projected result is exactly the r-dimensional MLP result embedded as (U,BV). This includes a nonzero active nonlinearity; projection need not turn the whole PDE linear. For the counterexample r=1 and f_r=0, so all nonlinear terms vanish.

## Validation
- 12 collapse/nonzero-active-generator checks in counterexample_rescue.py passed.
- 17 checks in validate_rescue.py passed (9 independent literal Eq.(6) comparisons; 8 conservative rank-r commutation comparisons).
- Maximum literal-form difference: 1.78e-15.
- IR nonlinear generator and IR-minus-projected-terminal maximum: exactly 0 in the main sweep.
- Separate Gaussian-moment validation: 65,536 replicas at each N in {1,4,27,256}; empirical MSE within 0.56% of the exact formula in all four settings.
- All main state estimates were finite. Bootstrap intervals are descriptive, per configuration; no omnibus significance claim or multiplicity-adjusted hypothesis test is made.

## Cost and limitations
The empirical recursive comparison deliberately traverses the full tree for IR too. The combined vectorized runtime is not a per-method wall-clock comparison. After proving the zero correction identity, a pruned implementation needs N terminal samples; its cost is O(dN) when d-dimensional Brownian vectors are generated, or O(N+d) with coordinate elimination and a dense output vector. The O(d epsilon^-2) bound belongs to this optimized implementation, NOT to the unpruned tree at arbitrary (n,m).

For the diagonal choice n=m and n>=2 minimal with n^n >= (4-2/pi)/epsilon^2, overshoot of the target N can cost a logarithmic factor. An arbitrary integer N (choose n=1,m=N) avoids that rounding factor.

The source lower bound is not contradicted: this is a different, symmetry-informed scheme. The target itself is intrinsically one-dimensional. Reduced MC is equally good, and an analytic Gaussian formula is even available. This is a controlled explanation of an algorithmic failure, not a claim that a previously intractable PDE has been solved or all MLP dimensionality problems disappear.

The previous draft's HJB/finance results are retained as historical exploratory results, not rerun in this update. Euclidean Lipschitz extension alone does not establish dimension-uniform MLP complexity in a maximum-norm theorem. The abstract unbiased-noise Picard inequality is not automatically a theorem for nested full-history MLP.

## Files
`counterexample_rescue.py`: full experiment, `validate_rescue.py`: independent validation, `per_replica_metrics.csv.gz`: all replica metrics, `main/summary.json`: aggregate results, `independent_validation.json`: validation. Plotting is provided separately. Full gradient vectors can be regenerated from the recorded streams and code; compact data retain the exact per-replica quantities used for all reported errors.