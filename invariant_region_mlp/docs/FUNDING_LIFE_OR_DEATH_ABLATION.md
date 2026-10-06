# Funding life-or-death ablation: certified geometry vs noisy-z suppression

## Question

Does the positive 100D nonlinear-funding result require the certified normalized-delta ellipsoid, or can the gain be explained by suppressing the noisy gradient channel?

We reuse the existing funding protocol:
- d = 100, T = 0.5, sigma = 0.2, mu = 0.06,
- R_l = 0.04, R_b = 0.06,
- reference value 21.299,
- Beta(1/2,1) random times,
- identical seed schedule and random-tree construction,
- settings (n,M) = (2,10), (3,8), (4,3),
- 100 paired root estimates per setting.

The driver is
\[
f(y,z)=-R_l y-\frac{\mu-R_l}{\sigma}\sum_i z_i
+(R_b-R_l)\left(\frac{1}{\sigma}\sum_i z_i-y\right)_+.
\]

We compare:
1. raw MLP,
2. certified Samplewise IR,
3. certified Batch-IR,
4. **z-suppressed driver** \(f(y,0)\),
5. fully zero driver \(f\equiv0\),
6. constant gradient scaling \(z\mapsto cz\),
7. deliberately too-tight radius factors \(aR(t)\), including \(a<1\).

## Main result

| (n,M) | Raw | Samplewise certified | Batch certified | **f(y,0)** | f=0 |
|---|---:|---:|---:|---:|---:|
| (2,10) | 1.52683 | 0.52147 | 0.45658 | **0.38039** | 0.57440 |
| (3,8)  | 1.10714 | 0.38394 | 0.30012 | **0.19302** | 0.46593 |
| (4,3)  | 3.20564 | 0.77079 | 0.68970 | **0.51383** | 0.60549 |

The z-suppressed driver beats certified Batch-IR at all three depths/budgets. On paired absolute error it wins 71/100, 62/100, and 62/100 roots respectively. The paired mean absolute-error differences (z=0 minus Batch-IR) are -0.0762, -0.1071, and -0.1759; descriptive paired t-test p-values are 1.58e-4, 1.05e-5, and 1.26e-4.

Importantly, setting the **entire** driver to zero is worse than setting only z to zero. Therefore the result is not simply “delete the nonlinear PDE term.” The value-dependent part of the funding driver is useful; the noisy z-channel is what hurts at these finite budgets.

## Too-tight radius sweep

MAE for radius factor \(a\) multiplying the model-derived normalized-delta radius:

### (n,M)=(2,10)

| a | Samplewise | Batch |
|---:|---:|---:|
|0|**0.38039**|**0.38039**|
|0.25|0.38997|0.38447|
|0.50|0.42445|0.40210|
|0.75|0.46710|0.42725|
|1.00 certified|0.52147|0.45658|

### (n,M)=(3,8)

| a | Samplewise | Batch |
|---:|---:|---:|
|0|**0.19302**|**0.19302**|
|0.25|0.21237|0.20166|
|0.50|0.25249|0.22935|
|0.75|0.30825|0.26122|
|1.00 certified|0.38394|0.30012|

### (n,M)=(4,3)

| a | Samplewise | Batch |
|---:|---:|---:|
|0|**0.51383**|**0.51383**|
|0.25|0.57043|0.55733|
|0.50|0.61597|0.59391|
|0.75|0.68524|0.63951|
|1.00 certified|0.77079|0.68970|

At every tested depth, tighter-than-certified contraction improves error monotonically all the way to a=0.

## Constant scaling sweep

The same conclusion appears without any certificate geometry. For \(z\mapsto cz\), the best tested c is zero at every depth.

Selected MAEs:

| (n,M) | c=0 | c=0.25/0.2 | c=0.5 | c=1 |
|---|---:|---:|---:|---:|
| (2,10) | **0.38039** | 0.53150 (c=.25) | 0.82625 | 1.52683 |
| (3,8) | **0.19302** | 0.29561 (c=.2) | 0.56064 | 1.10714 |
| (4,3) | **0.51383** | 0.78379 (c=.2) | 1.46406 | 3.20564 |

## Interpretation

The current funding experiment also fails the geometry-vs-shrinkage life-or-death test.

The certified ellipsoid is valid as a model-derived envelope, but at the tested finite budgets the best behavior is obtained by suppressing the estimated gradient contribution much more aggressively than certification requires. Batch-IR helps because its common factor shrinks z strongly, not because the certified boundary is statistically optimal.

Funding is nevertheless more informative than the Rosenbrock HJB case: fully setting f=0 is worse than f(y,0). This isolates the failure to the noisy gradient channel rather than the entire nonlinear/value feedback. That supports a revised mechanism story centered on **finite-sample gradient-noise rectification / nonlinear reuse bias**.

## Consequence for the paper

The present HJB + funding evidence does not support a headline claim that certified geometry is the source of the finite-budget accuracy gains. Both main positive benchmarks prefer stronger, uncertified suppression.

The method paper can still retain:
- correctness / inheritance results for valid projected generators,
- the exact subspace rescue theorem,
- negative controls and mechanism diagnostics.

But the empirical headline should be reconsidered unless a new active-driver benchmark is found where:
1. the nonlinear gradient contribution is genuinely needed,
2. the certificate is independently rigorous,
3. certified projection beats z=0, best constant scaling, and invalid tighter radii.

Otherwise the stronger paper direction is a diagnosis of finite-budget MLP failure through noisy nonlinear gradient reuse, with certified projection as one principled but not finite-budget-optimal regularizer.
