# HJB life-or-death ablation: certified geometry vs driver suppression

## Question

Does the headline Rosenbrock HJB result demonstrate a benefit from certified geometry, or is the gain mainly explained by suppressing a noisy nonlinear driver?

The existing headline protocol uses
- PDE: \(u_t+\Delta u-\|\nabla u\|^2=0\),
- \(z=\sqrt2\nabla u\), so \(f(z)=-\|z\|^2/2\),
- dimensions 100, 120, 140, 160,
- \(n=2, M=10\),
- 1000 interior unit-ball + 200 unit-sphere test points,
- 10 repetitions,
- certified radius \(R_0=\sqrt{15}\),
- final value never clipped.

We add the deliberately structure-violating baseline \(f\equiv0\), i.e. delete the nonlinear correction and retain only the root terminal Monte Carlo block. The rerun advances the random-number generator through the same unused child draws as the headline implementation, so the root terminal blocks remain aligned with the archived headline random-tree stream. We also sweep constant drivers \(f\equiv-\kappa\), \(\kappa\in\{0,0.25,\ldots,7.5\}\).

## Main result

| d | Raw | Samplewise \(R_0\) | Batch \(R_0\) | \(f\equiv0\) |
|---:|---:|---:|---:|---:|
|100|2.27995|0.62077|0.26530|**0.004231 ± 0.000060**|
|120|2.66715|0.62851|0.25586|**0.003562 ± 0.000074**|
|140|3.09462|0.62791|0.24556|**0.003167 ± 0.000071**|
|160|3.47939|0.64044|0.24287|**0.002752 ± 0.000065**|

The zero-driver baseline is approximately 63x, 72x, 78x, and 88x lower-error than Batch-IR at d=100,120,140,160 respectively.

The constant-driver sweep selects \(\kappa=0\) in all four dimensions. Thus the best constant correction in the tested range is the zero nonlinear correction.

## Small-radius diagnostic

A separate 100D diagnostic uses 360 fixed points (300 interior + 60 sphere) and 5 paired repetitions. It sweeps invalid radius factors \(a\in\{0,0.25,0.5,0.75,1\}\), with radius \(aR_0\).

| a | Samplewise | Batch |
|---:|---:|---:|
|0|**0.00426**|**0.00426**|
|0.25|0.05234|0.01479|
|0.50|0.19916|0.06450|
|0.75|0.39769|0.14785|
|1.00|0.61554|0.26431|

Error worsens monotonically as the radius is relaxed toward the certified radius. A constant scaling diagnostic \(z\mapsto cz\) shows the same pattern: the minimum is at \(c=0\), and error rises from 0.00426 at c=0 to 2.28338 at c=1.

## Generator-level diagnostic

On a quick 100D sample of 640 intermediate child states,
- mean true driver: -0.0543,
- mean absolute true driver: 0.0543,
- zero-driver generator MSE: **0.0141**,
- raw generator MSE: 3752.3.

For comparison, the existing headline study reports generator MSE 38.42 for Samplewise IR and 10.50 for Batch-IR. Thus the zero-driver alternative is not merely better at the final value; on this diagnostic it is dramatically closer to the true generator than either certified correction.

## Interpretation

This falsifies the current headline interpretation for this HJB protocol. The certified radius is valid, but it is far too loose to be statistically useful at n=2,M=10. The true nonlinear driver is very small on the sampled domain, while Monte Carlo gradient noise makes the raw quadratic driver enormous. Batch-IR helps mainly by shrinking that noise, but complete driver suppression helps much more.

Accordingly, the HJB experiment cannot presently support the claim that certified geometry is the source of the gain. It is better viewed as evidence for nonlinear noise rectification / driver bias in finite-budget MLP.

This does not falsify the general correctness theory for certified projection, nor the exact subspace rescue theorem. It does invalidate the current use of headline HJB as evidence that the certified radius gives a near-optimal finite-budget correction.

## Immediate consequence

Do not revise the manuscript around the current HJB headline. The next decisive step is to test an active-driver benchmark where:
1. the true nonlinear contribution is non-negligible,
2. the certificate is independently valid,
3. certified projection beats zero-driver, constant shrinkage, and invalid tighter radii.

Funding is the first existing candidate, but its certificate should be strengthened before it carries the main method claim.
