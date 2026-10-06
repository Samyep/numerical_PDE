# Funding high-budget rescue

## Outcome

**Funding verdict F-B - PARTIALLY RESCUED.** Samplewise certified IR improves Raw full-history MLP throughout the completed primary grid, including every 100-root final configuration. The improvement is therefore a real stabilization effect, not a single favorable low-budget point. However, the correct nonlinear-gradient estimator never beats the `z=0` wrong-gradient diagnostic at the representative low, medium, or highest sampled settings. Funding is not rescued as evidence that the correct z-channel is statistically useful at attainable work.

## Protocol and implementation

The benchmark is exactly the validated 100D Funding problem: `T=0.5`, `sigma=0.2`, `mu=0.06`, `R_l=0.04`, `R_b=0.06`, root `x_i=100`, and reference value `21.299`. The payoff and driver are

`g(x)=(max_i x_i-120)_+ - 2(max_i x_i-150)_+`,

`f(y,z)=-R_l y-((mu-R_l)/sigma) sum_i z_i +(R_b-R_l)(sum_i z_i/sigma-y)_+`.

The solver uses float64 geometric-Brownian transitions, Beta(1/2,1) random-time importance sampling, and the corrected terminal EBL weight `xi/sqrt(T-t)`. Samplewise IR projects `Delta_i=z_i/(sigma*x_i)` onto the model-derived envelope `||Delta||_2 <= exp(sigma^2(T-t)/2)` immediately before every generator reuse. The returned root state is never clipped. This report does not strengthen the theoretical status of the Funding envelope beyond the earlier life-or-death study.

## Reproduction gate

The old `(n,M)=(2,10)`, 100-root seed schedule was reproduced before new sweeps. Raw MAE was exactly `1.5268272487666519`; Samplewise IR was `0.5214655268205621`, differing from the archived value by `3.3e-16`. Raw/IR root terminal blocks and random-tree fingerprints are bitwise paired.

## Completed grid and primary results

Stage F1 covered 15 prescribed settings with 30 paired roots. Five higher-sampling probes and six ultra-high probes were then added adaptively. Ten representative configurations were extended to 100 paired roots, including `(3,48)`. Actual work is reported per root; depth and sampling allocation are not treated as interchangeable.

| n | M | roots | Raw MAE | IR MAE | Raw-IR paired gain | IR win frac. | f evals | samples | sec/root |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 10 | 100 | 1.5268 | 0.52147 | 1.0054 | 0.85 | 10 | 210 | 0.0008745 |
| 2 | 12 | 30 | 1.4485 | 0.43972 | 1.0088 | 0.967 | 12 | 300 | 0.001336 |
| 2 | 16 | 30 | 1.0699 | 0.34675 | 0.72315 | 0.9 | 16 | 528 | 0.002203 |
| 2 | 24 | 30 | 0.69918 | 0.26792 | 0.43126 | 0.8 | 24 | 1176 | 0.004988 |
| 2 | 32 | 100 | 0.61892 | 0.22943 | 0.38948 | 0.94 | 32 | 2080 | 0.008666 |
| 2 | 48 | 30 | 0.42318 | 0.13987 | 0.28331 | 0.933 | 48 | 4656 | 0.01456 |
| 2 | 64 | 100 | 0.37438 | 0.16632 | 0.20805 | 0.98 | 64 | 8256 | 0.03022 |
| 2 | 96 | 30 | 0.21212 | 0.10993 | 0.10219 | 0.9 | 96 | 18528 | 0.06356 |
| 3 | 8 | 100 | 1.1071 | 0.38394 | 0.72319 | 0.8 | 144 | 2248 | 0.009459 |
| 3 | 10 | 30 | 0.90188 | 0.2755 | 0.62638 | 0.967 | 220 | 4310 | 0.02012 |
| 3 | 12 | 30 | 0.69312 | 0.23652 | 0.4566 | 0.833 | 312 | 7356 | 0.03642 |
| 3 | 16 | 100 | 0.49832 | 0.18729 | 0.31103 | 0.83 | 544 | 17168 | 0.07319 |
| 3 | 20 | 30 | 0.4819 | 0.19751 | 0.28439 | 0.867 | 840 | 33220 | 0.1054 |
| 3 | 24 | 100 | 0.29469 | 0.13018 | 0.16451 | 0.88 | 1200 | 57048 | 0.2238 |
| 3 | 32 | 30 | 0.24848 | 0.15216 | 0.096324 | 0.8 | 2112 | 134176 | 0.4764 |
| 3 | 40 | 30 | 0.2182 | 0.14856 | 0.069636 | 0.7 | 3280 | 260840 | 0.9435 |
| 3 | 48 | 100 | 0.15391 | 0.12461 | 0.029299 | 0.62 | 4704 | 449328 | 1.873 |
| 4 | 3 | 100 | 3.2056 | 0.77079 | 2.4348 | 0.88 | 159 | 894 | 0.003755 |
| 4 | 4 | 30 | 1.9901 | 0.45875 | 1.5314 | 0.867 | 344 | 2612 | 0.01129 |
| 4 | 5 | 30 | 1.1218 | 0.38247 | 0.73937 | 0.833 | 635 | 6080 | 0.02493 |
| 4 | 6 | 30 | 1.2427 | 0.43812 | 0.80458 | 0.933 | 1056 | 12210 | 0.05195 |
| 4 | 8 | 100 | 0.56525 | 0.25923 | 0.30602 | 0.86 | 2384 | 37064 | 0.1405 |
| 4 | 10 | 30 | 0.28627 | 0.16038 | 0.12588 | 0.633 | 4520 | 88310 | 0.3209 |
| 4 | 12 | 30 | 0.21874 | 0.12437 | 0.094373 | 0.7 | 7656 | 180156 | 0.6211 |
| 5 | 2 | 30 | 6.7082 | 1.1246 | 5.5836 | 0.833 | 252 | 918 | 0.003587 |
| 5 | 3 | 100 | 2.3153 | 0.55016 | 1.7651 | 0.89 | 1032 | 5781 | 0.0216 |

## Suppression and shrinkage diagnostics

The constant shrinkage coefficient was selected on a separate 30-root validation seed across `(2,32)`, `(3,16)`, and `(3,24)`. Among `c in {0.25,0.5,0.75}`, `c=0.25` was best at all three validation settings and was frozen before the diagnostic test roots were examined. The deliberately invalid tighter envelope used factor `0.75`.

| n | M | method | roots | MAE | bias | win frac. vs Raw |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 10 | raw | 100 | 1.5268 | 1.4821 | - |
| 2 | 10 | samplewise | 100 | 0.52147 | 0.30359 | 0.85 |
| 2 | 10 | z_zero | 100 | 0.38039 | 0.042496 | 0.85 |
| 2 | 10 | f_zero | 100 | 0.5744 | 0.48513 | 0.87 |
| 2 | 10 | shrink_c0.25 | 100 | 0.5315 | 0.32229 | 0.85 |
| 2 | 10 | tight_a0.75 | 100 | 0.4671 | 0.2142 | 0.85 |
| 3 | 16 | raw | 100 | 0.49832 | 0.064341 | - |
| 3 | 16 | samplewise | 100 | 0.18729 | 0.035306 | 0.83 |
| 3 | 16 | z_zero | 100 | 0.073679 | 0.0073778 | 0.86 |
| 3 | 16 | f_zero | 100 | 0.43452 | 0.43452 | 0.52 |
| 3 | 16 | shrink_c0.25 | 100 | 0.15373 | -0.063657 | 0.88 |
| 3 | 16 | tight_a0.75 | 100 | 0.15561 | 0.0058702 | 0.83 |
| 3 | 24 | raw | 100 | 0.29469 | -0.088193 | - |
| 3 | 24 | samplewise | 100 | 0.13018 | -0.061961 | 0.88 |
| 3 | 24 | z_zero | 100 | 0.042338 | 0.0068269 | 0.98 |
| 3 | 24 | f_zero | 100 | 0.43847 | 0.43847 | 0.23 |
| 3 | 24 | shrink_c0.25 | 100 | 0.12358 | -0.10671 | 0.86 |
| 3 | 24 | tight_a0.75 | 100 | 0.11439 | -0.067052 | 0.88 |
| 3 | 48 | raw | 100 | 0.15391 | -0.10771 | - |
| 3 | 48 | samplewise | 100 | 0.12461 | -0.11199 | 0.62 |
| 3 | 48 | z_zero | 100 | 0.014596 | 0.0044813 | 0.98 |

At `(3,48)`, Raw and IR have entered the expected high-sampling convergence regime: their 100-root MAEs are 0.15391 and 0.12461, respectively. Yet `z=0` is 0.014596. Thus the biased approximation has reached a much lower practical error floor before the correct z-channel becomes worthwhile. `f=0` remains substantially worse where tested, confirming that the value-dependent driver channel is active.

## Interpretation and recommendation

Funding is a strong demonstration that certified projection stabilizes Raw recursive nonlinear reuse: IR wins the paired Raw comparison across the full work range and prevents deep/under-sampled degradation. It is not a clean positive active-gradient benchmark because `z=0` remains better even when Raw approaches IR. It should remain as a secondary mechanism/partial-rescue result, not be promoted beside active VB as evidence that certified geometry makes the correct gradient information useful.

No manuscript file was modified.

**FUNDING VERDICT F-B: PARTIALLY RESCUED.**
