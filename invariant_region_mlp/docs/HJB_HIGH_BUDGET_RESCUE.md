# Rosenbrock HJB high-budget rescue

## Outcome

**HJB verdict H-B - TRENDING BUT NOT REACHED.** Along the statistically sensible `n=2` high-sampling path, Raw and Samplewise certified IR errors decrease monotonically, generator MSE falls sharply, and Raw approaches IR. Nevertheless, even at `M=96` the nonlinear methods remain roughly two orders of magnitude above the `f=0` value-error floor. No attainable crossover was observed.

## Validated protocol

The PDE is `u_t + Delta u - ||grad u||^2 = 0`, with `z=sqrt(2) grad(u)` and `f(z)=-0.5||z||^2`. The terminal condition is `log((1+x^T A x)/2)`. The Rosenbrock coefficients exactly reuse the old JAX `PRNGKey(0/1)` construction; at d=100 this gives `trace(A)=304.8856089115143` and `lambda_max=6.383491595137865`. The non-oracle projected radius is `sqrt(15)`. Samplewise correction acts only on z immediately before generator reuse; final roots and values are never clipped. References use the stable scaled Hopf-Cole/Gauss-Laguerre implementation.

The implementation is float64, uses corrected `xi/sqrt(T-t)` terminal EBL normalization and the old uniform random-time protocol, shares complete random trees between Raw and IR, elides the algebraically zero level-0 generator term, and contains no SCaSML heuristic clipping.

## Reproduction gate

The complete old d=100 `(n,M)=(2,10)`, 1200-point, 10-repetition headline was reproduced exactly. Mean relative L2 is `2.2799487622119137` for Raw and `0.6207713181675171` for Samplewise IR, with zero difference from the archived means.

## Stage H1 and stopping-rule decision

H1 used d=100, 300 interior plus 60 boundary points, and 3 paired repetitions.

| n | M | method | reps | value relL2 | gradient relL2 | generator MSE | f evals | samples | sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 16 | raw | 3 | 1.3537 | 4530.4 | 945.43 | 5760 | 190080 | 0.71093 |
| 2 | 16 | samplewise | 3 | 0.54928 | 328.8 | 30.25 | 5760 | 190080 | 0.73731 |
| 2 | 16 | f_zero | 3 | 0.0033014 | 25.678 | - | 0 | 92160 | 0.35907 |
| 2 | 24 | raw | 3 | 0.90432 | 1529.8 | 455.45 | 8640 | 423360 | 1.6504 |
| 2 | 24 | samplewise | 3 | 0.46166 | 211.96 | 23.371 | 8640 | 423360 | 1.7036 |
| 2 | 24 | f_zero | 3 | 0.0029382 | 17.178 | - | 0 | 207360 | 0.81456 |
| 2 | 32 | raw | 3 | 0.66124 | 1116.6 | 138.07 | 11520 | 748800 | 2.9969 |
| 2 | 32 | samplewise | 3 | 0.39419 | 191.43 | 18.943 | 11520 | 748800 | 3.0059 |
| 2 | 32 | f_zero | 3 | 0.0027693 | 12.794 | - | 0 | 368640 | 1.4897 |
| 2 | 48 | raw | 3 | 0.44109 | 1844.7 | 77.637 | 17280 | 1676160 | 6.6295 |
| 2 | 48 | samplewise | 3 | 0.31185 | 245.36 | 14.772 | 17280 | 1676160 | 6.737 |
| 2 | 48 | f_zero | 3 | 0.0026863 | 8.5704 | - | 0 | 829440 | 3.2963 |
| 2 | 64 | raw | 3 | 0.32934 | 504.25 | 16.17 | 23040 | 2972160 | 11.478 |
| 2 | 64 | samplewise | 3 | 0.2584 | 153.36 | 8.4405 | 23040 | 2972160 | 11.549 |
| 2 | 64 | f_zero | 3 | 0.0026259 | 6.4263 | - | 0 | 1474560 | 5.73 |
| 3 | 8 | raw | 3 | 13025 | 1.6739e+08 | 995.53 | 51840 | 809280 | 3.2066 |
| 3 | 8 | samplewise | 3 | 0.83683 | 156.54 | 36.275 | 51840 | 809280 | 3.2161 |
| 3 | 8 | f_zero | 3 | 0.0030294 | 18.144 | - | 0 | 184320 | 0.70981 |
| 3 | 12 | raw | 3 | 13590 | 4.1718e+07 | 437.24 | 112320 | 2648160 | 10.37 |
| 3 | 12 | samplewise | 3 | 0.78245 | 137 | 30.144 | 112320 | 2648160 | 10.438 |
| 3 | 12 | f_zero | 3 | 0.0027119 | 9.9276 | - | 0 | 622080 | 2.3885 |
| 3 | 16 | raw | 3 | 4274.9 | 2.1142e+07 | 331.62 | 195840 | 6180480 | 24.757 |
| 3 | 16 | samplewise | 3 | 0.74097 | 86.549 | 27.16 | 195840 | 6180480 | 24.577 |
| 3 | 16 | f_zero | 3 | 0.0026327 | 6.4202 | - | 0 | 1474560 | 5.7363 |
| 3 | 24 | raw | 3 | 3524.9 | 1.1653e+07 | 184.27 | 432000 | 20537280 | 79.6 |
| 3 | 24 | samplewise | 3 | 0.66101 | 64.693 | 23.869 | 432000 | 20537280 | 79.883 |
| 3 | 24 | f_zero | 3 | 0.0025936 | 3.4941 | - | 0 | 4976640 | 19.122 |
| 4 | 4 | raw | 3 | 1.6326e+15 | 3.1943e+18 | 5093.2 | 123840 | 940320 | 3.7034 |
| 4 | 4 | samplewise | 3 | 0.91085 | 289.7 | 43.044 | 123840 | 940320 | 3.7681 |
| 4 | 4 | f_zero | 3 | 0.0033085 | 25.692 | - | 0 | 92160 | 0.35089 |
| 4 | 6 | raw | 3 | 8.1738e+12 | 1.391e+16 | 2399.6 | 380160 | 4395600 | 17.468 |
| 4 | 6 | samplewise | 3 | 0.88041 | 96.57 | 40.269 | 380160 | 4395600 | 17.469 |
| 4 | 6 | f_zero | 3 | 0.0027568 | 11.412 | - | 0 | 466560 | 1.7912 |
| 4 | 8 | raw | 3 | 4.2536e+13 | 1.3921e+17 | 1834 | 858240 | 13343040 | 48.093 |
| 4 | 8 | samplewise | 3 | 0.86344 | 64.928 | 40.705 | 858240 | 13343040 | 47.952 |
| 4 | 8 | f_zero | 3 | 0.0026267 | 6.4213 | - | 0 | 1474560 | 5.4409 |

Generator metrics use a fixed, deterministic cap of 512 child states per H1 method/repetition and 256 per H2 method/repetition; they diagnose the same early nonlinear-reuse locations under paired trees rather than pretending to enumerate every recursive child.

The continuation criterion was met only by the n=2 high-sampling path. From M=16 to M=64, Raw generator MSE fell from 945.43 to 16.17 (58.47x) and IR generator MSE fell from 30.25 to 8.4405 (3.584x). Value error also fell by more than 20 percent at successive high-budget points. In contrast, n=3 Raw reached errors from hundreds to tens of thousands, n=4 Raw reached `10^11-10^15`, and IR remained finite but worse than the n=2 path. These are deep/under-sampled failures, not evidence against the high-sampling trend.

## Stage H2

H2 therefore retained only n=2 with M in `{48,64,96}`, dimensions `{100,140,160}`, 500 interior plus 100 boundary points, and 5 paired repetitions. A matched M=96 `f=0` reference was also run.

| d | n | M | method | reps | value relL2 | std | generator MSE | samples | sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 100 | 2 | 48 | raw | 5 | 0.4412 | 0.0046002 | 79.591 | 2793600 | 11.459 |
| 100 | 2 | 48 | samplewise | 5 | 0.30762 | 0.0022052 | 18.027 | 2793600 | 11.337 |
| 100 | 2 | 64 | raw | 5 | 0.32985 | 0.0018868 | 38.415 | 4953600 | 20.343 |
| 100 | 2 | 64 | samplewise | 5 | 0.2558 | 0.0011081 | 13.575 | 4953600 | 20.259 |
| 100 | 2 | 96 | raw | 5 | 0.21805 | 0.00061958 | 13.795 | 11116800 | 44.061 |
| 100 | 2 | 96 | samplewise | 5 | 0.18978 | 0.0007688 | 7.8392 | 11116800 | 43.872 |
| 100 | 2 | 96 | f_zero | 5 | 0.0026359 | 8.7354e-06 | - | 5529600 | 22.034 |
| 140 | 2 | 48 | raw | 5 | 0.59355 | 0.0040176 | 142.44 | 2793600 | 15.703 |
| 140 | 2 | 48 | samplewise | 5 | 0.34738 | 0.0007718 | 19.046 | 2793600 | 15.759 |
| 140 | 2 | 64 | raw | 5 | 0.44372 | 0.0052011 | 89.368 | 4953600 | 28.423 |
| 140 | 2 | 64 | samplewise | 5 | 0.29508 | 0.0020018 | 15.41 | 4953600 | 28.402 |
| 140 | 2 | 96 | raw | 5 | 0.29188 | 0.0023266 | 14.807 | 11116800 | 61.832 |
| 140 | 2 | 96 | samplewise | 5 | 0.22549 | 0.0015005 | 7.9109 | 11116800 | 61.598 |
| 140 | 2 | 96 | f_zero | 5 | 0.0017802 | 1.2165e-05 | - | 5529600 | 30.139 |
| 160 | 2 | 48 | raw | 5 | 0.66597 | 0.0030042 | 467.83 | 2793600 | 17.82 |
| 160 | 2 | 48 | samplewise | 5 | 0.36171 | 0.0025925 | 30.18 | 2793600 | 17.971 |
| 160 | 2 | 64 | raw | 5 | 0.49929 | 0.002629 | 229.45 | 4953600 | 31.177 |
| 160 | 2 | 64 | samplewise | 5 | 0.30917 | 0.0016348 | 25.954 | 4953600 | 31.321 |
| 160 | 2 | 96 | raw | 5 | 0.33068 | 0.0011473 | 95.781 | 11116800 | 69.004 |
| 160 | 2 | 96 | samplewise | 5 | 0.24114 | 0.0005999 | 20.731 | 11116800 | 68.842 |
| 160 | 2 | 96 | f_zero | 5 | 0.0015042 | 3.2495e-06 | - | 5529600 | 31.899 |

At d=100, M=96, Raw and IR are approximately 0.218 and 0.190 while `f=0` is 0.00264. At d=160 they are approximately 0.331 and 0.241 while `f=0` is 0.00150. The paired IR improvement is highly stable, but it does not imply that the correct nonlinear correction is practically estimable. H3 was not run: scaling to 1200 points and 10 repetitions would confirm an already stable mean without plausibly closing a 60-160x error gap.

## Interpretation and recommendation

HJB should be presented as the weak-nonlinearity boundary/mechanism case. A valid certified radius materially stabilizes Raw MLP, and the high-sampling trend is real, but the nonlinear correction is so small relative to gradient-estimation noise that the biased `f=0` approximation remains vastly better. It should not remain a positive benchmark for practical certified-geometry accuracy.

No manuscript file was modified.

**HJB VERDICT H-B: TRENDING BUT NOT REACHED.**
