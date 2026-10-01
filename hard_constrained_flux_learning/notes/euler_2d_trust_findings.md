# 2D Euler: state-adaptive hard trust region

The first 2D Euler HCFL model is accurate on its training-like and moderate OOD
distributions, but an extreme pressure blast exposed a different failure mode:
the learned correction can be very inaccurate while still satisfying entropy
and remaining physically admissible.

For seed 0, the unguarded HCFL blast error was about 0.84--0.94 while coarse
HLLC was about 0.39. The safety limiters did not trigger because the state was
still admissible. This is an accuracy-under-distribution-shift problem, not a
physics-feasibility problem.

## Hard trust construction

Let both the learned candidate and HLLC base be projected into the same Tadmor
half-space:

[
F_	heta^H,qquad F_{m HLLC}^H.
]

Then use

[
F^{m trust}
=
F_{m HLLC}^H
+
	au(U)
left(
F_	heta^H-F_{m HLLC}^H
ight),
qquad 0le	au(U)le1.
]

Because the entropy-feasible set is convex, the trusted flux is still exactly
Tadmor feasible.

The trust coefficient is calibrated only from training trajectories. Define
the state score

[
z(U)
=
max_{m cells}
left[
rac{1}{d}
sum_j
left(
rac{P_j-mu_j}{sigma_j}
ight)^2
ight]^{1/2},
]

where (P=(ho,u,v,p)). Let (q_{0.99}) be the 99th percentile of this
score on training snapshots. At inference,

[
	au(U)
=
operatorname{clip}
left[
left(
rac{q_{0.99}}{z(U)}
ight)^2,
0.05,
1
ight].
]

No strong-test data are used to calibrate (q_{0.99}).

## Three-seed results

| regime | case | NRMSE | mean tau |
|---|---|---:|---:|
| moderate | ID | 0.00666 ± 0.00085 | 0.998 |
| moderate | OOD | 0.02205 ± 0.00594 | 0.776 |
| strong | blast | **0.37833 ± 0.00195** | 0.050 |
| strong | collision | **0.15758 ± 0.00023** | 0.247 |
| strong | quadrant | **0.07779 ± 0.00199** | 0.327 |

For comparison, coarse HLLC on the same strong tests gave approximately:

- blast: 0.3905;
- collision: 0.1636;
- quadrant: 0.0873.

Thus the state-adaptive trusted HCFL is slightly better than HLLC on all three
strong tests while retaining the learned advantage on moderate OOD. It does
sacrifice some moderate-OOD accuracy relative to the untrusted learned model
(about 0.0221 vs 0.0186).

## Interpretation

This layer addresses a failure mode that entropy/admissibility guarantees do
not cover: a learned flux can be physically legal but inaccurate far outside
its training support.

The state-adaptive trust layer gives the framework a clear hierarchy:

1. learned correction for accuracy inside support;
2. hard entropy projection for structural feasibility;
3. local admissibility limiting when positivity is threatened;
4. training-calibrated trust back toward a classical solver under severe OOD;
5. rare fully-discrete global entropy safeguard.

The trust coefficient is not claimed to be optimal. Its current value is as a
simple, auditable, training-calibrated fallback rule.
