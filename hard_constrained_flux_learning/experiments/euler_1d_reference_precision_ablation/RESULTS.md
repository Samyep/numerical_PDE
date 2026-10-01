# Seed-0 result: 512-cell versus 2048-cell training references

## Conclusion

The pre-registered screen is **negative**. A fourfold finer Rusanov + SSP-RK2
teacher changes the training targets materially, especially for random-extreme
trajectories, but produces only a small and inconsistent accuracy change in the
trained 64-cell HCFL solver. Do not spend compute on multi-seed confirmation of
this exact resolution-only change.

## Hypothesis and controlled change

**Hypothesis:** part of the remaining canonical Euler error is inherited from
the numerical diffusion and discretization error of the 512-cell reference
teacher, so training on otherwise identical 2048-cell references should improve
moderate OOD and several canonical rollouts.

Only the fine reference grid changed:

- `reference_512`: 512-cell Rusanov + SSP-RK2;
- `reference_2048`: 2048-cell Rusanov + SSP-RK2.

For a paired comparison, each initial condition was generated once at 2048
cells and conservatively averaged to create the 512-cell initial condition.
The resulting 64-cell initial snapshots are exactly equal. Both arms contain
the same 220 ordinary, 260 broad-random, and 100 random-extreme trajectories.

The direct-vector HLLC-HCFL architecture, 6,627 parameters, initialization,
normalization, minibatch sequence, Adam settings, 1,100 updates, hard Tadmor
projection, admissibility limiter, fully-discrete entropy safeguard, and
evaluation protocol were fixed. Both models were evaluated against the same
2048-cell references.

## How much did refinement change the targets?

Normalized RMSE between paired 512-cell and 2048-cell training references:

| training group | all snapshots | final snapshot |
|---|---:|---:|
| ordinary | 0.00432 | 0.00710 |
| broad random | 0.01184 | 0.01859 |
| random extreme | 0.02851 | 0.04466 |
| all 580 trajectories | 0.01449 | 0.02276 |

Initial-snapshot discrepancy is exactly zero for every group. The ablation
therefore changed evolved targets rather than coarse initial-condition
coverage. The largest effect is, as expected, on discontinuous extreme data.

## NRMSE comparison

Negative change means that training on the 2048-cell reference is better.

| evaluation | reference 512 | reference 2048 | relative change |
|---|---:|---:|---:|
| ordinary ID | **0.014086** | 0.014147 | +0.43% |
| broad random (in support) | **0.040894** | 0.040928 | +0.08% |
| moderate OOD, frequencies 4--6 | 0.007423 | **0.007315** | -1.46% |
| Sod | **0.021147** | 0.021894 | +3.53% |
| Lax | **0.088118** | 0.088325 | +0.23% |
| collision | 0.131690 | **0.127219** | -3.39% |
| strong pressure | **0.163986** | 0.164430 | +0.27% |
| near-vacuum expansion | 0.110692 | **0.105701** | -4.51% |
| canonical mean | 0.103127 | **0.101514** | -1.56% |

Only two of five canonical cases improve, the canonical mean gain is below the
pre-specified 5% threshold, and Sod regresses. Ordinary-ID performance is
effectively unchanged, and moderate OOD improves only slightly.

## Constraint and intervention metrics

All rollouts remain admissible and have zero measured Tadmor violations at the
declared `1e-5` reporting tolerance. The largest float64 Tadmor residual is
`1.87e-14` for the 512 arm and `1.73e-14` for the 2048 arm. The maximum measured
total-entropy change is negative for both arms.

The two difficult safety cases are:

| case / arm | min density | min pressure | local limiter rate | FD entropy rate |
|---|---:|---:|---:|---:|
| collision / 512 | 0.04382 | 1.0052e-5 | 0.149% | 3.175% |
| collision / 2048 | 0.03060 | 1.0053e-5 | 0.298% | 4.762% |
| near vacuum / 512 | 0.00981 | 1.0014e-5 | 1.488% | 23.810% |
| near vacuum / 2048 | 0.00890 | 1.0002e-5 | 1.488% | 26.984% |

The higher-precision arm improves error on both cases but invokes the global
entropy safeguard somewhat more often. Hard feasibility remains intact; it
does not explain the small accuracy differences.

## Compute

Reference generation for both paired training sets and the common evaluation
suite took 111.16 seconds on CPU. Total wall time was 126.35 seconds. Once the
data were restricted to 64 cells, training cost was unchanged: 5.59 seconds
for the 512 arm and 5.68 seconds for the 2048 arm.

## Interpretation and decision

This result rules out a simple explanation in which the current 512-cell
teacher resolution is the dominant cause of the remaining 1D Euler accuracy
gap. Refinement makes the target less diffusive and helps collision and
near-vacuum expansion, but the fixed one-step learner does not convert that
extra fidelity into a broad accuracy gain. The slightly higher logged training
loss for the 2048 arm is consistent with the sharper targets being harder to
fit, but a single stochastic loss sample is not sufficient to make a capacity
claim.

**Decision: stop the exact strategy of increasing only the Rusanov reference
grid from 512 to 2048; do not run more seeds.** This does not establish that
the 512-cell reference is publication quality. Converged high-order references
remain necessary for final benchmarking, but they should be treated as
evaluation/data-quality infrastructure rather than an expected standalone
accuracy fix. If reference quality is revisited as a learning intervention, a
higher-order teacher and/or training formulation capable of fitting sharper
targets is more informative than further Rusanov grid refinement alone.

Machine-readable outputs and both checkpoints are in [`results/`](results/).

