# Audited 64-cell Euler and SWE results

## Bottom line

The common four-cell HCFL method converged for all three independent seeds on
both systems.  In the final 64-cell test, it reduced rollout error relative to
the native first-order finite-volume method in all four aggregate groups while
completing every trajectory and passing the conservation, positivity,
interface-entropy, and fully-discrete entropy audits.  It did **not**
universally beat the stronger TVD MUSCL baseline: HCFL won the canonical Euler
group, but MUSCL was more accurate in the other three aggregate groups.

The defensible accuracy claim is therefore:

> Under the stated 64-cell protocol, HCFL adds substantial accuracy to a
> first-order shock-capturing baseline while retaining explicit physical and
> entropy safeguards.  It is not a universal replacement for every
> second-order method.

## Protocol

- Candidate resolution: 64 finite-volume cell averages.
- Final reference: independent 2048-cell HLLC + SSP-RK2 trajectories for
  Euler and HLL + SSP-RK2 trajectories for SWE.
- Reference initialization: piecewise-constant prolongation of the exact same
  64 initial averages, followed by conservative restriction to 64 cells.  The
  reference therefore receives no hidden subcell initial information.
- Random tests: 24 trajectories in each of three disjoint categories per
  system (72 trajectories per system).  Canonical tests: five Euler and four
  SWE problems.
- HCFL result: mean over seeds 0, 1, and 2; the reported uncertainty is the
  sample standard deviation over those seeds.
- Failed baseline trajectories remain in `all-finite` error and completion
  statistics; they are not silently removed.

The Euler training split contains 580 trajectories: 220 ordinary, 260 broad,
and 100 extreme.  Its disjoint validation split contains 136 trajectories: 44
ordinary, 52 broad, 20 extreme, and 20 randomized wave-coverage cases.  The
SWE split contains 580 training trajectories (220 ordinary, 160 broad, and 200
Froude-coverage) and 136 validation trajectories (44, 32, and 60).  Training
targets use a 512-cell strict Rusanov generator for Euler and a 1024-cell HLL
generator for SWE, conservatively restricted to 64 averages.  These training
references are distinct from the 2048-cell final scoring references.

## Convergence across independent seeds

| System | Seed 0 | Seed 1 | Seed 2 | Mean ± sample std | All converged |
|---|---:|---:|---:|---:|:---:|
| Euler | 0.03057 | 0.02415 | 0.02500 | 0.02657 ± 0.00349 | yes |
| SWE | 0.04355 | 0.03892 | 0.04244 | 0.04164 ± 0.00242 | yes |

Each entry is validation rollout NRMSE.  Runs stopped on a validation plateau
at the minimum learning rate before the 50,000-update cap; nonconverged
checkpoints are not included.

## Final accuracy at 64 cells

The table uses rollout NRMSE over all finite outputs.  HCFL entries are
mean ± sample standard deviation over three seeds.

| System / test group | Native FVM64 | TVD MUSCL64 | Learned flux64 | Residual FNO64 | HCFL64 |
|---|---:|---:|---:|---:|---:|
| Euler random (72) | 0.03467 | **0.02136** | 0.02967 | 0.04439 | 0.02559 ± 0.00120 |
| Euler canonical (5) | 0.14105 | 0.09432 | 0.12908 | 0.14766 | **0.08258 ± 0.00389** |
| SWE random (72) | 0.06967 | **0.04057** | 0.05685 | 0.05386 | 0.04207 ± 0.00189 |
| SWE canonical (4) | 0.06811 | **0.03158** | 0.05625 | 0.04761 | 0.03764 ± 0.00109 |

Relative to native first-order FVM64, HCFL lowers error by 26.2% (Euler
random), 41.5% (Euler canonical), 39.6% (SWE random), and 44.7% (SWE
canonical).  Relative to MUSCL64, HCFL is 12.4% better on Euler canonical,
but 19.8%, 3.7%, and 19.2% worse on Euler random, SWE random, and SWE
canonical, respectively.

HCFL is also more accurate than both same-data neural baselines in all four
groups.  This is evidence for the value of its numerical structure under this
matched protocol, not a claim that these compact controls reproduce every
published PINN, neural-operator, or neural-FV implementation.

## Completion and physical audit

| Method | Euler random | Euler canonical | SWE random | SWE canonical | Structural observation |
|---|---:|---:|---:|---:|---|
| HCFL64, every seed | 72/72 | 5/5 | 72/72 | 4/4 | Conservative FV update, hard Tadmor projection, positive accepted states, nonincreasing accepted-step entropy |
| Native FVM64 | 72/72 | 5/5 | 72/72 | 4/4 | Robust low-order numerical baseline |
| TVD MUSCL64 | 72/72 | 5/5 | 72/72 | 4/4 | Stronger accuracy baseline; no hard entropy audit in this implementation |
| Learned flux64 | 72/72 | 4/5 | 72/72 | 4/4 | FV conservation retained, but Euler near-vacuum pressure became negative |
| Residual FNO64 | 69/72 | 2/5 | 72/72 | 4/4 | No hard conservation, positivity, or entropy mechanism |

Across every HCFL aggregate row and seed:

- maximum relative periodic conservation drift was `1.02e-7`;
- maximum interface Tadmor residual was `1.72e-7`, below the audit tolerance
  `1e-5`;
- every maximum accepted-step total entropy change was nonpositive (the value
  closest to zero was `-5.40e-7`);
- the minimum Euler density, Euler pressure, and SWE depth remained positive:
  `1.44e-5`, `8.04e-4`, and `5.52e-2`, respectively; and
- all 24 HCFL aggregate entries (two systems × two test groups × two safety
  variants × three seeds) passed the fail-closed audit.

By contrast, the Euler FNO completed only 69/72 random trajectories and 2/5
canonical cases, reached negative density/pressure, and had up to 3.42%
relative conservation drift on the canonical set.  The SWE FNO completed its
rollouts but had up to 5.05% conservation drift and positive entropy changes.
The unconstrained learned-flux control retained FV conservation but failed the
Euler near-vacuum case with negative pressure.  These failures are counted in
the table rather than filtered out.

## Removing the low-order anchor

The requested `no F_low` variant uses the same trained network and hard Tadmor
proposal, but it never constructs or blends a low-order flux.  It accepts a
proposed update only if density/depth (and Euler pressure) stay positive and
the fully-discrete total entropy does not increase; otherwise it halves the
time step and retries.

| System / test group | Full safety stack | No `F_low` | Observation |
|---|---:|---:|---|
| Euler random | 0.02559 ± 0.00120 | 0.02541 ± 0.00133 | 1–2 total retry halvings per seed across the three random suites |
| Euler canonical | 0.08258 ± 0.00389 | 0.08258 ± 0.00389 | Numerically identical; no retries |
| SWE random | 0.04207 ± 0.00189 | 0.04207 ± 0.00189 | Numerically identical; no retries |
| SWE canonical | 0.03764 ± 0.00109 | 0.03764 ± 0.00109 | Numerically identical; no retries |

Thus `F_low` had almost no practical effect on these tests, and the simpler
accept/reject design is a strong empirical candidate.  The important caveat
is theoretical: every **accepted** no-`F_low` update still passes positivity
and fully-discrete entropy checks, but existence/termination of an accepted
time step is not supplied by the same constructive admissible low-order
endpoint used by the full method.  The full stack remains the safer version
for a theorem-level guarantee; the no-`F_low` result should be presented as an
empirical simplification unless that termination argument is proved.

## Reproducibility artifacts

- `benchmark64_summary.json`: complete machine-readable run summary.
- `benchmark64_metrics.csv`: per-case/per-method metrics.
- `benchmark64_aggregate.csv`: aggregate metrics used above.
- `benchmark64_audit.json`: fail-closed audit result.
- `benchmark64_accuracy.png`: accuracy and completion summary.
- `benchmark64_profiles.png`: representative solution profiles.
- `benchmark64_low_order_ablation.png`: full versus no-`F_low` comparison.
- `hcfl64_replicate_convergence.png`: independent-seed validation curves.
- `LITERATURE_COMPARISON.md`: fair-comparison boundary and primary references.
