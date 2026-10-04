# Entropy/OOD benchmark results

This is an exploratory benchmark, not an immutably preregistered confirmatory
study.  Four Euler initial conditions had appeared in earlier local work, and
the repository had no protocol commit preceding all generated results.  The
tables are useful evidence and fully reproducible, but they should not be used
alone for a paper-wide superiority claim.

## Bottom line

The experiments support a narrower and more defensible claim than “HCFL is
always the most accurate neural solver.”  HCFL is the only learned method in
this suite that combines full physical completion with the implemented hard
interface-entropy, positivity, conservation, and fully-discrete entropy
checks.  It is also substantially more accurate than its low-order finite
volume endpoint on every evaluated Euler and SWE stress split.

Vanilla FNO attains lower **raw NRMSE and normalized MAE** on the
in-distribution and radius-only OOD radial dam breaks, but it is not counted as
a successful physical solver:
its profiles ring under OOD shifts, it has large conservation residuals, and
38--60% of its saved time intervals violate the sampled entropy balance.  Its
strong-height and long-horizon errors also exceed HCFL by large margins.  The
paper may report FNO's low ID regression error, but must not call that result
physically admissible or an overall win.

The official data-free PINN checkpoints are the most accurate methods on their
own two native, per-instance Euler tasks.  They are not amortized solvers for a
held-out initial condition, and their sampled predictions do not exactly close
the finite-volume conservation or entropy balances.  Those facts make the
PINN result complementary rather than a same-task ranking.

In all tables, NMAE means mean absolute error after the same per-channel
training-standard-deviation normalization used by NRMSE.  Therefore NMAE and
NRMSE differ only in the absolute-versus-square aggregation, not in the tested
states or their scaling.  Euler uses conserved channels `(rho, rho*u, E)` and
2-D SWE uses primitive channels `(h, u, v)`.  Raw per-primitive-channel MAEs
are retained in the linked CSV files.

## 64-cell Euler stress suite

Five named periodic cases are rolled for 63 saved intervals, 4.2 times the
training horizon.  Four cases overlap an earlier local benchmark; the new
evidence is primarily the longer horizon plus the transonic case.  Every
failed trajectory remains in the completion denominator.  Errors and entropy
are reported only on physically completed trajectories.

| Method | Completion | Mean rollout NRMSE | Mean rollout NMAE | Worst completed NRMSE | Max entropy increase (completed only) | Max conservation drift |
|---|---:|---:|---:|---:|---:|---:|
| HLLC-64 | 5/5 | 0.12341 | 0.04153 | 0.20488 | -9.25e-3 | 1.38e-8 |
| MUSCL-HLLC-64 | 5/5 | 0.08424 | **0.02068** | 0.13581 | -2.10e-3 | 1.42e-8 |
| RoeNet-64 same-data adaptation | 3/5 | 0.11183 | 0.03042 | 0.24934 | -2.60e-2 | 4.92e-2 |
| Vanilla residual FNO-64 | 1/5 | 0.05549 | 0.02168 | 0.05549 | +4.52e-3 | 5.27e-2 |
| HCFL-64, 3 seeds | 15/15 | **0.06996 ± 0.00220** | 0.02355 ± 0.00045 | 0.15659 | -6.68e-3 | 1.12e-7 |

Across seeds, HCFL lowers mean NRMSE by 43.3% relative to HLLC-64 and by
16.9% relative to MUSCL-HLLC-64 while completing all cases.  NMAE gives a
more nuanced result: HCFL is 43.3% below HLLC but 13.9% above MUSCL.  Thus HCFL
reduces the larger localized errors that dominate RMSE, while MUSCL retains a
smaller average absolute error over the bulk of cells.  The FNO numbers are
conditional on its only completed case and must not be read as a 5-case mean.
RoeNet is an explicitly labelled architecture adaptation: it uses the official
64-wave learned decomposition but the matched 580/136 HCFL data split,
normalization, and a regularized inverse.  It is not an exact reproduction of
the paper's original 200-cell Sod experiment.  After an independent audit
found an effective checkpoint-selection tolerance bug, it was retrained with
strict completion-first selection; update 2000 is the corrected winner.

The readable final profiles are in
[`frozen64_profiles.png`](results/frozen64/frozen64_profiles.png).  Curves that
became nonphysical are omitted from their final-time panels and explicitly
marked as failed.  The unmasked diagnostic is retained as
[`frozen64_profiles_raw_failures.png`](results/frozen64/frozen64_profiles_raw_failures.png).

## Official data-free PINN checkpoints on native tasks

These numbers load the authors' official L-NN2 checkpoints directly.  An
eighth-order Gauss-Legendre rule converts their continuous prediction to 64
finite-volume cell averages.  This is not a same-task amortized/OOD comparison:
each checkpoint was optimized for its displayed initial condition.

| Case | Method | Rollout NRMSE | Rollout NMAE | Final NRMSE | Max relative conservation residual | Entropy-balance violating steps |
|---|---|---:|---:|---:|---:|---:|
| Sod | HLLC-64 | 0.01847 | 0.00476 | 0.02612 | 4.74e-9 | 0% |
| Sod | HCFL-64 | 0.01820 | 0.00405 | 0.02253 | 2.88e-8 | 0% |
| Sod | Official L-NN2 PINN | **0.00488** | **0.00090** | **0.00513** | 1.57e-3 | 38.1% |
| Lax | HLLC-64 | 0.08872 | 0.02152 | 0.14305 | 3.21e-8 | 0% |
| Lax | HCFL-64 | 0.07048 | 0.01544 | 0.12876 | 3.16e-8 | 0% |
| Lax | Official L-NN2 PINN | **0.01639** | **0.00376** | **0.01621** | 3.37e-3 | 31.7% |

The PINN is decisively more accurate on the tasks it was trained to solve.
HCFL's distinct advantage is an amortized conservative update with hard
forward safeguards, not superior fitting of one fixed Riemann problem.

## 2-D SWE radial dam-break benchmark

The task matches the PDEBench geometry: `[-2.5,2.5]^2`, `g=1`,
constant-extrapolation boundaries, and training radii uniformly distributed on
`[0.3,0.7]`.  The local reference is MC-HLL/SSP-RK2 on `128 x 128`,
conservatively averaged to `32 x 32`.  HCFL, FNO, and clawFNO receive the same
100 training trajectories and the same 20-trajectory validation split.

The FNO and clawFNO classes are imported from the official clawNO commit and
use the published radial-dam architecture (`modes=8`, `width=20`, 24-frame
one-shot output).  Because the data are locally generated, these are
**official-architecture adaptations**, not reproductions of the paper's
downloaded-data error table.

| Evaluation split | HLL-32 | HCFL-s6 | FNO (non-admissible) | clawFNO | Finite-positive count (HLL / HCFL / FNO / claw) |
|---|---:|---:|---:|---:|---:|
| ID radius, standard height | 0.37379 | 0.18158 | 0.02357 | — | 20/20 / 20/20 / 20/20 / 0/20 |
| Radius OOD | 0.35334 | 0.20751 | 0.12478 | — | 20/20 / 20/20 / 20/20 / 0/20 |
| Strong-height OOD | 0.60287 | **0.32998** | 0.62692 | — | 20/20 / 20/20 / 20/20 / 0/20 |
| Two-block / 2x horizon | 0.45789 | **0.22163** | 0.69303 | — | 10/10 / 10/10 / 10/10 / 0/10 |

The corresponding normalized MAE comparison is:

| Evaluation split | HLL-32 NMAE | HCFL-s6 NMAE | FNO NMAE (non-admissible) | clawFNO NMAE |
|---|---:|---:|---:|---:|
| ID radius, standard height | 0.16197 | 0.07146 | 0.00930 | — |
| Radius OOD | 0.15912 | 0.08145 | 0.04702 | — |
| Strong-height OOD | 0.27337 | **0.13553** | 0.23386 | — |
| Two-block / 2x horizon | 0.27767 | **0.11943** | 0.34833 | — |

Both error metrics use identical data and scaling across methods:
primitive-variable NRMSE/NMAE using training-split channel standard
deviations.  HCFL lowers HLL NRMSE by
51.4%, 41.3%, 45.3%, and 51.6% on the four rows.  FNO's raw NRMSE is 7.7 times
lower than HCFL in distribution and 1.66 times lower on radius OOD; these two
numbers measure regression fit only.  HCFL's raw NRMSE is 1.90 times lower on
strong-height OOD and 3.13 times lower at the doubled horizon.

NMAE reaches the same qualitative SWE conclusion: HCFL lowers HLL by 55.9%,
48.8%, 50.4%, and 57.0%.  FNO has lower NMAE on ID and radius OOD, while HCFL
has 1.73 times lower NMAE on strong-height OOD and 2.92 times lower NMAE at the
doubled horizon.  As with NRMSE, FNO's low conditional MAE does not override
its failed conservation/entropy/oscillation qualification.

FNO's numerical defects are not hidden by the error table:

| Split | Max conservation residual | Entropy-violating saved intervals | Mean centerline curvature / reference | Mean height-range overshoot |
|---|---:|---:|---:|---:|
| ID | 0.159 | 38.5% | 1.00 | 0.008 |
| Radius OOD | 1.52 | 51.9% | 1.16 | 0.062 |
| Strong-height OOD | 8.11 | 60.2% | 2.27 | 0.408 |
| 2x horizon | 8.19 | 38.1% | 0.71 | 0.000 |

The curvature ratio quantifies the visible centerline ringing: values above
one contain more discrete curvature than the reference.  The strong-height
case is more than twice the reference curvature and has a large overshoot;
the long-horizon case instead becomes too smooth while following the wrong
profile.  Thus FNO is **finite and positive here, but non-admissible**.  The
finite-positive counts above must not be called physical success rates.

clawFNO's conserved-variable training loss decreased, but no training epoch
produced even one fully physical validation trajectory (0/20).  Its diagnostic
epoch-258 checkpoint also fails every evaluation trajectory.  Encoding the
space-time continuity equation by a divergence-free potential does not imply
positivity or the SWE entropy inequality.  No conditional error or entropy is
reported for its failed cases.

### Exact HCFL forward-safety audit

The following diagnostics use every accepted internal flux update, not a
finite-difference estimate between saved frames.

| Split | Minimum depth | Max relative conservation closure | Max interface Tadmor residual | Max fully-discrete entropy balance | Entropy fallback batch-substeps |
|---|---:|---:|---:|---:|---:|
| ID | 0.5137 | 6.73e-8 | 5.77e-15 | -1.58e-4 | 1 / 960 |
| Radius OOD | 0.5387 | 4.69e-8 | 5.57e-15 | -1.82e-4 | 0 / 960 |
| Strong-height OOD | 0.3257 | 1.75e-7 | 5.85e-15 | -7.53e-4 | 0 / 1165 |
| 2x horizon | 0.5082 | 1.20e-7 | 5.55e-15 | 4.56e-7 | 93 / 960 |

The implemented fully-discrete numerical tolerance is `5e-7`; the worst accepted value
is `4.56e-7`.  The depth limiter never activates.  On the doubled horizon the
entropy blend is materially active, so the long-time guarantee is not achieved
by pretending the learned proposal is always feasible.  The HLL endpoint is a
deployment safety mechanism only: it is not used in the HCFL optimizer loss.

The main figures are
[`swe_radial_final_depth.png`](results/swe_radial/swe_radial_final_depth.png),
[`swe_radial_center_profiles.png`](results/swe_radial/swe_radial_center_profiles.png),
and [`swe_radial_aggregate.png`](results/swe_radial/swe_radial_aggregate.png).
They also expose a limitation: HCFL develops visible Cartesian/square
anisotropy on small-radius and long-horizon waves.  Hard entropy stability is
not rotational invariance.

## Reference and scope limitations

- A fixed radius-0.5 grid-refinement audit gives primitive NRMSE/NMAE
  0.0964/0.0305 for `64 -> 128` and 0.0343/0.0101 for `128 -> 256`, after
  restriction to 32 cells.  Thus the chosen 128-grid reference is materially
  converged relative to 64 and matches PDEBench's native resolution, but is
  not continuum-exact.
- The 2-D reference is a local MC-HLL implementation, not byte-identical
  PyClaw output.  This is why no published clawNO error number is claimed.
- The SWE experiment is flat-bed and strictly wet.  It supports no claim about
  dry beds, bathymetry, or well balancing.
- The Euler suite is periodic.  The native PINN comparison is transmissive.
  Those distinct boundary tasks are not pooled into a single average.
- This new SWE comparison is seed zero.  The Euler HCFL result has three
  independent seeds.  A paper-wide superiority claim would require operator
  and 2-D HCFL replicates as well.

## Reproducibility and provenance

- Exploratory design/implementation record: [`PROTOCOL.md`](PROTOCOL.md) and [`protocol.json`](protocol.json)
- Tuning history: [`TUNING_LEDGER.md`](TUNING_LEDGER.md)
- Euler raw metrics: [`frozen64_metrics.csv`](results/frozen64/frozen64_metrics.csv)
- PINN raw metrics: [`official_lnn2_native_metrics.csv`](results/official_pinn/official_lnn2_native_metrics.csv)
- SWE per-case metrics: [`swe_radial_case_metrics.csv`](results/swe_radial/swe_radial_case_metrics.csv)
- SWE machine-readable summary: [`swe_radial_summary.json`](results/swe_radial/swe_radial_summary.json)
- Reference refinement: [`reference_grid_audit.json`](results/swe_radial/reference_grid_audit.json)

Official source revisions used:

- RoeNet: `ef877957c1c0ddb16eac17006d75bdf7bd786d45`
- data-free PINN: `cebd5f0062903ac971ff7e18063ee546352d7127`
- clawNO: `1c549dbf1d06dc35a8a5df2b897c62aa2f9db186`
- PDEBench task definition: `4ff3e3a4aa1561721b5571fa3a048a0a463e0568`
