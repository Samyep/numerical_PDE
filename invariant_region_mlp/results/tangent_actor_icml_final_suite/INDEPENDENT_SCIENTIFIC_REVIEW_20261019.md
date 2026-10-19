# Independent scientific review — 2026-10-19

**Reviewed revision:** `b1d893a189c4d1df93e07399d9de944e2459a589`.  
**Disposition:** **100D stage conditionally complete; publication-ready 7/7 claim NOT approved; 175-cell dimension sweep NOT approved.**  
**Review gate:** keep `REVIEW_GATE_STATUS.json` at `pending_independent_review`; do not launch any job or change the frozen protocol in response to this review.

This document is a post-freeze independent cross-check. It does not replace or silently rewrite the frozen results, the artifact manifest, or the original paper-fidelity audit. The initial 16 canonical artifact checks can all pass while an important comparability error remains.

## P0 — an unpaired test bank invalidates the strict 7/7 assertion

I inspected the 29 rows in `MAIN_RESULTS.csv` at fixed seed 0 and compared `test_bank_sha256`, `states`, `method`, `dimension`, and `value_rmse` **within** each PDE.

| 100D PDE | Whole-function methods | Same bank SHA-256? | Same state count? | Interpretation |
|---|---:|---|---|---|
| HJB-Log | 4 | yes | yes, 960 | paired-bank comparison |
| OU/GMM | 4 | yes | yes, 960 | paired-bank comparison |
| Joint concave | 4 | yes | yes, 960 | paired-bank comparison |
| Joint convex | 4 | yes | yes, 960 | paired-bank comparison; Tangent-Actor tagged adapted |
| MLP failure | 4 | yes | yes, 960 | paired-bank comparison; all methods adapted |
| Allen–Cahn | 4 | yes | yes, 410 | paired-bank comparison; reference-refinement tier differs between source rows |
| **Nonconvex coupled** | **4** | **NO** | **NO** | **not a paired-bank comparison** |

For nonconvex coupled:
- Tangent-Actor: **960 states**, hash `85a9430c6f60762f5964a9de46122e8010cd1d96b67b286a6997a824d20c2180`, RMSE `0.0008041022870333`; historical `tangent_actor_accuracy_frontier/euler_convergence.csv` at 64 Euler steps and 2048 paths.
- DPI, PINN-HTE, SDGD-PINN: **512 states**, hash `3aef1dc9855657a66f1744f0cdfe3d12e80b009f52b158154a1da799e42044de`, corresponding RMSEs `0.0016991249965970503`, `0.003068255018974441`, and `0.00887679594351191`.
- The baseline bank points to `post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz`, whereas the frozen Tangent-Actor 960-state result has a different test-bank hash. A shared PDE and value metric do **not** guarantee a shared test population.

**Correct current wording:** Tangent-Actor records the smallest RMSE in **6/6 comparable 100D whole-function PDE banks**. A seventh nonconvex result is favorable numerically but is **provisional/unpaired** and must not be counted as a strict head-to-head win. The original 7/7 numerical observation is not fabricated, but the published-table interpretation requires correction.

**Required remediation:** Evaluate the frozen Tangent-Actor checkpoint on the 512-state baseline bank, or all three frozen baseline checkpoints on the 960-state Tangent-Actor bank (preferably both for robustness). Do not retrain or select checkpoints using test truth. Update `MAIN_RESULTS.csv`, `100D_RESULTS.csv`, `ACCURACY_COMPUTE_FRONTIER.csv`, Figure 2, the LaTeX table, report, and the two manifests only after the same-bank evaluation and full artifact regeneration; verify value references are computed identically. Rerun the complete tests and executable artifact audit. The required source checkpoints/raw workspace files have not been independently accessed through this review.

The current artifact auditor only checks that `test_bank_sha256` is *non-null*, not equal across methods. Run the separately added `verify_paired_banks.py` guard, which intentionally exits nonzero on the current data.

## Fidelity and training convergence

- DPI: **PASS** on 100D OU/GMM native paper metric, based on approximate Figure-4 readings rather than tabulated paper numbers.
- Deep BSDE: **PASS**, with correct scope restricted to the single-initial-state task. The published HJB and Allen–Cahn reference relative errors are 0.17% and 0.30%, while the seed-0 reproductions are approximately 0.2245% and 0.7057%. The historical three-seed Allen reproduction average is approximately 0.3204%; this historical evidence is not cherry-picked for the ranking.
- PINN-HTE: **CONDITIONAL**; native OU/GMM seed-0 value and gradient errors exceed the approximate paper bars. The terminal end of the plotted native convergence curve visibly deteriorates after approximately 18,000 updates. Completion of 20,000 updates does not establish stable convergence.
- SDGD-PINN: **CONDITIONAL**; native 250D HJB-Log relative L2 `0.00907227` versus the reported five-seed `0.006147` (approximately **47.6%** larger), with 204.39 versus 210 minutes. All 10,000 source-budget epochs were run with exact resume. Do not choose the minimum-error intermediate epoch post hoc.

The original algorithms and source budgets are documented, but a method-fidelity `PASS` on its native setting does not transfer the same numerical quality guarantee to custom PDEs.

## Native-point and compute caveats

The whole-function primary ranking excludes Deep BSDE by design. The native single-point table tells a less uniformly favorable story and **must be disclosed**:
- On joint concave, Deep BSDE native absolute error `3.812e-4` is smaller than Tangent-Actor's frozen 64-step native error `5.596e-4`, with substantially less compute.
- On joint convex, PINN-HTE and SDGD-PINN native absolute errors (about `2.22e-4` and `2.27e-4`) are smaller than Tangent-Actor's frozen 64-step native absolute error (about `3.28e-4`).
- On MLP failure, the DPI native point error (about `0.00384`) is smaller than Tangent-Actor's (about `0.0104`).
- On OU/GMM, SDGD-PINN native point error (about `0.191`) is smaller than Tangent-Actor's (about `0.254`), though SDGD consumes more training compute on that transfer.

The full-cost 100D points indicate **accuracy–cost tradeoffs**, not uniform compute superiority. For example, HJB-Log Tangent-Actor total time is about 1093 s versus about 304 s for PINN-HTE; joint concave about 2044 s versus about 1244 s for DPI. Tangent-Actor's 100D frozen-model MC inference also takes tens to hundreds of seconds on these full held-out banks, versus milliseconds for trained network baselines. Do not make a blanket speedup or deployment-latency claim.

## Effective dimension and numerical precision

Three joint/nonconvex PDEs have a known **effective dimension of two** and the MLP-failure PDE has effective dimension one. HJB-Log and Allen–Cahn have radial structure; OU/GMM has analytically structured Gaussian-mixture data. Full ambient 100D inputs remain real, but this suite alone does not establish advantages on nonseparable, intrinsically high-dimensional solutions. Disclose these structural restrictions and distinguish ambient from effective dimension.

The 128-step joint concave MCV RMSE (`2.1347e-4`) is near its MC contribution (`1.72335e-4`) and its independent reference uncertainty proxy (`8.45e-5`). The joint convex MCV RMSE (`1.43167e-4`) is near the MC contribution (`1.45160e-4`) and reference proxy (`7.78e-5`). Beyond the first Actor update, policy-value improvements are not reliably resolved by these error floors. A 117x variance reduction is **not** a 117x RMSE reduction. Allen–Cahn carries different reference-uncertainty metadata across methods (`3.86e-5` vs `9.64e-6`): audit the numerical truth arrays, not only the held-out state hash, before asserting identical numerical reference fidelity.

## Figures and tables (visual inspection)

- **Figure 1** accurately shows the frozen low-budget MLP dimension-amplification experiment, but pending-cell crosses are drawn close to zero; avoid visual implication that unrun dimensions have low error. The only measured Tangent-Actor ambient dimension there is 100.
- **Figure 2** mixes PDEs with widely different RMSE scales and depicts fixed-source-budget points alongside multiple Euler-step refinements of one method. It is primarily an **accuracy–cost scatter**, not a universally measured Pareto frontier. The nonconvex 960/512-bank inconsistency currently affects its interpretation.
- **Figure 3** does **not** yet measure dimension scaling: only `d=100` has methods completed. The script places unexecuted crosses at `floor = 0.45 * min(finite RMSE)`, which is an arbitrary y position close to favorable accuracy; these must be moved outside data coordinates or the figure omitted from a submission until more dimensions finish.
- **Figure 4** supports the decisive first Actor update and variance reduction; crowded dual-axis legends should be cleaned for camera-ready size.
- **Figure 5** shows post-freeze high-accuracy refinements, which are not separate source-budget baseline wins. Reference and finite-MC floors should remain visible.
- **LaTeX main table** currently lists six dimensions despite only one completed dimension; use a 100D-only table until the gated sweep exists. **Native-point table** has a “mean ± sd” heading even though each primary cell is a single seed; rename to “absolute error (seed 0)” or add actual replication.

## Accounting, integrity, and review gate

- GitHub artifact metadata reports 35/210 completed method cells (27 newly completed, 8 audited reuse), 175 pending; 29 main whole-function rows include one extra 250D native SDGD anchor, 44 native point rows, and 41 frontier rows.
- The tracked `FINAL_ARTIFACT_AUDIT.json` reports **16/16** implemented checks PASS, but the missing paired-bank check above is a new scientific defect. Passing the existing tests/audit is not sufficient to clear this review.
- The suite's incremental Slurm ledger is `35.2136` GH200 allocated GPU hours; historical source experiments `21.9911` additional hours are listed separately. The user's approximately 5 h 38 min / 45-credit SDGD accounting covers a different subset and pricing unit; do not add credits directly to allocated physical GPU-hours.
- A future 175-cell single-seed sweep is estimated at **139.53 GPU hours**, range **81.67–284.82**; such a model-based extrapolation is not an empirical completed experiment.
- Raw HPC checkpoints, live Slurm queues, and the full 20-test command were **not independently rerun or inspected** in this GitHub-level review. Those points remain claims of the run's recorded manifests, not observations independently revalidated here.

## Closure / release decision

**Closed as archival evidence:** corrected seed-0 100D training runs, four native fidelity anchors, MLP low-budget mechanism curve, Actor iteration and MCV measurements, incremental compute ledger, prior 20-test/16-audit records.

**Open, blocking a submission-wide head-to-head claim:** P0 nonconvex same-bank evaluation and corresponding regeneration; additional reference-array parity check for Allen–Cahn. **Open, blocking full-grid and general-dimensional claims:** 175 gated cells; no dimension scaling curve yet; single-seed uncertainty; structured/effective-low-dimensional PDE coverage; native-fidelity conditional ratings.

The `REVIEW_GATE_STATUS.json` status **must remain** `pending_independent_review`. This review does **not** grant permission to run the expensive dimension sweep. Prioritize cheap frozen-checkpoint same-bank evaluation and accurate paper framing before considering new training.

## Additional publication-table completeness issue

The generated `tables/native_point_comparison.tex` contains nonconvex-coupled single-point results for Deep BSDE, DPI, PINN-HTE, and SDGD-PINN, but **no Tangent-Actor nonconvex-coupled native-point row**. This absence is not evidence of failure or a win; a table caption or an explicit missing-data entry is necessary before making claims of complete native-point coverage. The `joint_convex` Tangent-Actor primary row is labeled `adapted` (sign-transform path) even though the whole-function discussion informally calls it a final-method win; retain and explain that classification consistently.

The primary native-point table's single-seed `Abs. error (mean ± sd)` header should also be replaced as noted above. No native-point checkpoint or numerical estimate was manufactured to fill the missing row.
