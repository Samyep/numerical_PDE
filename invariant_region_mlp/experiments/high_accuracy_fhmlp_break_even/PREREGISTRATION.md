# Task: High-Accuracy Amortization Break-Even — Tangent-Actor vs Full-History MLP

Repository: `Samyep/numerical_PDE`

Create a new branch:

`ir-mlp-high-accuracy-fhmlp-break-even`

## 1. Research objective

We previously completed an audited amortization experiment on a 100D Joint Nonseparable semilinear PDE.

Under the relatively loose accuracy tolerance RMSE ≤ 0.06, Tangent-Actor and FH-MLP crossed at approximately 524,288–1,048,576 queries.

We now want to answer:

**When both methods must achieve substantially higher accuracy, does the amortized computational advantage of Tangent-Actor emerge at fewer queries?**

Our primary deliverable is a publication-quality plot of measured cumulative computation versus query count, ideally showing a clearly observed crossover.

However, no crossing may be manufactured through selective accuracy thresholds, timing extrapolation, altered baselines, or omitted measurements.

Do not run Deep BSDE, DPI, PINN-HTE, or SDGD-PINN.

Do not modify any previously audited experiment results.

## 2. Source experiments

Study and reuse the validated implementations from:

- `invariant_region_mlp/experiments/single_point_break_even_joint/`
- `invariant_region_mlp/experiments/single_point_break_even/`
- `invariant_region_mlp/experiments/single_point_32_ice/`
- `invariant_region_mlp/results/single_point_break_even_joint/`
- `invariant_region_mlp/docs/SINGLE_POINT_BREAK_EVEN_JOINT_REPORT.md`

The previous result at RMSE ≤ 0.06 must remain unchanged.

In particular, preserve the exact validated FH-MLP recursion and the previously audited Tangent-Actor implementation.

Reuse available frozen TA training artifacts when their source revisions, configurations and hashes match. If checkpoints must be regenerated, use the same frozen algorithm and reference-free selection rules, not configurations chosen using the new test errors.

## 3. PDE benchmark

Use the existing 100-dimensional Joint Nonseparable PDE, with horizon \(T=0.5\):

\[
u_t+\Delta u+f(u,z)=0,
\qquad z=\sqrt2\nabla u,
\]

where

\[
f(y,z)=y-y^3-\frac12\|z+yv\|^2.
\]

Use the existing terminal condition, frozen basis vectors and independently validated two-dimensional periodic reference solver.

Do not substitute Full-Rank Log-Sum-Exp HJB, Dense Fourier, or a simplified PDE.

All methods must solve exactly the same equation, initial-state distribution, and terminal condition.

## 4. High-accuracy requirements

Freeze the following targets before conducting new accuracy calibration:

| Target | RMSE tolerance | Role |
|---|---:|---|
| Primary | 0.010 | Main crossover experiment |
| Secondary | 0.005 | Stricter accuracy challenge |

A method qualifies only if its RMSE is below the prescribed threshold on the calibration bank and subsequently passes on the independent final bank.

Accuracy matching means satisfying the same error tolerance, not having numerically identical RMSE.

Do not relax either threshold after observing results.

Do not change thresholds to obtain an attractive crossover.

If FH-MLP cannot achieve RMSE ≤ 0.005 within the tested budget, report `TARGET_NOT_REACHED_WITHIN_EVALUATED_BUDGET`. The primary 0.010 experiment should still proceed if qualified.

## 5. Accuracy calibration

Create fresh query banks independent of the previously inspected evaluation banks:

- Calibration: 64 query states
- Final accuracy: 128 query states
- Timing: a reproducible truth-free stream allowing up to \(2^{20}\) queries

Use the original \(t=0,\ x\sim N(0,I_{100})\) query distribution.

Calibration and final reference data must be generated using the validated independent solver, with numerical reference uncertainty checked to be negligible relative to the target thresholds.

Freeze seeds, bank hashes and all configuration candidates before opening new final errors.

### Tangent-Actor

Use three fixed training seeds and the mature terminal-aware actor.

Evaluate two inference modes:

1. `TA-direct`: prediction from the learned value head.
2. `TA-FK-MCV`: the original Feynman–Kac estimator with its martingale control variate.

Keep the inference modes separate.

For the main figure, use `TA-direct` only if it passes the chosen accuracy threshold on both calibration and final banks.

For FK–MCV, evaluate prespecified path budgets of 32, 128 and 512. Select the cheapest accuracy-qualified path budget using calibration only. Report its accuracy and timing separately, even if it is slower.

Do not call a direct value-head measurement an FK–MCV measurement.

The final report must flag that the current manuscript describes FK–MCV as its final evaluation method, whereas a direct-value-head timing figure uses a distinct inference option.

### Full-History MLP

Preserve the audited unprojected, float64, full-history recursive implementation, including its original gradient estimator.

Extend the previous candidate grid to include:

- \(n=2\): \(M\in\{32,48,64,96,128,192,256,384,512\}\)
- \(n=3\): \(M\in\{24,32,48,64,96\}\)

These are predeclared candidate sets, not a mandate to run configurations that exceed the frozen safety or GPU budget.

Use five independent estimator replicates for calibration and five fresh replicates for final validation.

Before executing expensive configurations, conduct truth-free profiling to estimate runtime, generator calls, and peak memory. Skip infeasible configurations only under predeclared resource rules, and preserve their `NOT_EVALUATED_COMPUTE_CAP` or `OOM` statuses.

No gradient clipping, oracle projection, modified nonlinear driver, or changed stochastic recursion is allowed.

For each accuracy target:

1. Determine which configurations satisfy the calibration RMSE threshold.
2. Among qualified configurations, select the one with the lowest measured representative per-query cost.
3. Freeze the selected \((n,M)\) and its configuration hash.
4. Evaluate the frozen configuration on the independent final bank.
5. If final accuracy fails, report the failure; do not select a replacement using final errors.

Record the complete configuration-level accuracy table, including nonmonotonic outcomes.

## 6. Actual H100 computation comparison

Use NVIDIA H100 80GB GPUs only.

All formal numerical measurements must run on actual H100 hardware. Do not mix GH200, RTX GPUs, or CPU timing into the figure.

The total cost for Tangent-Actor is

\[
T_{\mathrm{TA}}(N)
=
T_{\mathrm{train}}
+
T_{\mathrm{selection}}
+
T_{\mathrm{load}}
+
T_{\mathrm{inference}}(N).
\]

The total cost for FH-MLP is

\[
T_{\mathrm{FH}}(N)
=
T_{\mathrm{recursive\ evaluation}}(N).
\]

TA training cost must be charged even when an existing checkpoint is reused.

Distinguish the compute cost represented by the model from the additional GPU allocation required to conduct this experiment.

Use the same query stream for both methods.

**Every reported cumulative-time point must come from executing that number of queries. Do not multiply a single-query measurement to create the formal timing curve.**

Conduct three independent timing replicates with synchronized CUDA timing boundaries.

Use identical query prefixes across methods and record actual elapsed active computation.

Do not include scheduler queue delays or offline reference construction in the method-active curves, but retain them in resource accounting.

### Timing grid

First run a small truth-free performance pilot to estimate the approximate crossing region. This pilot may guide the timing grid but must not change the accuracy thresholds or selected solver configurations.

Use an initial geometric query grid:

\[
N\in\{2^k:k=0,\ldots,20\}.
\]

Add several prespecified intermediate query counts around the pilot-estimated crossover, using a deterministic refinement rule frozen before formal timing.

The purpose is to resolve the crossing more precisely than the previous interval of 524,288–1,048,576 queries.

Ensure that each method runs at least three measured points on either side of an observed crossing whenever the compute cap permits.

For the main result, the crossing interval must be determined only from common actually executed query counts.

Also check whether the existing FH-MLP implementation is limited by avoidable Python or GPU-launch overhead. If a numerically equivalent batched execution is practical, audit it on the same queries and random samples, and report the faster validated FH-MLP implementation as an additional conservative comparison. Do not replace the original recursion with a different estimator.

## 7. Break-even definition

For an accuracy-qualified pair, define the observed break-even as the first measured common query count \(N\) satisfying

\[
T_{\mathrm{FH}}(N)\geq T_{\mathrm{TA}}(N).
\]

Report:

- Last measured query count before crossover
- First measured query count at or after crossover
- Mean and standard deviation of both methods' times
- Individual replicate crossover intervals
- Accuracy and selected hyperparameters
- Actual H100 GPU allocation

Do not claim an exact crossing point from interpolation.

If no crossing is observed, report the maximum executed query count and `NOT_REACHED_WITHIN_TESTED_BUDGET`.

Do not extend the experiment selectively after examining whether a result looks favorable.

## 8. Publication-quality crossover figure

This is a major deliverable.

Create a clean, attractive figure suitable for a two-column ICML paper.

### Main panel

Title:

**High-Accuracy Amortization on a 100D Joint Nonseparable PDE**

Use:

- X-axis: Number of query points, logarithmic base 2
- Y-axis: Measured cumulative H100 computation (minutes)
- Tangent-Actor (direct value head): blue solid line
- Full-History MLP: orange solid line
- Circular markers at actual measured counts
- Light uncertainty ribbons showing ±1 SD across timing replicates
- A subtle shaded interval identifying the observed crossover bracket
- A vertical dotted annotation for the first measured crossing

Display the accuracy tolerance clearly, for example:

`Common accuracy requirement: RMSE ≤ 1e-2`

Do not use an interpolated accuracy–compute frontier or synthetic runtime points.

### Crossover zoom

If an observed crossing exists, add an inset magnifying the region around the intersection.

The inset should make both curves and the crossing interval clearly distinguishable.

Use actual measurements in both the main panel and inset.

Annotate the measured crossover interval and the relative computation saved at a representative observed query count beyond the crossing.

If no crossover is observed, omit the crossing annotation rather than drawing a fabricated intersection.

### Additional results

If the stricter 0.005 experiment also produces a qualified measured crossover, generate an additional panel or supplementary figure.

Plot TA-FK-MCV separately. It need not cross FH-MLP; a negative result must remain visible.

Use consistent professional typography, restrained colors, readable tick labels, and no unnecessary decorations.

Export:

- PDF vector figure
- 300-dpi PNG
- SVG
- Figure source code
- Publication-ready LaTeX caption

Render and inspect the final PDF for clipping, overlapping annotations, and legibility at ICML column width.

## 9. Resource budget and execution gates

Set a hard cap of **30 allocated H100 GPU-hours** for this new study, including preflight, failed allocations, calibration, formal timing, and retries.

Before formal execution, estimate the expected cost of the complete eligible workload.

Reserve resources for at least three independent timing replicates.

Execution order:

1. Verify source implementations and mathematical correctness.
2. Prepare fresh banks and validate references.
3. Conduct truth-free H100 performance calibration.
4. Run FH-MLP accuracy calibration.
5. Freeze selected configurations.
6. Verify final accuracy qualification.
7. Run actual-prefix cumulative timing on the eligible methods.
8. Compute crossover statistics.
9. Generate publication figures, reports, and audits.

Do not use final accuracy results to change configurations.

Do not launch the expensive timing stage for an accuracy-ineligible solver.

If the cost gate fails, report the projected cost and completed pilot results. Never silently reduce mathematical fidelity to fit the budget.

## 10. Required artifacts

Create isolated directories:

`invariant_region_mlp/experiments/high_accuracy_fhmlp_break_even/`

`invariant_region_mlp/results/high_accuracy_fhmlp_break_even/`

Deliver:

- `PREREGISTRATION.md`
- Frozen protocol and input manifests
- Mathematical and source-fidelity audit
- Calibration results for every attempted \((n,M)\)
- Frozen selections and hashes
- Independent final accuracy results
- Actual cumulative timing CSV
- Per-replicate timing records
- Break-even summary CSV and JSON
- Resource accounting
- Artifact audit
- Publication-quality PDF/PNG/SVG figures
- A concise Markdown scientific report
- A ready-to-insert ICML figure caption and results paragraph

Preserve every numerical failure, OOM, time cap, and infrastructure retry in the audit.

The final scientific conclusion must distinguish between an observed crossover, an extrapolated possibility, and a failure to reach the prescribed accuracy.

## 11. Completion

Complete the feasible study under the frozen allocation cap, run correctness tests, validate the generated artifacts, and inspect the publication figures.

Commit and push the isolated branch.

Do not modify the existing ICML manuscript or previous audited experiment directories.

In the final response, provide:

1. Branch and exact commit SHA.
2. Accuracy-qualified FH-MLP \((n,M)\) at both thresholds.
3. TA-direct and TA-FK-MCV accuracy.
4. Measured crossing intervals, if any.
5. Cumulative times at and beyond the crossing.
6. Actual H100 GPU-hours consumed.
7. Exact paths of the PDF figures and raw timing CSVs.
8. Any limitations that affect scientific interpretation.

**The objective is to determine whether higher accuracy requirements move the genuine amortization crossover earlier, not merely to produce a visually attractive crossing.**
