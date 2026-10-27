# Task: Deep BSDE Hyperparameter Budget Scaling for the ICML Tangent-Actor Paper

Repository: `Samyep/numerical_PDE`

Create a new branch named `ir-mlp-deep-bsde-budget-scaling-3seed`.

Use `ir-mlp-compute-scaling-3seed` as the primary source for the existing compute-scaling experiment and its audited results. Consult `ir-mlp-single-point-break-even` for the audited query-wise Deep BSDE implementation. Preserve the existing experiments and artifacts.

## 1. Scientific objective

We want to address one specific reviewer concern:

**Some baseline methods were trained with less computational budget than Tangent-Actor. Can increasing their hyperparameter budgets enable them to reach comparable accuracy?**

We already have completed, audited C0/C1/C2 sweeps for DPI and PINN-HTE. We now need equivalent Deep BSDE experiments.

Do not construct a Pareto frontier or fit an accuracy–compute curve. Report actual outcomes at prespecified hyperparameter configurations.

**Only Deep BSDE training is authorized for this task.**

Do not retrain Tangent-Actor, DPI, PINN-HTE, SDGD-PINN, or FH-MLP.

## 2. Existing experiments and benchmarks

Read and preserve:

- `invariant_region_mlp/experiments/whole_function_compute_scaling_3seed/`
- `invariant_region_mlp/results/whole_function_compute_scaling_3seed/scaling_summary.csv`
- `invariant_region_mlp/experiments/original_pde_benchmark_comparison/deep_bsde_port.py`
- The native scalar-\(Y_0\) Deep BSDE implementation under `single_point_break_even_joint` on the later branch

Use the same benchmark definitions and parameters as the previous whole-function study:

- B6 Full-Rank HJB, 100D, \(T=0.25\)
- Dense Fourier coupled PDE, 100D, \(T=0.5\)

Do not create simplified substitutes for these PDEs.

Verify the exact driver, diffusion coefficient, terminal condition, reference solution, and parameter hashes before running.

## 3. Evaluation points

Use the original frozen 960-state whole-function test bank for each PDE.

From each bank, select two disjoint subsets:

- Calibration: 32 query points
- Final evaluation: 64 query points

Freeze the selected indices and hashes before formal experiments. Use stratified random sampling over the time coordinate, based only on query inputs, to ensure reasonable coverage of \([0,T)\).

Do not use reference values or observed errors to select points.

Use training seeds `0,1,2` for every formal configuration.

**Crucial correctness requirement:** The original Deep BSDE implementation typically solves for \(u(0,x_0)\), but our frozen test points are \((t_i,x_i)\) with varying \(t_i\).

For every query, the adapted Deep BSDE solver must:

1. Initialize the forward process at \(X_{t_i}=x_i\).
2. Integrate over the remaining horizon \(T-t_i\).
3. Set the Euler timestep to \((T-t_i)/N\).
4. Use the original PDE's diffusion coefficient, drift, terminal condition, and driver.
5. Return its native learned scalar \(Y_{t_i}\) as the prediction.

Do not accidentally solve \(u(0,x_i)\) and compare it against \(u(t_i,x_i)\).

Check the exact sign convention and \(z=\sigma^\top\nabla u\). In these benchmarks, \(\sigma=\sqrt2 I\).

For B6, verify \(f(u,z)=-\frac12\|z\|^2\).

For Dense Fourier, verify \(f(u,z)=\frac12\tanh(u)\|z\|^2\).

Preserve the source-faithful Deep BSDE structure: scalar initial value, learned initial \(Z\), time-indexed subnetworks, terminal loss, and optimizer implementation. Do not introduce an amortized decoder.

Add unit tests for nonzero initial times and variable horizons.

## 4. Frozen Deep BSDE hyperparameter sweeps

Define three axes:

| Sweep | Hyperparameter | Level 1 | Level 2 | Level 3 |
|---|---|---:|---:|---:|
| C0 | Optimizer updates | 1,000 | 2,000 | 4,000 |
| C1 | Brownian paths per optimizer update (batch size) | 32 | 64 | 128 |
| C2 | Euler time-discretization steps | 10 | 20 | 40 |

The common reference configuration is:

`updates=2000, batch_size=64, time_steps=20`

For C0, keep C1 and C2 fixed at the reference configuration.

For C1, keep C0 and C2 fixed.

For C2, keep C0 and C1 fixed.

All other model and optimizer settings must remain unchanged across the sweeps. Start from the source-faithful Deep BSDE HJB profile: float64, hidden widths `(d+10,d+10)`, learning rate `1e-2`, clipped terminal loss, and the original time-slice network structure. Record and justify any benchmark-specific adaptation before examining test errors.

C2 is specifically a **time-discretization budget** for Deep BSDE, not an optimizer-step budget. Make that distinction explicit in the table.

For C0, train one continuous 4,000-update trajectory per query and seed, saving genuine 1,000/2,000/4,000 checkpoints. C1 and C2 require independent training when their configurations differ.

Never treat independent restarts as checkpoints of the same trajectory.

## 5. Joint scaling experiments

Additionally evaluate three simultaneous C0/C1/C2 configurations:

| Joint configuration | Updates | Batch size | Euler steps |
|---|---:|---:|---:|
| J1 | 1,000 | 32 | 10 |
| J2 | 2,000 | 64 | 20 |
| J3 | 4,000 | 128 | 40 |

Use the same 64 final queries and three seeds.

When multiple configurations are mathematically and computationally identical, reuse the exact checkpoint instead of retraining. In particular, the middle C1/C2 settings and J2 coincide with the common reference configuration.

Verify configuration signatures, query-input hashes, seed identity, and checkpoint equivalence before reusing any artifact.

Do not use joint scaling results to replace or retrospectively select individual-sweep results.

## 6. Compute preflight and execution controls

Deep BSDE must train a fresh model for every distinct query state.

Before formal execution, perform a small truth-free H100 timing and memory preflight. Estimate total GPU allocation for all planned configurations, seeds, and queries.

Set a documented maximum total allocation of 120 H100 GPU-hours for this new study, or a stricter limit if cluster resources require it.

If the projected complete experiment exceeds the available budget, stop before submitting the formal grid. Report the cost breakdown and feasibility blocker; do not silently remove configurations, reduce the 64-point test set, or shorten the official schedule.

Run genuine numerical jobs on NVIDIA H100. Local CPU tests are allowed, but they are not formal numerical results.

Implement resumable checkpoints, deterministic seeds, atomic output writing, scheduler accounting, and explicit handling of numerical failures.

If any query produces nonfinite internal trajectories, loss, gradients, or predictions, preserve the failure record. Do not silently drop failed queries and report RMSE on the remaining finite subset as if it represented all 64 points.

## 7. Accuracy and time metrics

For each formal configuration and training seed, compute RMSE over all 64 final query points:

\[
\mathrm{RMSE}_s
=
\sqrt{\frac1{64}\sum_{i=1}^{64}
(\widehat u_s(t_i,x_i)-u^\star(t_i,x_i))^2}.
\]

Report the arithmetic mean and sample standard deviation across the three seed-level RMSE values.

The final test reference values must not influence training, checkpoint selection, hyperparameter decisions, or stopping.

Record:

- RMSE mean ± seed SD
- RelMAE and relative L2
- Number of successful and failed queries
- Mean training seconds per query
- Total active training GPU seconds required for the full 64-query workload, computed separately per seed and then averaged
- Peak GPU memory
- Actual Slurm GPU allocation, including failed attempts and retries

Use cumulative 64-query training time as the main displayed computation cost for Deep BSDE. Do not display a single-query time as though it were the cost of obtaining a whole-function solution.

## 8. Required ICML tables

Produce a publication-ready LaTeX table in the same style as the existing three-seed compute-scaling table.

Keep DPI and PINN-HTE's historical results unchanged.

Remove Tangent-Actor rows.

Add Deep BSDE rows for C0, C1, and C2, with exactly three configurations per sweep and per PDE.

Required columns:

| Method | Sweep | B6 Level 1 | B6 Level 2 | B6 Level 3 | Fourier Level 1 | Fourier Level 2 | Fourier Level 3 |
|---|---|---|---|---|---|---|---|

For every Deep BSDE cell, show:

`10^3 × RMSE (mean ± SD)`

with its actual hyperparameter value and total 64-query training minutes underneath.

For the existing DPI and PINN-HTE cells, preserve the original parameter values, training times, and RMSE values.

**Important comparability limitation:** Existing DPI/PINN results use the 960-point frozen evaluation bank, while Deep BSDE uses 64 points. Do not present these as identical evaluation protocols.

If the original trained DPI/PINN checkpoints or per-point predictions are available, compute their accuracy on the identical 64-point subset **without any retraining**, and create a supplementary matched-evaluation table. Otherwise explicitly distinguish the evaluation sets in the main table caption and do not claim a direct paired accuracy win from that table alone.

Produce a second, smaller table for Deep BSDE joint configurations J1/J2/J3.

No new TA results or TA rows.

## 9. Deliverables

Create a new isolated experiment directory for Deep BSDE budget scaling.

Provide:

- Frozen preregistration and configuration manifest
- Code and correctness tests
- Calibration/final query index manifests with hashes
- Raw per-query training and failure records
- Per-seed accuracy and timing CSVs
- Aggregated C0/C1/C2 and joint-scaling CSVs
- Publication-ready `.tex` tables
- Rendered PDF table previews
- Resource accounting and artifact audit
- A concise Markdown report describing accuracy trends, actual computation spent, numerical instability, and interpretation limitations

The conclusion must distinguish between:

- Not reaching a target within the evaluated budget
- Numerical instability
- A genuine observed improvement
- Evidence of plateauing or nonmonotonic behavior

Never claim that Deep BSDE cannot converge at arbitrarily large compute budgets.

## 10. Completion

First audit the mathematical adaptation and run unit/smoke tests. Then perform timing preflight and freeze the experiment grid. Execute formal jobs only if the allocation gate passes.

Complete all permitted jobs, aggregate the results, render and inspect the tables, and verify artifacts and hashes.

Commit and push the new branch without modifying or overwriting the previous audited experiment branches.

In the final report, give the actual commit SHA, completed/expected job counts, all failed configurations, total H100 GPU-hours, and the paths of the final tables.

The primary objective is not to make Deep BSDE look worse. It is to give it larger, well-defined computational budgets and establish whether it actually closes the observed accuracy gap.
