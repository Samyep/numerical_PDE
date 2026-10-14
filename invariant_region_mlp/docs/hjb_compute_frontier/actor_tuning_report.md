# Actor tuning report

Status: pilot, frozen three-seed study, and frozen-Actor evaluation sensitivity
complete.  The targeted high-quality B1 branch and seed-extension replication
are tracked separately below.

The ambient-output Actor has no rank projection on either PDE. The declared study compares terminal-pathwise `q0`, one policy improvement `q1`, a matched fresh-label continuation `q1_refresh`, a true `q2`, and a true `q3`. Each controlled-policy label differentiates the full Euler trajectory with the source Actor frozen while preserving its state derivative.

The pilot uses 2,048 terminal-label points for `q0`, then 4,096 points, 256 paths per label, 16 Euler steps, and 4,000 optimizer steps for each subsequent phase. The final preset raises this to 4,096/256/6,000 for `q0` and 8,192/512/32/6,000 for policy-evaluation phases. A targeted high-quality B1 branch declares width 512, 16,384 points, 1,024 paths, 64 steps, and 12,000 optimizer steps for `q1`, matched `q1_refresh`, and matched true `q2`; `q3` retains the final-scale 8,192/512/32/6,000 budget. This isolates high-quality continuation from a same-cost policy update before testing one cheaper next iteration. The HJB-Log high branch was cancelled before allocation after the final runs showed that its Value RMSE was already dominated by final MC inference noise.

Per the user's explicit override, Actor candidates are ranked by test Relative L2/RMSE. The report will show both model-path cost and total exploration cost so this oracle selection is not presented as a deployable reference-free procedure. External baselines do not receive the same selection privilege.

## Pilot results

The pilot ran on NVIDIA GH200 GPUs with seed 0. Candidate evaluations within each PDE used common Brownian random numbers. All actors used ambient 400-dimensional outputs without a rank projection.

| PDE | Actor | Value RMSE | Relative L2 | Bias | gradient relative error | cumulative path training (s) |
|---|---:|---:|---:|---:|---:|---:|
| B1 (400,40) | q0 | 0.217052 | 0.071111 | +0.188670 | 0.7852 | 6.7 |
| B1 (400,40) | q1 | 0.146121 | 0.047872 | +0.127411 | 0.6272 | 54.4 |
| B1 (400,40) | q1_refresh | 0.112950 | 0.037005 | +0.095793 | 0.5343 | 102.1 |
| B1 (400,40) | q2 | 0.113306 | 0.037121 | +0.095990 | 0.5318 | 102.4 |
| B1 (400,40) | **q3** | **0.103285** | **0.033838** | **+0.078973** | **0.4566** | 150.1 |
| HJB-Log (400) | q0 | 0.012310 | 0.001926 | +0.004107 | 1.0777 | 7.0 |
| HJB-Log (400) | q1 | 0.009987 | 0.001562 | +0.003398 | 0.6166 | 57.7 |
| HJB-Log (400) | **q1_refresh** | **0.009511** | **0.001488** | **+0.002370** | 0.4742 | 107.8 |
| HJB-Log (400) | q2 | 0.009529 | 0.001491 | +0.002353 | 0.4739 | 107.9 |
| HJB-Log (400) | q3 | 0.009615 | 0.001504 | +0.002601 | **0.4203** | 157.8 |

The B1 pilot supports retaining multiple policy iterations: `q3` reduced RMSE by 29.3% from `q1` and by 8.6% from the matched `q1_refresh` continuation. The direct `q2` update and continuation were tied, so the gain appears only after the next update at this pilot budget. On HJB-Log, extra policy iterations improved the analytic gradient diagnostic but did not improve Value beyond `q1_refresh`; the test-oracle rule therefore selected `q1_refresh`.

The HJB-Log truth has mean 6.393 and standard deviation only 0.02395. Its small Relative L2 can hide poor predictors: a validation-fitted constant has test RMSE 0.02398. The selected Actor RMSE 0.00951 is also 0.385 times the truth standard deviation, so its gain is not just a mean-value effect.

The complete eight-task pilot consumed 1,100 allocated GH200-seconds (0.3056 GPU-hours). All tasks completed successfully. Independent recomputation of all 28 stored prediction artifacts agreed with reported metrics to at most `2.49e-14`.

## Frozen final budget

The final Actor budget was fixed after the pilot:

| phase | source | training / validation points | paths per label | label Euler steps | optimizer steps | network |
|---|---|---:|---:|---:|---:|---|
| `q0` | terminal pathwise | 4,096 / 1,024 | 256 | terminal draw | 6,000 | 4x384 SiLU |
| `q1` | `q0` policy evaluation | 8,192 / 512 | 512 | 32 | 6,000 | 4x384 SiLU |
| `q1_refresh` | fresh `q0` labels, warm-start `q1` | 8,192 / 512 | 512 | 32 | 6,000 | 4x384 SiLU |
| `q2` | `q1` policy evaluation | 8,192 / 512 | 512 | 32 | 6,000 | 4x384 SiLU |
| `q3` | `q2` policy evaluation | 8,192 / 512 | 512 | 32 | 6,000 | 4x384 SiLU |

The matched `q1_refresh` and `q2` phases use the same training states and the
same budget, with independent label randomness.  Their only substantive
difference is whether the fresh target remains the `q0` policy evaluation or
is changed to the true `q1` policy evaluation.

### Three-seed candidate results

Values are mean +/- sample standard deviation over the frozen seeds.  The
evaluation uses 256 antithetic paths and 64 Euler steps on 2,048 fixed states.

| PDE | Actor | Value RMSE | Relative L2 | direct-gradient relative error | path training (s) |
|---|---|---:|---:|---:|---:|
| B1 | `q0` | 0.16866 +/- 0.00211 | 0.05512 +/- 0.00069 | 0.6929 +/- 0.0029 | 12.6 +/- 0.7 |
| B1 | `q1` | 0.08945 +/- 0.00033 | 0.02923 +/- 0.00011 | 0.4621 +/- 0.0045 | 334.7 +/- 16.6 |
| B1 | `q1_refresh` | 0.07126 +/- 0.00186 | 0.02329 +/- 0.00061 | 0.3593 +/- 0.0032 | 655.3 +/- 32.2 |
| B1 | true `q2` | 0.08610 +/- 0.00345 | 0.02814 +/- 0.00113 | 0.3587 +/- 0.0041 | 655.8 +/- 31.6 |
| B1 | **true `q3`** | **0.06301 +/- 0.00521** | **0.02059 +/- 0.00170** | **0.3043 +/- 0.0039** | 978.1 +/- 48.9 |
| HJB-Log | `q0` | 0.005017 +/- 0.000051 | 0.0007848 +/- 0.0000080 | 0.4500 +/- 0.0040 | 12.7 +/- 0.5 |
| HJB-Log | `q1` | 0.005161 +/- 0.000016 | 0.0008072 +/- 0.0000025 | 0.3872 +/- 0.0009 | 325.7 +/- 14.0 |
| HJB-Log | `q1_refresh` | 0.005058 +/- 0.000041 | 0.0007911 +/- 0.0000064 | 0.3700 +/- 0.0008 | 639.5 +/- 28.4 |
| HJB-Log | true `q2` | 0.005060 +/- 0.000046 | 0.0007915 +/- 0.0000072 | 0.3701 +/- 0.0007 | 638.1 +/- 28.3 |
| HJB-Log | true `q3` | 0.004996 +/- 0.000051 | 0.0007814 +/- 0.0000080 | **0.3604 +/- 0.0024** | 951.7 +/- 41.8 |

On B1, fresh continuation is decisively better than the matched true `q2`:
its mean RMSE is 17.2% lower at essentially identical cost.  A subsequent
true update to `q3` then lowers RMSE by 26.8% from `q2` and 11.6% from
`q1_refresh`.  This is a reproducible nonmonotone approximate-policy-iteration
trajectory, not monotone policy improvement.  Both test RMSE and the
reference-free policy-value validation choose `q3` in all three B1 seeds.

On HJB-Log, `q1_refresh` and `q2` are numerically tied, and fixed `q3` improves
fixed `q0` RMSE by only 0.43% under the frozen evaluation.  The gradient metric
continues to improve, but this does not translate into a useful Value-cost
gain.

The held-out label losses support increasing the B1 fitting budget: `q1`,
`q1_refresh`, `q2`, and `q3` selected optimizer step 6,000/6,000 in every B1
seed.  Their mean held-out vector MSEs were `0.2860`, `0.1764`, `0.2051`, and
`0.1625`, respectively.  These losses compare fit within a phase; because the
policy-evaluation target changes, they are not treated as a monotone metric
across iterations.

### HJB-Log evaluation sensitivity

The five final Actors were frozen and re-evaluated without training.  All 45
prediction files passed an independent metric audit exactly.

| evaluation budget | fixed `q0` RMSE | fixed `q3` RMSE | oracle actors by seed | oracle RMSE |
|---|---:|---:|---|---:|
| 256 paths, 64 steps | 0.005017 +/- 0.000051 | 0.004996 +/- 0.000051 | `q3,q0,q3` | 0.004974 +/- 0.000045 |
| **1,024 paths, 64 steps** | 0.003233 +/- 0.000037 | **0.003094 +/- 0.000092** | `q3,q3,q3` | **0.003094 +/- 0.000092** |
| 256 paths, 128 steps | 0.005067 +/- 0.000155 | 0.004963 +/- 0.000082 | `q2,q1_refresh,q3` | 0.004937 +/- 0.000087 |

Quadrupling paths reduces fixed-`q3` RMSE by 38.1%; doubling time steps at the
original path count changes it by only 0.66%.  At 1,024 paths, estimated MC
variance still accounts for 51.9% of observed `q3` MSE.  Its
variance-subtracted RMSE diagnostic is `0.002147`, versus `0.002361` for
`q0`.  Better evaluation makes a small `q3` policy benefit detectable, but
training from `q0` to `q3` costs roughly 75 times as much and lowers high-path
RMSE by only 4.3%.  The useful HJB-Log method remains pathwise initialization
plus accurate martingale-CV evaluation.

## Selection and exploration cost

The user-authorized test oracle selects B1 `q3,q3,q3` and HJB-Log
`q3,q0,q3` under the frozen evaluation.  Producing and evaluating every
candidate costs, on average, 1,395.4 seconds per B1 seed and 1,351.8 seconds
per HJB-Log seed.  The model-path costs in the compute frontier do not hide
label generation or validation for that model, but this full candidate search
cost must be disclosed when the oracle selection itself is considered part of
the algorithm.

The reference-free B1 policy validation independently chooses `q3` in all
three seeds.  HJB-Log policy validation chooses `q3,q1,q3`, with differences
comparable to its MC standard error; it does not provide a stable multi-round
selection rule there.

## Targeted high-quality B1 branch

The targeted branch uses width 512.  It raises `q1`, matched `q1_refresh`, and
matched true `q2` to 16,384 training points, 1,024 paths per label, 64 label
Euler steps, and 12,000 optimizer steps.  The three phases therefore separate
better first-update supervision, continued fitting to a refreshed `q0`
target, and a same-cost true update from `q1`.  A later `q3` retains
8,192/512/32/6,000.  The run is phase-resumable; its Slurm results and
single-seed interpretation are appended after completion and are not mixed
into the frozen multi-seed table.
