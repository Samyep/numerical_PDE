# OU/GMM HJB — one-axis FH-MLP time/accuracy curve

## Fixed scientific question

How many **measured RTX 5070 Ti GPU seconds** does raw Full-History MLP require
to reach the same whole-function Value RMSE as frozen Tangent-Actor + MCV?

The comparison is the *exact same* 100D OU/GMM HJB, with 960 held-out
states and reference from the October 19 final-suite artifacts:

- Tangent-Actor frozen Value RMSE: **0.2197627563608737**.
- Tangent-Actor frozen total time: **1675.9736826 GH200 GPU seconds**
  = 1539.5905155 training/selection + 136.3831671 inference.
  This historical time is provenance only and is not plotted on the RTX time axis.
- Full-History depth is **n=4**, beta-time sampler alpha=0.5,
  MC bases **M = 4, 8, 12, 16, 24, 32, 48, 64** in ascending order.
- Seed 0, fixed 960 states. No FH depth search or Actor retraining.
  Stop after the first completed threshold hit, or once the
  precommitted two-hour total numerical-compute budget expires.

This is a *slice* through the FH compute frontier, not a guarantee of
optimal runtime over every possible FH depth.

## Why a new FH code path is mathematically necessary

The PDE contains the OU drift x dot grad(u):

    u_t + 0.5 Delta u + x dot grad u - 0.5||grad u||² - d = 0.

The old Brownian-only FH control explicitly excludes OU/GMM.
The new implementation uses the **correct OU process**:

    X_(t+h) = exp(h) x + sqrt((exp(2h)-1)/2) Z,

the OU BEL gradient score exp(h)Z/sqrt((exp(2h)-1)/2),
and the known terminal gradient times exp(h).

To make FH strong rather than artificially noisy, the *known spatially
constant* -d term is integrated analytically as -d*(T-t), with zero
gradient. No nonterminal exact PDE solution, learned Actor, or oracle
reduced coordinates enter the FH recurrence. The algorithm otherwise
uses the original multilevel telescoping recursion.

## Required frozen original HPC inputs

The original source tree must still contain:

    invariant_region_mlp/results/tangent_actor_icml_final_suite/MAIN_RESULTS.csv
    invariant_region_mlp/results/tangent_actor_icml_final_suite/inputs/ou_gmm_d100.npz
    invariant_region_mlp/results/tangent_actor_icml_final_suite/paired_predictions/main_ou_gmm_d100_tangent_actor_s0_n64p0_p1024_mcv.npz

The runner verifies test-bank hashes, actor predictions/target RMSE,
and the exact evolving-GMM truth numerically. It will fail rather than
fabricate or silently regenerate missing inputs.

## Run locally on the RTX 5070 Ti

Use a separate worktree and run the fail-closed audit first:

    $SOURCE = (Get-Location).Path
    python -m pytest -q invariant_region_mlp/experiments/ou_gmm_fh_frontier/test_runner.py
    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.audit --source-root "$SOURCE"
    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner --source-root "$SOURCE" --all --dry-run

Only after the audit passes, run one local GPU process:

    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner --source-root "$SOURCE" --all --device cuda

The runner caps charged numerical FH work to 2h cumulatively, with at
most 1h per M. The status JSON
records OUT_OF_TIME/FAILED and partial progress without an RMSE.
Do not run another major GPU process concurrently.

A single already-declared M may also be run in an allocated GPU session:

    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner --source-root "$SOURCE" --M 16 --device cuda

## Produce the requested figure (only from real runs)

    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.plot_frontier --source-root "$SOURCE" --out-dir invariant_region_mlp/results/ou_gmm_fh_frontier

It creates an actual measured runtime–RMSE curve as PNG/PDF, a CSV,
and a JSON crossing report, with the Actor target RMSE as a horizontal
line. The historical 1676-second GH200 total is text provenance only,
not a same-hardware timing marker.

The horizontal axis is *real 960-point RTX 5070 Ti inference time*; do not
multiply an already-batched measurement by 960 again.

If no completed FH point reaches 0.21976276, report
NOT_REACHED_WITHIN_TESTED_BUDGET; do not interpolate, extrapolate the
RMSE curve, or claim FH cannot achieve the accuracy at larger M/n.
If FH does reach the threshold, report its first observed fixed-grid
time and confirm with independent FH seeds before a paper headline.

Status: a measured curve exists only after the immutable-artifact audit
passes and the local RTX sweep completes. Missing original artifacts are
reported as BLOCKED_MISSING_ARTIFACTS; they are never regenerated silently.
