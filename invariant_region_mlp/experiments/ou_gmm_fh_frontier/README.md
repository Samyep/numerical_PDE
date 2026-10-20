# OU/GMM HJB — one-axis FH-MLP time/accuracy curve

## Fixed scientific question

How many **measured GH200 GPU seconds** does raw Full-History MLP require
to reach the same whole-function Value RMSE as frozen Tangent-Actor + MCV?

The comparison is the *exact same* 100D OU/GMM HJB, with 960 held-out
states and reference from the October 19 final-suite artifacts:

- Tangent-Actor frozen Value RMSE: **0.2197627563608737**.
- Tangent-Actor frozen total time: **1675.9736826 GPU seconds**
  = 1539.5905155 training/selection + 136.3831671 inference.
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

## Run on the GH200 cluster

Use a separate worktree to preserve the ongoing whole-function sweep:

    cd /work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
    git fetch origin ir-mlp-ou-gmm-fh-time-accuracy
    git worktree add ../source_ou_gmm_fh origin/ir-mlp-ou-gmm-fh-time-accuracy
    cd ../source_ou_gmm_fh
    source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
    conda activate base
    SOURCE=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
    python -m pytest -q invariant_region_mlp/experiments/ou_gmm_fh_frontier/test_runner.py
    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner --source-root "$SOURCE" --all --dry-run

When scheduling permits **one GH200 GPU job**, run:

    sbatch --export=ALL,JOB_SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source_ou_gmm_fh invariant_region_mlp/experiments/ou_gmm_fh_frontier/ou_gmm_fh_frontier.slurm

This requests 2h20 of Slurm wall time but caps *completed numerical
FH work* to 2h cumulatively, with at most 1h per M. The status JSON
records OUT_OF_TIME/FAILED and partial progress without an RMSE.
Check Slurm sacct for actual allocated GPU-hours, including failed
allocations. Do not submit concurrently if account rules forbid it.

A single already-declared M may also be run in an allocated GPU session:

    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner --source-root "$SOURCE" --M 16 --device cuda

## Produce the requested figure (only from real runs)

    python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.plot_frontier --source-root "$SOURCE" --out-dir invariant_region_mlp/results/ou_gmm_fh_frontier

It creates an actual measured runtime–RMSE curve as PNG/PDF, a CSV,
and a JSON crossing report, with the Actor target RMSE as a horizontal
line and its measured 1676-second total time as a marker.

The horizontal axis is *real 960-point GH200 inference time*; do not
multiply an already-batched measurement by 960 again.

If no completed FH point reaches 0.21976276, report
NOT_REACHED_IN_MEASURED_GRID; do not interpolate, extrapolate the
RMSE curve, or claim FH cannot achieve the accuracy at larger M/n.
If FH does reach the threshold, report its first observed fixed-grid
time and confirm with independent FH seeds before a paper headline.

Status: code and execution protocol are on GitHub. The live GH200 job
has NOT been submitted by ChatGPT; no new measured FH OU/GMM curve exists
until your cluster executes the supplied script.
