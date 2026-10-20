#!/bin/bash
#SBATCH --job-name=ou_gmm_fh_curve
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=120G
#SBATCH --time=02:20:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/ou_gmm_fh_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/ou_gmm_fh_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "$JOB_SOURCE_ROOT"
if [ -z "$JOB_SOURCE_ROOT" ]; then
    echo "Set JOB_SOURCE_ROOT to the separate OU-FH Git worktree" >&2
    exit 2
fi

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$JOB_SOURCE_ROOT"
export PYTHONPATH="$JOB_SOURCE_ROOT"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "source_checkout=$SOURCE_ROOT"
echo "ou_fh_code_checkout=$JOB_SOURCE_ROOT"
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -m pytest -q invariant_region_mlp/experiments/ou_gmm_fh_frontier/test_runner.py

# One GPU, frozen n=4, one-dimensional M sweep.
# Do not run concurrently with an existing whole-function GPU job if the
# account/partition permits only a single queued/running job.
srun --network=no_vni python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.runner \
  --source-root "$SOURCE_ROOT" \
  --out-dir "$JOB_SOURCE_ROOT/invariant_region_mlp/results/ou_gmm_fh_frontier" \
  --all --device cuda --max-cell-seconds 3600 --max-sweep-seconds 7200

python -m invariant_region_mlp.experiments.ou_gmm_fh_frontier.plot_frontier \
  --source-root "$SOURCE_ROOT" \
  --out-dir "$JOB_SOURCE_ROOT/invariant_region_mlp/results/ou_gmm_fh_frontier"

# The runner's inference seconds do not replace Slurm GPU-hour accounting.
echo "Check sacct for physical allocated GPU-hours including timed-out attempts."
