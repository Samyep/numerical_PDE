#!/bin/bash
#SBATCH --job-name=sp32_v3
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%A_%a.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%A_%a.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

manifest=${BUNDLE_MANIFEST:?submit with BUNDLE_MANIFEST}
index=${SLURM_ARRAY_TASK_ID:?submit as a task array}
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -m invariant_region_mlp.experiments.single_point_32_ice.consolidated_dispatch \
  run-bundle --manifest "$manifest" --index "$index" --device cuda
