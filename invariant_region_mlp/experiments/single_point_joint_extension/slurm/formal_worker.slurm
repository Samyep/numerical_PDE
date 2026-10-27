#!/bin/bash
#SBATCH --job-name=joint_formal
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%A_%a.out
#SBATCH --error=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%A_%a.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_joint_extension/slurm/common_env.sh"

python -m invariant_region_mlp.experiments.single_point_joint_extension.formal \
  run --index "${SLURM_ARRAY_TASK_ID:?array index required}" --device cuda
