#!/bin/bash
#SBATCH --job-name=dbbs-formal
#SBATCH --partition=ice-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --array=0-41%16
#SBATCH --output=invariant_region_mlp/results/deep_bsde_budget_scaling_3seed/slurm/formal-%A_%a.out
#SBATCH --error=invariant_region_mlp/results/deep_bsde_budget_scaling_3seed/slurm/formal-%A_%a.err
set -euo pipefail
script_dir="${DBBS_SCRIPT_DIR:-${SLURM_SUBMIT_DIR}/invariant_region_mlp/experiments/deep_bsde_budget_scaling_3seed/slurm}"
source "$script_dir/common_env.sh"
python -m invariant_region_mlp.experiments.deep_bsde_budget_scaling_3seed.worker \
  --index "${SLURM_ARRAY_TASK_ID:?array task is required}" --device cuda

