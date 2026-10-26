#!/bin/bash
#SBATCH --partition=coe-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=00:45:00
set -euo pipefail
script_dir="${WCS_SCRIPT_DIR:-${SLURM_SUBMIT_DIR}/invariant_region_mlp/experiments/whole_function_compute_scaling_3seed}"
source "${script_dir}/slurm_env.sh"
python -m invariant_region_mlp.experiments.whole_function_compute_scaling_3seed.run_study \
  evaluate --index "${SLURM_ARRAY_TASK_ID:?array index required}" --device cuda

