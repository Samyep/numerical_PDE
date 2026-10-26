#!/bin/bash
#SBATCH --partition=coe-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:10:00
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
source "${script_dir}/slurm_env.sh"
phase="${1:?phase is required}"
python -m invariant_region_mlp.experiments.whole_function_compute_scaling_3seed.run_study \
  train --phase "${phase}" --index "${SLURM_ARRAY_TASK_ID:?array index required}" --device cuda

