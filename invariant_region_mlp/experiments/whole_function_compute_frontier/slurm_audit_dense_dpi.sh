#!/bin/bash
#SBATCH --job-name=wcf_dpi_audit
#SBATCH --partition=coe-gpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=00:30:00
set -euo pipefail
source "${WCF_SCRIPT_DIR:-${SLURM_SUBMIT_DIR}/invariant_region_mlp/experiments/whole_function_compute_frontier}/slurm_env.sh"

python -m invariant_region_mlp.experiments.whole_function_compute_frontier.audit_dense_dpi_callback
