#!/bin/bash
#SBATCH --partition=ice-cpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --job-name=wcf_gate
set -euo pipefail
source "${WCF_SCRIPT_DIR:-${SLURM_SUBMIT_DIR}/invariant_region_mlp/experiments/whole_function_compute_frontier}/slurm_env.sh"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.run_frontier mark-ready
