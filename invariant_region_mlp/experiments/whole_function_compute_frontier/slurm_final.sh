#!/bin/bash
#SBATCH --partition=ice-cpu
#SBATCH --account=isye
#SBATCH --qos=coe-ice
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --job-name=wcf_final
set -euo pipefail
source "$(dirname "$0")/slurm_env.sh"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.collect_slurm_accounting --start "${WCF_SLURM_START:?missing WCF_SLURM_START}"
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.aggregate
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.plot_frontier
python -m invariant_region_mlp.experiments.whole_function_compute_frontier.report
