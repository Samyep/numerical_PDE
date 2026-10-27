#!/bin/bash
#SBATCH --job-name=spbe_prepare
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_break_even/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_break_even/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from experiment checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_break_even/slurm/common_env.sh"

python -m pytest -q \
  invariant_region_mlp/experiments/single_point_break_even/test_break_even.py
python -m invariant_region_mlp.experiments.single_point_break_even.prepare

