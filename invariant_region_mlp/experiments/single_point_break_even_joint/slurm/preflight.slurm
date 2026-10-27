#!/bin/bash
#SBATCH --job-name=spbej_preflight
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=invariant_region_mlp/results/single_point_break_even_joint/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_break_even_joint/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from experiment checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_break_even_joint/slurm/common_env.sh"
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -m pytest -q invariant_region_mlp/experiments/single_point_break_even_joint/test_break_even_joint.py
python -m invariant_region_mlp.experiments.single_point_break_even_joint.preflight

