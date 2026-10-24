#!/bin/bash
#SBATCH --job-name=joint_gate
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_joint_extension/slurm/common_env.sh"

nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -m pytest -q invariant_region_mlp/experiments/single_point_joint_extension/test_joint_extension.py
python -m invariant_region_mlp.experiments.single_point_joint_extension.preflight --device cuda
python -m invariant_region_mlp.experiments.single_point_joint_extension.smoke --device cuda
