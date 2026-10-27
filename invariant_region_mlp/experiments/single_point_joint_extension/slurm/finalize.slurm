#!/bin/bash
#SBATCH --job-name=joint_final
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_joint_extension/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_joint_extension/slurm/common_env.sh"

pending=$(python -m invariant_region_mlp.experiments.single_point_joint_extension.formal pending)
if [[ -n "$pending" ]]; then
  echo "formal bundles still pending: $pending" >&2
  exit 3
fi
python -m invariant_region_mlp.experiments.single_point_joint_extension.aggregate
