#!/bin/bash
#SBATCH --job-name=sp32_audit_v3
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

manifest=${BUNDLE_MANIFEST:?submit with BUNDLE_MANIFEST}
python -m invariant_region_mlp.experiments.single_point_32_ice.consolidated_dispatch \
  audit --manifest "$manifest"
