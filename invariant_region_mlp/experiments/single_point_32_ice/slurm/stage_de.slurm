#!/bin/bash
#SBATCH --job-name=sp32_de
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

python -m invariant_region_mlp.experiments.single_point_32_ice.aggregate
python -m pytest -q invariant_region_mlp/experiments/single_point_32_ice
