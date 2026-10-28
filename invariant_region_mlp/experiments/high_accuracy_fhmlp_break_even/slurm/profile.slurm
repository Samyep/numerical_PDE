#!/bin/bash
#SBATCH --job-name=hafh_profile
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=invariant_region_mlp/results/high_accuracy_fhmlp_break_even/slurm/%x_%A_%a.out
#SBATCH --error=invariant_region_mlp/results/high_accuracy_fhmlp_break_even/slurm/%x_%A_%a.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from experiment checkout root}}
source "$SOURCE_ROOT/invariant_region_mlp/experiments/high_accuracy_fhmlp_break_even/slurm/common_env.sh"
manifest=${TASK_MANIFEST:?TASK_MANIFEST is required}
index=${SLURM_ARRAY_TASK_ID:?submit as an array}
expected_sha=${CODE_SHA:?CODE_SHA is required}
actual_sha=$(git -C "$SOURCE_ROOT" rev-parse HEAD)
[[ "$actual_sha" == "$expected_sha" ]] || { echo "checkout moved: expected $expected_sha, found $actual_sha" >&2; exit 3; }
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -m invariant_region_mlp.experiments.high_accuracy_fhmlp_break_even.dispatch run-task --manifest "$manifest" --index "$index" --device cuda

