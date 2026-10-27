#!/bin/bash
set -euo pipefail

DBBS_SLURM_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DBBS_REPO_ROOT=${DBBS_REPO_ROOT:-$(cd "$DBBS_SLURM_DIR/../../../.." && pwd)}
DBBS_STORAGE_PARENT=$(cd "$DBBS_REPO_ROOT/.." && pwd)
DBBS_ORIGINAL_REPO=${DBBS_ORIGINAL_REPO:-$DBBS_STORAGE_PARENT/numerical_PDE}
DBBS_VENV=${DBBS_VENV:-$DBBS_ORIGINAL_REPO/.venv-main-highdim}
export DBBS_HISTORICAL_RESULTS=${DBBS_HISTORICAL_RESULTS:-$DBBS_ORIGINAL_REPO/invariant_region_mlp/results/main_highdim_expansion}
export DBBS_SCALING_RESULTS=${DBBS_SCALING_RESULTS:-$DBBS_STORAGE_PARENT/ir-mlp-compute-scaling-3seed-worktree/invariant_region_mlp/results/whole_function_compute_scaling_3seed}
source "$DBBS_VENV/bin/activate"
export PYTHONPATH="$DBBS_REPO_ROOT:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export WANDB_MODE=offline
cd "$DBBS_REPO_ROOT"

