#!/bin/bash
set -euo pipefail

WCS_REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
export WCF_HISTORICAL_RESULTS="${WCF_HISTORICAL_RESULTS:-/storage/ice1/4/8/ywang4280/numerical_PDE/invariant_region_mlp/results/main_highdim_expansion}"
export WCF_DPI_ROOT="${WCF_DPI_ROOT:-/storage/ice1/4/8/ywang4280/external/DeepPicardIteration}"
WCS_VENV="${WCS_VENV:-/storage/ice1/4/8/ywang4280/numerical_PDE/.venv-main-highdim}"
source "${WCS_VENV}/bin/activate"
export PYTHONPATH="${WCS_REPO_ROOT}:${WCF_DPI_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export WANDB_MODE=offline
cd "${WCS_REPO_ROOT}"

