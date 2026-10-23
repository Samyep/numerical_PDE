#!/bin/bash
# Shared, location-independent environment setup for the expansion suite.
# Callers may override any path/module through the documented environment
# variables without editing a tracked file.

EXPANSION_SLURM_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=${SOURCE_ROOT:-$(cd "$EXPANSION_SLURM_DIR/../../../.." && pwd)}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-$(cd "$SOURCE_ROOT/.." && pwd)/external}
DPI_ROOT=${DPI_ROOT:-$EXTERNAL_ROOT/DeepPicardIteration}
EXPANSION_VENV=${EXPANSION_VENV:-$SOURCE_ROOT/.venv-main-highdim}
PYTHON_MODULE=${PYTHON_MODULE:-python/3.11.9}
CUDA_MODULE=${CUDA_MODULE:-cuda/13.0.1}

if [[ ! -x "$EXPANSION_VENV/bin/python" ]]; then
  echo "missing expansion Python environment: $EXPANSION_VENV" >&2
  echo "set EXPANSION_VENV or run slurm/setup_ice.sh" >&2
  return 2 2>/dev/null || exit 2
fi
if [[ ! -d "$DPI_ROOT/.git" ]]; then
  echo "missing pinned DeepPicardIteration checkout: $DPI_ROOT" >&2
  return 2 2>/dev/null || exit 2
fi

# PACE's login environment can retain dependencies of the default modules.
# Force the purge; Lmod may still print harmless dependency-reload notices.
module --force purge
module load "$PYTHON_MODULE"
module load "$CUDA_MODULE"
source "$EXPANSION_VENV/bin/activate"

cd "$SOURCE_ROOT"
export SOURCE_ROOT DPI_ROOT EXPANSION_VENV EXTERNAL_ROOT
export PYTHONPATH="$DPI_ROOT:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
# Do not inherit DeltaAI's HPE-Slingshot `hsn` assumption. These jobs use one
# GPU each, so NCCL may auto-detect ICE's fabric. A task-specific override is
# still available when explicitly needed.
if [[ -n "${EXPANSION_NCCL_SOCKET_IFNAME:-}" ]]; then
  export NCCL_SOCKET_IFNAME="$EXPANSION_NCCL_SOCKET_IFNAME"
else
  unset NCCL_SOCKET_IFNAME
fi
export PYTORCH_CUDA_ALLOC_CONF=${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}
export WANDB_MODE=${WANDB_MODE:-offline}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-8}}
