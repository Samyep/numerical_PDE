#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=${SOURCE_ROOT:-$(cd "$SCRIPT_DIR/../../../.." && pwd)}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-$(cd "$SOURCE_ROOT/.." && pwd)/external}
DPI_ROOT=${DPI_ROOT:-$EXTERNAL_ROOT/DeepPicardIteration}
EXPANSION_VENV=${EXPANSION_VENV:-$SOURCE_ROOT/.venv-main-highdim}
PYTHON_MODULE=${PYTHON_MODULE:-python/3.11.9}
CUDA_MODULE=${CUDA_MODULE:-cuda/13.0.1}
DPI_COMMIT=e8616c2f531c9f95299b246a13c1180b8f02822a

module --force purge
module load "$PYTHON_MODULE"
module load "$CUDA_MODULE"

mkdir -p "$EXTERNAL_ROOT"
if [[ ! -d "$DPI_ROOT/.git" ]]; then
  git clone https://github.com/DeepOptimalControl/DeepPicardIteration.git "$DPI_ROOT"
fi
git -C "$DPI_ROOT" fetch --prune origin
git -C "$DPI_ROOT" checkout --detach "$DPI_COMMIT"
test "$(git -C "$DPI_ROOT" rev-parse HEAD)" = "$DPI_COMMIT"

if [[ ! -x "$EXPANSION_VENV/bin/python" ]]; then
  python -m venv "$EXPANSION_VENV"
  "$EXPANSION_VENV/bin/python" -m pip install --upgrade pip wheel setuptools
  "$EXPANSION_VENV/bin/python" -m pip install torch
fi

"$EXPANSION_VENV/bin/python" -m pip install \
  'numpy<2.5' scipy pandas matplotlib pytest \
  'lightning==2.5.0.post0' 'pytorch-lightning==2.5.0.post0' \
  'h5py==3.13.0' 'yacs==0.1.8' \
  'typer==0.15.2' 'wandb==0.19.7' 'tensorboardX==2.6.2.2'
"$EXPANSION_VENV/bin/python" -m pip install --no-deps -e "$DPI_ROOT"

mkdir -p "$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion"
"$EXPANSION_VENV/bin/python" -m pip freeze \
  > "$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/pip_freeze.txt"
echo "environment ready: $EXPANSION_VENV"
echo "DeepPicardIteration: $(git -C "$DPI_ROOT" rev-parse HEAD)"
