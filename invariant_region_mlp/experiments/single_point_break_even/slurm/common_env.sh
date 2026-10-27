#!/bin/bash
set -euo pipefail

BREAK_EVEN_SLURM_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=${SOURCE_ROOT:-$(cd "$BREAK_EVEN_SLURM_DIR/../../../.." && pwd)}
STORAGE_PARENT=$(cd "$SOURCE_ROOT/.." && pwd)
ORIGINAL_ICE_REPO=${ORIGINAL_ICE_REPO:-$STORAGE_PARENT/numerical_PDE}
EXPANSION_VENV=${EXPANSION_VENV:-$ORIGINAL_ICE_REPO/.venv-main-highdim}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-$STORAGE_PARENT/external}
export SOURCE_ROOT ORIGINAL_ICE_REPO EXPANSION_VENV EXTERNAL_ROOT

source "$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/common_env.sh"

