#!/bin/bash
set -euo pipefail

HIGH_ACC_SLURM_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=${SOURCE_ROOT:-$(cd "$HIGH_ACC_SLURM_DIR/../../../.." && pwd)}
STORAGE_PARENT=$(cd "$SOURCE_ROOT/.." && pwd)
ORIGINAL_ICE_REPO=${ORIGINAL_ICE_REPO:-$STORAGE_PARENT/numerical_PDE}
EXPANSION_VENV=${EXPANSION_VENV:-$ORIGINAL_ICE_REPO/.venv-main-highdim}
EXTERNAL_ROOT=${EXTERNAL_ROOT:-$STORAGE_PARENT/external}
IRMLP_JOINT_SOURCE_RESULTS=${IRMLP_JOINT_SOURCE_RESULTS:-$STORAGE_PARENT/ir-mlp-single-point-break-even-worktree/invariant_region_mlp/results/single_point_break_even_joint}
export SOURCE_ROOT ORIGINAL_ICE_REPO EXPANSION_VENV EXTERNAL_ROOT IRMLP_JOINT_SOURCE_RESULTS

source "$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/common_env.sh"

