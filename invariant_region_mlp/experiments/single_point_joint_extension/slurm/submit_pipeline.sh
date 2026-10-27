#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"
mkdir -p "$SOURCE_ROOT/invariant_region_mlp/results/single_point_joint_extension/slurm"

active=$(squeue -h -u "$USER" -o '%i|%j|%T' | grep -E '\|joint_(gate|formal|final)' || true)
if [[ -n "$active" ]]; then
  echo "refusing duplicate joint pipeline:" >&2
  echo "$active" >&2
  exit 6
fi

gate=$(sbatch --parsable --export="ALL,SOURCE_ROOT=$SOURCE_ROOT" "$SCRIPT_DIR/preflight_smoke.slurm")
gate=${gate%%;*}
formal=$(sbatch --parsable --dependency="afterok:$gate" --array="0-24%4" \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT" "$SCRIPT_DIR/formal_worker.slurm")
formal=${formal%%;*}
finalize=$(sbatch --parsable --dependency="afterok:$formal" \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT" "$SCRIPT_DIR/finalize.slurm")
finalize=${finalize%%;*}
echo "gate=$gate formal=$formal finalize=$finalize"
