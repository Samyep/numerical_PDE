#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/single_point_32_ice
for manifest in "$RESULTS_DIR"/dispatch_bundles/stage_c_d{20,50,100,200,400}.json; do
  [[ -f "$manifest" ]] || {
    echo "missing amended bundle manifest: $manifest" >&2
    exit 3
  }
done

active=$(squeue -h -u "$USER" -o '%i|%j|%T' | grep -E '\|sp32_(v2|audit_v2|launch_v2|de_v2)' || true)
[[ -z "$active" ]] || {
  echo "refusing a duplicate amended pipeline launch:" >&2
  echo "$active" >&2
  exit 6
}

start_dimension=${SP32_START_DIMENSION:-20}
case "$start_dimension" in
  20|50|100|200|400) ;;
  *) echo "invalid start dimension: $start_dimension" >&2; exit 2 ;;
esac
code_sha=$(git -C "$SOURCE_ROOT" rev-parse HEAD)
launcher=$(sbatch --parsable \
  --job-name="sp32_launch_v2_d${start_dimension}" \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,SP32_DIMENSION=$start_dimension,SP32_CODE_SHA=$code_sha" \
  "$SCRIPT_DIR/stage_c_bundle_launcher.slurm")
launcher=${launcher%%;*}
python -m invariant_region_mlp.experiments.single_point_32_ice.amended_dispatch \
  record-pipeline --job-id "$launcher" --dimension "$start_dimension"
echo "amended_stage_c_launcher=$launcher code_sha=$code_sha"
