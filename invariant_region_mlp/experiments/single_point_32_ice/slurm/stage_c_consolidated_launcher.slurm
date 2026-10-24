#!/bin/bash
#SBATCH --job-name=sp32_launch_v3
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

dimension=${SP32_DIMENSION:?missing SP32_DIMENSION}
code_sha=${SP32_CODE_SHA:?missing SP32_CODE_SHA}
case "$dimension" in
  100|200|400) ;;
  *) echo "V3 accepts only dimensions 100, 200, or 400: $dimension" >&2; exit 2 ;;
esac

actual_sha=$(git -C "$SOURCE_ROOT" rev-parse HEAD)
[[ "$actual_sha" == "$code_sha" ]] || {
  echo "checkout moved during V3 run: expected $code_sha, found $actual_sha" >&2
  exit 3
}

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/single_point_32_ice
manifest=$RESULTS_DIR/dispatch_bundles_v3/stage_c_d${dimension}.json
[[ -f "$manifest" ]] || {
  echo "missing V3 bundle manifest: $manifest" >&2
  exit 4
}

pending=$(python -m invariant_region_mlp.experiments.single_point_32_ice.consolidated_dispatch \
  pending-spec --manifest "$manifest")

worker=$SCRIPT_DIR/stage_c_consolidated_worker.slurm
auditor=$SCRIPT_DIR/audit_consolidated_stage.slurm
launcher=$SCRIPT_DIR/stage_c_consolidated_launcher.slurm
exclude_nodes=${SP32_EXCLUDE_NODES:-}
if [[ -n "$pending" ]]; then
  worker_options=(
    --parsable
    --array="${pending}%4"
    --job-name="sp32_v3_d${dimension}"
    --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,BUNDLE_MANIFEST=$manifest"
  )
  if [[ -n "$exclude_nodes" ]]; then
    worker_options+=(--exclude="$exclude_nodes")
  fi
  workers=$(sbatch "${worker_options[@]}" "$worker")
  workers=${workers%%;*}
  audit_dependency="--dependency=afterok:$workers"
else
  workers=
  audit_dependency=
fi

# shellcheck disable=SC2086
audit_job=$(sbatch --parsable $audit_dependency \
  --job-name="sp32_audit_v3_d${dimension}" \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,BUNDLE_MANIFEST=$manifest" "$auditor")
audit_job=${audit_job%%;*}

case "$dimension" in
  100) next_dimension=200 ;;
  200) next_dimension=400 ;;
  400) next_dimension= ;;
esac
if [[ -n "$next_dimension" ]]; then
  continuation=$(sbatch --parsable --dependency="afterok:$audit_job" \
    --job-name="sp32_launch_v3_d${next_dimension}" \
    --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,SP32_DIMENSION=$next_dimension,SP32_CODE_SHA=$code_sha" \
    "$launcher")
else
  continuation=$(sbatch --parsable --dependency="afterok:$audit_job" \
    --job-name=sp32_de_v3 --export="ALL,SOURCE_ROOT=$SOURCE_ROOT" \
    "$SCRIPT_DIR/stage_de.slurm")
fi
continuation=${continuation%%;*}

mkdir -p "$RESULTS_DIR/scheduler_events"
event_file=$RESULTS_DIR/scheduler_events/consolidated_d${dimension}.json
python - "$event_file" "$dimension" "$pending" "${workers:-}" \
  "$audit_job" "$continuation" "$code_sha" "$exclude_nodes" <<'PY'
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

(
    output,
    dimension,
    pending,
    workers,
    audit_job,
    continuation,
    code_sha,
    exclude_nodes,
) = sys.argv[1:]
payload = {
    "protocol_revision": "single_point_32_ice_v2_ta_seed0_bundled",
    "dispatch_revision": "single_point_32_ice_dispatch_v3_fh_by_pde",
    "recorded_utc": datetime.now(timezone.utc).isoformat(),
    "launcher_job_id": os.environ.get("SLURM_JOB_ID"),
    "dimension": int(dimension),
    "pending_bundle_indices": [int(v) for v in pending.split(",") if v],
    "worker_array_job_id": workers or None,
    "audit_job_id": audit_job,
    "continuation_job_id": continuation,
    "code_sha": code_sha,
    "excluded_nodes": [value for value in exclude_nodes.split(",") if value],
}
path = Path(output)
temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, path)
PY

echo "dimension=$dimension pending=${pending:-none} workers=${workers:-none} audit=$audit_job continuation=$continuation excluded=${exclude_nodes:-none}"
