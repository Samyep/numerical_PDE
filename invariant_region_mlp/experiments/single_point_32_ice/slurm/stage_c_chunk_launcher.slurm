#!/bin/bash
#SBATCH --job-name=sp32_launch
#SBATCH --account=isye
#SBATCH --partition=ice-cpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/single_point_32_ice/slurm/%x_%j.err

# Submit one bounded Stage-C array and its successor.  PACE's coe-ice QOS
# permits at most 50 submitted jobs per user, and every array element counts
# toward that limit.  A short dependency-chained launcher therefore keeps at
# most one <=40-element array in Slurm while preserving the frozen task order.
set -euo pipefail

SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the experiment checkout root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/single_point_32_ice/slurm
source "$SCRIPT_DIR/common_env.sh"

dimension=${SP32_DIMENSION:?missing SP32_DIMENSION}
start=${SP32_START_INDEX:?missing SP32_START_INDEX}
chunk_size=${SP32_CHUNK_SIZE:-40}
code_sha=${SP32_CODE_SHA:?missing SP32_CODE_SHA}

case "$dimension" in
  20|50|100|200|400) ;;
  *) echo "invalid Stage-C dimension: $dimension" >&2; exit 2 ;;
esac
[[ "$start" =~ ^[0-9]+$ ]] || {
  echo "invalid chunk start: $start" >&2
  exit 2
}
[[ "$chunk_size" =~ ^[0-9]+$ ]] && (( chunk_size >= 1 && chunk_size <= 40 )) || {
  echo "chunk size must be an integer in [1, 40]" >&2
  exit 2
}

actual_sha=$(git -C "$SOURCE_ROOT" rev-parse HEAD)
[[ "$actual_sha" == "$code_sha" ]] || {
  echo "checkout moved during the formal run: expected $code_sha, found $actual_sha" >&2
  exit 3
}

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/single_point_32_ice
manifest=$RESULTS_DIR/dispatch_queues/stage_c_d${dimension}.json
[[ -f "$manifest" ]] || {
  echo "missing frozen task manifest: $manifest" >&2
  exit 4
}

read -r task_count manifest_stage manifest_dimension < <(
  python - "$manifest" <<'PY'
import json
import sys

payload = json.load(open(sys.argv[1], encoding="utf-8"))
print(payload["task_count"], payload["stage"], payload["dimension"])
PY
)
[[ "$manifest_stage" == C && "$manifest_dimension" == "$dimension" ]] || {
  echo "manifest metadata does not match Stage C d=$dimension" >&2
  exit 4
}
(( start < task_count )) || {
  echo "chunk start $start is outside task count $task_count" >&2
  exit 4
}

end=$((start + chunk_size - 1))
(( end < task_count )) || end=$((task_count - 1))
chunk_index=$((start / chunk_size))
chunk_tag=$(printf '%03d' "$chunk_index")
worker=$SCRIPT_DIR/stage_c_worker.slurm
auditor=$SCRIPT_DIR/audit_stage.slurm
launcher=$SCRIPT_DIR/stage_c_chunk_launcher.slurm

workers=$(sbatch --parsable --array="${start}-${end}%4" \
  --job-name="sp32_c_d${dimension}" \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,TASK_MANIFEST=$manifest" "$worker")
workers=${workers%%;*}

next_kind=launcher
next_job=
if (( end + 1 < task_count )); then
  next_start=$((end + 1))
  next_job=$(sbatch --parsable --dependency="afterok:$workers" \
    --job-name="sp32_launch_d${dimension}_k$(printf '%03d' $((chunk_index + 1)))" \
    --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,SP32_DIMENSION=$dimension,SP32_START_INDEX=$next_start,SP32_CHUNK_SIZE=$chunk_size,SP32_CODE_SHA=$code_sha" \
    "$launcher")
  next_job=${next_job%%;*}
else
  audit_job=$(sbatch --parsable --dependency="afterok:$workers" \
    --job-name="sp32_audit_d${dimension}" \
    --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,TASK_MANIFEST=$manifest" "$auditor")
  audit_job=${audit_job%%;*}
  next_kind=audit
  next_job=$audit_job

  case "$dimension" in
    20) next_dimension=50 ;;
    50) next_dimension=100 ;;
    100) next_dimension=200 ;;
    200) next_dimension=400 ;;
    400) next_dimension= ;;
  esac
  if [[ -n "$next_dimension" ]]; then
    continuation=$(sbatch --parsable --dependency="afterok:$audit_job" \
      --job-name="sp32_launch_d${next_dimension}_k000" \
      --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,SP32_DIMENSION=$next_dimension,SP32_START_INDEX=0,SP32_CHUNK_SIZE=$chunk_size,SP32_CODE_SHA=$code_sha" \
      "$launcher")
    continuation=${continuation%%;*}
  else
    continuation=$(sbatch --parsable --dependency="afterok:$audit_job" \
      --job-name=sp32_de --export="ALL,SOURCE_ROOT=$SOURCE_ROOT" \
      "$SCRIPT_DIR/stage_de.slurm")
    continuation=${continuation%%;*}
  fi
fi

mkdir -p "$RESULTS_DIR/scheduler_events"
event_file=$RESULTS_DIR/scheduler_events/d${dimension}_chunk_${chunk_tag}.json
python - "$event_file" "$dimension" "$start" "$end" "$task_count" \
  "$workers" "$next_kind" "$next_job" "${audit_job:-}" "${continuation:-}" \
  "$code_sha" <<'PY'
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

(
    output, dimension, start, end, task_count, workers, next_kind,
    next_job, audit_job, continuation, code_sha,
) = sys.argv[1:]
payload = {
    "revision": "single_point_32_ice_v1",
    "recorded_utc": datetime.now(timezone.utc).isoformat(),
    "launcher_job_id": os.environ.get("SLURM_JOB_ID"),
    "dimension": int(dimension),
    "start_index": int(start),
    "end_index": int(end),
    "task_count": int(task_count),
    "worker_array_job_id": workers,
    "next_kind": next_kind,
    "next_job_id": next_job,
    "audit_job_id": audit_job or None,
    "continuation_job_id": continuation or None,
    "code_sha": code_sha,
}
path = Path(output)
temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
os.replace(temporary, path)
PY

echo "dimension=$dimension chunk=$chunk_tag indices=$start-$end workers=$workers next=$next_job continuation=${continuation:-none}"
