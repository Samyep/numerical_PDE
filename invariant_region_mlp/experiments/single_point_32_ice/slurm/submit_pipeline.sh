#!/bin/bash
# Submit Stage B and a self-advancing, chunked Stage-C dependency chain.
# Each stage uses at most four independent single-H100 writers.  Stage C is
# chunked because each Slurm array element counts against coe-ice's 50-job
# MaxSubmitJobsPerUser limit.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/single_point_32_ice
QUEUE_DIR=$RESULTS_DIR/dispatch_queues
WORKERS=${SP32_WORKERS:-4}
CHUNK_SIZE=${SP32_CHUNK_SIZE:-40}
[[ "$WORKERS" == 4 ]] || {
  echo "the frozen queues require exactly four workers" >&2
  exit 2
}
[[ "$CHUNK_SIZE" =~ ^[0-9]+$ ]] && (( CHUNK_SIZE >= 1 && CHUNK_SIZE <= 40 )) || {
  echo "SP32_CHUNK_SIZE must be an integer in [1, 40]" >&2
  exit 2
}
mkdir -p "$RESULTS_DIR/slurm"

python - <<'PY'
import json
from invariant_region_mlp.experiments.single_point_32_ice.protocol import require_stage_a
from invariant_region_mlp.experiments.single_point_32_ice.protocol import RESULTS_DIR
require_stage_a()
completion = RESULTS_DIR / "stage_a_complete.json"
if not completion.is_file():
    raise SystemExit(f"missing Stage-A completion marker: {completion}")
if json.loads(completion.read_text(encoding="utf-8")).get("status") != "COMPLETE":
    raise SystemExit("Stage-A CUDA smoke suite is not complete")
PY

active=$(squeue -h -u "$USER" -o '%i|%j|%T' | grep -E '\|sp32_(b|c|audit|launch|de)' || true)
[[ -z "$active" ]] || {
  echo "refusing a duplicate single-point pipeline launch:" >&2
  echo "$active" >&2
  exit 6
}

for manifest in "$QUEUE_DIR"/stage_b_d100.json "$QUEUE_DIR"/stage_c_d{20,50,100,200,400}.json; do
  [[ -f "$manifest" ]] || {
    echo "missing frozen task manifest: $manifest" >&2
    exit 3
  }
done

worker_b=$SCRIPT_DIR/stage_b_worker.slurm
auditor=$SCRIPT_DIR/audit_stage.slurm
launcher=$SCRIPT_DIR/stage_c_chunk_launcher.slurm
jobs=()

manifest=$QUEUE_DIR/stage_b_d100.json
stage_b_job=$(sbatch --parsable --array=0-43%4 \
  --export="ALL,TASK_MANIFEST=$manifest" "$worker_b")
stage_b_job=${stage_b_job%%;*}
jobs+=("stage_b_workers=$stage_b_job")
prior=$(sbatch --parsable --dependency="afterok:$stage_b_job" \
  --job-name=sp32_audit_b --export="ALL,TASK_MANIFEST=$manifest" "$auditor")
prior=${prior%%;*}
jobs+=("stage_b_audit=$prior")

code_sha=$(git -C "$SOURCE_ROOT" rev-parse HEAD)
controller=$(sbatch --parsable --dependency="afterok:$prior" \
  --job-name=sp32_launch_d20_k000 \
  --export="ALL,SOURCE_ROOT=$SOURCE_ROOT,SP32_DIMENSION=20,SP32_START_INDEX=0,SP32_CHUNK_SIZE=$CHUNK_SIZE,SP32_CODE_SHA=$code_sha" \
  "$launcher")
controller=${controller%%;*}
jobs+=("stage_c_controller=$controller")

python -m invariant_region_mlp.experiments.single_point_32_ice.orchestration \
  record-pipeline "${jobs[@]}"
printf '%s\n' "${jobs[@]}"
