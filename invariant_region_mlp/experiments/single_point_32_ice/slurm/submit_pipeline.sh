#!/bin/bash
# Submit Stage B and all five Stage-C dimensions as a strict dependency chain.
# Each stage uses four independent single-H100 writers.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/single_point_32_ice
QUEUE_DIR=$RESULTS_DIR/dispatch_queues
WORKERS=${SP32_WORKERS:-4}
[[ "$WORKERS" == 4 ]] || {
  echo "the frozen queues require exactly four workers" >&2
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

active=$(squeue -h -u "$USER" -o '%i|%j|%T' | grep -E '\|sp32_(b|c|audit|de)' || true)
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
worker_c=$SCRIPT_DIR/stage_c_worker.slurm
auditor=$SCRIPT_DIR/audit_stage.slurm
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

for dimension in 20 50 100 200 400; do
  manifest=$QUEUE_DIR/stage_c_d${dimension}.json
  workers=$(sbatch --parsable --dependency="afterok:$prior" --array=0-635%4 \
    --job-name="sp32_c_d${dimension}" \
    --export="ALL,TASK_MANIFEST=$manifest" "$worker_c")
  workers=${workers%%;*}
  jobs+=("stage_c_d${dimension}_workers=$workers")
  prior=$(sbatch --parsable --dependency="afterok:$workers" \
    --job-name="sp32_audit_d${dimension}" \
    --export="ALL,TASK_MANIFEST=$manifest" "$auditor")
  prior=${prior%%;*}
  jobs+=("stage_c_d${dimension}_audit=$prior")
done

final_job=$(sbatch --parsable --dependency="afterok:$prior" \
  --job-name=sp32_de "$SCRIPT_DIR/stage_de.slurm")
final_job=${final_job%%;*}
jobs+=("stage_d_e=$final_job")

python -m invariant_region_mlp.experiments.single_point_32_ice.orchestration \
  record-pipeline "${jobs[@]}"
printf '%s\n' "${jobs[@]}"
