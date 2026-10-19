#!/bin/bash
# Submit the approved Stage 2 and Stage 3 packs after a completed Stage-1 job.
# This controller runs on the login node and never performs GPU computation.
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULTS_DIR="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite"
PACK_SCRIPT="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_four_gpu_pack.slurm"
QUEUE_DIR="$RESULTS_DIR/whole_function_dispatch_queues"
PREDECESSOR_JOB=${1:?usage: auto_whole_function_controller.sh STAGE1_D50_JOB_ID}

source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"

timestamp() { date -u +%FT%TZ; }

wait_for_job() {
  local job=$1 state=""
  while squeue -h -j "$job" | grep -q .; do
    echo "$(timestamp) waiting job=$job $(squeue -h -j "$job" -o '%T|%M|%S|%R')"
    sleep 30
  done
  for _ in $(seq 1 10); do
    state=$(sacct -X -n -j "$job" -o State | awk 'NF {print $1; exit}')
    [[ -n "$state" ]] && break
    sleep 10
  done
  echo "$(timestamp) finished job=$job state=${state:-UNKNOWN}"
  case "$state" in
    COMPLETED*|FAILED*|TIMEOUT*) ;;
    *) echo "unexpected terminal Slurm state for job=$job: ${state:-UNKNOWN}" >&2; return 3 ;;
  esac
}

refresh_and_audit() {
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.aggregate_results
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.verify_paired_banks
  python - <<'PY'
import json
from pathlib import Path
audit = json.loads(Path(
    "invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_AUDIT.json"
).read_text())
if audit["status"] != "PASS":
    raise SystemExit("whole-function audit failed")
if not audit["checks"]["new_compute_within_120_gpu_hours"]:
    raise SystemExit("global 120 GH200-hour cap failed")
print({"status_counts": audit["status_counts"], "new_gpu_hours": audit["new_gpu_hours"]})
PY
}

require_stage1_complete() {
  python - <<'PY'
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status import METHODS, cell_budget
statuses = {
    dimension: {
        method: cell_budget(method, "joint_concave", dimension)["status"]
        for method in METHODS
    }
    for dimension in (20, 50)
}
print({"stage1_cumulative_statuses": statuses})
if any(status != "COMPLETE" for row in statuses.values() for status in row.values()):
    raise SystemExit("Stage 1 did not complete all eight pilot cells")
PY
}

remaining_gpu_seconds() {
  python - <<'PY'
import json
from pathlib import Path
p = Path("invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_JOB_ACCOUNTING.json")
d = json.loads(p.read_text())
print(max(0, int(120 * 3600 - float(d["allocated_gpu_seconds"]))))
PY
}

format_duration() {
  python - "$1" <<'PY'
import sys
seconds = int(sys.argv[1])
hours, rest = divmod(seconds, 3600)
minutes, seconds = divmod(rest, 60)
print(f"{hours:02d}:{minutes:02d}:{seconds:02d}")
PY
}

submit_stage() {
  local stage=$1 desired_wall_seconds=$2
  if squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'; then
    echo "duplicate protection: another icml_whole job is active" >&2
    return 6
  fi
  mkdir -p "$QUEUE_DIR"
  local task_file="$QUEUE_DIR/${stage}_$(date -u +%Y%m%dT%H%M%SZ).csv"
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_queue \
    --stage "$stage" --slots 4 --output "$task_file" >&2
  local task_sha
  task_sha=$(sha256sum "$task_file" | awk '{print $1}')
  local remaining max_wall wall duration
  remaining=$(remaining_gpu_seconds)
  max_wall=$((remaining / 4))
  wall=$desired_wall_seconds
  (( wall > max_wall )) && wall=$max_wall
  # Keep a small unreserved tail rather than risk a rounding race at 120 hours.
  (( wall > 60 )) && wall=$((wall - 60))
  if (( wall < 1800 )); then
    echo "insufficient global budget for a four-GPU stage: remaining_gpu_seconds=$remaining" >&2
    return 7
  fi
  duration=$(format_duration "$wall")
  local export_vars="ALL,TASK_FILE=$task_file,TASK_FILE_SHA256=$task_sha"
  sbatch --test-only --partition=ghx4 --time="$duration" --job-name="icml_whole_${stage}" \
    --export="$export_vars" "$PACK_SCRIPT" >&2
  local job
  job=$(sbatch --parsable --partition=ghx4 --time="$duration" --job-name="icml_whole_${stage}" \
    --export="$export_vars" "$PACK_SCRIPT")
  echo "$job" > "$RESULTS_DIR/${stage^^}_JOB_ID.txt"
  echo "$(timestamp) submitted stage=$stage job=$job duration=$duration task_sha=$task_sha" >&2
  printf '%s\n' "$job"
}

wait_for_job "$PREDECESSOR_JOB"
refresh_and_audit
require_stage1_complete
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite

stage2_job=$(submit_stage stage2 57600)
wait_for_job "$stage2_job"
refresh_and_audit
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite

# The frozen estimate requires about 45.4 GPU-hours for 400D.  Twelve hours
# on four GPUs supplies 48 hours while the reservation guard clamps this to
# the exact remaining global allowance if earlier stages ran long.
stage3_job=$(submit_stage stage3 43200)
wait_for_job "$stage3_job"
refresh_and_audit
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite
echo "$(timestamp) whole-function staged execution finished stage2=$stage2_job stage3=$stage3_job"
