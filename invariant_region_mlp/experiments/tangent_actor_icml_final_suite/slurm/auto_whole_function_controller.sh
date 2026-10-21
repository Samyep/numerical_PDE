#!/bin/bash
# Submit the approved Stage 2 and Stage 3 packs after a completed Stage-1 job.
# This controller runs on the login node and never performs GPU computation.
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULTS_DIR="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite"
PACK_SCRIPT="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_four_gpu_wave.slurm"
QUEUE_DIR="$RESULTS_DIR/whole_function_dispatch_queues"
PREDECESSOR_JOB=${1:?usage: auto_whole_function_controller.sh STAGE1_D50_JOB_ID}
WORKERS=4

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

pending_count() {
  local stage=$1
  python - "$stage" <<'PY'
import sys
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts import build_results_frame
stage = sys.argv[1]
dims = {"stage2": {20, 50, 200}, "stage3": {400}}[stage]
frame, _ = build_results_frame()
print(int((frame["dimension"].isin(dims) & (frame["status"] == "PENDING")).sum()))
PY
}

submit_wave() {
  local stage=$1 wave=$2
  if squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'; then
    echo "duplicate protection: another icml_whole job is active" >&2
    return 6
  fi
  mkdir -p "$QUEUE_DIR"
  local task_file="$QUEUE_DIR/${stage}_wave${wave}_$(date -u +%Y%m%dT%H%M%SZ).csv"
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_queue \
    --stage "$stage" --slots "$WORKERS" --output "$task_file" >&2
  local task_sha
  task_sha=$(sha256sum "$task_file" | awk '{print $1}')
  local remaining max_wall wall duration
  remaining=$(remaining_gpu_seconds)
  max_wall=$((remaining / WORKERS))
  # Give every dispatched cell its full six-hour scientific cap, plus ten
  # minutes for startup and teardown.  This avoids repeatedly queueing short
  # backfill waves and, critically, lets nonresumable DPI/PINN cells satisfy
  # their frozen high-runtime guard.  submit_wave subtracts a 60-second
  # suite-budget tail below, so start one minute above 06:10:00.
  wall=$((6 * 3600 + 11 * 60))
  (( wall > max_wall )) && wall=$max_wall
  # Keep a small unreserved tail rather than risk a rounding race at 120 hours.
  (( wall > 60 )) && wall=$((wall - 60))
  if (( wall < 900 )); then
    echo "insufficient global budget for a four-GPU wave: remaining_gpu_seconds=$remaining" >&2
    return 7
  fi
  duration=$(format_duration "$wall")
  local export_vars="ALL,TASK_FILE=$task_file,TASK_FILE_SHA256=$task_sha,WORKERS=$WORKERS"
  local job_name="icml_whole_${stage}longw${wave}"
  local resources=(--partition=ghx4 --nodes=2 --ntasks=4 --ntasks-per-node=2
    --gpus-per-node=2 --cpus-per-task=12 --mem=220G --time="$duration")
  sbatch --test-only "${resources[@]}" --job-name="$job_name" \
    --export="$export_vars" "$PACK_SCRIPT" >&2
  local job
  job=$(sbatch --parsable "${resources[@]}" --job-name="$job_name" \
    --export="$export_vars" "$PACK_SCRIPT")
  echo "$job" > "$RESULTS_DIR/${stage^^}_WAVE${wave}_JOB_ID.txt"
  echo "$(timestamp) submitted stage=$stage wave=$wave job=$job duration=$duration task_sha=$task_sha" >&2
  printf '%s\n' "$job"
}

next_wave_number() {
  local stage=$1 path base number maximum=0
  shopt -s nullglob
  for path in "$RESULTS_DIR/${stage^^}_WAVE"*_JOB_ID.txt; do
    base=$(basename "$path")
    if [[ "$base" =~ ^${stage^^}_WAVE([0-9]+)_JOB_ID\.txt$ ]]; then
      number=${BASH_REMATCH[1]}
      (( number > maximum )) && maximum=$number
    fi
  done
  shopt -u nullglob
  echo $((maximum + 1))
}

run_waves() {
  local stage=$1 max_waves=$2
  local wave first iteration before after job
  first=$(next_wave_number "$stage")
  for ((iteration=0; iteration<max_waves; iteration++)); do
    wave=$((first + iteration))
    before=$(pending_count "$stage")
    if (( before == 0 )); then
      echo "$(timestamp) stage=$stage complete before wave=$wave"
      break
    fi
    if ! job=$(submit_wave "$stage" "$wave"); then
      echo "$(timestamp) stage=$stage stopped before wave=$wave due to queue/budget guard" >&2
      break
    fi
    wait_for_job "$job"
    refresh_and_audit
    after=$(pending_count "$stage")
    echo "$(timestamp) stage=$stage wave=$wave pending_before=$before pending_after=$after job=$job"
  done
}

wait_for_job "$PREDECESSOR_JOB"
refresh_and_audit
require_stage1_complete
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite

run_waves stage2 16
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite

run_waves stage3 12
refresh_and_audit
pytest -q invariant_region_mlp/experiments/tangent_actor_icml_final_suite
echo "$(timestamp) whole-function staged execution finished"
