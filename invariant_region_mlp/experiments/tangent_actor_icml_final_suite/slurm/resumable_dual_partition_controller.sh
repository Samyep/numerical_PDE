#!/bin/bash
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULTS_DIR="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite"
QUEUE_DIR="$RESULTS_DIR/whole_function_dispatch_queues"
PACK_SCRIPT="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_four_gpu_wave.slurm"
REGULAR_PREFERENCE_SECONDS=7200
WAVE=1

source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"

timestamp() { date -Is; }

refresh_and_audit() {
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.aggregate_results
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.verify_paired_banks
}

pending_count() {
  python - <<'PY'
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts import build_results_frame
frame, _ = build_results_frame()
print(int((frame["status"] == "PENDING").sum()))
PY
}

queue_state() {
  squeue -h -j "$1" -o '%T' 2>/dev/null | head -1
}

regular_due_within_window() {
  local job=$1 estimate estimate_epoch now
  estimate=$(squeue -h -j "$job" --start -o '%S' 2>/dev/null | head -1)
  [[ -n "$estimate" && "$estimate" != "N/A" && "$estimate" != "Unknown" ]] || return 1
  estimate_epoch=$(date -d "$estimate" +%s 2>/dev/null) || return 1
  now=$(date +%s)
  (( estimate_epoch - now <= REGULAR_PREFERENCE_SECONDS ))
}

regular_walltime() {
  local gpus=$1
  python - "$gpus" <<'PY'
import json
import sys
from pathlib import Path

gpus = int(sys.argv[1])
payload = json.loads(Path(
    "invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_JOB_ACCOUNTING.json"
).read_text())
remaining = 120 * 3600 - float(payload["allocated_gpu_seconds"])
seconds = min(6 * 3600, int(remaining // gpus) - 120)
seconds = max(30 * 60, seconds)
seconds -= seconds % 60
hours, remainder = divmod(seconds, 3600)
minutes = remainder // 60
print(f"{hours:02d}:{minutes:02d}:00")
PY
}

submit_pair() {
  local pending workers nodes tasks_per_node gpus_per_node requested_gpus
  local task_sha regular_time export_vars int_name reg_name
  pending=$(pending_count)
  (( pending > 0 )) || return 3
  workers=$pending
  (( workers > 4 )) && workers=4
  if (( workers <= 2 )); then
    nodes=1
    tasks_per_node=$workers
    gpus_per_node=$workers
  else
    workers=4
    nodes=2
    tasks_per_node=2
    gpus_per_node=2
  fi
  requested_gpus=$((nodes * gpus_per_node))
  TASK_FILE="$QUEUE_DIR/resumable_remaining_wave${WAVE}_$(date -u +%Y%m%dT%H%M%SZ).csv"
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_queue \
    --stage remaining --slots "$workers" --output "$TASK_FILE" >&2
  task_sha=$(sha256sum "$TASK_FILE" | awk '{print $1}')
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting >&2
  regular_time=$(regular_walltime "$requested_gpus")
  export_vars="ALL,TASK_FILE=$TASK_FILE,TASK_FILE_SHA256=$task_sha,WORKERS=$workers"
  int_name="icml_whole_rem_intw${WAVE}"
  reg_name="icml_whole_rem_regw${WAVE}"

  if squeue -h -u "$USER" -o '%j' | grep -Eq "^(${int_name}|${reg_name})$"; then
    echo "$(timestamp) DUPLICATE_GUARD wave=$WAVE" >&2
    return 4
  fi

  sbatch --test-only --partition=ghx4-interactive --nodes="$nodes" --ntasks="$workers" \
    --ntasks-per-node="$tasks_per_node" --gpus-per-node="$gpus_per_node" \
    --cpus-per-task=12 --mem=220G --time=02:00:00 --job-name="$int_name" \
    --export="$export_vars" "$PACK_SCRIPT" >&2
  sbatch --test-only --partition=ghx4 --nodes="$nodes" --ntasks="$workers" \
    --ntasks-per-node="$tasks_per_node" --gpus-per-node="$gpus_per_node" \
    --cpus-per-task=12 --mem=220G --time="$regular_time" --job-name="$reg_name" \
    --export="$export_vars" "$PACK_SCRIPT" >&2

  REGULAR_JOB=$(sbatch --parsable --partition=ghx4 --nodes="$nodes" --ntasks="$workers" \
    --ntasks-per-node="$tasks_per_node" --gpus-per-node="$gpus_per_node" \
    --cpus-per-task=12 --mem=220G --time="$regular_time" --job-name="$reg_name" \
    --export="$export_vars" "$PACK_SCRIPT")
  INTERACTIVE_JOB=$(sbatch --parsable --partition=ghx4-interactive --nodes="$nodes" --ntasks="$workers" \
    --ntasks-per-node="$tasks_per_node" --gpus-per-node="$gpus_per_node" \
    --cpus-per-task=12 --mem=220G --time=02:00:00 --job-name="$int_name" \
    --export="$export_vars" "$PACK_SCRIPT")
  printf '%s\n' "$REGULAR_JOB" > "$RESULTS_DIR/REMAINING_REGULAR_WAVE${WAVE}_JOB_ID.txt"
  printf '%s\n' "$INTERACTIVE_JOB" > "$RESULTS_DIR/REMAINING_INTERACTIVE_WAVE${WAVE}_JOB_ID.txt"
  echo "$(timestamp) PAIR_SUBMITTED wave=$WAVE regular=$REGULAR_JOB interactive=$INTERACTIVE_JOB workers=$workers regular_time=$regular_time task_sha=$task_sha"
}

race_pair() {
  local interactive_state regular_state
  while true; do
    interactive_state=$(queue_state "$INTERACTIVE_JOB")
    regular_state=$(queue_state "$REGULAR_JOB")
    echo "$(timestamp) RACE wave=$WAVE interactive=$INTERACTIVE_JOB:${interactive_state:-TERMINAL} regular=$REGULAR_JOB:${regular_state:-TERMINAL}"

    if [[ "$interactive_state" == "RUNNING" ]]; then
      [[ -z "$regular_state" ]] || scancel "$REGULAR_JOB"
      WINNER_JOB=$INTERACTIVE_JOB
      echo "$(timestamp) RACE_WINNER wave=$WAVE job=$WINNER_JOB kind=interactive"
      return 0
    fi
    if [[ "$regular_state" == "RUNNING" ]]; then
      [[ -z "$interactive_state" ]] || scancel "$INTERACTIVE_JOB"
      WINNER_JOB=$REGULAR_JOB
      echo "$(timestamp) RACE_WINNER wave=$WAVE job=$WINNER_JOB kind=regular"
      return 0
    fi
    if [[ "$regular_state" == "PENDING" && "$interactive_state" == "PENDING" ]] \
        && regular_due_within_window "$REGULAR_JOB"; then
      scancel "$INTERACTIVE_JOB"
      WINNER_JOB=$REGULAR_JOB
      echo "$(timestamp) RACE_PREFERENCE wave=$WAVE job=$WINNER_JOB kind=regular reason=start_within_2h action=cancel_interactive"
      return 0
    fi
    if [[ -z "$interactive_state" && -n "$regular_state" ]]; then
      WINNER_JOB=$REGULAR_JOB
      return 0
    fi
    if [[ -z "$regular_state" && -n "$interactive_state" ]]; then
      WINNER_JOB=$INTERACTIVE_JOB
      return 0
    fi
    if [[ -z "$regular_state" && -z "$interactive_state" ]]; then
      echo "$(timestamp) RACE_FAILED wave=$WAVE reason=both_jobs_left_queue" >&2
      return 5
    fi
    sleep 15
  done
}

wait_for_winner() {
  local state=""
  while squeue -h -j "$WINNER_JOB" | grep -q .; do
    state=$(queue_state "$WINNER_JOB")
    echo "$(timestamp) WINNER_ACTIVE wave=$WAVE job=$WINNER_JOB state=${state:-UNKNOWN}"
    sleep 30
  done
  for _ in $(seq 1 12); do
    state=$(sacct -X -n -j "$WINNER_JOB" -o State | awk 'NF {print $1; exit}')
    [[ -n "$state" ]] && break
    sleep 10
  done
  echo "$(timestamp) WINNER_DONE wave=$WAVE job=$WINNER_JOB state=${state:-UNKNOWN}"
  [[ "$state" == COMPLETED* ]]
}

echo "$(timestamp) CONTROLLER_START pending=$(pending_count) regular_preference_seconds=$REGULAR_PREFERENCE_SECONDS pending_timeout=none"
while (( $(pending_count) > 0 )); do
  submit_pair
  race_pair
  wait_for_winner || {
    echo "$(timestamp) CONTROLLER_STOP reason=winner_not_completed job=$WINNER_JOB" >&2
    exit 6
  }
  refresh_and_audit
  WAVE=$((WAVE + 1))
done
echo "$(timestamp) CONTROLLER_DONE pending=0"
