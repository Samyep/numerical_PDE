#!/bin/bash
# Login-node controller for the frozen two-PDE expansion. Formal cells are
# separate one-GPU Slurm jobs; stages do not overlap and duplicate writers are
# rejected before submission.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"

RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion
QUEUE_DIR=$RESULTS_DIR/dispatch_queues
STAGE_A=$SCRIPT_DIR/stage_a.slurm
CELL_SCRIPT=$SCRIPT_DIR/single_cell.slurm
WORKERS=${WORKERS:-8}
CHECK_OLD_SUITE=${CHECK_OLD_SUITE:-0}
ACCOUNT=${EXPANSION_ACCOUNT:-isye}
PARTITION=${EXPANSION_PARTITION:-ice-gpu}
QOS=${EXPANSION_QOS:-coe-ice}
GPU_TYPE=${EXPANSION_GPU_TYPE:-h100}
CPUS_PER_CELL=${CPUS_PER_CELL:-8}
MEMORY_PER_CELL=${MEMORY_PER_CELL:-96G}
ANCHOR_WALL=${ANCHOR_WALL:-06:10:00}
ACTOR_WALL=${ACTOR_WALL:-02:00:00}
BASELINE_WALL=${BASELINE_WALL:-06:10:00}

[[ "$WORKERS" =~ ^[1-9][0-9]*$ ]]
mkdir -p "$QUEUE_DIR" "$RESULTS_DIR/slurm"

exec 9>"$RESULTS_DIR/controller.lock"
if ! flock -n 9; then
  echo "another expansion controller holds $RESULTS_DIR/controller.lock" >&2
  exit 5
fi
printf '%s\n' "$$" > "$RESULTS_DIR/controller.pid"

timestamp() { date -u +%FT%TZ; }

old_controller_active() {
  ps -u "$USER" -o args= | grep -E \
    'tangent_actor_icml_final_suite/slurm/auto_.*controller\.sh' \
    | grep -v grep >/dev/null
}

old_job_active() {
  squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'
}

wait_for_old_suite() {
  if [[ "$CHECK_OLD_SUITE" != 1 ]]; then
    if old_controller_active || old_job_active; then
      echo "old-suite controller/job is active on an independent-machine run" >&2
      return 4
    fi
    echo "$(timestamp) independent_machine_old_suite_check=idle"
    return 0
  fi
  while old_controller_active || old_job_active; do
    echo "$(timestamp) waiting_for_original_suite jobs=$(squeue -h -u "$USER" -o '%i:%j:%T' | grep ':icml_whole' | paste -sd, - || true)"
    sleep 60
  done
  echo "$(timestamp) original_suite_idle"
}

expansion_jobs() {
  squeue -h -u "$USER" -o '%i|%j|%T' | grep -E '\|icml_exp_' || true
}

assert_no_expansion_jobs() {
  local active
  active=$(expansion_jobs)
  [[ -z "$active" ]] || {
    echo "refusing duplicate expansion launch; active jobs:" >&2
    echo "$active" >&2
    return 6
  }
}

wait_for_jobs() {
  local ids=("$@") joined state job active
  joined=$(IFS=,; echo "${ids[*]}")
  while true; do
    active=$(squeue -h -j "$joined" -o '%i|%j|%T|%M|%S|%R' 2>/dev/null || true)
    [[ -z "$active" ]] && break
    echo "$(timestamp) waiting jobs=$joined"
    echo "$active"
    sleep 45
  done
  for job in "${ids[@]}"; do
    state=""
    for _ in $(seq 1 15); do
      state=$(sacct -X -n -P -j "$job" -o State | sed '/^$/d' | head -n1 | cut -d'|' -f1)
      [[ -n "$state" ]] && break
      sleep 4
    done
    echo "$(timestamp) terminal job=$job state=${state:-UNKNOWN}"
  done
}

refresh() {
  python -m invariant_region_mlp.experiments.main_highdim_expansion.artifacts
  python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting
}

global_remaining_seconds() {
  python - <<'PY'
from invariant_region_mlp.experiments.main_highdim_expansion.accounting import collect
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import GLOBAL_GPU_SECONDS_CAP
used = float(collect().get("allocated_gpu_seconds") or 0.0)
print(max(0, int(GLOBAL_GPU_SECONDS_CAP - used)))
PY
}

duration_seconds() {
  python - "$1" <<'PY'
import sys
parts = [int(value) for value in sys.argv[1].split(":")]
if len(parts) != 3:
    raise SystemExit("duration must be HH:MM:SS")
print(parts[0] * 3600 + parts[1] * 60 + parts[2])
PY
}

stage_runnable_count() {
  local stage=$1
  python - "$stage" <<'PY'
import sys
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import DIMENSIONS, METHODS, PDES
from invariant_region_mlp.experiments.main_highdim_expansion.queue import STAGES
from invariant_region_mlp.experiments.main_highdim_expansion.status import cell_budget
stage = sys.argv[1]
print(sum(
    STAGES[stage](method, pde, dimension)
    and cell_budget(method, pde, dimension)["status"] in {"PENDING", "CHECKPOINTED"}
    for pde in PDES for dimension in DIMENSIONS for method in METHODS
))
PY
}

mark_stage_budget_truncated() {
  local stage=$1 reason=$2 now
  now=$(timestamp)
  while IFS='|' read -r method pde dimension; do
    [[ -n "$method" ]] || continue
    python -m invariant_region_mlp.experiments.main_highdim_expansion.status record \
      --method "$method" --pde "$pde" --dimension "$dimension" \
      --status BUDGET_TRUNCATED --job-id controller --task-id "$stage" \
      --started-utc "$now" --finished-utc "$now" --actual-gpu-seconds 0 \
      --exit-code 0 --notes "$reason"
  done < <(python - "$stage" <<'PY'
import sys
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import DIMENSIONS, METHODS, PDES
from invariant_region_mlp.experiments.main_highdim_expansion.queue import STAGES
from invariant_region_mlp.experiments.main_highdim_expansion.status import cell_budget
stage = sys.argv[1]
for pde in PDES:
    for dimension in DIMENSIONS:
        for method in METHODS:
            if STAGES[stage](method, pde, dimension) and cell_budget(method, pde, dimension)["status"] in {"PENDING", "CHECKPOINTED"}:
                print(method, pde, dimension, sep="|")
PY
  )
}

submit_stage_a() {
  assert_no_expansion_jobs
  sbatch --test-only "$STAGE_A" >&2
  local job state
  job=$(sbatch --parsable "$STAGE_A")
  echo "$job" > "$RESULTS_DIR/STAGE_A_JOB_ID.txt"
  echo "$(timestamp) submitted stage_a job=$job"
  wait_for_jobs "$job"
  state=$(sacct -X -n -P -j "$job" -o State | sed '/^$/d' | head -n1 | cut -d'|' -f1)
  [[ "$state" == COMPLETED* ]] || return 7
  [[ -f "$RESULTS_DIR/stage_a_complete.json" ]] || return 8
  python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
  refresh
}

method_abbreviation() {
  case "$1" in
    tangent_actor) echo ta ;;
    dpi) echo dpi ;;
    pinn_hte) echo phte ;;
    sdgd_pinn) echo sdgd ;;
    *) return 2 ;;
  esac
}

pde_abbreviation() {
  case "$1" in
    b6_fullrank) echo b6 ;;
    dense_fourier) echo df ;;
    *) return 2 ;;
  esac
}

submit_cell() {
  local method=$1 pde=$2 dimension=$3 wall=$4 cap=$5 stage=$6
  local method_short pde_short job_name export_vars job active
  method_short=$(method_abbreviation "$method")
  pde_short=$(pde_abbreviation "$pde")
  job_name="icml_exp_${method_short}_${pde_short}_d${dimension}"
  active=$(squeue -h -u "$USER" -o '%j' | grep -Fx "$job_name" || true)
  [[ -z "$active" ]] || {
    echo "duplicate writer job already active for $method/$pde/d$dimension" >&2
    return 6
  }
  export_vars="ALL,METHOD=$method,PDE=$pde,DIMENSION=$dimension,CELL_CAP_SECONDS=$cap"
  local resources=(
    --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS"
    --nodes=1 --ntasks=1 --cpus-per-task="$CPUS_PER_CELL"
    --mem="$MEMORY_PER_CELL" --gres="gpu:$GPU_TYPE:1" --time="$wall"
    --job-name="$job_name" --export="$export_vars"
  )
  sbatch --test-only "${resources[@]}" "$CELL_SCRIPT" >&2
  job=$(sbatch --parsable "${resources[@]}" "$CELL_SCRIPT")
  echo "$(timestamp) submitted stage=$stage job=$job cell=$method/$pde/d$dimension wall=$wall cap=$cap" >&2
  printf '%s\n' "$job"
}

run_stage_waves() {
  local stage=$1 wall=$2 max_waves=$3
  local wall_seconds cap_seconds wave before after remaining capacity limit index
  local task_file row slot sequence method pde dimension estimate high_estimate extra job
  local -a rows jobs
  wall_seconds=$(duration_seconds "$wall")
  cap_seconds=$((wall_seconds - 240))
  (( cap_seconds > 21420 )) && cap_seconds=21420
  (( cap_seconds >= 300 )) || { echo "wall time too short: $wall" >&2; return 2; }

  for ((wave=1; wave<=max_waves; wave++)); do
    before=$(stage_runnable_count "$stage")
    if (( before == 0 )); then
      echo "$(timestamp) stage=$stage has_no_runnable_cells"
      return 0
    fi
    assert_no_expansion_jobs
    task_file="$QUEUE_DIR/${stage}_wave${wave}_$(date -u +%Y%m%dT%H%M%SZ).csv"
    python -m invariant_region_mlp.experiments.main_highdim_expansion.queue \
      --stage "$stage" --slots "$WORKERS" --output "$task_file" >&2
    mapfile -t rows < <(tail -n +2 "$task_file")
    remaining=$(global_remaining_seconds)
    capacity=$((remaining / wall_seconds))
    if (( capacity == 0 )); then
      echo "$(timestamp) global_budget_cannot_reserve stage=$stage remaining_seconds=$remaining wall_seconds=$wall_seconds"
      mark_stage_budget_truncated "$stage" "global 120-GPU-hour cap cannot reserve another full $wall cell allocation"
      return 0
    fi
    limit=${#rows[@]}
    (( limit > WORKERS )) && limit=$WORKERS
    (( limit > capacity )) && limit=$capacity
    jobs=()
    for ((index=0; index<limit; index++)); do
      row=${rows[$index]}
      IFS=',' read -r slot sequence method pde dimension estimate high_estimate extra <<< "$row"
      [[ -z "${extra:-}" ]]
      job=$(submit_cell "$method" "$pde" "$dimension" "$wall" "$cap_seconds" "$stage")
      jobs+=("$job")
    done
    wait_for_jobs "${jobs[@]}"
    refresh
    after=$(stage_runnable_count "$stage")
    echo "$(timestamp) stage=$stage wave=$wave runnable_before=$before runnable_after=$after jobs=$(IFS=,; echo "${jobs[*]}")"
  done

  after=$(stage_runnable_count "$stage")
  if (( after > 0 )); then
    mark_stage_budget_truncated "$stage" "controller exhausted the preregistered maximum wave count"
  fi
}

finalize() {
  local stage branch push_log
  for stage in anchor actor_scaling baseline_scaling; do
    if (( $(stage_runnable_count "$stage") > 0 )); then
      mark_stage_budget_truncated "$stage" "formal dispatch ended before another cap-compliant allocation was available"
    fi
  done
  python -m invariant_region_mlp.experiments.main_highdim_expansion.artifacts
  python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting
  python -m pytest -q invariant_region_mlp/experiments/main_highdim_expansion
  python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
  python -m invariant_region_mlp.experiments.main_highdim_expansion.artifacts
  git add .gitignore invariant_region_mlp/experiments/main_highdim_expansion \
    invariant_region_mlp/results/main_highdim_expansion \
    invariant_region_mlp/docs/ICML_MAIN_HIGHDIM_EXPANSION_REPORT.md
  if ! git diff --cached --quiet; then
    git commit -m "Record B6 and Dense Fourier high-dimensional benchmark results"
  fi
  branch=$(git branch --show-current)
  push_log=$RESULTS_DIR/slurm/final_push.log
  if ! git push -u origin "$branch" >"$push_log" 2>&1; then
    echo "$(timestamp) push_failed branch=$branch log=$push_log" >&2
    cat "$push_log" >&2
    return 9
  fi
  echo "$(timestamp) expansion_controller_finished branch=$branch git=$(git rev-parse HEAD)"
}

echo "$(timestamp) expansion_controller_start pid=$$ git=$(git rev-parse HEAD) branch=$(git branch --show-current) workers=$WORKERS"
wait_for_old_suite
assert_no_expansion_jobs
submit_stage_a

# Stage B: both PDEs x all four methods at d=100.
run_stage_waves anchor "$ANCHOR_WALL" 2

# Stage C: the eight remaining Tangent-Actor cells. The short allocations are
# exact-resume safe and allow scheduler backfill.
run_stage_waves actor_scaling "$ACTOR_WALL" 8

# Stage D: the remaining 24 fixed baseline transfers. Every DPI/PINN attempt
# receives the full allocation because those source runners cannot resume.
run_stage_waves baseline_scaling "$BASELINE_WALL" 12

finalize
