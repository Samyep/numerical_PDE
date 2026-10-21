#!/bin/bash
# Login-node controller.  It waits for the original seven-PDE suite to become
# fully idle, then dispatches the expansion in the preregistered stage order.
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/icml_two_highdim_expansion
OLD_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULTS_DIR=$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion
QUEUE_DIR=$RESULTS_DIR/dispatch_queues
STAGE_A=$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/stage_a.slurm
WAVE_SCRIPT=$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/four_gpu_wave.slurm
WORKERS=4

source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$QUEUE_DIR" "$RESULTS_DIR/slurm"

timestamp() { date -u +%FT%TZ; }

old_controller_active() {
  ps -u "$USER" -o args= | grep -E \
    'invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/auto_.*controller\.sh' \
    | grep -v grep >/dev/null
}

old_job_active() {
  squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'
}

wait_for_old_suite() {
  while old_controller_active || old_job_active; do
    echo "$(timestamp) waiting_for_original_suite jobs=$(squeue -h -u "$USER" -o '%i:%j:%T' | grep ':icml_whole' | paste -sd, - || true)"
    sleep 60
  done
  echo "$(timestamp) original_suite_idle"
  if [[ -f "$OLD_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_AUDIT.json" ]]; then
    (cd "$OLD_ROOT" && python - <<'PY'
import json
from pathlib import Path
p=Path("invariant_region_mlp/results/tangent_actor_icml_final_suite/WHOLE_FUNCTION_AUDIT.json")
d=json.loads(p.read_text())
print({"old_suite_audit": d.get("status"), "status_counts": d.get("status_counts"), "gpu_hours": d.get("new_gpu_hours")})
if d.get("status") != "PASS":
    raise SystemExit("old suite artifact audit is not PASS")
PY
    )
  fi
}

wait_for_job() {
  local job=$1 state=""
  while squeue -h -j "$job" | grep -q .; do
    echo "$(timestamp) waiting job=$job $(squeue -h -j "$job" -o '%T|%M|%S|%R')"
    sleep 45
  done
  for _ in $(seq 1 15); do
    state=$(sacct -X -n -P -j "$job" -o State | sed '/^$/d' | head -n1 | cut -d'|' -f1)
    [[ -n "$state" ]] && break
    sleep 4
  done
  echo "$(timestamp) terminal job=$job state=${state:-UNKNOWN}"
  case "$state" in COMPLETED*|FAILED*|TIMEOUT*|OUT_OF_MEMORY*|CANCELLED*) ;; *) return 3;; esac
  LAST_JOB_STATE=$state
}

no_expansion_job_active() {
  ! squeue -h -u "$USER" -o '%j' | grep -q '^icml_exp_'
}

refresh() {
  python -m invariant_region_mlp.experiments.main_highdim_expansion.artifacts
  python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting
}

stage_runnable_count() {
  local stage=$1
  python - "$stage" <<'PY'
import sys
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import DIMENSIONS, METHODS, PDES
from invariant_region_mlp.experiments.main_highdim_expansion.queue import STAGES
from invariant_region_mlp.experiments.main_highdim_expansion.status import cell_budget
stage=sys.argv[1]
print(sum(
    STAGES[stage](m,p,d) and cell_budget(m,p,d)["status"] in {"PENDING","CHECKPOINTED"}
    for p in PDES for d in DIMENSIONS for m in METHODS
))
PY
}

submit_stage_a() {
  no_expansion_job_active || { echo "duplicate expansion job active" >&2; return 6; }
  sbatch --test-only "$STAGE_A" >&2
  local job
  job=$(sbatch --parsable "$STAGE_A")
  echo "$job" > "$RESULTS_DIR/STAGE_A_JOB_ID.txt"
  echo "$(timestamp) submitted stage_a job=$job"
  wait_for_job "$job"
  [[ "$LAST_JOB_STATE" == COMPLETED* ]] || return 7
  [[ -f "$RESULTS_DIR/stage_a_complete.json" ]] || return 8
  python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
  refresh
}

submit_wave() {
  local stage=$1 wave=$2 partition=$3 duration=$4
  no_expansion_job_active || { echo "duplicate expansion job active" >&2; return 6; }
  local task_file="$QUEUE_DIR/${stage}_wave${wave}_$(date -u +%Y%m%dT%H%M%SZ).csv"
  python -m invariant_region_mlp.experiments.main_highdim_expansion.queue \
    --stage "$stage" --slots "$WORKERS" --output "$task_file" >&2
  local task_sha export_vars job_name
  task_sha=$(sha256sum "$task_file" | awk '{print $1}')
  export_vars="ALL,TASK_FILE=$task_file,TASK_FILE_SHA256=$task_sha,WORKERS=$WORKERS"
  job_name="icml_exp_${stage}w${wave}"
  local resources=(--partition="$partition" --nodes=2 --ntasks=4 --ntasks-per-node=2
    --gpus-per-node=2 --cpus-per-task=12 --mem=220G --time="$duration"
    --job-name="$job_name" --export="$export_vars")
  sbatch --test-only "${resources[@]}" "$WAVE_SCRIPT" >&2
  local job
  job=$(sbatch --parsable "${resources[@]}" "$WAVE_SCRIPT")
  echo "$job" > "$RESULTS_DIR/${stage^^}_WAVE${wave}_JOB_ID.txt"
  echo "$(timestamp) submitted stage=$stage wave=$wave job=$job partition=$partition duration=$duration sha=$task_sha" >&2
  printf '%s\n' "$job"
}

run_stage_waves() {
  local stage=$1 partition=$2 duration=$3 max_waves=$4
  local wave before after job
  for ((wave=1; wave<=max_waves; wave++)); do
    before=$(stage_runnable_count "$stage")
    if (( before == 0 )); then
      echo "$(timestamp) stage=$stage has_no_runnable_cells"
      return 0
    fi
    job=$(submit_wave "$stage" "$wave" "$partition" "$duration") || {
      echo "$(timestamp) stage=$stage submission_stopped wave=$wave" >&2
      return 0
    }
    wait_for_job "$job"
    refresh
    after=$(stage_runnable_count "$stage")
    echo "$(timestamp) stage=$stage wave=$wave runnable_before=$before runnable_after=$after state=$LAST_JOB_STATE"
    if (( after >= before )) && [[ "$LAST_JOB_STATE" != COMPLETED* ]]; then
      echo "$(timestamp) stage=$stage made_no_progress_after_noncompleted_job; stopping stage" >&2
      return 0
    fi
  done
}

finalize() {
  refresh
  pytest -q invariant_region_mlp/experiments/main_highdim_expansion
  git add .gitignore invariant_region_mlp/experiments/main_highdim_expansion \
    invariant_region_mlp/results/main_highdim_expansion \
    invariant_region_mlp/docs/ICML_MAIN_HIGHDIM_EXPANSION_REPORT.md
  if ! git diff --cached --quiet; then
    git commit -m "Record B6 and Dense Fourier high-dimensional benchmark results"
  fi
  git push -u origin ir-mlp-icml-two-highdim-expansion
  echo "$(timestamp) expansion_controller_finished"
}

echo "$(timestamp) expansion_controller_start pid=$$ git=$(git rev-parse HEAD)"
wait_for_old_suite
submit_stage_a

# Stage B: all four methods at the 100D formal anchor.  A 6h10 regular
# allocation accommodates source methods that cannot resume mid-training while
# preserving the cumulative six-hour cell cap and a ten-minute teardown tail.
run_stage_waves anchor ghx4 06:10:00 2

# Stage C: the eight remaining Tangent-Actor dimension cells.  These runners
# checkpoint at Actor-update/evaluation boundaries and can use fast backfill.
run_stage_waves actor_scaling ghx4-interactive 02:00:00 8

# Stage D: fixed source-faithful baselines.  Long SDGD cells checkpoint every
# 100 epochs; DPI/PINN start only when their frozen runtime high bound fits.
run_stage_waves baseline_scaling ghx4 06:10:00 12

finalize
