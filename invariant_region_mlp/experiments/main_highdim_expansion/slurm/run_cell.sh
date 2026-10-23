#!/bin/bash
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"
CELL_CAP_SECONDS=${CELL_CAP_SECONDS:-21420}
TEMPORARY_FAILURE_EXIT_CODE=${TEMPORARY_FAILURE_EXIT_CODE:-75}
method=${1:?usage: run_cell.sh METHOD PDE DIMENSION}
pde=${2:?usage: run_cell.sh METHOD PDE DIMENSION}
dimension=${3:?usage: run_cell.sh METHOD PDE DIMENSION}

case "$method" in tangent_actor|dpi|pinn_hte|sdgd_pinn) ;; *) exit 2;; esac
case "$pde" in b6_fullrank|dense_fourier) ;; *) exit 2;; esac
case "$dimension" in 20|50|100|200|400) ;; *) exit 2;; esac

job_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
task_id=${EXPANSION_CELL_TASK_ID:-${SLURM_ARRAY_TASK_ID:-0}}
tag="${method}__${pde}__d${dimension}"
cell_out=${CELL_OUT:-$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${job_id}_${task_id}_${tag}.out}
cell_err=${CELL_ERR:-$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${job_id}_${task_id}_${tag}.err}
mkdir -p "$(dirname "$cell_out")" "$(dirname "$cell_err")"
started_utc=$(date -u +%FT%TZ)
started_epoch=$(date +%s)

read -r prior_status remaining_seconds prior_seconds prior_attempts < <(
  python -m invariant_region_mlp.experiments.main_highdim_expansion.status budget \
    --method "$method" --pde "$pde" --dimension "$dimension"
)
case "$prior_status" in
  COMPLETE|OUT_OF_TIME|OUT_OF_MEMORY|FAILED|BUDGET_TRUNCATED)
    echo "SKIP cell=$tag status=$prior_status"
    exit 0 ;;
esac
if (( remaining_seconds <= 0 )); then
  echo "SKIP cell=$tag cumulative_six_hour_cap" >&2
  exit 3
fi
(( CELL_CAP_SECONDS > remaining_seconds )) && CELL_CAP_SECONDS=$remaining_seconds

record_status() {
  python -m invariant_region_mlp.experiments.main_highdim_expansion.status record \
    --method "$method" --pde "$pde" --dimension "$dimension" --status "$1" \
    --job-id "$job_id" --task-id "$task_id" --started-utc "$started_utc" \
    --finished-utc "$2" --actual-gpu-seconds "$3" --exit-code "$4" \
    --stdout "$cell_out" --stderr "$cell_err" --notes "$5"
}

python -m invariant_region_mlp.experiments.main_highdim_expansion.infrastructure \
  probe --context "cell:$tag" 2> >(tee -a "$cell_err" >&2)
probe_code=$?
if [[ $probe_code -ne 0 ]]; then
  finished_utc=$(date -u +%FT%TZ)
  seconds=$(( $(date +%s) - started_epoch ))
  if [[ $probe_code -eq $TEMPORARY_FAILURE_EXIT_CODE ]]; then
    status=CHECKPOINTED
    notes="transient ICE accelerator preflight failure; allocation charged and retry allowed"
  else
    status=FAILED
    notes="non-transient accelerator/environment preflight failure"
  fi
  record_status "$status" "$finished_utc" "$seconds" "$probe_code" "$notes"
  echo "END cell=$tag status=$status exit=$probe_code gpu_seconds=$seconds utc=$finished_utc"
  exit "$probe_code"
fi

echo "BEGIN cell=$tag utc=$started_utc git=$(git rev-parse HEAD) host=$(hostname) cap=$CELL_CAP_SECONDS"
python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
gate=$?
if [[ $gate -ne 0 ]]; then
  record_status FAILED "$(date -u +%FT%TZ)" 0 "$gate" "review gate rejected formal execution"
  exit "$gate"
fi
record_status RUNNING "" 0 0 "fixed seed 0; prior_gpu_seconds=$prior_seconds"

timeout --signal=TERM --kill-after=120 "$CELL_CAP_SECONDS" \
  python -m invariant_region_mlp.experiments.main_highdim_expansion.run_cell \
    --method "$method" --pde "$pde" --dimension "$dimension" \
    >"$cell_out" 2>"$cell_err"
code=$?
job_exit_code=$code
finished_utc=$(date -u +%FT%TZ)
seconds=$(( $(date +%s) - started_epoch ))
if [[ $code -eq 0 ]]; then
  status=COMPLETE; notes="runner completed under frozen protocol"
elif [[ $code -eq 124 ]]; then
  if (( prior_seconds + seconds >= 21420 )); then
    status=OUT_OF_TIME; notes="cumulative six-GPU-hour cap reached"
  else
    status=CHECKPOINTED; notes="allocation boundary; TA/SDGD resume exactly, DPI/PINN repeated work remains charged"
  fi
elif [[ $code -eq 137 ]] || grep -Eiq 'CUDA out of memory|OutOfMemoryError|oom-kill|oom_kill' "$cell_err"; then
  status=OUT_OF_MEMORY; notes="memory failure detected"
elif python -m invariant_region_mlp.experiments.main_highdim_expansion.infrastructure \
  classify --stderr "$cell_err" --context "cell:$tag" --record; then
  status=CHECKPOINTED
  notes="transient ICE accelerator failure; allocation charged and retry allowed"
  job_exit_code=$TEMPORARY_FAILURE_EXIT_CODE
else
  status=FAILED; notes="runner failed; inspect immutable cell logs"
fi
record_status "$status" "$finished_utc" "$seconds" "$code" "$notes"
echo "END cell=$tag status=$status exit=$code job_exit=$job_exit_code gpu_seconds=$seconds utc=$finished_utc"
exit "$job_exit_code"
