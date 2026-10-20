#!/bin/bash
set -uo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/icml_two_highdim_expansion
DPI_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/external/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps
CELL_CAP_SECONDS=${CELL_CAP_SECONDS:-21420}
method=${1:?usage: run_cell.sh METHOD PDE DIMENSION}
pde=${2:?usage: run_cell.sh METHOD PDE DIMENSION}
dimension=${3:?usage: run_cell.sh METHOD PDE DIMENSION}

case "$method" in tangent_actor|dpi|pinn_hte|sdgd_pinn) ;; *) exit 2;; esac
case "$pde" in b6_fullrank|dense_fourier) ;; *) exit 2;; esac
case "$dimension" in 20|50|100|200|400) ;; *) exit 2;; esac

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$DPI_DEPS:$DPI_ROOT:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export WANDB_MODE=offline
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

job_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}
task_id=${EXPANSION_CELL_TASK_ID:-${SLURM_ARRAY_TASK_ID:-0}}
tag="${method}__${pde}__d${dimension}"
cell_out=${CELL_OUT:-$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${job_id}_${task_id}_${tag}.out}
cell_err=${CELL_ERR:-$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${job_id}_${task_id}_${tag}.err}
started_utc=$(date -u +%FT%TZ)
started_epoch=$(date +%s)

read -r prior_status remaining_seconds prior_seconds prior_attempts < <(
  python -m invariant_region_mlp.experiments.main_highdim_expansion.status budget \
    --method "$method" --pde "$pde" --dimension "$dimension"
)
case "$prior_status" in
  COMPLETE|OUT_OF_TIME|OUT_OF_MEMORY|FAILED)
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
finished_utc=$(date -u +%FT%TZ)
seconds=$(( $(date +%s) - started_epoch ))
if [[ $code -eq 0 ]]; then
  status=COMPLETE; notes="runner completed under frozen protocol"
elif [[ $code -eq 124 ]]; then
  if (( prior_seconds + seconds >= 21420 )); then
    status=OUT_OF_TIME; notes="cumulative six-GH200-hour cap reached"
  else
    status=CHECKPOINTED; notes="allocation boundary; TA/SDGD resume exactly, DPI/PINN repeated work remains charged"
  fi
elif [[ $code -eq 137 ]] || grep -Eiq 'CUDA out of memory|OutOfMemoryError|oom-kill|oom_kill' "$cell_err"; then
  status=OUT_OF_MEMORY; notes="memory failure detected"
else
  status=FAILED; notes="runner failed; inspect immutable cell logs"
fi
record_status "$status" "$finished_utc" "$seconds" "$code" "$notes"
echo "END cell=$tag status=$status exit=$code gpu_seconds=$seconds utc=$finished_utc"
exit "$code"
