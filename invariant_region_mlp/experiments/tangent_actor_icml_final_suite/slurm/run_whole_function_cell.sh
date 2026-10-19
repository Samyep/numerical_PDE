#!/bin/bash
set -uo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DPI_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/external/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps
CELL_CAP_SECONDS=${CELL_CAP_SECONDS:-21600}
method=${1:?usage: run_whole_function_cell.sh METHOD PDE DIMENSION}
pde=${2:?usage: run_whole_function_cell.sh METHOD PDE DIMENSION}
dimension=${3:?usage: run_whole_function_cell.sh METHOD PDE DIMENSION}

case "$method" in tangent_actor|dpi|pinn_hte|sdgd_pinn) ;; *) echo "invalid method=$method" >&2; exit 2;; esac
case "$pde" in hjb_log|ou_gmm|joint_concave|joint_convex|nonconvex_coupled|mlp_failure|allen_cahn) ;; *) echo "invalid pde=$pde" >&2; exit 2;; esac
case "$dimension" in 20|50|100|200|400) ;; *) echo "invalid dimension=$dimension" >&2; exit 2;; esac

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
task_id=${WHOLE_CELL_TASK_ID:-${SLURM_ARRAY_TASK_ID:-0}}
tag="${method}__${pde}__d${dimension}"
cell_out=${CELL_OUT:-$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${job_id}_${task_id}_${tag}.out}
cell_err=${CELL_ERR:-$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${job_id}_${task_id}_${tag}.err}
started_utc=$(date -u +%FT%TZ)
started_epoch=$(date +%s)

record_status() {
  local status=$1 finished_utc=$2 seconds=$3 code=$4 notes=$5
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status record \
    --method "$method" --pde "$pde" --dimension "$dimension" --status "$status" \
    --job-id "$job_id" --task-id "$task_id" --started-utc "$started_utc" \
    --finished-utc "$finished_utc" --actual-gpu-seconds "$seconds" --exit-code "$code" \
    --stdout "$cell_out" --stderr "$cell_err" --notes "$notes"
}

echo "BEGIN cell=$tag utc=$started_utc git=$(git rev-parse HEAD) host=$(hostname) cap_seconds=$CELL_CAP_SECONDS"
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.review_gate
gate_status=$?
if [[ $gate_status -ne 0 ]]; then
  record_status FAILED "$(date -u +%FT%TZ)" 0 "$gate_status" "review gate rejected execution"
  exit "$gate_status"
fi
record_status RUNNING "" 0 0 "fixed seed 0; cumulative six GPU-hour cell cap"

command=()
case "$method" in
  tangent_actor)
    case "$pde" in
      hjb_log|ou_gmm|allen_cahn)
        command=(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor
          --family "$pde" --dimension "$dimension" --seed 0
          --evaluation-paths 1024 --evaluation-steps 64,128) ;;
      joint_concave|joint_convex|nonconvex_coupled)
        command=(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor
          --family "$pde" --dimension "$dimension" --seed 0
          --evaluation-paths 1024 --evaluation-steps 64,128) ;;
      mlp_failure)
        command=(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_mlp_failure_actor
          --dimension "$dimension" --seed 0 --evaluation-paths 1024 --evaluation-steps 64,128) ;;
    esac ;;
  dpi|pinn_hte)
    command=(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.official_dpi_transfer
      --method "$method" --pde "$pde" --dimension "$dimension" --seed 0 --mode official --profile auto) ;;
  sdgd_pinn)
    command=(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.sdgd_transfer
      --family "$pde" --dimension "$dimension" --seed 0 --profile auto) ;;
esac

timeout --signal=TERM --kill-after=120 "$CELL_CAP_SECONDS" "${command[@]}" >"$cell_out" 2>"$cell_err"
code=$?
finished_utc=$(date -u +%FT%TZ)
finished_epoch=$(date +%s)
seconds=$((finished_epoch - started_epoch))
if [[ $code -eq 0 ]]; then
  status=COMPLETE; notes="runner returned complete under frozen protocol"
elif [[ $code -eq 124 ]]; then
  status=OUT_OF_TIME; notes="current allocation ended before completion; cumulative cap remains six GPU hours"
elif [[ $code -eq 137 ]] || grep -Eiq 'CUDA out of memory|OutOfMemoryError|oom-kill|oom_kill' "$cell_err"; then
  status=OUT_OF_MEMORY; notes="memory failure detected in exit code or stderr"
else
  status=FAILED; notes="runner failed; inspect immutable per-cell stderr"
fi
record_status "$status" "$finished_utc" "$seconds" "$code" "$notes"
echo "END cell=$tag status=$status exit_code=$code gpu_seconds=$seconds utc=$finished_utc"
exit "$code"
