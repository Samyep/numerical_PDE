#!/bin/bash
#SBATCH --job-name=icml_exp_wave
#SBATCH --account=isye
#SBATCH --partition=ice-gpu
#SBATCH --qos=coe-ice
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:h100:4
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.out
#SBATCH --error=invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=${SOURCE_ROOT:-${SLURM_SUBMIT_DIR:?submit from the repository root}}
SCRIPT_DIR=$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm
source "$SCRIPT_DIR/common_env.sh"
TASK_FILE=${TASK_FILE:?submit with immutable TASK_FILE}
TASK_FILE_SHA256=${TASK_FILE_SHA256:?submit with TASK_FILE_SHA256}
WORKERS=${WORKERS:-4}
[[ "$(sha256sum "$TASK_FILE" | awk '{print $1}')" == "$TASK_FILE_SHA256" ]]
[[ "$WORKERS" =~ ^[1-4]$ ]]

python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate
python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,uuid,name,memory.total --format=csv,noheader || true

worker="$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/run_queue_worker.sh"
pids=()
for ((slot=0; slot<WORKERS; slot++)); do
  slot_out="$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${SLURM_JOB_ID}_slot${slot}.out"
  slot_err="$SOURCE_ROOT/invariant_region_mlp/results/main_highdim_expansion/slurm/${SLURM_JOB_ID}_slot${slot}.err"
  srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --mem=64G --gpus-per-task=1 --gpu-bind=single:1 \
    --output="$slot_out" --error="$slot_err" "$worker" "$TASK_FILE" "$slot" &
  pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
  if wait "$pid"; then
    code=0
  else
    code=$?
  fi
  [[ $code -eq 0 ]] && continue
  if [[ $code -eq 75 ]]; then
    status=75
  elif [[ $status -eq 0 ]]; then
    status=$code
  fi
done
exit "$status"
