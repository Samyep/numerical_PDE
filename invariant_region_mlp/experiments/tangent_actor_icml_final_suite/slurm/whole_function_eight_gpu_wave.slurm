#!/bin/bash
#SBATCH --job-name=icml_whole_wave
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=440G
#SBATCH --time=02:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TASK_FILE=${TASK_FILE:?submit with an immutable queue CSV}
TASK_FILE_SHA256=${TASK_FILE_SHA256:?submit with TASK_FILE_SHA256}
WORKERS=${WORKERS:-8}
actual_sha=$(sha256sum "$TASK_FILE" | awk '{print $1}')
[[ "$actual_sha" == "$TASK_FILE_SHA256" ]] || { echo "task file changed after submission" >&2; exit 2; }
[[ "$WORKERS" =~ ^[1-8]$ ]] || { echo "WORKERS must be in [1,8]" >&2; exit 2; }

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.review_gate
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting \
  --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print("torch.cuda.is_available=", torch.cuda.is_available(), "count=", torch.cuda.device_count())'

worker="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_queue_worker.sh"
pids=()
for ((slot=0; slot<WORKERS; slot++)); do
  slot_out="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_queue_slot${slot}.out"
  slot_err="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_ID}_queue_slot${slot}.err"
  srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --mem=110G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
    --output="$slot_out" --error="$slot_err" "$worker" "$TASK_FILE" "$slot" &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
exit "$status"
