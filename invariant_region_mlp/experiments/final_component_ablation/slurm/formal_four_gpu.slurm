#!/bin/bash
#SBATCH --job-name=ta_final_ablation
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=220G
#SBATCH --time=06:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/final_component_ablation/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/final_component_ablation/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TASK_FILE="$SOURCE_ROOT/invariant_region_mlp/experiments/final_component_ablation/formal_queue.csv"
TASK_HASH_FILE="$TASK_FILE.sha256"
WORKER="$SOURCE_ROOT/invariant_region_mlp/experiments/final_component_ablation/slurm/run_queue_worker.sh"

cd "$SOURCE_ROOT"
sha256sum --check "$TASK_HASH_FILE"
module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn

echo "git_head=$(git rev-parse HEAD)"
echo "job_id=$SLURM_JOB_ID nodes=$SLURM_JOB_NUM_NODES tasks=$SLURM_NTASKS"
python -V
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print("torch", torch.__version__, "cuda", torch.cuda.is_available(), "count", torch.cuda.device_count())'
python -m invariant_region_mlp.experiments.final_component_ablation.preflight

worker_count=${SLURM_NTASKS:?SLURM_NTASKS is required}
pids=()
for ((slot = 0; slot < worker_count; slot++)); do
  slot_out="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation/slurm/${SLURM_JOB_ID}_slot${slot}.out"
  slot_err="$SOURCE_ROOT/invariant_region_mlp/results/final_component_ablation/slurm/${SLURM_JOB_ID}_slot${slot}.err"
  srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    --mem=50G --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
    --output="$slot_out" --error="$slot_err" \
    "$WORKER" "$TASK_FILE" "$slot" "$worker_count" &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
exit "$status"
