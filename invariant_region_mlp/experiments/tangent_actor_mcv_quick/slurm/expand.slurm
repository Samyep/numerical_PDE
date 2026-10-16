#!/bin/bash
#SBATCH --job-name=ta_mcv_expand
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=384G
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=00:30:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_mcv_quick/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_mcv_quick/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python -V
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'visible_devices', torch.cuda.device_count())"

BENCHMARKS=(nonseparable_d100_k2 nonseparable_d100_k2 joint_concave_d400_k2 joint_concave_d400_k2)
SEEDS=(0 1 0 1)
CHECKPOINTS=(a3 all all all)
TAGS=(expand_a3 expand_seed1 expand_d400_s0 expand_d400_s1)

for TASK in 0 1 2 3; do
  TASK_OUT=$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_mcv_quick/slurm/expand_${SLURM_JOB_ID}_${TASK}.out
  TASK_ERR=$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_mcv_quick/slurm/expand_${SLURM_JOB_ID}_${TASK}.err
  srun --network=no_vni --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=12 --mem=90G \
    --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
    python -m invariant_region_mlp.experiments.tangent_actor_mcv_quick.evaluate \
      --benchmark "${BENCHMARKS[$TASK]}" --seed "${SEEDS[$TASK]}" \
      --checkpoints "${CHECKPOINTS[$TASK]}" --run-tag "${TAGS[$TASK]}" \
      >"$TASK_OUT" 2>"$TASK_ERR" &
done

wait
