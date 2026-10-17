#!/bin/bash
#SBATCH --job-name=ta_frontier_d400_lr
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RANK=${RANK:?submit with --export=ALL,RANK=8|16|32}
TRAINING_SEED=${TRAINING_SEED:-0}

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

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.run_low_rank_training \
  --family primary --d 400 --lam 0.2 --seed "$TRAINING_SEED" --rank "$RANK"
srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.evaluate_low_rank \
  --family primary --d 400 --lam 0.2 --seed "$TRAINING_SEED" --rank "$RANK"
