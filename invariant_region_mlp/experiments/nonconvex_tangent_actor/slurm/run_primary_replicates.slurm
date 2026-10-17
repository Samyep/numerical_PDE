#!/bin/bash
#SBATCH --job-name=nc_ta_replicates
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:30:00
#SBATCH --array=1-2
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TRAINING_SEED=$SLURM_ARRAY_TASK_ID

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

srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.run_training \
  --family primary --d 100 --lam 0.2 --seed "$TRAINING_SEED" --budget main --estimator mcv
srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.evaluate_postfreeze \
  --family primary --d 100 --lam 0.2 --seed "$TRAINING_SEED" --budget main --estimator mcv
