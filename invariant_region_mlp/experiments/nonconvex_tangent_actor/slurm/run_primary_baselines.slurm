#!/bin/bash
#SBATCH --job-name=nc_ta_baselines
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%j.err

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

for configuration in 1:64 1:256 2:8 2:16 2:32 2:64 3:4 3:8 3:12 3:16 3:24 3:32 3:48 3:64; do
  n=${configuration%%:*}
  M=${configuration##*:}
  srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.full_history_baseline \
    --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed 0 --n "$n" --M "$M"
done

for baseline_seed in 1 2; do
  for configuration in 3:32 3:64; do
    n=${configuration%%:*}
    M=${configuration##*:}
    srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.full_history_baseline \
      --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed "$baseline_seed" --n "$n" --M "$M"
  done
done

srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.pinn_baseline \
  --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed 0
