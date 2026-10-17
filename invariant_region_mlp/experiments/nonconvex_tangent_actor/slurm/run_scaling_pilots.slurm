#!/bin/bash
#SBATCH --job-name=nc_ta_scaling
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:15:00
#SBATCH --array=0-5%2
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/nonconvex_tangent_actor/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source

case "$SLURM_ARRAY_TASK_ID" in
  0) FAMILY=primary; DIMENSION=400; LAMBDA=0.2; ETA=0.0; ESTIMATOR=mcv ;;
  1) FAMILY=primary; DIMENSION=100; LAMBDA=0.0; ETA=0.0; ESTIMATOR=mcv ;;
  2) FAMILY=primary; DIMENSION=100; LAMBDA=0.05; ETA=0.0; ESTIMATOR=mcv ;;
  3) FAMILY=primary; DIMENSION=100; LAMBDA=0.5; ETA=0.0; ESTIMATOR=mcv ;;
  4) FAMILY=coupled; DIMENSION=100; LAMBDA=0.2; ETA=1.0; ESTIMATOR=mcv ;;
  5) FAMILY=primary; DIMENSION=100; LAMBDA=0.2; ETA=0.0; ESTIMATOR=plain ;;
  *) echo "unexpected array index: $SLURM_ARRAY_TASK_ID" >&2; exit 2 ;;
esac

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
  --family "$FAMILY" --d "$DIMENSION" --lam "$LAMBDA" --eta "$ETA" --seed 0 --budget pilot --estimator "$ESTIMATOR"
srun --network=no_vni python -m invariant_region_mlp.experiments.nonconvex_tangent_actor.evaluate_postfreeze \
  --family "$FAMILY" --d "$DIMENSION" --lam "$LAMBDA" --eta "$ETA" --seed 0 --budget pilot --estimator "$ESTIMATOR"
