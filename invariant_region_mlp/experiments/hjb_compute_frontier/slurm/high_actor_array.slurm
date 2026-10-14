#!/bin/bash
#SBATCH --job-name=hjbfront_actorhi
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=04:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/hjb_compute_frontier/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/hjb_compute_frontier/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
PYTHON_BIN=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python
EQUATIONS=(b1 hjb_log)
EQUATION=${EQUATIONS[$SLURM_ARRAY_TASK_ID]}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT/invariant_region_mlp/experiments${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python -V
git rev-parse HEAD
printf 'high-quality Actor equation=%s seed=0\n' "$EQUATION"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'device', torch.cuda.get_device_name(0))"

srun "$PYTHON_BIN" -u -m hjb_compute_frontier.run_experiment actor \
  --equation "$EQUATION" --seed 0 --profile high
