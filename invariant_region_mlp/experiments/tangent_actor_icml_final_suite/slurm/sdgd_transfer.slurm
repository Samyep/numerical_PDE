#!/bin/bash
#SBATCH --job-name=icml_sdgd_transfer
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "${FAMILY:?set FAMILY to a frozen suite PDE key}"
DIMENSION=${DIMENSION:-100}
TRAINING_SEED=${TRAINING_SEED:-0}
PROFILE=${PROFILE:-auto}
if [[ "$TRAINING_SEED" != "0" ]]; then
  echo "Single-seed protocol permits only TRAINING_SEED=0" >&2
  exit 2
fi

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.sdgd_transfer \
  --family "$FAMILY" --dimension "$DIMENSION" --seed "$TRAINING_SEED" --profile "$PROFILE"
