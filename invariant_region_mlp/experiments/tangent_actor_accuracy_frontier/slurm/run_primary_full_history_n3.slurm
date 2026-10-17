#!/bin/bash
#SBATCH --job-name=ta_frontier_fh_n3
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --array=0-5%3
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
CONFIGURATIONS=("384:0" "384:1" "384:2" "512:0" "512:1" "512:2")
CONFIGURATION=${CONFIGURATIONS[$SLURM_ARRAY_TASK_ID]}
M=${CONFIGURATION%%:*}
BASELINE_SEED=${CONFIGURATION##*:}

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
test -f invariant_region_mlp/results/tangent_actor_accuracy_frontier/references/primary_refinement_audit.json

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
  --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed "$BASELINE_SEED" \
  --n 3 --M "$M" --role endpoint

