#!/bin/bash
#SBATCH --job-name=ta_frontier_stagearep
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_accuracy_frontier/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TRAINING_SEED=${TRAINING_SEED:?submit with --export=ALL,TRAINING_SEED=1|2}
if [[ "$TRAINING_SEED" != 1 && "$TRAINING_SEED" != 2 ]]; then
  echo "replication seed must be 1 or 2" >&2
  exit 2
fi

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

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.precision_evaluation \
  --training-seed "$TRAINING_SEED" --paths 2048 --levels 64 128 256 512 \
  --fine-steps 512 --parent-batch 32 --path-chunk 128 --bootstrap-replicates 4000

for M in 384 512; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
    --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed "$TRAINING_SEED" \
    --n 3 --M "$M" --role endpoint
done
