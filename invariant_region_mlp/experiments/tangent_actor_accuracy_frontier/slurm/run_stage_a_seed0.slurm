#!/bin/bash
#SBATCH --job-name=ta_frontier_stagea0
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
test -f invariant_region_mlp/results/nonconvex_tangent_actor/checkpoints/primary_d100_k2_lam0p2_eta0_main_mcv_s0_mcv_a4.pt

# Required resource profile precedes every larger n=4 decision.
srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
  --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed 0 \
  --n 4 --M 2 --role profile

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
  --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed 0 \
  --n 3 --M 384 --role endpoint

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.precision_evaluation \
  --training-seed 0 --paths 2048 --levels 64 128 256 512 \
  --fine-steps 512 --parent-batch 32 --path-chunk 128 --bootstrap-replicates 4000

srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
  --family primary --d 100 --lam 0.2 --actor-seed 0 --baseline-seed 0 \
  --n 3 --M 512 --role endpoint
