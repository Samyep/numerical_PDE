#!/bin/bash
#SBATCH --job-name=ta_frontier_cfh
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
python -c "import json; p=json.load(open('invariant_region_mlp/results/tangent_actor_accuracy_frontier/COUPLED_FH_FRONTIER_FREEZE.json')); assert not p['ground_truth_used_for_grid_selection']"
test -f invariant_region_mlp/results/tangent_actor_accuracy_frontier/predictions/post_freeze_coupled_d100_k2_lam0p2_eta1_main_mcv_s0.npz

for M in 256 384 512; do
  ROLE=frontier
  if [[ "$M" == 512 ]]; then ROLE=endpoint; fi
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
    --family coupled --d 100 --lam 0.2 --eta 1.0 --actor-seed 0 --baseline-seed 0 \
    --n 3 --M "$M" --role "$ROLE"
done

for M in 16 32 64 96; do
  ROLE=frontier
  if [[ "$M" == 96 ]]; then ROLE=endpoint; fi
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
    --family coupled --d 100 --lam 0.2 --eta 1.0 --actor-seed 0 --baseline-seed 0 \
    --n 4 --M "$M" --role "$ROLE"
done

for BASELINE_SEED in 1 2; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
    --family coupled --d 100 --lam 0.2 --eta 1.0 --actor-seed 0 --baseline-seed "$BASELINE_SEED" \
    --n 3 --M 512 --role endpoint
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_accuracy_frontier.full_history_frontier \
    --family coupled --d 100 --lam 0.2 --eta 1.0 --actor-seed 0 --baseline-seed "$BASELINE_SEED" \
    --n 4 --M 96 --role endpoint
done
