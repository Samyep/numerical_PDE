#!/bin/bash
#SBATCH --job-name=icml_joint_actor
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
FAMILY=${FAMILY:?submit with FAMILY=joint_concave|joint_convex|nonconvex_coupled}
DIMENSION=${DIMENSION:-100}
TRAINING_SEED=${TRAINING_SEED:-0}
EVALUATION_PATHS=${EVALUATION_PATHS:-1024}
EVALUATION_STEPS=${EVALUATION_STEPS:-64,128}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "git_head=$(git rev-parse HEAD)"
echo "host=$(hostname) family=$FAMILY dimension=$DIMENSION seed=$TRAINING_SEED"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'visible_devices', torch.cuda.device_count())"

RUN_ARGS=(
  --family "$FAMILY" --dimension "$DIMENSION" --seed "$TRAINING_SEED"
  --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
)
if [[ "${EVALUATE_ONLY:-0}" == "1" ]]; then
  RUN_ARGS+=(--evaluate-only)
fi
srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor "${RUN_ARGS[@]}"
