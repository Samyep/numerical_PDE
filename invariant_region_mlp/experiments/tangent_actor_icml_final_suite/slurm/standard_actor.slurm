#!/bin/bash
#SBATCH --job-name=icml_std_actor
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "${FAMILY:?set FAMILY to hjb_log, ou_gmm, or allen_cahn}"
DIMENSION=${DIMENSION:-100}
TRAINING_SEED=${TRAINING_SEED:-0}
if [[ "$TRAINING_SEED" != "0" ]]; then
  echo "Single-seed protocol permits only TRAINING_SEED=0" >&2
  exit 2
fi
EVALUATION_PATHS=${EVALUATION_PATHS:-1024}
EVALUATION_STEPS=${EVALUATION_STEPS:-64,128}

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
RUN_ARGS=(
  --family "$FAMILY" --dimension "$DIMENSION" --seed "$TRAINING_SEED"
  --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
)
if [[ "${EVALUATE_ONLY:-0}" == "1" ]]; then
  RUN_ARGS+=(--evaluate-only)
fi
srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor "${RUN_ARGS[@]}"
