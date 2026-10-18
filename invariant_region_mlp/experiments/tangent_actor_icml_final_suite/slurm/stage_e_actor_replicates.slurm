#!/bin/bash
#SBATCH --job-name=icml_actor_reps
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=12:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DIMENSION=${DIMENSION:-100}
TRAINING_SEEDS=${TRAINING_SEEDS:-"0"}
FAMILIES=${FAMILIES:-"joint_concave joint_convex mlp_failure"}
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

echo "git_head=$(git rev-parse HEAD)"
echo "dimension=$DIMENSION seeds=$TRAINING_SEEDS families=$FAMILIES"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

for seed in $TRAINING_SEEDS; do
  if [[ "$seed" != "0" ]]; then
    echo "Single-seed protocol permits only training seed 0; got $seed" >&2
    exit 2
  fi
  for family in $FAMILIES; do
    echo "BEGIN family=$family dimension=$DIMENSION seed=$seed utc=$(date -u +%FT%TZ)"
    case "$family" in
      joint_concave|joint_convex|nonconvex_coupled)
        srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor \
          --family "$family" --dimension "$DIMENSION" --seed "$seed" \
          --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
        ;;
      mlp_failure)
        srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_mlp_failure_actor \
          --dimension "$DIMENSION" --seed "$seed" \
          --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
        ;;
      hjb_log|ou_gmm|allen_cahn)
        srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor \
          --family "$family" --dimension "$DIMENSION" --seed "$seed" \
          --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
        ;;
      *)
        echo "Unsupported FAMILY=$family" >&2
        exit 2
        ;;
    esac
    echo "END family=$family dimension=$DIMENSION seed=$seed utc=$(date -u +%FT%TZ)"
  done
done
