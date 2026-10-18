#!/bin/bash
#SBATCH --job-name=icml_stagec_actor
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=2-00:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TRAINING_SEED=${TRAINING_SEED:-0}
SMOOTH_DIMS=${SMOOTH_DIMS:-"10 400"}
MLP_DIMS=${MLP_DIMS:-"20 1000"}
FAMILIES=${FAMILIES:-"hjb_log ou_gmm joint_concave joint_convex nonconvex_coupled allen_cahn"}
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
echo "seed=$TRAINING_SEED smooth_dims=$SMOOTH_DIMS mlp_dims=$MLP_DIMS families=$FAMILIES"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.review_gate
if [[ "$TRAINING_SEED" != "0" ]]; then
  echo "Single-seed protocol permits only TRAINING_SEED=0" >&2
  exit 2
fi

for dimension in $SMOOTH_DIMS; do
  for family in $FAMILIES; do
    echo "BEGIN family=$family dimension=$dimension seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
    case "$family" in
      hjb_log|ou_gmm|allen_cahn)
        srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_standard_actor \
          --family "$family" --dimension "$dimension" --seed "$TRAINING_SEED" \
          --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
        ;;
      joint_concave|joint_convex|nonconvex_coupled)
        srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_joint_actor \
          --family "$family" --dimension "$dimension" --seed "$TRAINING_SEED" \
          --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
        ;;
      *)
        echo "Unsupported FAMILY=$family" >&2
        exit 2
        ;;
    esac
    echo "END family=$family dimension=$dimension seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
  done
done

for dimension in $MLP_DIMS; do
  echo "BEGIN family=mlp_failure dimension=$dimension seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
  srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.run_mlp_failure_actor \
    --dimension "$dimension" --seed "$TRAINING_SEED" \
    --evaluation-paths "$EVALUATION_PATHS" --evaluation-steps "$EVALUATION_STEPS"
  echo "END family=mlp_failure dimension=$dimension seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
done
