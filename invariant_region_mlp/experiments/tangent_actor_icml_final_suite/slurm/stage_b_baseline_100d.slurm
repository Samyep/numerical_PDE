#!/bin/bash
#SBATCH --job-name=icml_stageb_base
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DPI_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/external/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps
BASELINE=${BASELINE:?submit with BASELINE=deep_bsde|dpi|pinn_hte|sdgd_pinn}
DIMENSION=${DIMENSION:-100}
TRAINING_SEED=${TRAINING_SEED:-0}

case "$BASELINE" in
  deep_bsde)
    DEFAULT_FAMILIES="ou_gmm joint_concave joint_convex nonconvex_coupled mlp_failure"
    ;;
  dpi|pinn_hte)
    DEFAULT_FAMILIES="hjb_log joint_concave joint_convex nonconvex_coupled mlp_failure allen_cahn"
    ;;
  sdgd_pinn)
    DEFAULT_FAMILIES="ou_gmm joint_concave joint_convex nonconvex_coupled mlp_failure allen_cahn"
    ;;
  *)
    echo "Unsupported BASELINE=$BASELINE" >&2
    exit 2
    ;;
esac
FAMILIES=${FAMILIES:-$DEFAULT_FAMILIES}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
if [[ "$BASELINE" == "dpi" || "$BASELINE" == "pinn_hte" ]]; then
  export PYTHONPATH="$DPI_DEPS:$DPI_ROOT:$PYTHONPATH"
fi
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export WANDB_MODE=offline
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "git_head=$(git rev-parse HEAD)"
echo "baseline=$BASELINE dimension=$DIMENSION seed=$TRAINING_SEED families=$FAMILIES"
if [[ "$BASELINE" == "dpi" || "$BASELINE" == "pinn_hte" ]]; then
  echo "upstream_dpi_head=$(git -C "$DPI_ROOT" rev-parse HEAD)"
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

for family in $FAMILIES; do
  echo "BEGIN baseline=$BASELINE family=$family dimension=$DIMENSION seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
  case "$BASELINE" in
    deep_bsde)
      srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.deep_bsde_transfer \
        --family "$family" --dimension "$DIMENSION" --seed "$TRAINING_SEED"
      ;;
    dpi|pinn_hte)
      srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.official_dpi_transfer \
        --method "$BASELINE" --pde "$family" --dimension "$DIMENSION" \
        --seed "$TRAINING_SEED" --mode official
      ;;
    sdgd_pinn)
      srun --network=no_vni python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.sdgd_transfer \
        --family "$family" --dimension "$DIMENSION" --seed "$TRAINING_SEED"
      ;;
  esac
  echo "END baseline=$BASELINE family=$family dimension=$DIMENSION seed=$TRAINING_SEED utc=$(date -u +%FT%TZ)"
done
