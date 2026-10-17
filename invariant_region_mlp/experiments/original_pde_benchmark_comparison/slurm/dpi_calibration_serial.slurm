#!/bin/bash
#SBATCH --job-name=orig_dpi_cal_s
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DPI_ROOT=/work/nvme/bgvf/ywang151/original_pde_source_audit/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps
DPI_INPUTS=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/official_work/dpi_inputs
METHODS=(dpi pinn_hte d_dbsde)

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$DPI_INPUTS"
export PYTHONPATH="$DPI_DEPS:$DPI_ROOT:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python -V
git -C "$SOURCE_ROOT" rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

for METHOD in "${METHODS[@]}"; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.official_dpi_runner \
    --method "$METHOD" --seed 0 --calibration
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.evaluate_official_dpi \
    --method "$METHOD" --seed 0 --calibration
done
