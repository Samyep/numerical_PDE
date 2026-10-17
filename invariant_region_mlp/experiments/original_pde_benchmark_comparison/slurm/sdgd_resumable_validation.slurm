#!/bin/bash
#SBATCH --job-name=orig_sdgd_resume_val
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=00:30:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
SDGD_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/sdgd_deps
module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SDGD_DEPS:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

python -V
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import openpyxl, torch; print('openpyxl', openpyxl.__version__, 'torch', torch.__version__, 'cuda', torch.cuda.is_available())"
srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.official_sdgd_resumable \
  --seed 0 --start 0 --end 100 --total-epochs 100 --tag equivalence100
srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.official_sdgd_resumable \
  --seed 0 --tag equivalence100 --compare-validation
