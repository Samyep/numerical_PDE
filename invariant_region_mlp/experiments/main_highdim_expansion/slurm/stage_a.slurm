#!/bin/bash
#SBATCH --job-name=icml_exp_stagea
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --mem=220G
#SBATCH --time=01:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/icml_two_highdim_expansion/invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/icml_two_highdim_expansion/invariant_region_mlp/results/main_highdim_expansion/slurm/%x_%j.err

set -euo pipefail
SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/icml_two_highdim_expansion
DPI_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/external/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$DPI_DEPS:$DPI_ROOT:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export WANDB_MODE=offline
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

python -m invariant_region_mlp.experiments.main_highdim_expansion.accounting --reserve-job-id "$SLURM_JOB_ID"
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv,noheader
python -c 'import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))'

srun --ntasks=1 --gpus-per-task=1 python -m invariant_region_mlp.experiments.main_highdim_expansion.calibrate --device cuda
python -m invariant_region_mlp.experiments.main_highdim_expansion.prepare
python -m invariant_region_mlp.experiments.main_highdim_expansion.reuse_audit
pytest -q invariant_region_mlp/experiments/main_highdim_expansion
srun --ntasks=1 --gpus-per-task=1 python -m invariant_region_mlp.experiments.main_highdim_expansion.validate --device cuda

for method in dpi pinn_hte; do
  for pde in b6_fullrank dense_fourier; do
    srun --ntasks=1 --gpus-per-task=1 python -m invariant_region_mlp.experiments.main_highdim_expansion.baseline_transfer \
      --method "$method" --pde "$pde" --dimension 100 --seed 0 --mode smoke --profile auto
  done
done
for pde in b6_fullrank dense_fourier; do
  srun --ntasks=1 --gpus-per-task=1 python -m invariant_region_mlp.experiments.main_highdim_expansion.sdgd_transfer \
    --family "$pde" --dimension 100 --seed 0 --device cuda --smoke --profile auto
done

srun --ntasks=1 --gpus-per-task=1 python -m invariant_region_mlp.experiments.main_highdim_expansion.characterize \
  --device cuda --dense-heat-pairs 4096
python -m invariant_region_mlp.experiments.main_highdim_expansion.runtime_profile
python -m invariant_region_mlp.experiments.main_highdim_expansion.review_gate

python - <<'PY'
from datetime import datetime, timezone
from invariant_region_mlp.experiments.main_highdim_expansion.protocol import RESULTS_DIR, write_json_atomic
write_json_atomic(RESULTS_DIR / "stage_a_complete.json", {
    "status": "COMPLETE",
    "finished_utc": datetime.now(timezone.utc).isoformat(),
    "job_id": __import__("os").environ["SLURM_JOB_ID"],
})
PY
