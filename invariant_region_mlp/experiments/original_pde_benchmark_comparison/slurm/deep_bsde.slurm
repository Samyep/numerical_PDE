#!/bin/bash
#SBATCH --job-name=orig_deep_bsde
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=384G
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

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

run_task() {
  local task=$1
  local benchmark=$2
  local seed=$3
  local out=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/deep_bsde_${SLURM_JOB_ID}_${task}.out
  local err=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/deep_bsde_${SLURM_JOB_ID}_${task}.err
  srun --network=no_vni --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=12 --mem=90G \
    --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
    python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.deep_bsde_port \
      --benchmark "$benchmark" --seed "$seed" >"$out" 2>"$err" &
}

run_task 0 hjb_log_d100 0
run_task 1 hjb_log_d100 1
run_task 2 hjb_log_d100 2
run_task 3 allen_cahn_d100 0
wait
run_task 4 allen_cahn_d100 1
run_task 5 allen_cahn_d100 2
wait

