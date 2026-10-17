#!/bin/bash
#SBATCH --job-name=orig_actor_broad
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=384G
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=01:30:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "${BENCHMARK:?submit with --export=ALL,BENCHMARK=<name>}"
TRIALS=(pilot192 wide384 deep_elu quality128)

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

for TASK in 0 1 2 3; do
  TRIAL=${TRIALS[$TASK]}
  TASK_OUT=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/actor_${BENCHMARK}_${SLURM_JOB_ID}_${TASK}.out
  TASK_ERR=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/actor_${BENCHMARK}_${SLURM_JOB_ID}_${TASK}.err
  srun --network=no_vni --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task=12 --mem=90G \
    --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
    python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor \
      --benchmark "$BENCHMARK" --trial "$TRIAL" --seed 0 --iterations 1 \
      >"$TASK_OUT" 2>"$TASK_ERR" &
done

wait

