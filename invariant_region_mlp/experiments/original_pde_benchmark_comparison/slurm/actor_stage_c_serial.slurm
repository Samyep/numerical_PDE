#!/bin/bash
#SBATCH --job-name=orig_actor_stage_c
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --exclude=gh097
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=02:00:00
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/original_pde_benchmark_comparison/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
: "${BENCHMARK:?set BENCHMARK to hjb_log_d100, dpi_ou_gmm_d100, or allen_cahn_d100}"
DECISION=$SOURCE_ROOT/invariant_region_mlp/results/original_pde_benchmark_comparison/stage_b_selection_${BENCHMARK}.json

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

mapfile -t SELECTED < <(python -c 'import json,sys; print(*json.load(open(sys.argv[1]))["selected_for_continuation"], sep="\n")' "$DECISION")
if [[ ${#SELECTED[@]} -ne 2 ]]; then
  echo "expected two selected trials in $DECISION" >&2
  exit 1
fi

python -V
git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

for TRIAL in "${SELECTED[@]}"; do
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor \
    --benchmark "$BENCHMARK" --trial "$TRIAL" --seed 0 --iterations 3
done

srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor \
  --benchmark "$BENCHMARK" --trial data8192 --seed 0 --iterations 3
srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor \
  --benchmark "$BENCHMARK" --trial focused --seed 0 --iterations 4
if [[ "$BENCHMARK" != "allen_cahn_d100" ]]; then
  srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.run_actor \
    --benchmark "$BENCHMARK" --trial focused_plain --seed 0 --iterations 3
fi

srun --network=no_vni python -m invariant_region_mlp.experiments.original_pde_benchmark_comparison.rank_actor_checkpoints \
  --benchmark "$BENCHMARK" --seed 0 --paths 512 --steps 64 --replicates 3
