#!/bin/bash
#SBATCH --job-name=icml_baseline_pack
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=4
#SBATCH --mem=440G
#SBATCH --time=02:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%j.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DPI_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/external/DeepPicardIteration
DPI_DEPS=/work/nvme/bgvf/ywang151/original_pde_source_audit/dpi_deps
: "${PACK_TASKS:?set PACK_TASKS to four comma-separated GPU slots; join sequential tasks in a slot with +}"

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$DPI_DEPS:$DPI_ROOT:$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NCCL_SOCKET_IFNAME=hsn
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export WANDB_MODE=offline
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

git rev-parse HEAD
git -C "$DPI_ROOT" rev-parse HEAD
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv,noheader
python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'devices', torch.cuda.device_count())"

IFS=',' read -r -a slots <<< "$PACK_TASKS"
if [[ ${#slots[@]} -ne 4 ]]; then
  echo "PACK_TASKS must contain exactly four GPU slots; received ${#slots[@]}" >&2
  exit 2
fi

run_task() {
  local spec=$1
  local method family dimension profile extra
  IFS=':' read -r method family dimension profile extra <<< "$spec"
  if [[ -n "${extra:-}" || -z "${profile:-}" || ! "$dimension" =~ ^[0-9]+$ ]]; then
    echo "Invalid task specification: $spec" >&2
    return 2
  fi
  case "$family" in
    hjb_log|ou_gmm|joint_concave|joint_convex|nonconvex_coupled|mlp_failure|allen_cahn) ;;
    *) echo "Unknown PDE family: $family" >&2; return 2 ;;
  esac

  local tag="${method}_${family}_d${dimension}"
  local task_out="invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${tag}.out"
  local task_err="invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${tag}.err"
  echo "BEGIN task=$spec utc=$(date -u +%FT%TZ) stdout=$task_out"

  case "$method" in
    dpi|pinn_hte)
      srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
        --output="$task_out" --error="$task_err" \
        python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.official_dpi_transfer \
        --method "$method" --pde "$family" --dimension "$dimension" \
        --seed 0 --mode official --profile "$profile"
      ;;
    sdgd)
      srun --exclusive --exact --nodes=1 --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
        --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=single:1 --network=no_vni \
        --output="$task_out" --error="$task_err" \
        python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.sdgd_transfer \
        --family "$family" --dimension "$dimension" --seed 0 --profile "$profile"
      ;;
    *)
      echo "Unknown method: $method" >&2
      return 2
      ;;
  esac
  echo "END task=$spec utc=$(date -u +%FT%TZ)"
}

run_slot() {
  local chain=$1
  local -a slot_tasks
  IFS='+' read -r -a slot_tasks <<< "$chain"
  for task in "${slot_tasks[@]}"; do
    run_task "$task"
  done
}

pids=()
for slot in "${slots[@]}"; do
  run_slot "$slot" &
  pids+=("$!")
done

status=0
for index in "${!pids[@]}"; do
  if ! wait "${pids[$index]}"; then
    echo "FAILED slot=${slots[$index]}" >&2
    status=1
  fi
done
exit "$status"
