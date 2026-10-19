#!/bin/bash
#SBATCH --job-name=icml_whole_cell
#SBATCH --account=bgvf-dtai-gh
#SBATCH --partition=ghx4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=verbose,closest
#SBATCH --time=06:00:00
#SBATCH --exclude=gh097
#SBATCH --output=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%A_%a.out
#SBATCH --error=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source/invariant_region_mlp/results/tangent_actor_icml_final_suite/slurm/%x_%A_%a.err

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TASK_FILE=${TASK_FILE:?submit with an immutable three-column method,pde,dimension task file}
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASK_FILE")
IFS=',' read -r method pde dimension extra <<< "$line"
if [[ -z "${method:-}" || -z "${pde:-}" || ! "${dimension:-}" =~ ^[0-9]+$ || -n "${extra:-}" ]]; then
  echo "Invalid or missing task line for array index $SLURM_ARRAY_TASK_ID: $line" >&2
  exit 2
fi
exec env CELL_CAP_SECONDS=21420 \
  "$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh" \
  "$method" "$pde" "$dimension"
