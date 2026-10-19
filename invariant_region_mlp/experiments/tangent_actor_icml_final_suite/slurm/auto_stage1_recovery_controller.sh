#!/bin/bash
# Inspect the d50 pilot, resume only exact-checkpoint methods when necessary,
# then hand the verified Stage 1 result to the Stage 2/3 controller.
set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULTS_DIR="$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite"
SINGLE="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_stage1_resume.slurm"
DUAL="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/whole_function_stage1_dual_resume.slurm"
FULL="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/auto_whole_function_controller.sh"
PREDECESSOR=${1:?usage: auto_stage1_recovery_controller.sh D50_JOB_ID}

source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"

wait_for_job() {
  local job=$1 state=""
  while squeue -h -j "$job" | grep -q .; do
    echo "$(date -u +%FT%TZ) recovery waiting job=$job $(squeue -h -j "$job" -o '%T|%M|%S|%R')"
    sleep 30
  done
  for _ in $(seq 1 10); do
    state=$(sacct -X -n -j "$job" -o State | awk 'NF {print $1; exit}')
    [[ -n "$state" ]] && break
    sleep 10
  done
  case "$state" in COMPLETED*|FAILED*|TIMEOUT*) ;; *) return 3 ;; esac
}

refresh() {
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.aggregate_results
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_job_accounting
  python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_artifacts
}

stage1_statuses() {
  python - <<'PY'
from invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status import METHODS, cell_budget
print(" ".join(f"{method}={cell_budget(method, 'joint_concave', 50)['status']}" for method in METHODS))
PY
}

wait_for_job "$PREDECESSOR"
refresh
read -ra pairs <<< "$(stage1_statuses)"
declare -A status
for pair in "${pairs[@]}"; do status[${pair%%=*}]=${pair#*=}; done
incomplete=()
for method in tangent_actor dpi pinn_hte sdgd_pinn; do
  [[ "${status[$method]}" == COMPLETE ]] || incomplete+=("$method")
done
echo "$(date -u +%FT%TZ) d50 statuses: ${pairs[*]}"

trigger_job=$PREDECESSOR
if (( ${#incomplete[@]} > 0 )); then
  for method in "${incomplete[@]}"; do
    if [[ "${status[$method]}" != CHECKPOINTED || "$method" == dpi || "$method" == pinn_hte ]]; then
      echo "non-resumable or terminal d50 failure: method=$method status=${status[$method]}" >&2
      exit 4
    fi
  done
  if squeue -h -u "$USER" -o '%j' | grep -q '^icml_whole'; then
    echo "duplicate suite job active" >&2
    exit 6
  fi
  if (( ${#incomplete[@]} == 2 )); then
    export_vars="ALL,DIMENSION=50,METHOD_A=${incomplete[0]},METHOD_B=${incomplete[1]},CAP_A_SECONDS=3300,CAP_B_SECONDS=3300"
    sbatch --test-only --partition=ghx4-interactive --time=01:00:00 --export="$export_vars" "$DUAL"
    trigger_job=$(sbatch --parsable --partition=ghx4-interactive --time=01:00:00 --export="$export_vars" "$DUAL")
  elif (( ${#incomplete[@]} == 1 )); then
    export_vars="ALL,METHOD=${incomplete[0]},PDE=joint_concave,DIMENSION=50,RESUME_CAP_SECONDS=3300"
    sbatch --test-only --partition=ghx4-interactive --time=01:00:00 --export="$export_vars" "$SINGLE"
    trigger_job=$(sbatch --parsable --partition=ghx4-interactive --time=01:00:00 --export="$export_vars" "$SINGLE")
  else
    echo "unexpected number of checkpointed d50 cells: ${#incomplete[@]}" >&2
    exit 5
  fi
  echo "$trigger_job" > "$RESULTS_DIR/STAGE1_D50_RESUME_JOB_ID.txt"
  echo "$(date -u +%FT%TZ) submitted d50 checkpoint resume job=$trigger_job methods=${incomplete[*]}"
  wait_for_job "$trigger_job"
  refresh
  for method in tangent_actor dpi pinn_hte sdgd_pinn; do
    read -r final_status _ < <(python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status budget \
      --method "$method" --pde joint_concave --dimension 50)
    [[ "$final_status" == COMPLETE ]] || { echo "d50 remains incomplete: $method=$final_status" >&2; exit 7; }
  done
fi

exec /bin/bash "$FULL" "$trigger_job"
