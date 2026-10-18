#!/bin/bash
# Continue the corrected 100D/fidelity-anchor work under the interactive
# partition's one-submitted-job-per-user limit. Run this lightweight poller on
# the login host after submitting the first four-GPU pack.

set -euo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
RESULT_ROOT=$SOURCE_ROOT/invariant_region_mlp/results/tangent_actor_icml_final_suite
RAW_DIR=$RESULT_ROOT/raw
PACK_SCRIPT=$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/parallel_baseline_pack.slurm
PYTHON=/sw/user/python/miniforge3-pytorch-2.10.0/bin/python
PREDECESSOR_JOB_ID=${1:?usage: auto_stage_b_controller.sh PREDECESSOR_JOB_ID}
MAX_BASELINE_RETRIES=${MAX_BASELINE_RETRIES:-2}
MAX_SDGD_CYCLES=${MAX_SDGD_CYCLES:-3}

BASELINE_PACK='pinn_hte:joint_convex:100:auto+dpi:nonconvex_coupled:100:auto,pinn_hte:nonconvex_coupled:100:auto+dpi:joint_concave:100:auto+dpi:allen_cahn:100:auto,pinn_hte:mlp_failure:100:auto+dpi:joint_convex:100:auto+dpi:hjb_log:100:auto,pinn_hte:allen_cahn:100:auto+dpi:mlp_failure:100:auto+sdgd:hjb_log:250:hjb_adversarial'
SDGD_PACK='sdgd:hjb_log:250:hjb_adversarial+sdgd:allen_cahn:100:auto,sdgd:hjb_log:100:auto+sdgd:joint_concave:100:auto,sdgd:ou_gmm:100:auto+sdgd:joint_convex:100:auto,sdgd:nonconvex_coupled:100:auto+sdgd:mlp_failure:100:auto'

BASELINE_RESULTS=(
  pinn_hte_hjb_log_d100_semilinear_paper_official_s0.json
  pinn_hte_joint_concave_d100_semilinear_paper_official_s0.json
  pinn_hte_joint_convex_d100_semilinear_paper_official_s0.json
  pinn_hte_nonconvex_coupled_d100_semilinear_paper_official_s0.json
  pinn_hte_mlp_failure_d100_semilinear_paper_official_s0.json
  pinn_hte_allen_cahn_d100_semilinear_paper_official_s0.json
  dpi_hjb_log_d100_semilinear_paper_official_s0.json
  dpi_joint_concave_d100_semilinear_paper_official_s0.json
  dpi_joint_convex_d100_semilinear_paper_official_s0.json
  dpi_nonconvex_coupled_d100_semilinear_paper_official_s0.json
  dpi_mlp_failure_d100_semilinear_paper_official_s0.json
  dpi_allen_cahn_d100_semilinear_paper_official_s0.json
)

SDGD_RESULTS=(
  sdgd_pinn_hjb_log_d250_hjb_adversarial_official_s0.json
  sdgd_pinn_hjb_log_d100_hjb_adversarial_official_s0.json
  sdgd_pinn_ou_gmm_d100_hjb_adversarial_official_s0.json
  sdgd_pinn_joint_concave_d100_semilinear_sde_official_s0.json
  sdgd_pinn_joint_convex_d100_semilinear_sde_official_s0.json
  sdgd_pinn_nonconvex_coupled_d100_semilinear_sde_official_s0.json
  sdgd_pinn_mlp_failure_d100_semilinear_sde_official_s0.json
  sdgd_pinn_allen_cahn_d100_semilinear_sde_official_s0.json
)

timestamp() { date -u +%FT%TZ; }

wait_for_job() {
  local job_id=$1
  while squeue -h -j "$job_id" | grep -q .; do
    sleep 30
  done
  local state=""
  for _ in $(seq 1 10); do
    state=$(sacct -n -X -j "$job_id" --format=State -P | sed '/^$/d' | head -n1)
    [[ -n "$state" ]] && break
    sleep 2
  done
  echo "$(timestamp) job_id=$job_id state=$state" >&2
  case "$state" in
    COMPLETED*|TIMEOUT*) ;;
    *) echo "Refusing continuation after unexpected state $state" >&2; return 3 ;;
  esac
}

results_complete() {
  "$PYTHON" - "$RAW_DIR" "$@" <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
missing = []
for name in sys.argv[2:]:
    path = root / name
    try:
        complete = json.loads(path.read_text()).get("status") == "complete"
    except (FileNotFoundError, json.JSONDecodeError):
        complete = False
    if not complete:
        missing.append(name)
if missing:
    print("incomplete=" + ",".join(missing), file=sys.stderr)
    raise SystemExit(1)
PY
}

wait_for_interactive_slot() {
  while squeue -h -u "$USER" -p ghx4-interactive | grep -q .; do
    sleep 30
  done
}

submit_pack() {
  local job_name=$1
  local pack_tasks=$2
  wait_for_interactive_slot
  if squeue -h -u "$USER" --name="$job_name" | grep -q .; then
    echo "Duplicate job name already queued: $job_name" >&2
    return 4
  fi
  env PACK_TASKS="$pack_tasks" sbatch --test-only \
    --job-name="$job_name" --nodes=2 --ntasks=4 --gpus-per-node=2 --mem=220G \
    "$PACK_SCRIPT" >&2
  env PACK_TASKS="$pack_tasks" sbatch --parsable \
    --job-name="$job_name" --nodes=2 --ntasks=4 --gpus-per-node=2 --mem=220G \
    "$PACK_SCRIPT"
}

echo "$(timestamp) controller_start predecessor=$PREDECESSOR_JOB_ID"
wait_for_job "$PREDECESSOR_JOB_ID"

baseline_retry=0
while ! results_complete "${BASELINE_RESULTS[@]}"; do
  if (( baseline_retry >= MAX_BASELINE_RETRIES )); then
    echo "Baseline retry cap reached" >&2
    exit 5
  fi
  baseline_retry=$((baseline_retry + 1))
  job_id=$(submit_pack "icml_base_retry${baseline_retry}" "$BASELINE_PACK")
  echo "$(timestamp) submitted_baseline_retry=$baseline_retry job_id=$job_id"
  wait_for_job "$job_id"
done
echo "$(timestamp) corrected_baselines_complete"

sdgd_cycle=0
while ! results_complete "${SDGD_RESULTS[@]}"; do
  if (( sdgd_cycle >= MAX_SDGD_CYCLES )); then
    echo "SDGD cycle cap reached" >&2
    exit 6
  fi
  sdgd_cycle=$((sdgd_cycle + 1))
  job_id=$(submit_pack "icml_sdgd_pack${sdgd_cycle}" "$SDGD_PACK")
  echo "$(timestamp) submitted_sdgd_cycle=$sdgd_cycle job_id=$job_id"
  wait_for_job "$job_id"
done
echo "$(timestamp) corrected_stage_b_and_sdgd_anchor_complete"
