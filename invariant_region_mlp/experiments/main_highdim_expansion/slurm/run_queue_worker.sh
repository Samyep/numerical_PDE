#!/bin/bash
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/common_env.sh"
TASK_FILE=${1:?usage: run_queue_worker.sh TASK_FILE SLOT}
SLOT=${2:?usage: run_queue_worker.sh TASK_FILE SLOT}
HELPER="$SOURCE_ROOT/invariant_region_mlp/experiments/main_highdim_expansion/slurm/run_cell.sh"

allocation_seconds_left() {
  python - "$SLURM_JOB_ID" <<'PY'
from datetime import datetime
import re, subprocess, sys
info=subprocess.check_output(["scontrol","show","job","-o",sys.argv[1]],text=True)
match=re.search(r"\bEndTime=(\S+)",info)
if not match or match.group(1) in {"Unknown","N/A"}: raise SystemExit(4)
print(max(0,int(datetime.fromisoformat(match.group(1)).timestamp()-datetime.now().timestamp())))
PY
}

tail -n +2 "$TASK_FILE" | while IFS=',' read -r slot sequence method pde dimension estimate high_estimate extra; do
  [[ "$slot" == "$SLOT" ]] || continue
  [[ -z "${extra:-}" ]] || { echo "invalid queue row" >&2; continue; }
  read -r prior_status cell_remaining prior_used prior_attempts < <(
    python -m invariant_region_mlp.experiments.main_highdim_expansion.status budget \
      --method "$method" --pde "$pde" --dimension "$dimension"
  )
  case "$prior_status" in COMPLETE|OUT_OF_TIME|OUT_OF_MEMORY|FAILED|BUDGET_TRUNCATED) continue;; esac
  allocation_left=$(allocation_seconds_left) || exit 4
  available=$((allocation_left - 240))
  (( available >= 300 )) || break
  cap=$cell_remaining
  (( cap > available )) && cap=$available
  if [[ "$method" == "dpi" || "$method" == "pinn_hte" ]]; then
    required=$(python - "$high_estimate" <<'PY'
import math,sys
print(min(21600,math.ceil(float(sys.argv[1])*3600)))
PY
)
    if (( required > available || required > cell_remaining )); then
      echo "QUEUE_DEFER slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} required=$required available=$available remaining=$cell_remaining"
      continue
    fi
  fi
  export EXPANSION_CELL_TASK_ID="${SLOT}_${sequence}"
  export CELL_CAP_SECONDS=$cap
  echo "QUEUE_BEGIN slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} cap=$cap estimate=$estimate high=$high_estimate"
  "$HELPER" "$method" "$pde" "$dimension"
  code=$?
  echo "QUEUE_END slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} exit=$code"
done
