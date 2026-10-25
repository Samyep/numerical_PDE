#!/bin/bash
set -uo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
TASK_FILE=${1:?usage: run_whole_function_queue_worker.sh TASK_FILE SLOT}
SLOT=${2:?usage: run_whole_function_queue_worker.sh TASK_FILE SLOT}
HELPER="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh"

allocation_seconds_left() {
  python - "$SLURM_JOB_ID" <<'PY'
from datetime import datetime
import re
import subprocess
import sys
info = subprocess.check_output(["scontrol", "show", "job", "-o", sys.argv[1]], text=True)
match = re.search(r"\bEndTime=(\S+)", info)
if not match or match.group(1) in {"Unknown", "N/A"}:
    raise SystemExit("Slurm did not expose a finite allocation end time")
end = datetime.fromisoformat(match.group(1)).timestamp()
print(max(0, int(end - datetime.now().timestamp())))
PY
}

module load default
source /sw/user/python/miniforge3-pytorch-2.10.0/etc/profile.d/conda.sh
conda activate base
cd "$SOURCE_ROOT"
export PYTHONPATH="$SOURCE_ROOT${PYTHONPATH:+:$PYTHONPATH}"

tail -n +2 "$TASK_FILE" | while IFS=',' read -r slot sequence method pde dimension estimate high_estimate extra; do
  [[ "$slot" == "$SLOT" ]] || continue
  if [[ -n "${extra:-}" ]]; then
    echo "invalid queue row: slot=$slot sequence=$sequence" >&2
    continue
  fi
  read -r prior_status cell_remaining prior_used prior_attempts < <(
    python -m invariant_region_mlp.experiments.tangent_actor_icml_final_suite.whole_function_status budget \
      --method "$method" --pde "$pde" --dimension "$dimension"
  )
  case "$prior_status" in
    COMPLETE|OUT_OF_TIME|OUT_OF_MEMORY|FAILED)
      echo "QUEUE_SKIP slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} status=$prior_status"
      continue ;;
  esac
  allocation_left=$(allocation_seconds_left) || exit 4
  available=$((allocation_left - 240))
  if (( available < 300 )); then
    echo "QUEUE_STOP slot=$SLOT reason=allocation_margin seconds_left=$allocation_left"
    break
  fi
  cap=$cell_remaining
  (( cap > available )) && cap=$available
  # PINN-HTE remains a single nonresumable training call.  DPI is resumable at
  # exact Picard-round boundaries: its wrapper atomically commits the model,
  # RNG stream, and dataset-size state after every completed round.
  if [[ "$method" == "pinn_hte" ]]; then
    required=$(python - "$high_estimate" <<'PY'
import math
import sys
print(min(21600, math.ceil(float(sys.argv[1]) * 3600) + 60))
PY
)
    if [[ "$prior_status" == "CHECKPOINTED" || $required -gt $available ]]; then
      echo "QUEUE_DEFER slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} reason=nonresumable required=$required available=$available"
      continue
    fi
  fi
  export WHOLE_CELL_TASK_ID="${SLOT}_${sequence}"
  export CELL_CAP_SECONDS=$cap
  echo "QUEUE_BEGIN slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} cap=$cap estimate_hours=$estimate high_hours=$high_estimate"
  "$HELPER" "$method" "$pde" "$dimension"
  code=$?
  echo "QUEUE_END slot=$SLOT seq=$sequence cell=${method}/${pde}/d${dimension} exit=$code"
done
