#!/bin/bash
set -uo pipefail

SOURCE_ROOT=/work/nvme/bgvf/ywang151/main_b1_sixmethods_20261012/source
DIMENSION=${1:?usage: run_stage1_fast_lane.sh DIMENSION}
HELPER="$SOURCE_ROOT/invariant_region_mlp/experiments/tangent_actor_icml_final_suite/slurm/run_whole_function_cell.sh"
case "$DIMENSION" in 20|50) ;; *) echo "invalid Stage-1 dimension" >&2; exit 2;; esac

# Frozen from the pre-run upper runtime estimates. Nonresumable DPI/PINN run
# first; the resumable Actor receives the remaining conservative allowance.
methods=(dpi pinn_hte tangent_actor)
caps=(1600 850 2500)
status=0
for index in 0 1 2; do
  method=${methods[$index]}
  cap=${caps[$index]}
  echo "FAST_LANE_BEGIN method=$method dimension=$DIMENSION cap_seconds=$cap"
  if env WHOLE_CELL_TASK_ID="fast_${index}" CELL_CAP_SECONDS="$cap" \
    "$HELPER" "$method" joint_concave "$DIMENSION"; then
    code=0
  else
    code=$?
    status=1
  fi
  echo "FAST_LANE_END method=$method dimension=$DIMENSION exit=$code"
done
exit "$status"
