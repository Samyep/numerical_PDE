#!/bin/bash
set -euo pipefail

REGULAR_JOB=${1:?usage: race_existing_jobs.sh REGULAR_JOB INTERACTIVE_JOB}
INTERACTIVE_JOB=${2:?usage: race_existing_jobs.sh REGULAR_JOB INTERACTIVE_JOB}
REGULAR_PREFERENCE_SECONDS=${REGULAR_PREFERENCE_SECONDS:-7200}

timestamp() { date -Is; }

queue_state() {
  squeue -h -j "$1" -o '%T' 2>/dev/null | head -1
}

is_active() {
  [[ "$1" == "RUNNING" || "$1" == "CONFIGURING" || "$1" == "COMPLETING" ]]
}

regular_due_within_window() {
  local estimate estimate_epoch now
  estimate=$(squeue -h -j "$REGULAR_JOB" --start -o '%S' 2>/dev/null | head -1)
  [[ -n "$estimate" && "$estimate" != "N/A" && "$estimate" != "Unknown" ]] || return 1
  estimate_epoch=$(date -d "$estimate" +%s 2>/dev/null) || return 1
  now=$(date +%s)
  (( estimate_epoch - now <= REGULAR_PREFERENCE_SECONDS ))
}

echo "$(timestamp) RACE_START regular=$REGULAR_JOB interactive=$INTERACTIVE_JOB regular_preference_seconds=$REGULAR_PREFERENCE_SECONDS"
while true; do
  regular_state=$(queue_state "$REGULAR_JOB")
  interactive_state=$(queue_state "$INTERACTIVE_JOB")
  echo "$(timestamp) RACE interactive=$INTERACTIVE_JOB:${interactive_state:-TERMINAL} regular=$REGULAR_JOB:${regular_state:-TERMINAL}"

  if is_active "$interactive_state"; then
    [[ -z "$regular_state" ]] || scancel "$REGULAR_JOB"
    echo "$(timestamp) RACE_WINNER job=$INTERACTIVE_JOB kind=interactive"
    exit 0
  fi
  if is_active "$regular_state"; then
    [[ -z "$interactive_state" ]] || scancel "$INTERACTIVE_JOB"
    echo "$(timestamp) RACE_WINNER job=$REGULAR_JOB kind=regular"
    exit 0
  fi
  if [[ "$regular_state" == "PENDING" && "$interactive_state" == "PENDING" ]] \
      && regular_due_within_window; then
    scancel "$INTERACTIVE_JOB"
    echo "$(timestamp) RACE_PREFERENCE job=$REGULAR_JOB kind=regular reason=start_within_2h action=cancel_interactive"
    exit 0
  fi
  if [[ -z "$interactive_state" && -n "$regular_state" ]]; then
    echo "$(timestamp) RACE_REMAINDER job=$REGULAR_JOB kind=regular"
    exit 0
  fi
  if [[ -z "$regular_state" && -n "$interactive_state" ]]; then
    echo "$(timestamp) RACE_REMAINDER job=$INTERACTIVE_JOB kind=interactive"
    exit 0
  fi
  if [[ -z "$regular_state" && -z "$interactive_state" ]]; then
    echo "$(timestamp) RACE_FAILED reason=both_jobs_left_queue" >&2
    exit 5
  fi
  sleep 15
done
