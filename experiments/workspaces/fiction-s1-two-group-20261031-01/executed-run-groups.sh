#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/home/steven/mjx/fiction/experiments/workspaces/fiction-s1-two-group-20261031-01
TASK_BIN=/tmp/fiction-s1-two-group-build/tension_check
TASK_MESH=/home/steven/mjx/fiction/simulator/output/validation/s1-20261030-01/meshes/M2.msh
for group in control iteration300 penalty1e10; do
    case "$group" in
        control) iterations=100; stiffness=1e8;;
        iteration300) iterations=300; stiffness=1e8;;
        penalty1e10) iterations=100; stiffness=1e10;;
    esac
    destination="$TASK_ROOT/results/$group"
    mkdir -p "$destination"
    set +e
    /usr/bin/time -v timeout 120s "$TASK_BIN" --mesh "$TASK_MESH" --csv "$destination/metrics.csv" --final-max-iter "$iterations" --final-dbc-stiffness "$stiffness" --trace-dir "$destination" > "$destination/stdout.txt" 2> "$destination/time.txt"
    status=$?
    set -e
    echo "$status" > "$destination/exitcode"
    if [ "$status" -ne 0 ] && [ "$status" -ne 1 ]; then
        echo "Unexpected error in $group: $status" >&2
        exit "$status"
    fi
    echo "$group exit=$status"
done
