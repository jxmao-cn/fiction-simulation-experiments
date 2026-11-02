#!/usr/bin/env bash
set -euo pipefail
TASK_ROOT=/home/steven/mjx/fiction/experiments/workspaces/fiction-s1-diagnostic-20261031-01
MESH_ROOT=/home/steven/mjx/fiction/simulator/output/validation/s1-20261030-01/meshes
# 复现输出另存，避免覆盖此次冻结证据；该脚本不会修改原工程。
RUN_DIR="$TASK_ROOT/reproductions/$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$RUN_DIR"
for variant in baseline variant-a; do
    mkdir -p "$RUN_DIR/$variant"
    cmake -S "$TASK_ROOT/$variant/cpp" -B "$TASK_ROOT/$variant/build" -DCMAKE_BUILD_TYPE=Release > "$RUN_DIR/$variant/configure.log" 2>&1
    cmake --build "$TASK_ROOT/$variant/build" -j3 > "$RUN_DIR/$variant/build.log" 2>&1
    ctest --test-dir "$TASK_ROOT/$variant/build" --output-on-failure > "$RUN_DIR/$variant/ctest.log" 2>&1
done
run_s1() {
    local variant="$1" mesh="$2" status
    set +e
    /usr/bin/time -v timeout 120s "$TASK_ROOT/$variant/build/tension_check" --mesh "$MESH_ROOT/$mesh.msh" --csv "$RUN_DIR/$variant/$mesh.csv" > "$RUN_DIR/$variant/$mesh.stdout.txt" 2> "$RUN_DIR/$variant/$mesh.time.txt"
    status=$?
    set -e
    echo "$status" > "$RUN_DIR/$variant/$mesh.exitcode"
    return "$status"
}
# baseline 预期 exit 1；捕获它，避免 set -e 使后续对照无法执行。
if run_s1 baseline M2; then
    echo 'baseline unexpectedly passed; inspect results before continuing' >&2
    exit 1
else
    status=$?
    if [ "$status" -ne 1 ]; then exit "$status"; fi
fi
run_s1 variant-a M1
if run_s1 variant-a M2; then
    echo 'variant-a M2 passed; this script still does not start M3 automatically'
else
    status=$?
    if [ "$status" -ne 1 ]; then exit "$status"; fi
    echo 'variant-a M2 failed acceptance; stopped without M3'
fi
echo "Results: $RUN_DIR"
