#!/usr/bin/env bash
# v0.3 停止判据候选确认：combined 在 S2 M3/M6 + S6a 1/4 体上的重验收
#   性质：候选确认（isolated copy），不构成协议通过；生产 fiction/simulator/cpp 不动。
#   判据：--stop-mode combined --grad-tol 1e-6 --actual-step-tol 1e-12 --actual-step-count 2
#   4 个 case 顺序执行（并发墙钟不作性能结论）；M6 单 case 上限 600 s。
# 用法: ./run_combined.sh
set -uo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

RUN_ID="s2s6a-combined-20261031-01"
WORK="/home/steven/mjx/fiction/experiments/workspaces/fiction-s2s6a-combined-20261031-01"
BASE_COPY="/home/steven/mjx/fiction/experiments/workspaces/fiction-stopcriterion-20261031-01"
ARCHIVE="/home/steven/mjx/fiction/simulator/output/validation"
BUILD="/tmp/fiction-s2s6a-combined-build"
S2_M3_MESH="${ARCHIVE}/s2-cantilever-20261031-01/meshes/M3.msh"
S2_M6_MESH="${ARCHIVE}/s2-extension-20261031-01/meshes/M6.msh"
S6A_MESH="/home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh"
OUT="${ARCHIVE}/${RUN_ID}"
S2_BIN="${BUILD}/cantilever_check"
S6A_BIN="${BUILD}/s6a_multibody_check"

S2_GRAD_SCALE="24.525"
S6A_BODY1_SCALE="5086.259120998311"
S6A_BODY4_SCALE="20345.036483993244"
COMMON_ARGS=(--stop-mode combined --grad-tol 1e-6 --actual-step-tol 1e-12 --actual-step-count 2)

mkdir -p "${OUT}"
cd "${WORK}" || exit 1

# ---------------- 现场记录 ----------------
{
  echo "run_id=${RUN_ID}"
  echo "date=$(date -Is)"
  echo "host=$(hostname)"
  echo "work=${WORK}"
  echo "base_copy=${BASE_COPY}"
  echo "build=${BUILD}"
  echo "build_type=Release; CHOLMOD=ON (SuiteSparse /usr); FICTION_USE_VIEWER=ON"
  echo "s2_m3_mesh=${S2_M3_MESH}"
  echo "s2_m6_mesh=${S2_M6_MESH}"
  echo "s6a_mesh=${S6A_MESH}"
  echo "s2_grad_scale=${S2_GRAD_SCALE}"
  echo "s6a_body1_grad_scale=${S6A_BODY1_SCALE}"
  echo "s6a_body4_grad_scale=${S6A_BODY4_SCALE}"
  echo "stop_args=--stop-mode combined --grad-tol 1e-6 --actual-step-tol 1e-12 --actual-step-count 2"
  echo "thread_env=OMP/OPENBLAS/MKL/BLIS/VECLIB=1"
  echo "order=M3 -> s6a/body1 -> s6a/body4 -> M6(600s cap)；顺序执行，墙钟只作粗参照"
} > "${OUT}/env.txt"
sha256sum "${S2_M3_MESH}" "${S2_M6_MESH}" "${S6A_MESH}" > "${OUT}/asset_sha256.txt"
sha256sum "${S2_BIN}" "${S6A_BIN}" > "${OUT}/binary_sha256.txt"
cp /tmp/fiction-s2s6a-combined-configure.log "${OUT}/configure.log"
cp /tmp/fiction-s2s6a-combined-build.log "${OUT}/build.log"
cp /tmp/fiction-s2s6a-combined-ctest.log "${OUT}/ctest.log"
sha256sum cpp/src/time_integrator.cpp cpp/include/time_integrator.h \
         cpp/examples/cantilever_check.cpp cpp/examples/s6a_multibody_check.cpp \
         cpp/include/scene/multibody_scene.h cpp/CMakeLists.txt \
         > "${OUT}/source_files_sha256.txt"
tar czf "${OUT}/source.tar.gz" --exclude='*.msh' cpp run_combined.sh tools 2>/dev/null
sha256sum "${OUT}/source.tar.gz" > "${OUT}/source_hash.txt"
# 实现 diff：相对本候选的基（stopcriterion 隔离副本），只含源码
( cd /home/steven/mjx && diff -ruN \
    --exclude='*.msh' --exclude='build' \
    "${BASE_COPY}/cpp" "${WORK}/cpp" > "${OUT}/source-change.diff" )
echo "source-change.diff lines=$(wc -l < "${OUT}/source-change.diff")"

run_case() {
  local case_dir="$1"; shift
  mkdir -p "${case_dir}"
  { echo "cd ${PWD} && $*"; } > "${case_dir}/command.txt"
  /usr/bin/time -v "$@" > "${case_dir}/run.log" 2> "${case_dir}/time_rss.txt"
  echo $? > "${case_dir}/exit_code.txt"
}

# ---------------- 1) S2 M3 静态（cantilever_check 硬编码 tol=1e-10, maxIter=100） ----------------
d="${OUT}/s2/M3"
run_case "${d}" "${S2_BIN}" --mesh "${S2_M3_MESH}" \
  --csv "${d}/M3.csv" --state "${d}/M3-state.csv" --trace-dir "${d}/trace" \
  "${COMMON_ARGS[@]}" --grad-scale "${S2_GRAD_SCALE}"
echo "M3 exit=$(cat "${d}/exit_code.txt")"

# ---------------- 2) S6a 1 体（tol=1e-8, maxIter=100, 10 步） ----------------
d="${OUT}/s6a/body1"
run_case "${d}" "${S6A_BIN}" --mesh "${S6A_MESH}" --out "${d}" --only-body1 \
  --run-id "${RUN_ID}" "${COMMON_ARGS[@]}" --grad-scale "${S6A_BODY1_SCALE}"
echo "s6a/body1 exit=$(cat "${d}/exit_code.txt")"

# ---------------- 3) S6a 4 体（tol=1e-8, maxIter=100, 10 步；含块对角 + 4 次独立单体） ----------------
d="${OUT}/s6a/body4"
run_case "${d}" "${S6A_BIN}" --mesh "${S6A_MESH}" --out "${d}" \
  --run-id "${RUN_ID}" "${COMMON_ARGS[@]}" --grad-scale "${S6A_BODY4_SCALE}"
echo "s6a/body4 exit=$(cat "${d}/exit_code.txt")"
# 程序自身的 summary.json 是原始证据，改名保留；汇总脚本另写 case summary.json
if [ -f "${d}/summary.json" ]; then mv "${d}/summary.json" "${d}/program_summary.json"; fi

# ---------------- 4) S2 M6 静态（同一 cantilever_check，硬编码 tol=1e-10；上限 600 s） ----------------
d="${OUT}/s2/M6"
mkdir -p "${d}"
{ echo "cd ${PWD} && timeout 600 ${S2_BIN} --mesh ${S2_M6_MESH} --csv ${d}/M6.csv --state ${d}/M6-state.csv --trace-dir ${d}/trace ${COMMON_ARGS[*]} --grad-scale ${S2_GRAD_SCALE}"; } > "${d}/command.txt"
/usr/bin/time -v timeout 600 "${S2_BIN}" --mesh "${S2_M6_MESH}" \
  --csv "${d}/M6.csv" --state "${d}/M6-state.csv" --trace-dir "${d}/trace" \
  "${COMMON_ARGS[@]}" --grad-scale "${S2_GRAD_SCALE}" \
  > "${d}/run.log" 2> "${d}/time_rss.txt"
echo $? > "${d}/exit_code.txt"
echo "M6 exit=$(cat "${d}/exit_code.txt")"

# ---------------- 逐 case config.json + 汇总 + 预登记判定 ----------------
python3 "${WORK}/tools/summarize_combined.py" --out "${OUT}"
echo "summarize exit=$?"
exit 0
