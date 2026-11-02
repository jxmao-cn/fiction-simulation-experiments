#!/usr/bin/env bash
# 停止判据诊断运行脚本（隔离副本，不改生产）
#   病例：S2 M3 静态 + S6a 1 体 10 步
#   判据：step(基线) / grad / actual / combined，共 8 个 run，顺序执行（并发墙钟不作性能结论）
# 用法: ./run_stopcriterion.sh
set -uo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

RUN_ID="stopcriterion-20261031-01"
WORK="/home/steven/mjx/fiction/experiments/workspaces/fiction-stopcriterion-20261031-01"
BUILD="/tmp/fiction-stopcriterion-build"
S2_MESH="/home/steven/mjx/fiction/simulator/output/validation/s2-cantilever-20261031-01/meshes/M3.msh"
S6A_MESH="/home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh"
OUT="/home/steven/mjx/fiction/simulator/output/validation/${RUN_ID}"
S2_BIN="${BUILD}/cantilever_check"
S6A_BIN="${BUILD}/s6a_multibody_check"

# 力尺度（预登记）：S2 总重 24.525 N；S6a 1 体总重 = 518.4769746175648 * 9.81 N
S2_GRAD_SCALE="24.525"
S6A_GRAD_SCALE="$(python3 -c 'print(repr(518.4769746175648*9.81))')"

mkdir -p "${OUT}"
cd "${WORK}" || exit 1

# ---------------- 现场记录 ----------------
{
  echo "run_id=${RUN_ID}"
  echo "date=$(date -Is)"
  echo "host=$(hostname)"
  echo "work=${WORK}"
  echo "build=${BUILD}"
  echo "build_type=Release; CHOLMOD=ON (SuiteSparse /usr); FICTION_USE_VIEWER=ON"
  echo "s2_mesh=${S2_MESH}"
  echo "s2_grad_scale=${S2_GRAD_SCALE}"
  echo "s6a_mesh=${S6A_MESH}"
  echo "s6a_grad_scale=${S6A_GRAD_SCALE}"
  echo "thread_env=OMP/OPENBLAS/MKL/BLIS/VECLIB=1"
  echo "note=8 个 run 顺序执行；墙钟只作粗参照，不作性能结论"
} > "${OUT}/env.txt"
sha256sum "${S2_MESH}" "${S6A_MESH}" > "${OUT}/asset_sha256.txt"
sha256sum "${S2_BIN}" "${S6A_BIN}" > "${OUT}/binary_sha256.txt"
sha256sum cpp/src/time_integrator.cpp cpp/include/time_integrator.h \
         cpp/examples/cantilever_check.cpp cpp/examples/s6a_multibody_check.cpp \
         cpp/include/scene/multibody_scene.h cpp/CMakeLists.txt \
         > "${OUT}/source_files_sha256.txt"
tar czf "${OUT}/source.tar.gz" --exclude='*.msh' cpp run_stopcriterion.sh tools 2>/dev/null
sha256sum "${OUT}/source.tar.gz" > "${OUT}/source_hash.txt"

# ---------------- S2 M3 静态：4 个判据 ----------------
run_case() {
  local case_dir="$1"; shift
  mkdir -p "${case_dir}"
  { echo "cd ${PWD} && $*"; } > "${case_dir}/command.txt"
  /usr/bin/time -v "$@" > "${case_dir}/run.log" 2> "${case_dir}/time_rss.txt"
  echo $? > "${case_dir}/exit_code.txt"
}

S2_DIR="${OUT}/s2/M3"
for mode in step grad actual combined; do
  d="${S2_DIR}/${mode}"
  args=(--mesh "${S2_MESH}" --csv "${d}/M3.csv" --state "${d}/M3-state.csv"
        --stop-mode "${mode}" --grad-tol 1e-7 --grad-scale "${S2_GRAD_SCALE}"
        --actual-step-tol 1e-12 --actual-step-count 2)
  run_case "${d}" "${S2_BIN}" "${args[@]}"
done

# ---------------- S6a 1 体 10 步：4 个判据 ----------------
S6A_DIR="${OUT}/s6a/body1"
for mode in step grad actual combined; do
  d="${S6A_DIR}/${mode}"
  args=(--mesh "${S6A_MESH}" --out "${d}" --only-body1
        --stop-mode "${mode}" --grad-tol 1e-7 --grad-scale "${S6A_GRAD_SCALE}"
        --actual-step-tol 1e-12 --actual-step-count 2)
  run_case "${d}" "${S6A_BIN}" "${args[@]}"
done

# ---------------- 逐 case config.json + 汇总 + 判定 ----------------
python3 "${WORK}/tools/summarize_stopcriterion.py" \
  --out "${OUT}" --s2-grad-scale "${S2_GRAD_SCALE}" --s6a-grad-scale "${S6A_GRAD_SCALE}"
echo "summarize exit=$?"
exit 0
