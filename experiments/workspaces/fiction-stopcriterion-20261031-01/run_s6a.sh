#!/usr/bin/env bash
# S6a 运行脚本：只读资产 + 1 体/4 体全局场景 + 4 次独立单体对照 + 停滞诊断
# 用法: ./run_s6a.sh
set -uo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

RUN_ID="s6a-multibody-20261031-01"
WORK="/home/steven/mjx/fiction/experiments/workspaces/fiction-s6a-20261031-01"
BUILD="/tmp/fiction-s6a-build"
MESH="/home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh"
OUT="/home/steven/mjx/fiction/simulator/output/validation/${RUN_ID}"
BIN="${BUILD}/s6a_multibody_check"

mkdir -p "${OUT}"
cd "${WORK}" || exit 1

# --- 现场记录 ---
{
  echo "run_id=${RUN_ID}"
  echo "date=$(date -Is)"
  echo "host=$(hostname)"
  echo "mesh=${MESH}"
  echo "binary=${BIN}"
  echo "build_type=Release (cmake -DCMAKE_BUILD_TYPE=Release, CHOLMOD path)"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS} MKL_NUM_THREADS=${MKL_NUM_THREADS} BLIS_NUM_THREADS=${BLIS_NUM_THREADS}"
  echo "note=线程环境变量已设置，未验证 CHOLMOD/BLAS 内部是否全部遵守"
} > "${OUT}/env.txt"

sha256sum "${MESH}" > "${OUT}/asset_sha256.txt"
sha256sum "${BIN}" > "${OUT}/binary_sha256.txt"
tar czf "${OUT}/source.tar.gz" cpp tools
sha256sum "${OUT}/source.tar.gz" > "${OUT}/source_hash.txt"
sha256sum cpp/include/scene/multibody_scene.h \
         cpp/examples/s6a_multibody_check.cpp \
         cpp/src/time_integrator.cpp \
         cpp/CMakeLists.txt > "${OUT}/source_files_sha256.txt"

cat > "${OUT}/config.json" <<'JSON'
{
  "run_id": "s6a-multibody-20261031-01",
  "stage": "S6a",
  "scope": "无接触多体全局场景组织；1 体和 4 体",
  "asset": "Archive/IPC/input/tetMeshes/sphere1K.msh",
  "instances": [1, 4],
  "layout": "4 体沿 x/y 成 2x2，间距 1.5 m（> 1.2 x 直径 1.0 m），纯平移，不旋转、不缩放",
  "material": {"model": "NeoHookean tet (P1)", "youngs_pa": 100000.0, "poisson": 0.3, "density_kg_m3": 1000.0},
  "gravity_m_s2": [0.0, 0.0, -9.81],
  "dbc": {"region": "每体局部 z >= 0.45 的顶部 cap", "stiffness_initial": 1e8, "target": "参考位置"},
  "time": {"scheme": "implicit Euler", "dt_s": 0.005, "steps": 10},
  "newton": {"tol": 1e-8, "max_iter": 100},
  "contacts": false,
  "gates": {
    "per_step_converged": true,
    "finite": true,
    "min_J_gt": 0.3,
    "dbc_error_lt_m": 1e-6,
    "block_diagonal_cross_body_entries": 0,
    "independent_compare_max_position_diff_lt_m": 1e-8
  }
}
JSON

echo "cd ${WORK} && ${BIN} --mesh ${MESH} --out ${OUT}" > "${OUT}/command.txt"

# --- 主跑（无 trace，计时干净） ---
/usr/bin/time -v "${BIN}" --mesh "${MESH}" --out "${OUT}" \
  > "${OUT}/run.log" 2> "${OUT}/time_rss.txt"
echo $? > "${OUT}/exit_code.txt"

# --- 诊断 1：1 体开 trace 重跑（时间戳逐迭代轨迹；注意 trace 文件每步被截断，仅留末步） ---
mkdir -p "${OUT}/diagnose"
/usr/bin/time -v "${BIN}" --mesh "${MESH}" --out "${OUT}/diagnose" \
  --trace-dir "${OUT}/diagnose" \
  > "${OUT}/diagnose/run.log" 2> "${OUT}/diagnose/time_rss.txt"
echo $? > "${OUT}/diagnose/exit_code.txt"

# --- 诊断 2：第 6 步（1 体停滞步）迭代预算扫描 ---
echo "sweep of --max-iter at step 6 (1 body, diagnostic only)" > "${OUT}/diagnose/step6_sweep.txt"
for mi in 2 3 5 10 100; do
  "${BIN}" --mesh "${MESH}" --out "${OUT}/diagnose/step6_maxiter_${mi}" \
    --steps 6 --max-iter "${mi}" --only-body1 \
    > "${OUT}/diagnose/step6_maxiter_${mi}.log" 2>&1
  echo "max_iter=${mi} exit=$?" >> "${OUT}/diagnose/step6_sweep.txt"
done

# --- 诊断 3：4 体场景 maxIter=3 对照（看停滞步是否只增加迭代不改变状态） ---
"${BIN}" --mesh "${MESH}" --out "${OUT}/diagnose4b_maxiter3" --max-iter 3 \
  > "${OUT}/diagnose4b_maxiter3.log" 2>&1
echo $? > "${OUT}/diagnose4b_maxiter3_exit_code.txt"

exit 0
