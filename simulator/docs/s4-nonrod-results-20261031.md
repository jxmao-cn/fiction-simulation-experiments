# S4 非杆几何验证结果（torus，2026-10-31）

- run-id：`s4-nonrod-20261031-01`
- 预登记：`Archive/fiction/docs/s3-s6a-execution-design-20261031.md` 第 3 节
- 隔离工作区：`/home/steven/mjx/fiction/experiments/workspaces/fiction-s4-20261031-01/`（复制自 `Archive/fiction/cpp` 与 `Archive/fiction/tools`，生产目录未改动）
- 构建目录：`/tmp/fiction-s4-build`（Release；CHOLMOD 求解路径）
- 原始结果根：`Archive/fiction/output/validation/s4-nonrod-20261031-01/`
- 精选图/短片：`personal-homepage/outputs/physics-cpp-blog/validation/s4-nonrod-20261031-01/`

## 1. 结论摘要

同一套 Neo-Hookean tet FEM + 隐式欧拉 / projected Newton 框架成功加载并求解了一个**非杆、带孔的体网格**（torus，200 节点 / 491 tet）：静态满载收敛，2 s 动态 400 步全部收敛并通过全部安全门槛（有限、`min J > 0.3`、DBC 误差 `< 1e-6`）。本阶段检查非杆体网格的求解状态与安全门槛。

| 项目 | 结果 |
|---|---|
| 资产 | `Archive/IPC/input/tetMeshes/torus.msh`（首选资产，未启用备选 `cow_head.msh`） |
| 规模 | 200 节点 / 491 四面体 / 400 表面三角面 / 191 自由节点 / 9 固定节点 |
| 静态满载 | 收敛，13 次 Newton，`min J = 0.871477`，最大位移 0.196671 m，墙钟 0.049 s |
| 动态 2 s | 400/400 步完成，全部 `converged=1`，`min J` 全程最低 0.765382，最大位移 0.382756 m |
| 资源 | 进程墙钟 5.49 s，峰值 RSS 408 764 kB（约 399 MiB） |
| 退出码 | 0 |
| 现有回归 | CTest 10/10 通过（物理批次之前运行） |

## 2. 预登记参数符合性

| 参数 | 预登记值 | 实际值 |
|---|---|---|
| 材料 | Neo-Hookean tet FEM | 同 |
| `E` | 1.0e6 Pa | 1.0e6 Pa |
| `nu` | 0.3 | 0.3 |
| `rho` | 1000 kg/m³ | 1000 kg/m³ |
| 重力 | (0, 0, -9.81) m/s² | 同 |
| 固定区域 | `y >= 0.45` 的节点 | 9 个节点，实际范围 `x∈[-0.154509, 0.154509]`、`y∈[0.457365, 0.5]`、`z∈[-0.058778, 0.058779]` |
| 静态 | `dt=1`、`staticSolve=true`、`tol=1e-8`、`maxIter=100`、DBC 初始 `k=1e8` | 同 |
| 动态 | `dt=0.005`、400 步（2 s）、`tol=1e-4`、`maxIter=50` | 同 |
| 接触 | 无 | 无 |

动态初态采用**参考形状 + 零速度 + 全重力**（与 S3 的动态约定一致），因此 2 s 内记录的是受突然加载后的瞬态摆动，而不是静态平衡附近的线性化响应；这一点在 `summary.txt` 中记为 `dynamic_initial_condition=reference_shape_zero_velocity`。

## 3. 资产审计（只读）

审计脚本：`Archive/fiction/output/validation/s4-nonrod-20261031-01/audit_torus.py`，结果 `asset_audit.json`。审计在求解前完成，通过后才进入构建与求解。

| 检查项 | 结果 |
|---|---|
| 节点 / tet | 200 / 491（与预登记一致） |
| 包围盒 | `[-0.5, -0.5, -0.095106]` – `[0.5, 0.5, 0.095106]`（与预登记一致） |
| 体积（有向求和 = 绝对值和） | 0.0726544475741816 m³ |
| 质量（`rho·V`） | 72.6544475741815 kg |
| 参考定向 | 491 个 tet 全部 `det > 0`；最小有向体积 `3.46891e-05`，无零/负体积 |
| 表面闭合 | 400 三角面、600 条无向边；每条无向边恰好出现 2 次，每条有向边恰好出现 1 次；Euler 示性数 `V-E+F = 0`（与 torus 拓扑一致） |
| 节点标签 / 坐标 | 唯一、全部有限 |
| 备选资产 | 未使用 `cow_head.msh`（torus 加载与求解均未失败） |

## 4. 静态满载结果

| 指标 | 数值 |
|---|---|
| `converged` / `maxIterHit` / `lineSearchFailed` | 1 / 0 / 0 |
| Newton 迭代 | 13 |
| `min J` / `max J` | 0.871477 / 1.116342 |
| `min` / `max` 主伸长 | 0.858509 / 1.141073 |
| 最大位移（相对参考形状） | 0.196671 m |
| DBC 实测最大误差 | 0（受约束自由度经罚项 + 消去后精确满足；`report.dbcMaxError` 同为 0） |
| 弹性能 / 重力能 / 总势能 | 28.831438 / -57.983894 / -29.152456 J |
| 墙钟 | 0.049 s |
| 门槛（收敛、有限、`min J>0.3`、DBC 误差 `<1e-6`） | 通过 |

静态解相对参考形状的最大位移约 0.197 m（环面直径 1 m），对应软材料（`E=1e6`）在自重下的明显下垂，`min J` 仍为 0.87，未接近反转。

## 5. 动态 2 s 结果

逐 step 数据：`dynamic_steps.csv`（400 行，列含 `converged / iterations / min_J / max_J / min_stretch / max_stretch / dbc_max_error / measured_dbc_error / max_disp / max_speed / elastic_energy / gravity_energy / total_potential / wall_seconds / solver / gate_pass`）。复核脚本 `check_dynamic.py` 输出 `dynamic_summary.json`。

| 指标 | 数值 |
|---|---|
| 完成步数 | 400 / 400（`failStep=-1`） |
| 每步 `converged` | 全部为 1 |
| `maxIterHit` / `lineSearchFailed` | 全程 0 |
| Newton 迭代 | 每步 1–2 次，均值 1.94；377 步用到 2 次 |
| `min J` 全程最低 | 0.765382（step 72，t=0.36 s） |
| `max J` 全程最高 | 1.213588 |
| `min` / `max` 主伸长 | 0.750236 / 1.260537 |
| 最大位移 | 0.382756 m（step 73，t=0.365 s） |
| 末步最大位移 | 0.217549 m |
| 最大节点速度 | 1.746426 m/s（全程）；末步 1.106036 m/s |
| DBC 实测最大误差 | 0（全程） |
| 弹性能范围 | 0.003112 – 103.652546 J |
| 总势能范围 | -28.903825 – -0.152243 J |
| 单步墙钟 | 中位约 0.00926 s，最大 0.00994 s |
| 动态总墙钟 | 5.130 s（含 81 帧软光栅化渲染） |
| 求解器路径 | CHOLMOD（全部步） |
| 门槛（每步收敛、有限、`min J>0.3`、DBC 误差 `<1e-6`） | 全部通过 |

末步最大节点速度 1.106 m/s，说明 2 s 时系统**仍在运动，未静止**。动态最低 `min J` 出现在 t≈0.36 s 的瞬态峰值附近，与最大位移同步，呈现突然满载后的瞬态摆动。

## 6. 视觉产物

固定视角（`yaw=45°`、`pitch=25°`、相机距离 `1.7×bbox 对角线`，不环绕），地面网格仅渲染、无接触：

| 文件 | 物理时刻 | 说明 |
|---|---|---|
| `keyframe_t0000.png` | t = 0 s | 参考形状（初始） |
| `keyframe_t1000.png` | t = 1.0 s | 中间态 |
| `keyframe_t2000.png` | t = 2.0 s | 末态 |

对应 OBJ 关键帧也在同一目录（`keyframe_t0000.obj` 等）。

短片：`video.mp4`，960×720，25 fps，81 帧，**容器时长 3.24 s**，对应**物理时间 2.0 s**（每 5 个时间步一帧：`(81-1)×5×0.005 = 2.0 s`）。短片以3.24 s容器时间展示2.0 s物理形变序列。

精选副本已放到：
`personal-homepage/outputs/physics-cpp-blog/validation/s4-nonrod-20261031-01/`（3 张关键帧 + `video.mp4`）。

## 7. 门槛判定

| 门槛 | 结果 |
|---|---|
| 现有 CTest 回归（先于物理批次） | 10/10 通过 |
| 静态收敛、有限 | 通过 |
| 静态 `min J > 0.3` | 通过（0.871477） |
| 静态 DBC 误差 `< 1e-6` | 通过（0） |
| 动态每步收敛 | 通过（400/400） |
| 动态每步有限 | 通过 |
| 动态每步 `min J > 0.3` | 通过（最低 0.765382） |
| 动态每步 DBC 误差 `< 1e-6` | 通过（0） |
| 失败即停 / 不调参重跑 | 未触发任何失败，未做二次调参运行 |

## 8. 实验配置

- 本阶段检查求解状态、安全门槛与形变序列。
- 网格很粗（200 节点 / 491 tet），环面截面仅由少量四面体构成，局部应力/弯曲刚度受离散化影响。
- 固定区域只按 `y >= 0.45` 选取了 9 个节点（顶部一小段弧），因此整体存在明显的悬挂/摆动态。
- 动态只跑了 2 s，末态未静止（最大速度 1.1 m/s）。
- 无接触；地面网格用于渲染，求解中没有地面接触能量。
- 未验证第三方库（CHOLMOD/BLAS 等）是否完全遵守线程数环境变量；线程数已按要求设为 1。

## 9. 证据清单（run-id 目录内）

| 文件 | 内容 |
|---|---|
| `summary.txt` | 配置、资产规模、静态/动态关键数字、视频与关键帧路径 |
| `static_steps.csv` | 静态单行指标（`static_solve=1`） |
| `dynamic_steps.csv` | 400 行动态逐步指标 |
| `static_state.csv` / `final_state.csv` | 静态/末态逐节点参考坐标、当前位置、位移、DBC 标记 |
| `asset_audit.json` + `audit_torus.py` | 资产只读审计结果与脚本 |
| `dynamic_summary.json` + `check_dynamic.py` | 逐行门槛复核与极值汇总 |
| `logs/command.txt` | 实际运行命令（含环境变量） |
| `logs/run.log` | 程序标准输出 |
| `logs/time.log` | `/usr/bin/time -v`（墙钟、峰值 RSS） |
| `logs/exit_code.txt` | 退出码 0 |
| `logs/ctest.log` | 物理批次前的 CTest 结果（10/10） |
| `logs/build.log` / `logs/cmake.log` | Release 构建日志 |
| `logs/sha256.txt` | 源码、CMakeLists、二进制哈希 |
| `source_snapshot/` | `s4_nonrod_check.cpp` 与注册后的 `CMakeLists.txt` 快照 |
| `frames/` | 81 帧软光栅化 PNG |
| `keyframe_t0000/1000/2000.{png,obj}` | 3 张固定视角关键帧及 OBJ |
| `video.mp4` | 3.24 s 容器 / 2.0 s 物理 短片 |

关键哈希：

```text
s4_nonrod_check.cpp  ae46b3fb482c458cb46ebd59b942418b0d0729e6ce012dd9b1ce1e8ad75110ab
CMakeLists.txt       abba584945630d311681c703d5bb43d4eb8b9698f0dd6e34bfb15167dfc6f5ad
s4_nonrod_check      b43a2196fb52ff28521fa1dca8a3ce6bb9301173947a194da360a4022d3d1f52
```

## 10. 复现

```bash
cd /home/steven/mjx/fiction/experiments/workspaces/fiction-s4-20261031-01
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S cpp -B /tmp/fiction-s4-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/fiction-s4-build -j16
ctest --test-dir /tmp/fiction-s4-build --output-on-failure
/tmp/fiction-s4-build/s4_nonrod_check \
  --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/torus.msh \
  --out /home/steven/mjx/fiction/simulator/output/validation/s4-nonrod-20261031-01 \
  --steps 400 --export-every 5 --fps 25 --video-size 960 720
```

静态与动态在同一次进程内顺序执行：静态未过门槛时程序写 `FAILURE.txt` 并以非零码退出，不进入动态；动态某步失败时保留已完成 CSV、末态与现场并立即停止。
