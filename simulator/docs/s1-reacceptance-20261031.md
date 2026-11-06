# S1 主工程合并与重新验收

日期：2026-10-31。状态：**已合并到工作区，并按明确记录的罚参数日程完成 S1，三档网格共 12/12 载荷通过。**

## 1. 合并内容与范围

用户明确授权合并并重新验收。以合并前源码快照为基线，只修改四个既有文件，并新增一个回归程序：

| 文件 | 内容 |
|---|---|
| `cpp/src/time_integrator.cpp` | 增大罚刚度后立即刷新方向、满足节点集合和报告；合入默认关闭的逐轮/逐节点诊断，节点 CSV 17 位精度 |
| `cpp/include/time_integrator.h` | 可选 trace 路径与标签 |
| `cpp/examples/tension_check.cpp` | 最后载荷的刚度/迭代 CLI 参数，初始和末态记录；帮助信息说明验收命令 |
| `cpp/CMakeLists.txt` | 注册新回归测试 |
| `cpp/examples/dbc_penalty_refresh_check.cpp` | 一维最小反例：增刚度后应改变旧方向，并降低**同一新目标**的能量 |

`scope-check.json` 核实合并前快照中的其他源码未改变。旧工作区未提交修改均保留。未改变材料、通用 DBC 默认参数、增刚度门槛、线搜索、原验收容差和最大迭代数；未改硬 DBC 或罚参数跨载荷继承策略。

## 2. 验收采用的实际配置

运行前配置见 [S1 重新验收规格](s1-reacceptance-spec-20261031.md)。精确沿用已验证的两组实验高初始 k 配置：前三个载荷初始 k=1e8，最后 5% 载荷初始 k=1e10；每载荷仍自适应增大 k，maxIter=100、tol=1e-10。

runner 默认值为兼容历史仍保留 1e8；验收命令显式传入`--final-dbc-stiffness 1e10`，对应本次S1罚参数日程。

网格重新生成并自检，与原几何一致：长1 m、截面0.05×0.05 m，y/z 居中，体积0.0025 m³、质量2.5 kg；E=1e6 Pa、nu=0.3，静态，两端给定仿射位置，中部从上一载荷求解结果继续。

## 3. S1 结果

| 网格 | 节点 / tet | 载荷 1 / 1.01 / 1.02 / 1.05 的 Newton 次数 | 结果 |
|---|---:|---|---|
| M1 | 99 / 240 | 0 / 18 / 20 / 14 | 4/4 通过 |
| M2 | 525 / 1920 | 0 / 21 / 24 / 15 | 4/4 通过 |
| M3 | 3321 / 15360 | 0 / 23 / 24 / 16 | 4/4 通过 |

各档均满足下列逐载荷判据：

| 非零载荷中的最大误差 | M1 | M2 | M3 | 上限 |
|---|---:|---:|---:|---:|
| 端面位置（m，无穷范数） | 9.41e-11 | 9.41e-11 | 9.41e-11 | 1e-6 |
| 中段仿射（m） | 5.31e-11 | 5.22e-11 | 5.30e-11 | 1e-7 |
| 全部自由节点仿射（m） | 7.16e-11 | 7.93e-11 | 8.61e-11 | 1e-5 |
| 能量相对误差 | 2.79e-8 | 2.81e-8 | 2.67e-8 | 1e-5 |
| 反力相对误差 | 1.44e-8 | 1.45e-8 | 1.38e-8 | 1e-5 |
| 所有 tet 的 J 最大参考偏差 | 4.86e-10 | 1.33e-9 | 3.08e-9 | 1e-8 |

零载荷另通过原绝对能量、反力和位移条件；所有单元正体积，所有载荷无迭代超限或线搜索失败。S1检查指定均匀形变下的材料、装配与静态边界，三档网格分别完成相同载荷路径。

## 4. 回归与负对照

- Release、CHOLMOD：10/10 CTest 通过，含 2 个 GUI 占位项。
- LDLT-only 独立构建：S0/T0、报告检查、新 DBC 刷新检查、simulator 和 mesh_check 共 6/6 通过。没有用 LDLT 重跑三档 S1。
- 新 DBC 刷新反例在合并前核心库上返回失败（位置不前进），在合并后通过。负对照使用合并前头文件与核心库，避免混用新旧 Config 的 ABI；它验证回归用例能够捕获旧问题。

本轮完整 S1 三档（含输出、5% trace）的单次墙钟为 0.14/1.32/16.72 s，峰值 RSS 为 8944/19100/100260 KiB。这些为包含输出与trace的单次运行记录。

## 5. 证据、复现与回滚

run-id：[s1-merged-20261031-01](../output/validation/s1-merged-20261031-01)。保存 M1/M2/M3 CSV、日志、退出码、5%逐轮/逐节点 trace、自检网格、summary.json、两后端构建/CTest、新回归负对照、合并前后源码与哈希。

```bash
cd /home/steven/mjx/fiction/simulator
cmake -S cpp -B /tmp/fiction-s1-merged-20261031-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/fiction-s1-merged-20261031-build -j 3
ctest --test-dir /tmp/fiction-s1-merged-20261031-build --output-on-failure
# 示例复现：输出新文件，避免覆盖冻结验收结果
/tmp/fiction-s1-merged-20261031-build/tension_check \
  --mesh output/validation/s1-merged-20261031-01/meshes/M3.msh \
  --final-dbc-stiffness 1e10 --final-max-iter 100 --csv /tmp/s1-M3-reproduce.csv
```

仅本轮代码的逆补丁为 `merge-only.patch`，不是包含用户旧修改的整个 git diff。已在独立 `/tmp/fiction-s1-merged-rollback-check` 实际逆向应用，四个既有文件与合并前逐字节一致，新回归文件按预期删除；主工程保留合并状态。需要撤销时先检查当前没有后续改动，再执行：

```bash
cd /home/steven/mjx/fiction/simulator
patch -R --dry-run -p1 < output/validation/s1-merged-20261031-01/merge-only.patch
patch -R -p1 < output/validation/s1-merged-20261031-01/merge-only.patch
```

逆补丁撤销本轮代码与 CMake，保留实验数据和说明文档；回滚后重新构建，不使用旧二进制代表回滚源码。不要 git reset 或将全目录快照覆盖到已有新修改的工作区。

S1 到此完成并停止，后续 S2/动态展示/规模与多体仍需按原阶段继续。
