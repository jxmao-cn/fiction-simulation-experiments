# 停止判据诊断结果（2026-10-31，run-id `stopcriterion-20261031-01`）

日期：2026-10-31。性质：**隔离副本内的诊断**。生产 `Archive/fiction/cpp`、协议 v0.2 与既有 run-id 均未改动（见 §8 源码与运行记录）。预登记见 [stopcriterion-diagnostic-spec-20261031.md](stopcriterion-diagnostic-spec-20261031.md)。

## 1. 一句话结论

- **actual 与 combined 在本诊断的两个病例上都通过三项门槛**（物理解不变、假阴性消除、不提前误停）。
- **grad 只在 S2 M3 上通过，在 S6a 1 体上不通过**：它把末态改变了 $3.2\times10^{-6}$ m（门槛 $10^{-8}$ m）。原因是"相对梯度容差"未约束 $H^{-1}$ 的条件数：软体上同一个相对梯度对应的位移残差比刚性梁大 3–4 个量级，动态步还要再被 $1/h$ 放大。
- 两种病例中，原 `maxIter=100` 的停滞步在 actual/combined 下全部变成 `converged=1` 且 `stopReason != maxIter`，**假阴性被消除**；没有出现"梯度仍然很大就停"的预登记 false positive。
- 推荐候选：**combined**（与 actual 同等准确，成本更低：S6a 36 次迭代 / 159 次线搜索，对比 actual 的 55 / 150）。combined 通过了本诊断的三项门槛。

## 2. 方法

### 2.1 隔离副本与构建

- 副本：`/home/steven/mjx/fiction/experiments/workspaces/fiction-stopcriterion-20261031-01/`，以 `fiction-s6a-20261031-01` 为基（已含 `scene/multibody_scene.h` 与 `s6a_multibody_check.cpp`），并把 `Archive/fiction/cpp/examples/cantilever_check.cpp` 原样复制进来（`cpp/CMakeLists.txt` 中该目标已注册，无需新增）。
- 构建：`cmake -S cpp -B /tmp/fiction-stopcriterion-build -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_CHOLMOD=ON`，CHOLMOD 生效（`/usr` SuiteSparse）。`fiction_core` 无新增编译警告（`-Wall -Wextra -Wpedantic`）。
- CTest：**10/10 通过**（`ctest.log`），回归未失败。

### 2.2 判据实现（只改副本的 `TimeIntegrator`）

`Config` 新增：`stopMode`（`step|grad|actual|combined`，默认 `step`）、`gradTolRel`（默认 `1e-7`）、`gradScale`（默认 `1`，调用方传总重）、`actualStepTol`（默认 `1e-12` m）、`actualStepCount`（默认 `2`）。`StepReport` 新增 `stopReason`（`step|grad|actual|combined|maxIter|lineSearch|linearSolve`）。

停止判定放在 Newton 循环顶部（与原判据同一位置，保证"判据成立时不消耗 `maxIter`"的语义不变）：

$$
\text{grad:}\quad \frac{\lVert g_{\mathrm{free}}\rVert_\infty}{g_{\mathrm{scale}}} \le 10^{-7},\qquad
\text{actual:}\quad \text{连续 } 2 \text{ 次 } \lVert \Delta x\rVert_\infty \le 10^{-12}\,\mathrm{m},
$$

$$
\text{combined:}\quad \frac{\lVert p\rVert_\infty}{h}\le \mathrm{tol}
\ \ \text{或}\ \
\left(\text{actual 条件}\ \wedge\ \frac{\lVert g_{\mathrm{free}}\rVert_\infty}{g_{\mathrm{scale}}} \le 10^{-6}\right),
$$

三者都额外要求 DBC 满足。自由梯度 $\lVert g_{\mathrm{free}}\rVert_\infty$ 抽成私有 helper `computeFreeGradientInf`（总梯度减 DBC 罚项、只统计非 DBC 节点），与末态诊断**完全同一算法**，末态报告也改走该 helper，停止判据与报告共用这一实现。

`converged` 现在要求：生效判据在停止状态成立 **且** DBC 满足 **且** `maxIterHit/lineSearchFailed/linearSolve` 均未触发、全部量有限。`maxIterHit`、`lineSearchFailed` 与 DBC 罚刚度加倍逻辑保持不变。

### 2.3 CLI 与力尺度

`cantilever_check` 与 `s6a_multibody_check` 新增 `--stop-mode`、`--grad-tol`、`--grad-scale`、`--actual-step-tol`、`--actual-step-count`；缺省 `--stop-mode step`，且 `step` 模式下 `gradScale` 不参与判定，因此默认行为与改动前一致（§3 用逐位对比证明）。

- S2 M3：`gradScale = 24.525` N（总重 $2.5\,\mathrm{kg}\times 9.81$）。
- S6a 1 体：`gradScale = 518.4769746175648 \times 9.81 = 5086.259120998311` N。

### 2.4 运行

8 个 run，逐 case 独立目录，**顺序执行**（同一时刻只有一个进程），线程变量 `OMP/OPENBLAS/MKL/BLIS/VECLIB=1`。驱动脚本：副本内 `run_stopcriterion.sh`；汇总脚本：`tools/summarize_stopcriterion.py`（生成 `config.json`、逐 case `summary.json`、`cases.csv`、`judgment.json`）。

## 3. 基线复现（先证明改动没动数值）

| 病例 | 迭代 | 线搜索 | stopReason | 与既有 run 末态逐节点最大差 |
|---|---:|---:|---|---:|
| S2 M3（step, tol=1e-10） | 100 | 3470 | `maxIter` | **0.0 m**（`M3-state.csv` 逐位相同） |
| S6a 1 体（step, tol=1e-8） | 317 | 7445 | 6/8/9 步 `maxIter` | **0.0 m**（`body1_state.csv` 逐位相同） |

此外 S2 基线复现出参考行的全部诊断量（方向率 `1.0545646946513488e-09`、自由梯度 `1.2264105903048161e-07` N、尖端 `7.407508555371e-3 m`），S6a 逐步迭代数 `[2,2,2,2,3,100,3,100,100,3]` 与既有 run 完全一致。**stop 模式 = 原算法，逐位复现。**

## 4. 八个 case 的结果

### 4.1 S2 M3 静态（`cantilever_check`，tol=1e-10，maxIter=100）

| 判据 | stopReason | 迭代 | 线搜索 | converged | 停止时 $g_{\mathrm{free}}/g_{\mathrm{scale}}$ | 尖端位移 (mm) | 尖端差 (m) | 节点最大差 (m) | 墙钟 (s) | RSS (MiB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| step（基线） | `maxIter` | 100 | 3470 | 0 | 5.001e-09 | 7.407508555 | — | 0.0 | 15.28 | 97.1 |
| grad | `grad` | 3 | 0 | 1 | 5.711e-09 | 7.407508660 | 1.04e-10 | 1.51e-10 | 0.63 | 96.1 |
| actual | `actual` | 6 | 27 | 1 | 5.003e-09 | 7.407508556 | 7.12e-13 | 1.03e-12 | 1.08 | 96.9 |
| combined | `combined` | 6 | 27 | 1 | 5.003e-09 | 7.407508556 | 7.12e-13 | 1.03e-12 | 1.08 | 96.6 |

三者 DBC 均 81/81 满足，`minJ = 0.999511375` 与基线一致。变体下 `numerical_pass` 标志变为 1，覆盖收敛/有限/minJ/DBC；S2 的梁理论相对误差为 25.94%，与基线同值。

### 4.2 S6a 1 体（`s6a_multibody_check --only-body1`，10 步，tol=1e-8，maxIter=100）

| 判据 | 逐步 stopReason | 总迭代 | 总线搜索 | 全部步 converged | 停滞步(6/8/9) 迭代 | 停止时最大 $g_{\mathrm{free}}/g_{\mathrm{scale}}$ | 末态最大位移 (m) | 末态节点最大差 (m) | 墙钟 (s) | RSS (MiB) |
|---|---|---:|---:|---|---|---:|---:|---:|---:|---:|
| step（基线） | 6/8/9 = `maxIter` | 317 | 7445 | 否 | 100/100/100 | 1.05e-14 | 1.3471948779e-2 | 0.0 | 18.98 | 48.4 |
| grad | 全 `grad` | 10 | 0 | 是 | 1/1/1 | 4.07e-10 | 1.3471948151e-2 | **3.21e-06**（未过） | 1.22 | 47.9 |
| actual | 全 `actual` | 55 | 150 | 是 | 5/6/5 | 4.41e-15 | 1.3471948779e-2 | 4.95e-11 | 3.65 | 48.2 |
| combined | 全 `combined` | 36 | 159 | 是 | 5/5/6 | 1.03e-14 | 1.3471948780e-2 | 1.61e-10 | 2.59 | 47.5 |

逐步迭代数：

- baseline：`[2,2,2,2,3,100,3,100,100,3]`，线搜索 `[0,0,0,0,0,2441,0,2451,2553,0]`
- grad：`[1,1,1,1,1,1,1,1,1,1]`，线搜索全 0
- actual：`[4,5,6,5,6,5,6,6,5,7]`，线搜索 `[4,10,22,3,7,10,19,36,0,39]`
- combined：`[2,2,2,2,3,5,3,5,6,6]`，线搜索 `[0,0,0,0,0,47,0,41,44,27]`

## 5. 物理解差（门槛 $<10^{-8}$ m）

| 病例 | grad | actual | combined |
|---|---:|---:|---:|
| S2 节点最大差 / 尖端差 | 1.51e-10 / 1.04e-10 ✅ | 1.03e-12 / 7.12e-13 ✅ | 1.03e-12 / 7.12e-13 ✅ |
| S6a 末态节点最大差 | **3.21e-06 ❌** | 4.95e-11 ✅ | 1.61e-10 ✅ |

**grad 在 S6a 上的失败机制**：它每步只做 1 次 Newton 迭代就停，停止时残余 Newton 方向 $\lVert p\rVert_\infty = 1.66\times10^{-7}$ 增长到 $1.78\times10^{-6}$ m（方向率 $3.3\times10^{-5}$–$3.6\times10^{-4}$，乘 $h=0.005$）。每步遗留约 $10^{-7}$ m 的位移残差，又通过 $v=(x-x_n)/h$ 被放大 $1/h=200$ 倍带进下一步，10 步累积到 $3.2\times10^{-6}$ m。根本原因是**梯度判据不控制条件数**：位移误差是 $H^{-1}g$，而 S6a 的软体（$E=10^5$ Pa）在近奇异弯曲模态上的 $H$ 远小于 S2 的梁（$E=5.9\times10^8$ Pa），所以同一个相对梯度容差在 S2 对应 $1.2\times10^{-9}$ m 的残余方向（静态，$h_\mathrm{eff}=1$，没有速度放大），在 S6a 对应 $10^{-7}$–$10^{-6}$ m 的残余方向，并叠加 $1/h=200$ 的速度放大。

## 6. 假阴性消除与 false positive 检查

### 6.1 假阴性（原 `maxIter=100` 的步）

| 病例 | 基线停滞步 | grad | actual | combined |
|---|---|---|---|---|
| S2 M3（整步） | `maxIter`，100 迭代，`converged=0` | `grad`，3 迭代，`converged=1` ✅ | `actual`，6 迭代，`converged=1` ✅ | `combined`，6 迭代，`converged=1` ✅ |
| S6a 步 6 / 8 / 9 | `maxIter`，100/100/100，`converged=0` | 1/1/1，`converged=1` ✅ | 5/6/5，`converged=1` ✅ | 5/5/6，`converged=1` ✅ |

**三种判定都在两种病例上消除了假阴性**（`stopReason` 均不是 `maxIter`）。

停滞的机制在 S6a 基线里看得很清楚：第 6/8/9 步第 4 次迭代后 $x_\mathrm{trial}=x+\alpha p$ 在浮点上已等于 $x$，能量比较退化成"相等"，线搜索要连砍约 24 次才接受；100 次迭代、约 2400 次能量评估，位置一个 ulp 都没动。

### 6.2 提前误停（预登记定义：停止时 $g_{\mathrm{free}}/g_{\mathrm{scale}}\ge 10^{-6}$）

- 全部 8 个 run 的**每一个停止点**都满足 $g_{\mathrm{free}}/g_{\mathrm{scale}} < 10^{-6}$：S2 最大 5.71e-09，S6a 最大 4.07e-10。
- 因此**按预登记定义没有 false positive**。
- grad 在 S6a 上改变了物理解（§5），失败项为"物理解不变"门槛，机制是**判据对位移过松**导致提前停止。actual 与 combined 在两种病例上同时满足位移比较和相对自由梯度门槛。

## 7. 单次运行成本

- S2：基线 100 迭代 + 3470 线搜索 ≈ 15.3 s/次求解；actual、combined 6 + 27 ≈ 1.08 s；grad 3 + 0 ≈ 0.63 s。
- S6a 10 步：基线 317 + 7445 ≈ 19.0 s；combined 36 + 159 ≈ 2.59 s；actual 55 + 150 ≈ 3.65 s；grad 10 + 0 ≈ 1.22 s。
- RSS 48–99 MiB，四种判据无实质差别（矩阵与网格相同）。
- **8 个 run 顺序执行**，各记录一次墙钟；CPU 绑定、频率和温度采用系统当时设置。

## 8. 源码与运行记录

- `Archive/fiction/cpp/src/time_integrator.cpp`、`include/time_integrator.h`、`examples/cantilever_check.cpp` 的 mtime 为 14:57–15:07（本次会话在 17:30 之后），且与 `fiction-s6a-20261031-01` 快照逐文件 `diff` 相同；`fiction-s6a-20261031-01` 也未被写入（新建副本为独立目录）。
- 未 `git commit` / `git push`；未修改协议文档、既有 run-id 与生产源码。
- 未出现 NaN/Inf（对所有 CSV 的数值单元做了 `math.isfinite` 扫描，0 个非有限值），无 `lineSearch`/`linearSolve` 停止，无阈值放宽、无调参重跑。
- 病例为两个单一网格问题，S6a 使用 1 体 10 步。S6a 的"参考末态"本身含 3 个 maxIter 步（pilot 已证 maxIter 3 vs 100 末态差 $6.1\times10^{-16}$ m，可作参考）。`actualStepTol=1e-12` m 与坐标 ulp（约 $2.2\times10^{-16}$ m）差 4 个量级；若换到 $10^3$ m 量级的场景需重新核对该阈值。
- 变体改动集中在 `--stop-mode`，能量、装配、线搜索、DBC 沿用基线。

## 9. 推荐结论

1. **本诊断通过**：`actual`、`combined` —— 两种病例都满足"物理解不变 <1e-8 m""假阴性消除""停止时相对自由梯度 <1e-6"。
2. **本诊断不通过**：`grad` —— S2 通过、S6a 不通过（末态差 3.2e-6 m）。相对梯度容差不能单独作为停止判据，除非按问题刚度重新标定，而那属于另一次预登记。
3. **若只能选一个**：选 `combined`。它保留原 `step` 判据作为良态问题的快路径（S2 与 actual 同为 6 迭代；S6a 前 5 步与第 7 步仍是 2–3 迭代），只在停滞步启用"位移冻结 + 相对梯度已小"的出口，成本低于纯 actual（36 vs 55 迭代），准确度与 actual 同级（1.6e-10 vs 5.0e-11 m）。
4. **后续**：若要写进协议，需另立 v0.3 停止判据小节、在新 run-id 下覆盖 S2 全网格、S6a 4 体与后续批次重新验收；本文件记录候选判据在 S2 M3 + S6a 1 体上的诊断结果。

## 10. run-id 证据清单

根目录：`Archive/fiction/output/validation/stopcriterion-20261031-01/`

| 路径 | 内容 |
|---|---|
| `env.txt` | 主机、构建、网格、力尺度、线程变量、执行方式 |
| `ctest.log` | 10/10 通过 |
| `asset_sha256.txt` / `binary_sha256.txt` | 输入网格与两个可执行文件哈希 |
| `source_files_sha256.txt` / `source.tar.gz` / `source_hash.txt` | 副本源码快照与哈希 |
| `source-change.diff` | 副本相对 `fiction-s6a-20261031-01` 的全部改动（412 行） |
| `cases.csv` | 8 个 case 的一行汇总（含门槛判定列） |
| `judgment.json` | 逐 case/逐判据的门槛判定与对照基线差值 |
| `s2/M3/{step,grad,actual,combined}/` | 各自 `config.json`、`command.txt`、`M3.csv`、`M3-state.csv`、`summary.json`、`run.log`、`time_rss.txt`、`exit_code.txt` |
| `s6a/body1/{step,grad,actual,combined}/` | 各自 `config.json`、`command.txt`、`asset_audit.json`、`body1_steps.csv`、`body1_state.csv`、`summary.json`、`run.log`、`time_rss.txt`、`exit_code.txt` |

退出码说明：`s2/M3/step` 为 1、`s6a/body1/step` 为 2，与既有 run 一致（原判据下 `converged=0` 导致 `NUMERICAL_FAIL` / `S6A_BODY1_FAIL`）；三个变体在两种病例上均为 0。这些退出码对应本次比较的原停止状态与变体状态。

## 11. 与预登记规格的偏差

无。未改门槛、未调参数、未重跑；未出现需要保存现场的失败（无回归失败、无 NaN、无提前误停）。

§6.2 按预登记定义分别记录两个判断：grad 在 S6a 通过相对自由梯度门槛，未通过"物理解不变"门槛。
