# combined 停止判据候选确认结果（2026-10-31，run-id `s2s6a-combined-20261031-01`）

日期：2026-10-31。性质：**候选判据确认（重验收）**，不是协议通过。协议 v0.2、既有 run-id 与生产 `Archive/fiction/cpp` 均未改动（生产源码 mtime 仍为 14:57–15:07，本次运行在 17:38 之后，见 §9）。预登记门槛与 [protocol-v0.3 草案](protocol-v0.3-stopcriterion-draft-20261031.md) §4 一致，未放宽、未调参重跑。

## 1. 一句话结论

- `combined`（`--grad-tol 1e-6 --actual-step-tol 1e-12 --actual-step-count 2`）在 **S2 M3 / S2 M6 / S6a 1 体 / S6a 4 体** 四个 case 上**全部通过**预登记门槛：物理解不变、每步/整步 converged、块对角 0 跨体、4 体 vs 独立解一致、无提前误停。
- S2 的精度门槛没有变化：M3 相对梁理论误差仍 **25.94%**，M6 仍 **4.26%**；判据只消除代数停滞，不改善空间离散精度。
- 最大末态差出现在 M6：$3.21\times10^{-9}$ m（门槛 $10^{-8}$ m，余量约 3.1 倍）。M3 为 $1.27\times10^{-12}$ m，S6a 1 体 $1.91\times10^{-10}$ m，4 体 $1.82\times10^{-10}$ m。
- **本结果只能写成“候选判据在 S2/S6a 上重验收通过”，不能写成协议 v0.3 已生效。**

## 2. 方法

### 2.1 隔离副本与构建

- 源：`/home/steven/mjx/fiction/experiments/workspaces/fiction-stopcriterion-20261031-01/`（已含 `combined` 实现与 `s6a_multibody_check`）整目录复制为 `/home/steven/mjx/fiction/experiments/workspaces/fiction-s2s6a-combined-20261031-01/`。生产 `Archive/fiction/cpp`、`stopcriterion` 原副本、既有 run-id 均未写入。
- 构建：`cmake -S cpp -B /tmp/fiction-s2s6a-combined-build -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_CHOLMOD=ON`，CHOLMOD 生效（`/usr` SuiteSparse）。`fiction_core` 与两个 example 均无新增编译警告（`-Wall -Wextra -Wpedantic`，build.log 中 warning/error 计数 0）。
- CTest：**10/10 通过**（`ctest.log`），回归未失败。

### 2.2 相对 `stopcriterion` 副本的实现改动（`source-change.diff`，111 行）

只动两个 example，不动 `time_integrator` 的判据实现：

1. `s6a_multibody_check.cpp`：4 次独立单体求解原本**不传**停止判据（恒用默认 `step`），与全局 `combined` 求解不可比。现改为传播 `stopMode/gradTolRel/gradScale/actualStepTol/actualStepCount`，独立单体用**单体重** $5086.259$ N。
2. `s6a_multibody_check.cpp`：4 体 job 内部的 1 体诊断场景原本沿用 4 体总重作为 `gradScale`（相对梯度被系统性缩小 4 倍），现按**本场景总重**缩放（1 体 $5086.259$ N，4 体 $20345.036$ N）；`summary.json` 增加 `stop_criterion` 块显式记录三个尺度。
3. `s6a_multibody_check.cpp`：新增 `--run-id`，只写入 `summary.json` 溯源字段（原程序把 `s6a-multibody-20261031-01` 写死，会污染本 run-id 的证据）。
4. `cantilever_check.cpp`：资产自检的体积改用 **long double** 独立求和，门槛仍为 $10^{-12}$ 相对（详见 §9 偏差说明）。M6（414720 tet）的 double 顺序累加体积为 $0.0025000000000128$，相对偏差 $5.1\times10^{-12}$，会被原自检误判为几何错误；这是**输入资产自检**，不是验收门槛。求和方式与生成 M6 参考的 `s2-extension-20261031-01` 一致。

第 1、2 项是为了让“4 体 vs 独立单体”和“每场景相对梯度”这两个预登记检查真正成立；第 3 项是证据溯源；第 4 项是数值资产自检的舍入问题。**没有改动能量、装配、线搜索、DBC、时间积分与停止判据本身。**

### 2.3 运行配置

| case | 程序 | 网格 | tol / maxIter | gradScale [N] |
|---|---|---|---:|---:|
| S2 M3 | `cantilever_check`（硬编码 tol=1e-10、maxIter=100、静态） | `s2-cantilever-20261031-01/meshes/M3.msh`（3321 节点 / 15360 tet） | 1e-10 / 100 | 24.525（总重 $2.5\,\mathrm{kg}\times9.81$） |
| S2 M6 | 同上 | `s2-extension-20261031-01/meshes/M6.msh`（75625 节点 / 414720 tet） | 1e-10 / 100 | 24.525 |
| S6a 1 体 | `s6a_multibody_check --only-body1` | `IPC/input/tetMeshes/sphere1K.msh`（1760 节点 / 6851 tet） | 1e-8 / 100，10 步 | 5086.259120998311（$518.4769746175648\times9.81$） |
| S6a 4 体 | `s6a_multibody_check`（默认 4 体） | 同一资产纯平移 4 份（7040 节点 / 27404 tet） | 1e-8 / 100，10 步 | 20345.036483993244（4× 单体重） |

四个 case 统一 `--stop-mode combined --grad-tol 1e-6 --actual-step-tol 1e-12 --actual-step-count 2`。顺序执行（同一时刻一个进程），`OMP/OPENBLAS/MKL/BLIS/VECLIB=1`；M6 单 case 186.5 s 完成，仍设 `timeout 600` 上限。执行顺序 M3 → S6a 1 体 → S6a 4 体 → M6（把最慢、风险最高的放最后，避免超时损失其余证据）。

### 2.4 预登记判定

- 物理解不变：末态最大节点差 $<10^{-8}$ m（S6a 4 体另加 4 体 vs 4 次独立单体 $<10^{-8}$ m）。
- 每步 converged：S6a 1/4 体 10/10；S2 静态 `converged=1` 且 `stopReason=combined`。
- 块对角：4 体跨体 triplet / 组装后 entry / 跨体 tet 均为 0。
- 不提前误停：每个停止点 $\lVert g_{\mathrm{free}}\rVert_\infty/g_{\mathrm{scale}}<10^{-6}$。
- S2 精度：`|relative_error - 参考值| <= 1e-6`（即 25.94% / 4.26% 不变），且 S2 精度门槛本身不变。
- 全部量有限、退出码 0。

## 3. 四个 case 总表

| case | stopReason | 迭代 | 线搜索 | converged | 停止时 $g_{\mathrm{free}}/g_{\mathrm{scale}}$ | 末态最大节点差 vs v0.2 | 墙钟 (s) | 峰值 RSS (MiB) | exit |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| S2 M3 | `combined` | 6 | 27 | 1 | $5.00\times10^{-9}$ | $1.27\times10^{-12}$ m | 1.09 | 97.7 | 0 |
| S2 M6 | `combined` | 7 | 50 | 1 | $6.43\times10^{-10}$ | $3.21\times10^{-9}$ m | 186.50 | 2519.3 | 0 |
| S6a 1 体 | `combined`×10 | 36（总） | 159（总） | 10/10 | $1.02\times10^{-14}$（最大） | $1.91\times10^{-10}$ m | 2.68 | 48.0 | 0 |
| S6a 4 体 | `combined`×10 | 31（总） | 70（总） | 10/10 | $1.97\times10^{-15}$（最大） | $1.82\times10^{-10}$ m | 24.32 | 157.2 | 0 |

对照基准：S2 M3 → `s2-cantilever-20261031-01/M3-state.csv`（v0.2 strict，step 判据 `maxIter` 末态）；S2 M6 → `s2-extension-20261031-01/M6-tol1e8/state.csv`（v0.2 tol=1e-8）；S6a → `s6a-multibody-20261031-01/body{1,4}_state.csv`（v0.2 step 判据，含 `maxIter` 步）。

## 4. S2 静态细节

### 4.1 S2 M3

| 量 | 值 |
|---|---:|
| 迭代 / 线搜索 | 6 / 27 |
| `stopReason` | `combined` |
| 方向率 $\lVert p\rVert_\infty/h$ | $1.0555\times10^{-9}$（tol 的 10.6 倍，方向率已冻结） |
| 自由梯度 / 相对值 | $1.2269\times10^{-7}$ N / $5.00\times10^{-9}$ |
| 尖端位移 | $7.407508556\times10^{-3}$ m（参考 $7.407508555\times10^{-3}$ m） |
| 相对梁理论误差 | 25.9413%（参考 25.9413%，绝对差 $7.1\times10^{-11}$） |
| 末态最大节点差 vs v0.2 M3 | $1.27\times10^{-12}$ m |
| 与 `stopcriterion` combined 末态差 | **0.0 m（逐位相同）** |
| `min J` / DBC | 0.9995114 / 81/81 |

迭代轨迹（`s2/M3/trace/iterations.csv`）：方向率 $7.41\times10^{-3}\to3.13\times10^{-5}\to2.55\times10^{-8}\to1.20\times10^{-9}\to1.056\times10^{-9}\to1.056\times10^{-9}$；第 4 次起方向率冻结在 tol 之上，实际步长 $\lVert\Delta x\rVert_\infty$ 为 $1.5\times10^{-10}\to2.6\times10^{-13}\to2.6\times10^{-13}$，第 5、6 次连续满足 $10^{-12}$ 且相对梯度 $5\times10^{-9}<10^{-6}$，判据成立即停。这正是诊断里记录的浮点不动点，`combined` 在 6 次迭代解决，而不是耗满 100 次。

### 4.2 S2 M6（本 run 最严的 case）

| 量 | 值 |
|---|---:|
| 迭代 / 线搜索 | 7 / 50 |
| `stopReason` | `combined` |
| 方向率 | $1.6595\times10^{-9}$（tol 的 16.6 倍，已冻结） |
| 自由梯度 / 相对值 | $1.5780\times10^{-8}$ N / $6.43\times10^{-10}$ |
| 尖端位移 | $0.009576142682$ m（参考 $0.009576145248$ m，差 $2.6\times10^{-9}$ m） |
| 相对梁理论误差 | 4.25973%（参考 4.25970%，绝对差 $2.6\times10^{-7}$） |
| 末态最大节点差 vs v0.2 M6 (tol=1e-8) | $3.21\times10^{-9}$ m |
| `min J` / DBC | 0.9991480 / 625/625 |

迭代轨迹：第 3 次方向率 $4.06\times10^{-9}$（已低于 v0.2 的 tol=1e-8），第 4 次 $3.11\times10^{-9}$，第 5–7 次冻结在 $1.659\times10^{-9}$；实际步长 $1.0\times10^{-9}\to1.6\times10^{-9}\to6.3\times10^{-15}\to3.5\times10^{-18}$（后两次线搜索分别二分 18 与 29 次才接受），第 6、7 次连续满足位移冻结且相对梯度 $6.4\times10^{-10}<10^{-6}$，第 7 次停止。**M6 与 M3 是同一机制：方向率冻结在 tol 之上，位移落到坐标 ulp 附近。**

M6 的末态差 $3.21\times10^{-9}$ m 同时包含两件事：(a) 参考本身是 tol=1e-8 的平衡点，(b) 本 run 是 tol=1e-10 的平衡点；两者之差是**代数残差**，远小于空间离散误差（尖端 $4.26\times10^{-4}$ m 量级），也仍在 $10^{-8}$ m 门槛内。

## 5. S6a 细节

### 5.1 逐步结果

| step | 1 体 迭代 / LS / stopReason | 1 体 $g/g_s$ | 4 体 迭代 / LS / stopReason | 4 体 $g/g_s$ |
|---:|---|---:|---|---:|
| 1 | 2 / 0 / combined | 9.8e-17 | 2 / 0 / combined | 2.4e-17 |
| 2 | 2 / 0 / combined | 4.5e-16 | 2 / 0 / combined | 1.1e-16 |
| 3 | 2 / 0 / combined | 1.0e-15 | 2 / 0 / combined | 2.5e-16 |
| 4 | 2 / 0 / combined | 2.2e-15 | 2 / 0 / combined | 5.5e-16 |
| 5 | 3 / 0 / combined | 1.1e-18 | 3 / 0 / combined | 2.7e-19 |
| 6 | 5 / 47 / combined | 5.3e-15 | 3 / 1 / combined | 6.8e-16 |
| 7 | 3 / 0 / combined | 2.7e-18 | 6 / 28 / combined | 1.6e-15 |
| 8 | 5 / 41 / combined | 8.3e-15 | 3 / 0 / combined | 9.1e-19 |
| 9 | 6 / 44 / combined | 9.2e-15 | 5 / 41 / combined | 2.0e-15 |
| 10 | 6 / 27 / combined | 1.0e-14 | 3 / 0 / combined | 1.4e-18 |

- 1 体：总 36 迭代 / 159 线搜索，10/10 converged，末态最大位移 $1.3471948779\times10^{-2}$ m，末态差 $1.91\times10^{-10}$ m；与 `stopcriterion` 诊断的 combined 结果逐位一致（迭代 $[2,2,2,2,3,5,3,5,6,6]$、线搜索 159 全同）。
- 4 体：总 31 迭代 / 70 线搜索，10/10 converged，末态最大位移 $1.3471948779\times10^{-2}$ m，末态差 $1.82\times10^{-10}$ m。
- **4 体与独立单体的逐步迭代数分布不同**（如第 6 步 3 vs 5、第 7 步 6 vs 3），原因是两者线性求解的浮点路径不同，$\lVert\Delta x\rVert_\infty$ 在 $10^{-12}$ 量级两侧摆动，`actual` 条件成立的迭代号随之变化。这不是物理差异：逐步轨迹差从第 6 步起约 $5\times10^{-11}$–$2\times10^{-10}$ m，末态 $2.04\times10^{-10}$ m，远低于 $10^{-8}$ m。
- 两组的停止时相对梯度都 $\le1.02\times10^{-14}$（1 体）/ $1.97\times10^{-15}$（4 体），离 $10^{-6}$ 有 8 个量级以上；`min J` 最小 0.98878，无翻转 tet；全部量有限；退出码 0。

### 5.2 块对角与独立解对照（4 体）

| 检查 | 值 | 门槛 |
|---|---:|---|
| 末态跨体 Hessian triplet | 0 | = 0 |
| 末态组装后跨体 entry | 0 | = 0 |
| 跨体 tet | 0 | = 0 |
| 节点按 body 连续 | true | true |
| 组装非零元 | 772344 | — |
| 4 体 vs 4 次独立单体（逐步最大） | $2.04\times10^{-10}$ m | $<10^{-8}$ m |
| 4 体无接触：bbox 最小间隙 / 表面节点最小距离 | 0.5016 m / 0.5023 m | $>0$ |

4 次独立单体各自 10/10 converged，逐步迭代数均为 $[2,2,2,2,3,5,3,5,6,6]$（总 36），与 1 体 case 完全一致。

## 6. 假阴性消除与 false positive

- **假阴性**：S2 M3 原 `maxIter`（100 迭代，`converged=0`）→ 6 迭代 `combined`；S6a 1 体原本第 6/8/9 步 `maxIter` → 全部 `combined`；4 体同样 10/10 `combined`。S2 M6 本次没有重跑 step 判据（未在预登记矩阵内），但它停在方向率冻结于 16.6×tol 的状态，按 M3 已证实的同一机制，原 step 判据在该状态不会满足；这条是机制推断，不是本 run 的实测。
- **false positive**：按预登记定义（停止时 $g_{\mathrm{free}}/g_{\mathrm{scale}}\ge10^{-6}$），四个 case 的全部停止点最大值为 $5.00\times10^{-9}$（S2 M3），**没有出现提前误停**。
- 需要区分：S6a 4 体第 6 步只用 3 次迭代就停（独立求解要 5 次），但停止时相对梯度 $6.8\times10^{-16}$、末态差 $1.82\times10^{-10}$ m，属“已经到浮点不动点”，不是“梯度仍大就停”。

## 7. 成本（只作数量级参照）

| case | solver 墙钟 (s) | 进程墙钟 | 峰值 RSS |
|---|---:|---:|---:|
| S2 M3 | 1.09 | 1.09 | 97.7 MiB |
| S2 M6 | 186.6 | 186.5 | 2.46 GiB |
| S6a 1 体 | 2.68 | 2.68 | 48.0 MiB |
| S6a 4 体（含块对角 + 4 次独立单体） | 10.67（全局 4 体） + 10.44（4× 独立） + 0.37（块对角） | 24.32 | 157.2 MiB |

M6 的 2.46 GiB 与 186 s 来自 75625 节点 / 414720 tet 的 CHOLMOD 分解，与判据无关；作为对照，同一网格在 v0.2 step 判据、tol=1e-8 下是 3 迭代 / 92.8 s（参考 run），本 run 因 tol=1e-10 多做 4 次迭代、其中 2 次带深度线搜索（二分 18/29 次），所以更慢。**未做重复测量、CPU 绑定或频率控制，墙钟不是性能结论。**

## 8. 边界

- 覆盖范围：S2 两个网格（M3 粗、M6 细）、S6a 1 体与 4 体各 10 步；单一材料/单一资产；无接触、无长时动态、无大规模批次。
- `actualStepTol=1e-12` m 只对米级场景确认；坐标 ulp 约 $2.2\times10^{-16}$ m，差 4 个量级，换到 $10^3$ m 量级需按 $L_\mathrm{char}$ 重标定。
- `gradScale` 用“本场景总重”：S2 静态取 $24.525$ N（虽然静态梁的自由梯度只有 $10^{-7}$ N 量级，比值极小）；S6a 1 体 $5086.259$ N、4 体 $20345.036$ N、独立单体 $5086.259$ N。无重力/驱动场景（S3t）不在本次范围。
- S2 精度门槛不变：M3 仍 25.94%、M6 仍 4.26%，`combined` 只消除代数停滞，不改善离散误差；M3 的 25.94% 仍是**未通过** S2 精度门槛的状态，与 v0.2 一致。
- M6 的参考末态本身是 tol=1e-8 的平衡点，$3.21\times10^{-9}$ m 的差包含参考的代数残差；它是四个 case 里余量最小的一项（3.1 倍），若以后换更粗的参考容差需重算这条门槛。
- 本 run 不覆盖 `--stop-mode` 之外的行为（能量、装配、线搜索、DBC、时间积分均未改）。

## 9. 与执行计划的偏差（如实记录，现场保存在 `_attempt1/`、`_attempt2/`）

1. **attempt1（17:53）**：M6 在 0.2 s 内被 `cantilever_check` 的**资产自检**拒绝（`S2 ERROR: 体积或质量不符`，exit 2）。原因：414720 个 tet 的 double 顺序累加体积为 $0.0025000000000128$，相对偏差 $5.1\times10^{-12}$，触发 $10^{-12}$ 相对门槛。这是输入网格的浮点累加误差，不是物理/收敛失败。处理：改用 long double 独立求和（与 `s2-extension-20261031-01` 生成 M6 参考时的资产自检相同），**门槛保持 $10^{-12}$ 不变**；未改任何验收门槛、未调物理参数。
2. **attempt2（17:58）**：四个 case 全部通过，但发现 4 体 job 内部的 1 体诊断场景沿用了 4 体总重作为 `gradScale`（相对梯度被缩小 4 倍），且 `summary.json` 的尺度标签错误。处理：按“本场景总重”缩放并修正标签；4 体全局结果与 1 体 case 数值不变。
3. **最终 run（18:01–18:05）**：四个 case 使用同一二进制重新完整跑一遍，`binary_sha256.txt` 与所有 case 一一对应。attempt1/attempt2 的现场原样保留，未删除。

没有出现 NaN/Inf（四个 case 的所有 CSV 数值单元非有限值计数 0），没有 `lineSearch`/`linearSolve` 停止，没有超时，没有为通过而放宽验收门槛或调参重跑。未 `git commit` / `git push`。

## 10. run-id 证据清单

根目录：`Archive/fiction/output/validation/s2s6a-combined-20261031-01/`

| 路径 | 内容 |
|---|---|
| `env.txt` | 主机、构建、网格、三个 gradScale、线程变量、执行顺序 |
| `configure.log` / `build.log` / `ctest.log` | 配置（CHOLMOD=ON）、构建（0 warning）、CTest 10/10 |
| `asset_sha256.txt` / `binary_sha256.txt` | 三个输入资产与两个可执行文件哈希 |
| `source_files_sha256.txt` / `source.tar.gz` / `source_hash.txt` | 源码快照与哈希 |
| `source-change.diff` | 相对 `fiction-stopcriterion-20261031-01` 的全部源码改动（111 行，仅两个 example） |
| `cases.csv` | 四个 case 一行汇总（含门槛判定列） |
| `judgment.json` | 逐 case 预登记门槛、末态差、停止判据、false positive 结论 |
| `s2/M3/` | `config.json`、`command.txt`、`M3.csv`、`M3-state.csv`、`summary.json`、`run.log`、`time_rss.txt`、`exit_code.txt`、`trace/iterations.csv(.nodes.csv)` |
| `s2/M6/` | 同上（`M6.csv`、`M6-state.csv`、`trace/`） |
| `s6a/body1/` | `config.json`、`command.txt`、`asset_audit.json`、`body1_steps.csv`、`body1_state.csv`、`summary.json`、`run.log`、`time_rss.txt`、`exit_code.txt` |
| `s6a/body4/` | 同上（`body4_*`）＋ `block_diagonal.json`、`trajectory_compare.csv`、`independent_body{0..3}_steps.csv`、`program_summary.json`（程序原生汇总，未覆盖） |
| `_attempt1/`、`_attempt2/` | 偏差说明 §9 的两次前序现场（保留，不参与判定） |

源码副本：`/home/steven/mjx/fiction/experiments/workspaces/fiction-s2s6a-combined-20261031-01/`（含 `run_combined.sh`、`tools/summarize_combined.py`）。构建目录 `/tmp/fiction-s2s6a-combined-build`。

## 11. 结论

1. **候选判据 `combined` 在 S2 M3/M6 + S6a 1/4 体上重验收通过**：物理解不变（最大 $3.21\times10^{-9}$ m $<10^{-8}$ m）、每步/整步 converged、块对角 0 跨体、4 体 vs 独立解 $2.04\times10^{-10}$ m、无提前误停（最大相对梯度 $5.0\times10^{-9}<10^{-6}$）。
2. S2 精度门槛与结论不变（M3 25.94% 未过、M6 4.26% 通过）；本 run 不改变任何 S2 精度结论。
3. **这不等于协议 v0.3 生效**：`combined` 只是“候选判据在 S2/S6a 上重验收通过”。是否成为 v0.3 默认判据、是否接受 $1e-12$ m / 2 次 / $10^{-6}$ / 总重这套默认值，仍需按草案 §5 由用户批准；S3t、S3、接触（S6b）、非米级场景与后续批次不在本次范围。
