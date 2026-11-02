# 停止判据诊断（预登记，不改生产）

日期：2026-10-31。目的：S2 M3 与 S6a 1 体都出现“方向率冻结在 tol 之上、实际位置不再变化”的假阴性。本诊断在隔离副本中比较三种停止判据，确认哪一种能在**不改变物理解**的前提下消除假阴性。**这是诊断，不是重验收；生产 `Archive/fiction/cpp`、协议 v0.2 和已有 run-id 均不改动。**

## 固定病例

| 病例 | 来源 | 参数 | 基线 |
|---|---|---|---|
| S2 M3 静态 | `Archive/fiction/output/validation/s2-cantilever-20261031-01/meshes/M3.msh` | NH E=5.90e8、ν=0.3、ρ=1000、重力 -z、x=0 全端面固定、staticSolve、tol=1e-10、maxIter=100 | 该 run 的 M3 末态：尖端 7.40750851 mm，方向率 1.0546e-9，自由梯度 1.2264e-7 N |
| S6a 1 体 | `fiction-s6a-20261031-01` 的 `s6a_multibody_check --only-body1` | sphere1K、NH E=1e5、ν=0.3、ρ=1000、重力 -z、顶部 cap 固定、dt=0.005、10 步、tol=1e-8、maxIter=100 | 该 run 的 body1：第 6/8/9 步 maxIter=100；末态最大位移 1.3472e-2 m |

力尺度：S2 `gradScale = 24.525 N`（总重）；S6a 1 体 `gradScale = 518.4769746175648 × 9.81 = 5086.259 N`。

## 三种判据

在隔离副本 `fiction-stopcriterion-20261031-01` 的 `TimeIntegrator` 中增加配置项，不改变能量、装配、线搜索或 DBC：

1. **grad**：`freeGradientInf / gradScale <= 1e-7` 且 DBC 满足。
2. **actual**：连续 2 次迭代 `actualDxInf <= 1e-12 m` 且 DBC 满足。
3. **combined**：原判据 `||p||_inf/h <= tol` **或**（actual 条件 **且** `freeGradientInf / gradScale <= 1e-6`）。

每次停止在 `StepReport` 中写 `stopReason`（step|grad|actual|combined|maxIter|lineSearch|linearSolve）。`freeGradientInf` 的算法与现有末态诊断一致：总梯度减去 DBC 罚项，只统计非 DBC 节点。

## 运行

- 每个病例 × 每个判据（grad / actual / combined）+ 原判据基线，共 8 个 run；每个 run 独立目录，保存 config、命令、CSV、末态、summary、时间/RSS。
- S2 M3 用 `cantilever_check` 的副本（静态）；S6a 1 体用 `s6a_multibody_check --only-body1` 的副本（10 步）。
- 构建 Release/CHOLMOD，先跑 CTest；线程变量按既有约定设为 1。
- 不调参数、不改门槛；如果某判据提前停止导致物理解变化，如实报告为 false positive，不放宽比较阈值。

## 判定

对每个病例：

1. **物理解不变**：变体末态与基线末态的最大节点位置差 `< 1e-8 m`（S2 尖端差 `< 1e-8 m`）。
2. **假阴性消除**：原 maxIter=100 的步在变体下 `converged=1`，且 `stopReason` 不是 `maxIter`。
3. **不提前误停**：变体停止时的 `freeGradientInf / gradScale < 1e-6`；若 actual 判据在梯度仍大时停止，标为 false positive。
4. **成本**：记录迭代数、线搜索次数、墙钟；不把并发墙钟当性能结论。

推荐结论只能是“某判据在本诊断中通过/不通过”，不能直接写成协议通过；协议修订需另立版本并在新 run-id 下重新验收。

## 证据

run-id：`stopcriterion-20261031-01`。报告：`Archive/fiction/docs/stopcriterion-diagnostic-results-20261031.md`。原始数据根：`Archive/fiction/output/validation/stopcriterion-20261031-01/`。
