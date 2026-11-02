# 协议 v0.3 修订草案：Newton 停止判据（2026-10-31）

状态：**草案，待用户批准。** 协议 v0.2、已有 run-id 和生产 `Archive/fiction/cpp` 仍是当前权威；本文件不构成协议通过，也不改写任何既有结果。

## 1. 问题

S2 M3（tol=1e-10）和 S6a 1 体（tol=1e-8）都出现同一假阴性：方向率冻结在 tol 之上，实际位置变化降到远低于坐标 ulp，继续迭代不改变状态，但 `converged=0`。

| 证据 | S2 M3 | S6a 1 体 |
|---|---|---|
| 冻结方向率 | 1.05e-9（tol 的 10.5 倍） | 1.7e-8–3.4e-8（tol 的 1.7–3.4 倍） |
| 实际位置变化 | 4.1e-25 m | 约 1e-16 m 以下 |
| 自由梯度 | 1.23e-7 N | 2.7e-11–5.9e-11 N |
| 继续迭代的效果 | 1000 次与 100 次末态逐位相同 | maxIter 3 vs 100 差 6.1e-16 m |

机制：绝对步长范数判据 + 近奇异弯曲模态 + 能量比较分辨率，三者共同造成浮点不动点。不是材料、装配、边界或线搜索错误。

## 2. 提议的判据

在 `TimeIntegrator::Config` 增加：

```text
stopMode = combined            // step | grad | actual | combined
actualStepTol = 1e-12 m        // 米级场景；非米级需按 L_char 缩放并另验
actualStepCount = 2
gradTolRel = 1e-6
gradScale = 场景特征外力        // 重力加载场景取总重；驱动场景另定并记录
```

`combined` 的停止条件：

$$
\left(\frac{\lVert p\rVert_\infty}{h}\le \mathrm{tol}\right)
\ \ \text{或}\ \
\left(\text{连续 } k \text{ 次 }\lVert \Delta x\rVert_\infty\le \varepsilon_x
\ \wedge\ \frac{\lVert g_{\mathrm{free}}\rVert_\infty}{g_{\mathrm{scale}}}\le \varepsilon_g\right)
$$

并额外要求 DBC 满足、全部量有限、未触发 `maxIter`/`lineSearch`/`linearSolve` 失败。`StepReport` 必须记录 `stopReason`。

**不单独采用 `grad`**：诊断中 grad 在 S6a 上每步只 1 次迭代就停，末态差 3.2e-6 m（门槛 1e-8 m）；原因是位移误差是 $H^{-1}g$，软体近奇异模态把小的相对梯度放大成不可忽略的位移残差，动态步再被 $1/h$ 放大。`actual` 单独也可行（S6a 末态差 4.95e-11 m），但迭代数 55 vs combined 36；combined 保留原 step 快路径，只在停滞步启用出口。

**不改变**：能量、装配、线搜索、DBC 罚、时间积分、`tol` 在 step 快路径中的语义、所有精度门槛。

## 3. 适用范围与标定

- 已在 S2 M3（静态、米级）与 S6a 1 体（动态、米级）上诊断通过；这是有限证据，不是普遍证明。
- 非米级场景（如 Armadillo 151 网格单位）不能直接复用 `1e-12 m`；需按 $L_\mathrm{char}$ 缩放 `actualStepTol`，并在代表性 case 上重新确认。
- 无重力/驱动场景（如 S3t）的 `gradScale` 需用初始自由梯度范数或驱动力矩尺度，单独记录；本草案不自动覆盖。
- 接触阶段（S6b）另行设计；接触集每迭代变化，不能直接套用本判据。

## 4. 候选重验收（新 run-id，需批准后执行）

| 阶段 | 候选配置 | 对照 | 通过条件 |
|---|---|---|---|
| S2 | M3、M6 静态，`combined`，tol=1e-10（M3）/1e-8（M6） | v0.2 M3 strict、M6 tol=1e-8 末态 | 末态最大节点差 <1e-8 m；M3 误差仍 25.94%、M6 仍 4.26%；S2 精度门槛不变 |
| S6a | 1 体、4 体，`combined`，tol=1e-8，10 步 | v0.2 S6a 1/4 体末态 | 每步 converged；块对角 0 跨体；4 体 vs 独立解 <1e-8 m |
| S3t/S3 | 不在本候选内 | — | 两阶段主跑已全部 converged；S3 的 dt 判据问题另议 |

候选重验收使用 `stopcriterion-20261031-01` 的副本和 `combined` 实现，新建 run-id `s2s6a-combined-20261031-01`；保存候选协议、实现 diff、CTest、config、CSV、state、summary、时间/RSS。**结果只能写成“候选判据在 S2/S6a 上重验收通过/不通过”，不能写成协议 v0.3 已生效。**

## 5. 批准项

1. 是否接受 `combined` 作为 v0.3 的默认停止判据。
2. 是否接受 `actualStepTol=1e-12 m`、`actualStepCount=2`、`gradTolRel=1e-6`、`gradScale=总重` 的默认值。
3. 是否批准候选重验收矩阵（S2 M3/M6 + S6a 1/4 体）。
4. 是否允许在候选 run-id 中复用 `stopcriterion-20261031-01` 的隔离副本与 `combined` 实现。

任一项未批准，不进入候选重验收；生产代码与协议 v0.2 不动。
