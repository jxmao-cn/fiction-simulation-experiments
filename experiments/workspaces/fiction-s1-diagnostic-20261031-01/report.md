# S1 单因素受控诊断：增刚后重算方向（A）

日期：2026-10-31。执行范围：隔离副本 `baseline/` 与 `variant-a/`；原工程 `/home/steven/mjx/fiction/simulator` 全程只读。输出目录保留供审查，不合并回原工程。

## 结论

假说 A 在 M2 的指定失败载荷 1.02 上得到支持：原版复现 100 次迭代上限失败；仅在罚刚度翻倍后重算当前目标方向的 Variant-A 在 1.02 完整通过 S1 所有判据。M2 后续载荷 1.05 仍在 100 次迭代上限失败，因此 A 不是完整 M2/S1 修复，不能把这一结果表述为全载荷通过。按约定 1.05 失败后停止，没有运行 M3，也未叠加其他参数或算法变化。

Variant-A 的 M1 四个载荷均 PASS。两个副本的 CTest 均为 9/9 通过；其中两个 GUI 项是 build-only 占位，不代表交互式 GUI 测试。

## 实验设计与边界

A 唯一算法改动在 `cpp/src/time_integrator.cpp`：检测到方向率不高于 tol、DBC 尚未全满足而将 `dbc.stiffness` 翻倍后，立即重算当前罚目标能量、搜索方向 `p` 和 `satisfied`，同步报告的刚度和 solver；若新线性求解失败、方向非有限或能量非有限则退出。其他物理、静态模式、`tol=1e-10`、`maxIter=100`、初始 `k=1e8`、线搜索和 tension_check 的完整验收判据均未改变。

两套构建均为 Release、CHOLMOD enabled，使用 `-j3`。网格直接只读原记录中的 M1/M2；SHA-256 分别为 `049aac54366f09dfcac4a46cf1898c940611aa3dfe85cfa12f23666862303702` 与 `38877c77e6c7b4a174ea1804ef3dc57993bfff0df8c51de2e3d8cc12859bdb0d`。每次 S1 命令设 `timeout 120s`。此次工具编排没有保存单个进程的精确 elapsed-time 字段；exitcode 文件修改时间只能给完成时间，不能冒充运行耗时。所有命令、CSV、stdout 和 exit code 均已留存。

## 结果摘要

| 副本 / 网格 | load | iter | 收敛/停止情况 | DBC 满足 | 关键完整验收 | exit |
|---|---:|---:|---|---:|---|---:|
| baseline / M2 | 1.02 | 100 | `converged=0`, `maxIterHit=1`, `lineSearchFailed=0` | 1/50 | `dbcMaxError=7.785367e-6`；能量相对误差 `8.0306e-4`；反力相对误差 `4.1281e-4`；`pass=0` | 1 |
| Variant-A / M1 | 1.01 | 18 | converged | 18/18 | CSV `pass=1` | 0 |
| Variant-A / M1 | 1.02 | 20 | converged | 18/18 | CSV `pass=1` | 0 |
| Variant-A / M1 | 1.05 | 23 | converged | 18/18 | CSV `pass=1` | 0 |
| Variant-A / M2 | 1.01 | 21 | converged | 50/50 | CSV `pass=1` | 0 |
| Variant-A / M2 | 1.02 | 24 | converged | 50/50 | `dbcMaxError=9.0171e-11`；能量相对误差 `1.39146e-8`；反力相对误差 `7.14555e-9`；端点误差 `9.0171e-11`；全体 J 最大参考误差 `1.16208e-9`；CSV `pass=1` | 1* |
| Variant-A / M2 | 1.05 | 100 | `converged=0`, `maxIterHit=1`, `lineSearchFailed=0` | 1/50 | `dbcMaxError=4.8020e-6`；能量相对误差 `1.9572e-4`；反力相对误差 `9.9502e-5`；最大 J 参考误差 `4.0240e-5`；CSV `pass=0` | 1 |

`*` M2 程序采用“首个失败即停止”语义，因此在 1.05 失败后进程 exit 1。CSV 中此前 1.02 的 `pass=1` 仍是该载荷完整验收通过，不应误读为整次 M2 命令成功。

Variant-A 相比 baseline M2 1.02 在报告中显著改善：方向/边界满足（50/50），不再耗尽迭代，物理和几何验收误差均低于原阈值；但 1.05 仍显示高载荷阶段问题。这不是在完全相同的 1.01 状态上分叉的单轮实验：两实现经各自求解得到的 warm-start 状态并非逐字节相同。因此，现有观测支持“在本次加载路径上，增刚后刷新方向有益”，但不足以证明旧方向是原 1.02 失败的唯一根因，也不证明数值收敛或全局收敛性。

## 回归与回滚核验

- baseline CTest：9/9 passed。
- Variant-A CTest：9/9 passed。
- `fix-a.patch` 的逆向 dry-run 与实际逆向应用在独立 `rollback-check/` 副本完成；恢复后 `time_integrator.cpp` SHA-256 与 baseline 完全相同：`84ee9de8ed6118c748d5625b1ddbea7f4b632cb909ac7f8c22978553a9216dc6`。
- 原工程源码没有被改动，也不需要对原工程回滚。Variant-A 保留补丁以便审查；没有清理诊断目录。

## 已知限制

1. 仅验证指定的 M1/M2 网格和静态 tension_check 路径；M3 未执行，动态求解与其他材料/网格未覆盖。
2. 载荷 1.05 的失败并未逐迭代记录新旧方向导数、alpha、实际坐标增量或零位移接受情况，故不能据此确认/排除假说 B，也不能把失败归因到某个确定子机制。
3. CTest 中 GUI 项只验证目标构建；未启动窗口。完整 S1 验收以 CSV `pass` 列为准，不能只看 `converged`。
4. 当前 CLI 仅接受 `--mesh` 和 `--csv`，加载序列及物理参数来自 `tension_check.cpp` 编译常量；本次保持其规格原样。
5. 进程精确运行耗时没有被本轮 shell wrapper 记录；结果目录有 stdout/CSV/exit code，exitcode 文件时间戳是完成时间，不是耗时。

## 产物索引

- `fix-a.patch`：仅 A 的相对路径补丁。
- `results/baseline/M2.csv`、`M2.stdout.txt`、`M2.exitcode`：原版 M2 复现。
- `results/baseline/ctest.txt`、`ctest.exitcode`：原版回归。
- `results/variant-a/M1.csv`、`M1.stdout.txt`、`M1.exitcode`：Variant-A M1。
- `results/variant-a/M2.csv`、`M2.stdout.txt`、`M2.exitcode`：Variant-A M2。
- `results/variant-a/ctest.txt`、`ctest.exitcode`：Variant-A 回归。
- `reproduce.sh`：可复制的配置、构建与复现命令。
- `rollback.md`：审查后撤销补丁的方式与 SHA 校验。
- `comparison.json`：机器可读对照摘要。
