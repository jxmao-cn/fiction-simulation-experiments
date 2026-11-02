# S6a 单体：主工程 combined 判据 10 步重验收

日期：2026-11-02。结论：**本次单体10步重验收通过**。仅覆盖无接触sphere1K与本次固定配置；不改写原S6a失败记录，不代表四体/16体或完整S6a阶段通过。生产源码没有在本轮修改。

## 运行与核心来源

规格在运行前登记于 [s6a-single-main-spec-20261102.md](s6a-single-main-spec-20261102.md)。run-id为 [s6a-single-main-20261102-01](../output/validation/s6a-single-main-20261102-01/)。

测试harness复用既有隔离副本的场景构造辅助头文件和runner观察逻辑，加入失败即停、物理力残差/归一化量输出和每步trace归档。它通过CMake `add_subdirectory` 编译当前生产 `simulator/cpp`，链接**新构建的 `production/libfiction_core.a`**；生产include先于harness辅助include。保存的 `header-dependencies.txt` 确认使用生产 `time_integrator.h`，没有使用旧隔离核心；`link-command.txt` 确认实际静态库来源。源码、资产、二进制及核心库哈希均已记录。

Release/CHOLMOD独立构建 `/tmp/fiction-s6a-single-main-20261102-build`；生产子目录CTest **11/11通过**，含2项GUI占位。没有执行四体求解。

## 固定配置与每步结果

1760节点、6851四面体、顶部局部z>=0.45的63个节点固定；NH E=1e5 Pa、nu=0.3、rho=1000，重力-z9.81 m/s²。dt=0.005 s，10步共0.05 s；tol=1e-8，maxIter=100，DBC k=1e8。显式combined，gradTolRel=1e-6、gradScale=5086.259120998311 N、actualStepTol=1e-12 m、连续2次。

| 时间步 | Newton次数 | 物理力残差最大分量 N | 实际成立的停止出口 |
|---:|---:|---:|---|
| 1 | 2 | 1.9862e-8 | 原step快路径 |
| 2 | 2 | 9.2011e-8 | 原step快路径 |
| 3 | 2 | 2.0287e-7 | 原step快路径 |
| 4 | 2 | 4.4955e-7 | 原step快路径 |
| 5 | 3 | 2.1835e-10 | 原step快路径 |
| 6 | 5 | 1.0732e-6 | 连续小步+相对力残差备用出口 |
| 7 | 3 | 5.4690e-10 | 原step快路径 |
| 8 | 5 | 1.6923e-6 | 连续小步+相对力残差备用出口 |
| 9 | 6 | 1.8658e-6 | 连续小步+相对力残差备用出口 |
| 10 | 6 | 2.0853e-6 | 连续小步+相对力残差备用出口 |

10/10步均收敛、合计36次Newton；`maxIterHit`和`lineSearchFailed`均未触发。所有状态有限，min J最低0.98877877055，DBC误差全程0，63/63约束满足。最大物理残差为2.0853251e-6 N，除单体重后为4.0999191e-10，低于本批1e-6门槛。动态CSV的物理力残差等于raw IP梯度除dt²，已重新核对。

每行`stop_reason=combined`表示所选判据，不表示10步全由备用出口结束。根据方向率及逐Newton trace，第1–5和7步满足原step条件；第6、8、9、10步方向率仍高于tol，依靠连续2次实际最大分量位移<=1e-12 m且力残差已小的备用条件结束。这四步的末两次实际位移均可在独立归档trace中查到，不能用整个时间步位移替代。

## 完整末态比较

| 对照 | 最大节点欧氏位置差 | 门槛 |
|---|---:|---:|
| 原S6a单体（含旧maxIter步骤的末态） | 1.906284176e-10 m | <1e-8 m |
| corrected隔离副本单体 | 0 m | <1e-8 m |

完整末态与隔离修正结果逐节点相同，支持主工程接入后的行为一致。旧原始结果仅作数值一致性对照，并非独立物理真解。

## 证据与限制

运行目录保存spec/config、source.tar.gz、源码与资产hashes、二进制/核心库哈希、实际命令、编译/CTest日志、完整state.csv、10行steps.csv、summary.json、时间/RSS。每步iteration和nodes trace都保存为独立 `iterations.csv.step-N.csv` 及 `.nodes.csv`，避免原核心每步覆盖同一路径后只剩末步。

`state_change_inf`字段继承旧runner名称，其值实际是整个时间步节点欧氏位移的最大值；逐Newton最大分量位移在trace的`actual_dx_inf`。二者定义不同。本轮计时含trace和状态归档开销，未作重复性能测量。

本批执行失败即停，10步均通过所以未触发停止。范围仅一个同质物体、短时段、固定物理参数。它支持本次单体重验收通过，不支持异质体、接触、长时稳定或连续体精度结论。本小步完成后停止，交由用户指导下一项。
