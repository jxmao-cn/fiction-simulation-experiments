# 停止判据合入主工程与小规模确认

日期：2026-11-02。用户接受限定病例的停止判据修正，并授权进行下一步。本轮将可选判据和诊断接口合入生产工作区，默认仍为step；代码处于本地待提交状态，原S2/S6a验收状态保持原记录。

## 合入范围

共六个文件，完整范围见 [summary.json](../output/validation/stopcriterion-merge-20261102-01/summary.json)：

- `cpp/include/time_integrator.h`、`cpp/src/time_integrator.cpp`：可选step/grad/actual/combined模式、停止原因、raw IP梯度和物理力残差；动态力残差除以dt²，combined阈值读取配置。
- `cpp/examples/cantilever_check.cpp`：增加模式/残差阈值/力尺度/实际步长阈值/连续小步次数CLI，CSV同时记录停止原因、raw IP梯度、物理力残差[N]和配置值；材料、几何、验收容差保持原设置。
- `cpp/examples/stopcriterion_units_check.cpp`、`cpp/CMakeLists.txt`：加入实际求解器量纲与错误早停负例。
- `Python/time_integrator.py`：新增兼容原调用的关键词选项、动态残差换算、停止原因及可选报告；支持可选迭代上限与非有限失败。原默认调用保持step模式及原[x,v]返回接口。

不复制多体runner或其他场景实现。既有未提交文件均保留；快照逐文件核对未发现本轮六文件之外的修改。

## 验证

Release/CHOLMOD独立构建 `/tmp/fiction-stop-merge-20261102-build` 完成。CTest **11/11通过**，其中2项为GUI占位；另直接运行simulator全部既有场景通过，包含2D-in-3D回归。

新增C++负例验证：dt=0.005/0.01/0.02下物理残差均为1 N、raw随dt²改变；大残余力与极小实际步长不能从combined备用出口误停；配置阈值确实生效；四份无耦合复制不改变该测试的残差判断。

Python通过语法检查和实际step_forward控制逻辑测试：时间步换算、错误早停拒绝、阈值生效、四种模式NaN拒绝、非法配置拒绝、默认step初始满足时零迭代。该控制逻辑测试使用真实函数和mock能量/搜索适配器。旧Python保留独立2D接触/摩擦实现及`Python_slovement`导入约定。

| 主工程M3病例 | 迭代 | 停止原因 | 与原始M3完整末态最大节点差 |
|---|---:|---|---:|
| 默认step | 100 | maxIter，退出码1，预期复现 | 0 m |
| 显式combined，gradTolRel=1e-6 | 6 | combined，退出码0 | 1.274092324e-12 m |

combined的尖端挠度7.40750850694 mm，梁参考误差为25.94%。

## 复现命令

从fiction工作区根执行，输出应使用新目录，避免覆盖冻结结果：

```bash
cmake -S simulator/cpp -B /tmp/fiction-stop-merge-reproduce -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/fiction-stop-merge-reproduce -j 3
ctest --test-dir /tmp/fiction-stop-merge-reproduce --output-on-failure
/tmp/fiction-stop-merge-reproduce/cantilever_check \
  --mesh simulator/output/validation/s2-cantilever-20261031-01/meshes/M3.msh \
  --csv /tmp/fiction-stop-combined.csv --state /tmp/fiction-stop-combined-state.csv \
  --stop-mode combined --grad-tol 1e-6 --grad-scale 24.525 \
  --actual-step-tol 1e-12 --actual-step-count 2
```

去掉模式参数即可复现默认step，预期该原始M3病例退出码1，对应原始maxIter状态。

## 判据配置与证据

combined的step快路径按方向率判定，力残差门槛用于停滞备用出口。grad和actual保留为独立诊断模式，分别按相对梯度与实际小步判定。

各生产demo默认模式保持原设置。调用combined时须配置单位、特征力与物体归一化；异质体、接触及不同几何尺度需分别标定gradScale。物理确认使用上表两组M3病例。

证据根为 [stopcriterion-merge-20261102-01](../output/validation/stopcriterion-merge-20261102-01/)：保存本轮前源码快照、仅本轮merge-only.patch、合后六文件快照与哈希、求解二进制、构建/CTest/simulator日志、Python测试脚本与日志、两组M3 CSV/完整末态和独立差值summary。回滚应先确认没有后续修改，再针对merge-only.patch检查逆向应用；不要用git reset覆盖既有工作。

本小步完成后停止，等待用户指导下一项。
