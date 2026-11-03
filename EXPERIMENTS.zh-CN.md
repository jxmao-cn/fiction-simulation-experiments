# 实验索引

本索引将内部阶段编号对应到实验主题和代码快照。编号用于追溯记录；下表同时明确说明各阶段检查的内容。

| 阶段编号 | 实验内容 | 代码与记录 |
|---|---|---|
| S0 / T0 | Neo-Hookean 材料导数与单元装配；线性振子的隐式 Euler 时间更新 | [`material_validation.cpp`](simulator/cpp/examples/material_validation.cpp)、[`dynamic_reference_check.cpp`](simulator/cpp/examples/dynamic_reference_check.cpp)；相关记录见 [`simulator/docs/`](simulator/docs/) |
| S1 | 杆的均匀拉伸与 Dirichlet 边界条件（DBC）罚刚度诊断 | [`tension_check.cpp`](simulator/cpp/examples/tension_check.cpp)；隔离运行见 [S1 诊断工作区](experiments/workspaces/fiction-s1-diagnostic-20261031-01/) 和 [S1 双组工作区](experiments/workspaces/fiction-s1-two-group-20261031-01/) |
| S2 | 悬臂静态弯曲与空间网格加密 | [`cantilever_check.cpp`](simulator/cpp/examples/cantilever_check.cpp)；[S2 扩展工作区](experiments/workspaces/fiction-s2-extension-20261031-01/) |
| S3 | 悬臂动态弯曲与时间步加密 | [原 S3 工作区](experiments/workspaces/fiction-s3-20261031-01/)、[S3 时间步与网格工作区](experiments/workspaces/fiction-s3-dtmesh-20261102-01/)、[S3 代数精度工作区](experiments/workspaces/fiction-s3-algebra-main-20261102-01/) |
| S3t | 规定扭转运动与视觉样稿 | [S3/S3t 视觉工作区](experiments/workspaces/fiction-s3-s3t-visual-main-20261102-01/) |
| S4 | 非杆状 torus 几何的变形 | [S4 工作区](experiments/workspaces/fiction-s4-20261031-01/) |
| S5 | Armadillo 网格规模 pilot 与较大网格运行 | [S5 pilot 工作区](experiments/workspaces/fiction-s5-20261031-01/)、[54K 主运行工作区](experiments/workspaces/fiction-s5-54k-main-20261102-01/) |
| S6a | 无接触单体、四体和 16 体场景 | [原 S6a 工作区](experiments/workspaces/fiction-s6a-20261031-01/)、[单体](experiments/workspaces/fiction-s6a-single-main-20261102-01/)、[四体](experiments/workspaces/fiction-s6a-four-main-20261102-01/)和 [16 体](experiments/workspaces/fiction-s6a-sixteen-main-20261102-01/)工作区 |
| S6b | 隔离工作区中的接触探索 | [C0 工作区](experiments/workspaces/fiction-s6b-c0-20261102-01/)、[C1 工作区](experiments/workspaces/fiction-s6b-c1-20261102-01/)；仓库不包含外部 IPC Toolkit |
| 停止判据 | 动态物理力残差和求解器停止规则的隔离诊断 | [初始诊断](experiments/workspaces/fiction-stopcriterion-20261031-01/)、[修正后诊断](experiments/workspaces/fiction-stopcriterion-corrected-20261031-01/)和 [combined 判据对照](experiments/workspaces/fiction-s2s6a-combined-20261031-01/) |

[`experiments/workspaces/`](experiments/workspaces/) 中的快照只包含选定的程序文件、头文件、构建配置、脚本、轻量文本和小型测试网格，不包含各工作区的构建目录、二进制文件、第三方 checkout 或原始结果。历史脚本可能依赖本机路径和未随仓库分发的大型网格。完整原始数据、源码压缩快照和实验报告仍保存在原 fiction 本地工作区的 `simulator/output/validation/` 下；本仓库不替代这些原件。

## 未随仓库分发的大型网格

| 文件 | SHA-256 | 原机器上的历史输入位置 |
|---|---|---|
| Armadillo13K.msh | `03beaf3d4c9692a10664a3a90f1baa8c10d43d7ba47126437d420b68a36e6bce` | `fiction/assets/ipc/tetMeshes/Armadillo13K.msh` |
| Armadillo28K.msh | `73395043488c6c1a99117ce750e4a7ac1f2f3d8849ccd8f18be277a99c0e359b` | `fiction/assets/ipc/tetMeshes/Armadillo28K.msh` |
| Armadillo54K.msh | `e1faa13af8543e01bda3620ae83fe9a1e5ddff1280861d5896b850a84493d520` | `Archive/IPC/input/tetMeshes/Armadillo54K.msh` |

对应场景沿用原研究中的网格坐标和指定加载方向。运行时请核对 `simulator/docs/` 下相应的规格与实验脚本。
