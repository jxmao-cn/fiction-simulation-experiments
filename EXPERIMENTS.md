# 实验代码导航

| 阶段 | 内容 | 代码入口 |
|---|---|---|
| S0 / T0 | 材料导数、装配、线性振子 | `simulator/cpp/examples/material_validation.cpp`、`dynamic_reference_check.cpp` |
| S1 | 均匀拉伸、DBC 罚参数诊断 | `simulator/cpp/examples/tension_check.cpp`；`experiments/workspaces/fiction-s1-*` |
| S2 | 静态悬臂与网格加密 | `simulator/cpp/examples/cantilever_check.cpp`；`experiments/workspaces/fiction-s2-extension-20261031-01/` |
| S3 / S3t | 动态弯曲、扭转、时间步与视觉样稿 | `experiments/workspaces/fiction-s3-20261031-01/`、`fiction-s3-dtmesh-20261102-01/` 和 `fiction-s3-s3t-visual-main-20261102-01/` |
| S4 | torus 圆环 | `experiments/workspaces/fiction-s4-20261031-01/` |
| S5 | Armadillo 规模 pilot | `experiments/workspaces/fiction-s5-20261031-01/`、`fiction-s5-54k-main-20261102-01/` |
| S6a | 单体、四体、16体无接触场景 | `experiments/workspaces/fiction-s6a-20261031-01/` 和 `fiction-s6a-{single,four,sixteen}-main-20261102-01/` |
| 停止条件 | 动态物理力残差与停止规则的隔离诊断 | `experiments/workspaces/fiction-stopcriterion-20261031-01/`、`fiction-stopcriterion-corrected-20261031-01/`、`fiction-s2s6a-combined-20261031-01/` |

`experiments/workspaces/` 仅复制了程序、头文件、构建配置、脚本、轻量文本和小型测试网格，不包含其中的构建目录、二进制、第三方 checkout 或原始结果。旧脚本可能依赖本地路径和未包含的大型网格。完整数据、源码压缩快照和实验报告的历史原件仍位于原 fiction 工作区的 `simulator/output/validation/`；此仓库不覆盖它们。

## 未随仓库分发的大网格

| 文件 | SHA-256 | 机器上的历史输入位置 |
|---|---|---|
| Armadillo13K.msh | `03beaf3d4c9692a10664a3a90f1baa8c10d43d7ba47126437d420b68a36e6bce` | `fiction/assets/ipc/tetMeshes/Armadillo13K.msh` |
| Armadillo28K.msh | `73395043488c6c1a99117ce750e4a7ac1f2f3d8849ccd8f18be277a99c0e359b` | `fiction/assets/ipc/tetMeshes/Armadillo28K.msh` |
| Armadillo54K.msh | `e1faa13af8543e01bda3620ae83fe9a1e5ddff1280861d5896b850a84493d520` | `Archive/IPC/input/tetMeshes/Armadillo54K.msh` |

对应场景在原研究中沿用网格坐标和指定加载方向，运行时应核对各自 `simulator/docs/` 规格与实验脚本。
