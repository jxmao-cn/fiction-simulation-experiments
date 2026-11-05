# 实验代码与结果导航

本页区分“当前检出可直接运行的主工程”和“保留历史版本的实验工作区”。构建与小网格样例从 [README](README.md#快速开始) 开始；不要直接批量执行历史 runner。

## 阶段与入口

| 阶段 | 内容与代码入口 | 报告及结论范围 |
|---|---|---|
| S0 / T0 | `simulator/cpp/examples/material_validation.cpp`、`dynamic_reference_check.cpp` | [材料与积分检查](simulator/docs/stage0-validation-20261030.md)；相关小测试接入当前 CTest |
| S1 | `simulator/cpp/examples/tension_check.cpp`；`experiments/workspaces/fiction-s1-*` | [拉伸重验](simulator/docs/s1-reacceptance-20261031.md)；对应指定罚参数日程 |
| S2 | `simulator/cpp/examples/cantilever_check.cpp`；`fiction-s2-extension-20261031-01/` 等工作区 | [空间/时间离散报告](simulator/docs/s2-s3-refinement-report-20261102.md)；区分原 M3 记录、停止判据修正和后续网格加密 |
| S3 / S3t | `fiction-s3-20261031-01/`、`fiction-s3-dtmesh-20261102-01/`、`fiction-s3-s3t-visual-main-20261102-01/` | [时间步对照](simulator/docs/s3-dtmesh-results-20261102.md)、[补充时间步](simulator/docs/s3-dtmesh-000625-addendum-20261102.md)、[扭转](simulator/docs/s3t-twist-results-20261031.md)；原 0–2 s、2% 轨迹对照未通过，后续趋势不改判原记录 |
| S4 | `fiction-s4-20261031-01/` | [圆环实验](simulator/docs/s4-nonrod-results-20261031.md)；能力展示，不是解析精度验证 |
| S5 | `fiction-s5-20261031-01/`、`fiction-s5-54k-main-20261102-01/` | [13K/28K](simulator/docs/s5-scale-results-20261031.md)、[54K](simulator/docs/s5-54k-main-results-20261102.md)；各 10 步成本记录，不是严格 benchmark |
| S6a | `fiction-s6a-20261031-01/`、`fiction-s6a-{single,four,sixteen}-main-20261102-01/` | [单体](simulator/docs/s6a-single-main-results-20261102.md)、[四体](simulator/docs/s6a-four-main-results-20261102.md)、[16 体](simulator/docs/s6a-sixteen-main-results-20261102.md)；无接触短段与独立单体对照 |
| S6b（隔离探索） | `fiction-s6b-c0-20261102-01/`、`fiction-s6b-c1-20261102-01/` | C0 依赖/工具单元检查完成，C1 稳态间隙与反力门槛未通过，后续阶段搁置；此选集未包含 S6b 结果报告和外部依赖 checkout |
| 停止条件 | `fiction-stopcriterion-20261031-01/`、`fiction-stopcriterion-corrected-20261031-01/`、`fiction-s2s6a-combined-20261031-01/` | [修正审阅](simulator/docs/stopcriterion-corrected-review-20261102.md)、[主工程合入](simulator/docs/stopcriterion-merge-20261102.md)；默认 step 与可选 combined 不混为一个判据 |

表中仅写目录名的工作区都位于 [`experiments/workspaces/`](experiments/workspaces/) 下；通配符和花括号表示一组目录，不是可点击文件路径。

## 哪些可以直接运行？

- **主工程构建、CTest 与 README 小网格样例**：使用仓库内代码和输入，可直接运行；CHOLMOD、查看器和视频功能按需安装依赖。
- **结构化杆网格生成**：`simulator/tools/make_structured_bar.py` 只依赖 Python 标准库；先用 `--help` 查看参数，在自己的项目 Python 环境中运行。
- **阶段 runner**：它们是原执行现场的记录。部分脚本硬编码 `/home/steven/...` 或 `/tmp/fiction-...-build`，也可能假设结果目录和第三方源码已经存在。重跑前应检查路径、构建、输入和输出目录，使用新的 run-id，不能覆盖原始结果或直接沿用通过/失败结论。
- **历史 Python 参考实现**：保留旧包名 `Python_slovement`，目前不能从 `simulator/Python/` 直接按旧 README 命令运行；不作为本仓库的支持入口。
- **S6b 依赖探索**：IPC Toolkit、Catch2 等外部源码没有随选集分发，其 CMake 不能在没有这些依赖时直接配置成功。

## 未随仓库分发的内容

`experiments/workspaces/` 保留了程序、头文件、构建配置、脚本、轻量文本和部分小网格，不含构建目录、二进制、第三方 checkout 或完整原始结果。完整轨迹、日志、源码压缩快照、视频及部分图表仍在原本地工作区。历史报告中指向 `output/validation/`、`figures/` 的链接可能无法打开，且其中行号对应历史源码；这是证据分发范围，不能用新生成的文件冒充旧 run。

### 大网格

| 文件 | SHA-256 | 历史输入位置（非仓库内现存文件） |
|---|---|---|
| Armadillo13K.msh | `03beaf3d4c9692a10664a3a90f1baa8c10d43d7ba47126437d420b68a36e6bce` | `fiction/assets/ipc/tetMeshes/Armadillo13K.msh` |
| Armadillo28K.msh | `73395043488c6c1a99117ce750e4a7ac1f2f3d8849ccd8f18be277a99c0e359b` | `fiction/assets/ipc/tetMeshes/Armadillo28K.msh` |
| Armadillo54K.msh | `e1faa13af8543e01bda3620ae83fe9a1e5ddff1280861d5896b850a84493d520` | `Archive/IPC/input/tetMeshes/Armadillo54K.msh` |

这三项均不包含于本仓库。对应场景沿用模型坐标与指定加载方向，准备输入后仍需核对各阶段规格和脚本。现有来源记录只证明原实验输入的身份，不是完整下载/授权说明。
