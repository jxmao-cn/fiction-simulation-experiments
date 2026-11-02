# fiction-simulation-experiments

fiction 三维固体仿真器的源码与实验工作区，涵盖材料与时间积分检查、均匀拉伸、悬臂弯曲、扭转、非杆几何、规模和无接触多体。仓库当前为私有研究快照，便于查看已完成实验使用的实现与配置。

## 从这里开始

| 目录 | 内容 |
|---|---|
| `simulator/cpp/` | 当前 C++20 主工程：能量、隐式 Euler、Projected Newton、网格与场景、可执行样例和小型测试输入 |
| `simulator/Python/` | 历史 Python 参考实现；原先使用的导入名见其中说明 |
| `simulator/tools/` | 结构化杆网格生成、S2 运行与绘图等工具 |
| `experiments/workspaces/` | 按 run 分隔的实验源码和脚本；每个目录保留该实验的独立代码，不覆盖主工程 |
| `simulator/docs/` | S0/T0 至 S6 相关规格、报告与诊断记录 |
| `planning/` | 历史验证协议 |
| `assets/ipc/` | torus、sphere1K 小型网格及原始来源记录；其他大网格需单独准备 |
| `SOURCE_MANIFEST.json` | 本次导入各文件的原始位置和 SHA-256 |

本快照从本地 `fiction` 工作区整理，生产工程原 Git HEAD 为 `b3745128119eb5def111159dbd3310a9e31fcc6f`；快照也包含当时未提交或未跟踪的工作文件，因此不是该 Git commit 的原样检出。不同工作区的 `cpp/` 是各阶段的隔离版本；原始实验复现以相应 run 的源码快照、参数和输入为准。

## 构建主工程

需要 C++20、CMake 3.20+ 和 Eigen 3.4+；SuiteSparse/CHOLMOD 可选。无图形环境可关闭 viewer：

```bash
cmake -S simulator/cpp -B build -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_VIEWER=OFF
cmake --build build -j 3
ctest --test-dir build --output-on-failure
```

viewer、PNG 和视频功能依构建选项分别使用 GLFW/OpenGL、zlib 与 ffmpeg。不同阶段的 runner、命令及设计见 [实验索引](EXPERIMENTS.md)。脚本以普通文本文件导入，运行 shell 脚本时可使用 `bash path/to/script.sh`。历史脚本部分保留旧机器的绝对路径和构建位置；使用前应按自己的环境检查，不应覆盖原始记录。

## 输入与第三方内容

结构化方杆可由 `simulator/tools/make_structured_bar.py` 生成。torus 和 sphere1K 来自 IPC 项目，原始来源与哈希见 `assets/ipc/sources.json`；`simulator/cpp/tests/data/` 中的 rod 网格由 ElastoGen 示例转换。Armadillo13K、28K、54K 网格不包含在此源码快照中，文件哈希和本地原位置见 [实验索引](EXPERIMENTS.md)。

隔离实验中探索接触的源码保留在对应工作区；外部 IPC Toolkit、Catch2 源码与构建缓存不包含在仓库中。完整原始轨迹、日志、可执行文件、源码压缩快照、视频和博客草稿留在本地 fiction 项目。`simulator/docs/` 中的历史路径和 run-id 作为实验记录保留，部分链接需要回到原始工作区访问。

仓库尚未指定项目源码许可证。部分网格来自第三方，`assets/ipc/UPSTREAM_LICENSE` 是上游仓库许可副本；公开发布之前需逐项确认项目许可证与网格授权。
