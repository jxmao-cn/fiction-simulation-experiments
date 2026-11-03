# fiction-simulation-experiments

fiction 三维固体仿真器的源码与实验工作区，涵盖材料与时间积分检查、均匀拉伸、悬臂弯曲、扭转、非杆几何、规模实验和无接触多体场景。仓库当前为私有研究快照，便于查看已记录实验所用的实现与配置。

## 从这里开始

| 阅读或浏览 | 内容 |
|---|---|
| [实验文章（中文）](blog/fiction.md) · [English](blog/fiction-en.md) | 含插图的叙述性介绍；图像资源共用 [`blog/assets/`](blog/assets/) |
| [博客文件说明](blog/README.md) | 文章、图片资源、可编辑插图和来源记录的组织方式 |
| [实验索引](EXPERIMENTS.zh-CN.md) · [English](EXPERIMENTS.md) | 按阶段浏览代码与记录 |
| [`simulator/cpp/`](simulator/cpp/) | 当前 C++20 主工程：能量、隐式 Euler、Projected Newton、网格与场景、可执行样例和小型测试输入 |
| [`simulator/Python/`](simulator/Python/) | 历史 Python 参考实现；原先使用的导入名见其中说明 |
| [`simulator/tools/`](simulator/tools/) | 结构化杆网格生成、S2 运行与绘图等工具 |
| [`experiments/workspaces/`](experiments/workspaces/) | 按 run 分隔的实验源码和脚本；每个目录保留该实验的独立代码，不覆盖主工程 |
| [`simulator/docs/`](simulator/docs/) | S0/T0 至 S6 相关规格、报告与诊断记录 |
| [`planning/`](planning/) | 历史验证协议 |
| [`assets/ipc/`](assets/ipc/) | torus、sphere1K 小型网格及其来源记录；其他大网格需单独准备 |
| [`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) | 本次导入文件的原始位置和 SHA-256 |

本快照从本地 `fiction` 工作区整理。生产工程原 Git HEAD 为 `b3745128119eb5def111159dbd3310a9e31fcc6f`。快照也包含当时未提交或未跟踪的工作文件，因此不是该 Git commit 的原样检出。不同工作区中的 `cpp/` 是各阶段的隔离版本；复现原实验时，应以相应 run 的源码快照、参数和输入为准。

## 构建主工程

需要 C++20、CMake 3.20+ 和 Eigen 3.4+；SuiteSparse/CHOLMOD 可选。无图形环境时可以关闭 viewer：

```bash
cmake -S simulator/cpp -B build -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_VIEWER=OFF
cmake --build build -j 3
ctest --test-dir build --output-on-failure
```

viewer、PNG 和视频功能依构建选项分别使用 GLFW/OpenGL、zlib 与 ffmpeg。各阶段的 runner、命令及设计见[实验索引](EXPERIMENTS.zh-CN.md)。脚本以普通文本文件导入；运行 shell 脚本时可使用 `bash path/to/script.sh`。历史脚本部分保留旧机器的绝对路径和构建位置；使用前应按自己的环境检查，不要覆盖原始记录。

## 输入与第三方内容

结构化方杆可由 [`simulator/tools/make_structured_bar.py`](simulator/tools/make_structured_bar.py) 生成。torus 和 sphere1K 来自 IPC 项目，原始来源与哈希见 [`assets/ipc/sources.json`](assets/ipc/sources.json)。`simulator/cpp/tests/data/` 中的 rod 网格由 ElastoGen 示例转换。Armadillo13K、28K、54K 网格不包含在此源码快照中；文件哈希和本机历史输入位置见[实验索引](EXPERIMENTS.zh-CN.md)。

隔离实验中探索接触的源码保留在对应工作区；外部 IPC Toolkit、Catch2 源码与构建缓存不包含在仓库中。完整原始轨迹、日志、可执行文件和冻结源码快照留在本地 fiction 项目；选入文章的图和视频现已包含在 `blog/` 中。`simulator/docs/` 中的历史路径和 run-id 作为实验记录保留，部分链接需要回到原始工作区访问。

仓库尚未指定项目源码许可证。部分网格来自第三方，[`assets/ipc/UPSTREAM_LICENSE`](assets/ipc/UPSTREAM_LICENSE) 是上游仓库许可副本。本仓库没有在此作出项目许可证或网格授权方面的新决定。

## 核对文章文件

```bash
python tools/check_blog.py
```

检查中英文加粗、数学表达式、图号、段首格式、相对链接和随文章分发文件的哈希；不运行仿真。
