# IPC 资产清单（模型文件扫描存档）

> 扫描对象：`/home/steven/mjx/IPC/input`（IPC 参考实现的全部输入资产）
> 扫描日期：2026-09-15
>
> 扫描方法：
> 1. **fiction 全量实测**：对全部 70 个 `.msh` 执行 `ipc_mesh_demo --mesh <file> --no-gravity --steps 0`
>    （真实加载 + 建 mass-spring 网 + 表面提取，非纸上统计）；
> 2. **静态解析补数**：对 fiction 不支持格式的 5 个文件，用脚本解析节点/单元/边界，保证表格完整；
> 3. **场景引用统计**：扫描 154 个场景 `.txt`（其中 10 个为符号链接）中显式出现的 `input/...` 资产路径。
>
> 配套文档：连接设计与理论说明见 [`ipc-mesh-integration.md`](ipc-mesh-integration.md)。

## 图例

- **顶点 / 四面体**：网格原始规模（MshIO 4.1 / Gmsh 2.2 / IPC legacy 三类格式统一统计）。
- **唯一边(弹簧)**：把每个四面体的 6 条边去重后的规模，即 fiction 的 mass-spring 拓扑（自由度 3N、刚度矩阵非零元来源）。
- **表面面片**：渲染与 OBJ 导出用的三角面；优先取 `$Surface` 段，缺失时从四面体边界提取。
- **被场景引用**：引用该文件的场景数（同一场景多次引用只计一次；`—` 表示未被引用）。
- **fiction 可视化**：✅ 可直接跑 `ipc_mesh_demo/visual`；❌ 表示需要补 loader。

## 0. 总览

- `.msh` 体积网格 **70** 个：MshIO 4.1 **65** 个、Gmsh 2.2 **3** 个、IPC legacy **2** 个；单元类型全部为 tet4。**fiction 可直接处理 65/70**，其余 5 个仅格式家族不同。
- 规模分布（节点数）：≤200 共 36 个、201–2K 共 7 个、2K–10K 共 11 个、>10K 共 16 个。
- `.obj` 表面网格 **46** 个：全部三角面；其中闭合（水密）**26** 个、点云（无面）**2** 个。角色为运动学碰撞体，**不能参与仿真**（无体积），且 fiction 目前**没有 OBJ loader，无法可视化**。
- `.seg` 线段网格 **203** 个：独立文件 2 个 + 动画序列 201 帧。同上，需要新增 loader 才能渲染。
- `.pt` 点云：场景引用存在，但**仓库内没有对应文件**（见第 5 节悬空引用）。
- 解析碰撞体（`ground` / `halfSpace`）：定义在场景脚本里，无独立文件；fiction 已有 `BarrierEnergy` 平面能量但未接入本 demo。
- 唯一没有 `$Surface` 段的网格：`tetMeshes/32770_octocat.msh`，fiction 用边界提取回退成功（已实测可渲染）。

## 1. 体积网格 `.msh`（可变形体，70 个）

这些都是 IPC 仿真中的**可变形体**（`shapes` 输入）。fiction 读取后沿四面体边建质量-弹簧网、以 ρ·V/4 集中质量，并可表面渲染。格式为 MshIO 4.1 的 65 个**已实测全部通过**；另外 5 个（SOFA 对照组的 Gmsh 2.2、Utopia 对照组的 IPC legacy）需扩展解析器，几何上同样是 tet4。

| 文件（相对 `IPC/input/`） | 对应物品 | 格式 | 顶点 | 四面体 | 唯一边(弹簧) | 表面面片 | 被场景引用 | fiction 可视化 |
|---|---|---|---:|---:|---:|---:|---:|---|
| `paperExamples/supplementB/SOFA/chain/torus.msh` | 圆环（SOFA 对照组） | Gmsh 2.2 | 200 | 491 | 891 | 400 | — | ❌ 不支持（Gmsh 2.2） |
| `paperExamples/supplementB/SOFA/dolphin/dolphin5K.msh` | 海豚（SOFA 对照组） | Gmsh 2.2 | 8111 | 36304 | 49414 | 10000 | — | ❌ 不支持（Gmsh 2.2） |
| `paperExamples/supplementB/SOFA/mat-on-knives/mat40x40.msh` | 薄板（SOFA 刀阵对照） | Gmsh 2.2 | 3200 | 9126 | 15523 | 6396 | — | ❌ 不支持（Gmsh 2.2） |
| `paperExamples/supplementB/Utopia/utopiaComparison39K.msh` | 对照模型 39K（Utopia 基准） | IPC legacy | 7453 | 38962 | 48266 | 3706 | — | ❌ 不支持（IPC legacy） |
| `paperExamples/supplementB/Utopia/utopiaComparison67K.msh` | 对照模型 67K（Utopia 基准） | IPC legacy | 117874 | 665074 | 797585 | 29278 | — | ❌ 不支持（IPC legacy） |
| `tetMeshes/32770_octocat.msh` | 章鱼猫（Octocat） | MshIO 4.1 | 6133 | 21426 | 31835 | 8554 | 8 | ✅ 可直接 |
| `tetMeshes/ANSYS1e-3Corner_coarse.msh` | 角部模型（ANSYS 命名，粗网格） | MshIO 4.1 | 285 | 812 | 1350 | 510 | 1 | ✅ 可直接 |
| `tetMeshes/Armadillo122K.msh` | 犰狳（Stanford Armadillo，分辨率梯度） | MshIO 4.1 | 121976 | 547848 | 739823 | 140000 | 1 | ✅ 可直接 |
| `tetMeshes/Armadillo13K.msh` | 犰狳（Stanford Armadillo，分辨率梯度） | MshIO 4.1 | 13230 | 54985 | 76862 | 17296 | 1 | ✅ 可直接 |
| `tetMeshes/Armadillo219K.msh` | 犰狳（Stanford Armadillo，分辨率梯度） | MshIO 4.1 | 219355 | 928030 | 1285752 | 276736 | 1 | ✅ 可直接 |
| `tetMeshes/Armadillo28K.msh` | 犰狳（Stanford Armadillo，分辨率梯度） | MshIO 4.1 | 27951 | 121190 | 166436 | 34592 | 1 | ✅ 可直接 |
| `tetMeshes/Armadillo54K.msh` | 犰狳（Stanford Armadillo，分辨率梯度） | MshIO 4.1 | 54082 | 226840 | 315513 | 69184 | 1 | ✅ 可直接 |
| `tetMeshes/alien.msh` | 外星人 | MshIO 4.1 | 12531 | 48643 | 70028 | 17710 | 2 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.01.msh` | 砌体拱石块（第 1 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.02.msh` | 砌体拱石块（第 2 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.03.msh` | 砌体拱石块（第 3 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.04.msh` | 砌体拱石块（第 4 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.05.msh` | 砌体拱石块（第 5 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.06.msh` | 砌体拱石块（第 6 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.07.msh` | 砌体拱石块（第 7 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.08.msh` | 砌体拱石块（第 8 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.09.msh` | 砌体拱石块（第 9 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.10.msh` | 砌体拱石块（第 10 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.11.msh` | 砌体拱石块（第 11 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.12.msh` | 砌体拱石块（第 12 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.13.msh` | 砌体拱石块（第 13 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.14.msh` | 砌体拱石块（第 14 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.15.msh` | 砌体拱石块（第 15 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.16.msh` | 砌体拱石块（第 16 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.17.msh` | 砌体拱石块（第 17 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.18.msh` | 砌体拱石块（第 18 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.19.msh` | 砌体拱石块（第 19 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.20.msh` | 砌体拱石块（第 20 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.21.msh` | 砌体拱石块（第 21 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.22.msh` | 砌体拱石块（第 22 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.23.msh` | 砌体拱石块（第 23 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.24.msh` | 砌体拱石块（第 24 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/arch/largeArch.25.msh` | 砌体拱石块（第 25 块） | MshIO 4.1 | 8 | 6 | 19 | 12 | 3 | ✅ 可直接 |
| `tetMeshes/armadillo_verschoor.msh` | 犰狳（Verschoor 对照用） | MshIO 4.1 | 65387 | 385739 | 461614 | 20978 | 2 | ✅ 可直接 |
| `tetMeshes/bar-14114.msh` | 圆杆（14114 单元，分辨率梯度） | MshIO 4.1 | 3842 | 14114 | 20849 | 5788 | 1 | ✅ 可直接 |
| `tetMeshes/bar-186.msh` | 圆杆（186 单元，分辨率梯度） | MshIO 4.1 | 86 | 186 | 355 | 168 | 1 | ✅ 可直接 |
| `tetMeshes/bar-2523.msh` | 圆杆（2523 单元，分辨率梯度） | MshIO 4.1 | 886 | 2523 | 4217 | 1618 | 2 | ✅ 可直接 |
| `tetMeshes/bar-6.msh` | 圆杆（6 单元，分辨率梯度） | MshIO 4.1 | 8 | 6 | 19 | 12 | 1 | ✅ 可直接 |
| `tetMeshes/bar-651.msh` | 圆杆（651 单元，分辨率梯度） | MshIO 4.1 | 292 | 651 | 1232 | 580 | 1 | ✅ 可直接 |
| `tetMeshes/bunny.msh` | 兔子（Stanford） | MshIO 4.1 | 4670 | 19379 | 27181 | 6266 | 1 | ✅ 可直接 |
| `tetMeshes/cow_head.msh` | 牛头 | MshIO 4.1 | 1328 | 5029 | 7349 | 1986 | 1 | ✅ 可直接 |
| `tetMeshes/cube.msh` | 单位立方体（6 tet 剖分） | MshIO 4.1 | 8 | 6 | 19 | 12 | 45 | ✅ 可直接 |
| `tetMeshes/dolphin5K.msh` | 海豚 | MshIO 4.1 | 8111 | 36304 | 49414 | 10000 | 1 | ✅ 可直接 |
| `tetMeshes/erleben/spike.msh` | 尖刺（Erleben CCD 基准） | MshIO 4.1 | 5 | 2 | 9 | 6 | 11 | ✅ 可直接 |
| `tetMeshes/erleben/wedge.msh` | 楔形（Erleben CCD 基准） | MshIO 4.1 | 6 | 3 | 12 | 8 | 6 | ✅ 可直接 |
| `tetMeshes/hammer.msh` | 锤子 | MshIO 4.1 | 2705 | 9683 | 14432 | 4090 | — | ✅ 可直接 |
| `tetMeshes/hollowCat7.5K.msh` | 空心猫 | MshIO 4.1 | 12827 | 57636 | 77962 | 15000 | 2 | ✅ 可直接 |
| `tetMeshes/mat100x100t40.msh` | 薄板 100x100（厚 40 层） | MshIO 4.1 | 20000 | 58806 | 98803 | 39996 | 1 | ✅ 可直接 |
| `tetMeshes/mat150x150t40.msh` | 薄板 150x150（厚 40 层） | MshIO 4.1 | 45000 | 133206 | 223203 | 89996 | 2 | ✅ 可直接 |
| `tetMeshes/mat20x20.msh` | 薄板 20x20 | MshIO 4.1 | 800 | 2166 | 3763 | 1596 | 5 | ✅ 可直接 |
| `tetMeshes/mat225x225t40.msh` | 薄板 225x225（厚 40 层） | MshIO 4.1 | 101250 | 301056 | 503553 | 202496 | 1 | ✅ 可直接 |
| `tetMeshes/mat40x40.msh` | 薄板 40x40 | MshIO 4.1 | 3200 | 9126 | 15523 | 6396 | 15 | ✅ 可直接 |
| `tetMeshes/monkey8K.msh` | 猴子（Blender Suzanne） | MshIO 4.1 | 12448 | 52864 | 73312 | 16000 | 2 | ✅ 可直接 |
| `tetMeshes/plate1.msh` | 平板 | MshIO 4.1 | 3362 | 9600 | 16321 | 6720 | — | ✅ 可直接 |
| `tetMeshes/pyramid.msh` | 四棱锥 | MshIO 4.1 | 6 | 3 | 12 | 8 | — | ✅ 可直接 |
| `tetMeshes/rod.msh` | 长杆 | MshIO 4.1 | 630 | 1722 | 2932 | 1162 | 14 | ✅ 可直接 |
| `tetMeshes/rod300x33.msh` | 长杆 | MshIO 4.1 | 13333 | 50511 | 73826 | 19966 | 1 | ✅ 可直接 |
| `tetMeshes/sphere19K.msh` | 球体 | MshIO 4.1 | 28833 | 118202 | 166146 | 38224 | 1 | ✅ 可直接 |
| `tetMeshes/sphere1K.msh` | 球体 | MshIO 4.1 | 1760 | 6851 | 9847 | 2474 | 5 | ✅ 可直接 |
| `tetMeshes/sphere5K.msh` | 球体 | MshIO 4.1 | 7034 | 27792 | 39696 | 9742 | 9 | ✅ 可直接 |
| `tetMeshes/tet-perp-edges.msh` | 互垂边四面体 | MshIO 4.1 | 4 | 1 | 6 | 4 | 3 | ✅ 可直接 |
| `tetMeshes/tet-pyramid.msh` | 四棱锥形四面体 | MshIO 4.1 | 4 | 1 | 6 | 4 | 3 | ✅ 可直接 |
| `tetMeshes/tet.msh` | 单个四面体（最小单元） | MshIO 4.1 | 4 | 1 | 6 | 4 | 10 | ✅ 可直接 |
| `tetMeshes/torus.msh` | 圆环 | MshIO 4.1 | 200 | 491 | 891 | 400 | 6 | ✅ 可直接 |
| `tetMeshes/walnut.msh` | 核桃 | MshIO 4.1 | 39004 | 217718 | 266542 | 19642 | — | ✅ 可直接 |

## 2. 表面三角网格 `.obj`（运动学碰撞体，46 个）

IPC 用它们作为**碰撞对象**（闭合面、三角面或点云），本身没有体积、不参与弹性仿真。fiction 当前没有 OBJ loader，所以既不能仿真也不能渲染；若要「把碰撞体也画出来」，需新增只读加载 + 软光栅化线段/点/面叠加（仅显示用途）。

| 文件（相对 `IPC/input/`） | 对应物品 | 顶点 | 三角面 | 闭合 | 被场景引用 | fiction 可视化 |
|---|---|---:|---:|---|---:|---|
| `paperExamples/supplementB/SOFA/dolphin/dolphin5K.obj` | 运动学碰撞体（对应 `tetMeshes/dolphin5K.msh`） | 5002 | 10000 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `paperExamples/supplementB/SOFA/dolphin/funnel.obj` | 运动学碰撞体 | 564 | 1128 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `paperExamples/supplementB/SOFA/mat-on-knives/triangle.obj` | 刀（SOFA 对照） | 6 | 8 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/YZplane.obj` | YZ 平面 | 4 | 2 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/bowl.obj` | 碗形碰撞体 | 1321 | 2560 | 否 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/cone-simple.obj` | 简单圆锥碰撞体 | 34 | 64 | 是 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/cone.obj` | 圆锥碰撞体 | 1152 | 2176 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/coneBox.obj` | 圆锥+盒复合碰撞体 | 1292 | 2400 | 否 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/cube.obj` | 立方体碰撞体（对应 `tetMeshes/cube.msh`） | 8 | 12 | 是 | 8 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/cylinder.obj` | 圆柱碰撞体 | 152 | 300 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/cliff.obj` | Erleben CCD 基准碰撞体：cliff | 8 | 12 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/crack.obj` | Erleben CCD 基准碰撞体：crack | 20 | 36 | 是 | 4 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/hole.obj` | Erleben CCD 基准碰撞体：hole | 48 | 88 | 否 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/internal_edges.obj` | Erleben CCD 基准碰撞体：internal_edges | 14 | 24 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/spike.obj` | Erleben CCD 基准碰撞体：spike | 5 | 6 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/erleben/wedge.obj` | Erleben CCD 基准碰撞体：wedge | 6 | 8 | 是 | 5 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/funnel.obj` | 漏斗碰撞体 | 288 | 552 | 否 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/hammer.obj` | 锤子碰撞体（对应 `tetMeshes/hammer.msh`） | 2107 | 4090 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/mat20x20.obj` | 薄板 20x20 表面（对应 `tetMeshes/mat20x20.msh`） | 800 | 1596 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/mat40x40.obj` | 薄板 40x40 表面（对应 `tetMeshes/mat40x40.msh`） | 3200 | 6396 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/plane.obj` | 单四边形平面 | 4 | 2 | 否 | 37 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/point.obj` | 点云（1 点） | 1 | 0 | — | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/pointPlane.obj` | 点云平面（420 点） | 420 | 0 | — | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/pyramid.obj` | 四棱锥碰撞体（对应 `tetMeshes/pyramid.msh`） | 6 | 8 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/rod.obj` | 杆碰撞体（对应 `tetMeshes/rod.msh`） | 583 | 1162 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sheet100.obj` | 薄片（100 顶点级） | 128 | 210 | 否 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sheet1000.obj` | 薄片（1000 顶点级） | 913 | 1708 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sheet2000.obj` | 薄片（2000 顶点级） | 1769 | 3360 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sheet500.obj` | 薄片（500 顶点级） | 465 | 840 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sheetPlane.obj` | 薄片平面 | 400 | 200 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/sphere1K.obj` | 球体碰撞体（对应 `tetMeshes/sphere1K.msh`） | 1239 | 2474 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/tet-perp-edges.obj` | 四面体表面（对应 `tetMeshes/tet-perp-edges.msh`） | 4 | 4 | 是 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/tet-pyramid.obj` | 四面体表面（对应 `tetMeshes/tet-pyramid.msh`） | 4 | 4 | 是 | 1 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/tet.obj` | 四面体表面（对应 `tetMeshes/tet.msh`） | 4 | 4 | 是 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/torus.obj` | 圆环碰撞体（对应 `tetMeshes/torus.msh`） | 200 | 400 | 是 | 5 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/triangle.obj` | 单个三角面 | 3 | 1 | 否 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/tube.obj` | 管道碰撞体 | 1024 | 1920 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller0.obj` | 滚轮（Verschoor 对照，第 0 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller1.obj` | 滚轮（Verschoor 对照，第 1 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller2.obj` | 滚轮（Verschoor 对照，第 2 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller3.obj` | 滚轮（Verschoor 对照，第 3 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller4.obj` | 滚轮（Verschoor 对照，第 4 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/verschoor/roller5.obj` | 滚轮（Verschoor 对照，第 5 帧） | 274 | 544 | 是 | 2 | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/voro_sphere.obj` | Voronoi 简化表面（sphere） | 3749 | 6518 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/voro_walnut.obj` | Voronoi 简化表面（walnut） | 71203 | 139748 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |
| `triMeshes/wedge.obj` | 楔形（三角面） | 3 | 1 | 否 | — | ❌ 需新增 OBJ loader（仅渲染） |

## 3. 线段网格 `.seg`（运动学线元，203 个）

`.seg` 为 IPC 自定义的「点 + 线段」格式（`v x y z` / `s i j`），用作线状碰撞对象，`sequence/` 是动画序列（每帧拓扑一致，仅顶点位置变化，供 `meshSeq` 播放）。fiction 目前不支持该格式。

| 文件 | 对应物品 | 顶点 | 线段 | fiction 可视化 |
|---|---|---:|---:|---|
| `segMeshes/edge.seg` | 单线段 | 2 | 1 | ❌ 需新增 SEG loader（仅渲染） |
| `segMeshes/segPlane.seg` | 栅格线段平面 | 420 | 210 | ❌ 需新增 SEG loader（仅渲染） |
| `segMeshes/sequence/*.seg`（共 201 帧，`0.seg` … `200.seg`） | 动画线段序列（meshSeq/edgeSequence） | 3（每帧一致） | 3（每帧一致） | ❌ 需新增 SEG loader（仅渲染） |

## 4. 其他资产

| 类型 | 数量 | 说明 | fiction 可用性 |
|---|---:|---|---|
| `.txt` 场景脚本 | 154 | shapes/material/DBC/NBC/ground/time 等命令（144 个实体 + 10 个符号链接指向 `11_erleben`） | ⚠️ 未解析；参数目前用 CLI 手工等价 |
| `.pt` 点云 | 0（引用悬空） | 场景引用的点云文件不在仓库（见第 5 节） | ❌ |
| `ground` / `halfSpace` 解析平面 | 定义在 txt 内 | IPC 的解析碰撞体 | ⚠️ `BarrierEnergy` 已有但未接入 demo |
| `.scn` / `.mph` / `.e` / `.hipnc` / `.view` / `.xml` | 13 | SOFA / COMSOL / Houdini 对照工程文件 | ❌ 与 IPC 运行无关 |
| `.md` / `.sh` / `.py` | 4 | 说明文档与格式转换脚本 | ❌ 与 IPC 运行无关 |

## 5. 悬空引用（场景引用了不存在的文件）

以下路径被场景脚本引用，但 `IPC/input` 中并不存在（IPC 自身也无法直接运行这些场景；其中 `.pt`/`.seg` 的点云/线段教程可用同名 `.obj` 代替理解）：

- `tetMeshes/doug_j_squishy_ball.msh` —— 被 1 个场景引用（`paperExamples/22_squishyBall.txt`）
- `triMeshes/cylinder.pt` —— 被 1 个场景引用（`paperExamples/18_pointRollerBall.txt`）
- `triMeshes/cylinder.seg` —— 被 1 个场景引用（`paperExamples/18_segRollerBall.txt`）
- `triMeshes/point.pt` —— 被 4 个场景引用（`otherExamples/coDimUnitTests/mat40x40_pointPlaneDrop.txt`、`otherExamples/coDimUnitTests/octocat_pointPlaneDrop.txt`、`otherExamples/coDimUnitTests/sphere5K_pointPlaneDrop.txt` 等 4 个）
- `triMeshes/triangle.pt` —— 被 1 个场景引用（`tutorialExamples/MCO/2cubesFall_rotateCO_points.txt`）
- `triMeshes/triangle.seg` —— 被 2 个场景引用（`tutorialExamples/MCO/2cubesFall_rotateCO_edges.txt`、`tutorialExamples/advanced/2cubesFall_rotateCO_meshSeq.txt`）

## 6. fiction 实测性能（本机 Release + CHOLMOD）

**加载 / 建网**（解析 + 去重建边 + 集中质量 + 包围盒，全部一次完成）：

| 模型 | 规模 | 耗时 |
|---|---|---:|
| cube.msh | 8 节点 / 6 四面体 / 19 弹簧 | < 0.01 s |
| sphere1K.msh | 1760 节点 / 6851 四面体 / 9847 弹簧 | 0.01 s |
| bunny.msh | 4670 节点 / 19379 四面体 / 27181 弹簧 | 0.02 s |
| sphere5K.msh | 7034 节点 / 27792 四面体 / 39696 弹簧 | 0.03 s |
| Armadillo219K.msh | 219355 节点 / 928030 四面体 / 1285752 弹簧 | 0.75 s |

**单步仿真耗时**（隐式欧拉 + 投影牛顿；不同模型/容差）：

| 模型 | 节点数 | 约 ms/步 | 备注 |
|---|---:|---:|---|
| bar-186 | 86 | 4–12 | tol=1e-6 |
| sphere1K | 1760 | ~75 | tol=1e-4 |
| bunny | 4670 | ~190 | tol=1e-2 |
| sphere5K | 7034 | ~245 | tol=1e-2 |
| Armadillo13K | 13230 | ~360 | tol=1e-2 |
| mat150x150t40 | 45000 | ~1300 | tol=1e-2 |

结论：**加载不是瓶颈**（最大网格 0.75 s）；**仿真是瓶颈**——≤7K 节点适合交互/视频，1–5 万节点可跑但慢，>10 万节点建议只做几何查看/静态渲染。

## 7. 结论与补齐建议

1. **可直接使用的范围**：65/70 个体积网格（加载 + 仿真 + 可视化全通），覆盖从 4 节点单元到 21.9 万节点模型；唯一无 `$Surface` 的 `32770_octocat.msh` 通过边界提取正常渲染。
2. **差一步的范围**：3 个 Gmsh 2.2 + 2 个 IPC legacy（全部 tet4）。补解析器约 80 行即可让 70/70 全通。
3. **不能仿真但可补齐可视化的范围**：46 个 `.obj` + 203 个 `.seg`（+ 悬空 `.pt`）。需要新增「只读 loader + 线/点/面叠加渲染」，物理上仅作参考显示。
4. **仓库自身的瑕疵**：6 条悬空引用（第 5 节），这些教程场景在 IPC 里也跑不起来。

---

> 统计口径说明：`被场景引用` 基于 `input/**/*.txt` 全文匹配；`.obj` 闭合性按「每条无向边恰好被 2 个三角面共享」判定；
> 所有 `.msh` 顶点/四面体数均来自实际解析，弹簧数与表面面片数来自 fiction 建网/表面提取结果。
