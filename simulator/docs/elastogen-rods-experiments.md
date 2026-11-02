# Elastogen cantilever 实验复现（fiction NH-FEM 版）

> 记录用 fiction 复现 Elastogen 论文两张插图对应实验的设计、数据来源、参数与验证方法。
> 资产清单见 [`ipc-asset-inventory.md`](ipc-asset-inventory.md)，连接原理见 [`ipc-mesh-integration.md`](ipc-mesh-integration.md)。

## 1. 论文出处与场景定义

| 图 | 文献位置 | 场景 | 材料 | 要素 |
|---|---|---|---|---|
| `rods.pdf` | 正文 Fig. `rods-bend`(b)（`sec/5_exp.tex` Cantilever 行） | 16×3×3 节点悬臂方杆、一端夹持、重力弯曲 | Co-rotational / NH / StVK × ν=0.21/0.32/0.43 | E 固定、`Δt=0.001`、432 DoF |
| `rods_twist.pdf` | 补充 Fig. `rods-twist`(a)（`sec/X_suppl.tex`） | 同杆 + **NH 扭转轨迹**（驱动端随时间绕轴旋转） | NH | "loops"为 RNN-2 求解器 Jacobi 迭代数（非物理量） |

论文网格来自 `Elastogen/neuralProj/assets/rod.bgeo`（Houdini `.bgeo`，16×3×3 节点 / 15×2×2 体素 / dx=1，见 `sim_main.py` demo=test）。该数据集**不在 IPC**（IPC 只有近似几何的方杆 `bar-*.msh`）。

## 2. 数据获取：bgeo → tet 转换

`tools/bgeo_to_msh.py`（需 `elastogen` conda 环境读取 `pyHouGeoIO` 二进制扩展）：

```
conda run -n elastogen python tools/bgeo_to_msh.py \
    --input /home/steven/mjx/Archive/Elastogen/neuralProj/assets/rod.bgeo   --output cpp/tests/data/rod16x3x3.msh
conda run -n elastogen python tools/bgeo_to_msh.py \
    --input /home/steven/mjx/Archive/Elastogen/neuralProj/assets/rod_z.bgeo --output cpp/tests/data/rod_z16x3x3.msh
```

转换要点（细节见脚本头注释）：
- 六面体 → **Kuhn(Freudenthal) 6-tet 剖分**，主对角线方向按 (i+j+k) 奇偶翻转 → 相邻共享面三角对齐（保形网格）；
- 产物：144 顶点 / 360 tet / 631 唯一边 / 256 边界面 / 体积 60（=15·2·2 ✓，自检含逐体体积与边界闭合）；
- `rod_z16x3x3.msh` 轴向沿 −z（PDNet 的 twist 场景用）。

## 3. fiction 侧新增能力（本次实现）

| 能量/机制 | 文件 | 说明 |
|---|---|---|
| **Neo-Hookean tet FEM** | `cpp/include/energies/NeoHookeanTetEnergy.h` | P1 常应变四面体；与 IPC 的 `NeoHookeanEnergy` 同公式：E = μ/2(F:F−3) − μ ln J + λ/2 ln²J；12×12 局部块经 `makePSD`；λ/μ 由 E/ν 换算（μ=E/2(1+ν), λ=Eν/((1+ν)(1−2ν))） |
| `--material {spring\|nh} --poisson ν` | `cpp/include/scene/ipc_scene.h` | 弹性模型二选一装配（弹簧沿用 auto k = E·l̄） |
| **DBC 绝对目标轨迹** | `types.h`(DBCInfo.targetPositions) + `time_integrator.cpp` | 存在时优先于 velocity/limit 机制；驱动语义：目标在步末时刻到达、满足判据后自由度消去 |
| `--twist-axis <x|y|z> --twist-rate <deg/s> --drive-box ...` | `ipc_scene.h`（`refreshIpcDbcTargets`） | 驱动端节点目标 = R(ω·t)·x_rest，绕过 bbox 中心、沿轴的直线（与 PDNet 软 pin 追踪刚体旋转轨迹同义，我们用硬 DBC 实现） |

质量沿用集中质量 m_i = ρ·V/4（规则格网下与 PDNet 的 ρdx³/节点 完全一致）。

## 4. 数值验证（`mesh_check` 检验 5.5）

- 等参位形 E=0、grad=0；Hessian 对称 + PSD（投影生效，min_eig ≥ 0）；
- **旋转不变性**：刚体旋转位形零能量；
- 三配置 FD（拉伸 1.15× / 剪切 0.3 / 体积压缩 0.8-0.9）：val/addGrad/addHess 与中心差分一致（相对 <1e-6 / <1e-4）；
- PSD 域判定：拉伸位形下「投影后 == 精确海森」，压缩/剪切下投影版保持 PSD（与弹簧检验第 4 项同一方法论）；
- ν 单调性：体积压缩位形 E(ν=0.45) > E(ν=0.05)。

## 5. 复现命令与产出

```
# bend：NH 重力弯曲（rods.pdf 场景，10s 实时视频）
./cpp/build/ipc_mesh_demo --mesh cpp/tests/data/rod16x3x3.msh \
    --material nh --youngs 5e7 --poisson 0.32 \
    --pin-box -1 -1 -1 2 3 3 --gravity 0 0 -9.81 \
    --dt 0.002 --steps 5000 --export-every 20 --tol 1e-3 \
    --video output/videos/elastogen_cantilever/bend_rods_nh.mp4

# twist：NH 扭转（rods_twist 场景；10s 内 360°）
./cpp/build/ipc_mesh_demo --mesh cpp/tests/data/rod16x3x3.msh \
    --material nh --youngs 5e7 --poisson 0.32 \
    --pin-box -1 -1 -1 2 3 3 --drive-box 13 -1 -1 16 3 3 \
    --twist-axis x --twist-rate 36 --no-gravity \
    --dt 0.002 --steps 5000 --export-every 20 --tol 1e-3 \
    --video output/videos/elastogen_cantilever/twist_rods_nh.mp4

# 对照：同场景 mass-spring（无泊松比，弯曲形态定性相似）
./cpp/build/ipc_mesh_demo --mesh cpp/tests/data/rod16x3x3.msh \
    --material spring --youngs 5e7 \
    --pin-box -1 -1 -1 2 3 3 --gravity 0 0 -9.81 \
    --dt 0.002 --steps 5000 --export-every 20 --tol 1e-3 \
    --video output/videos/elastogen_cantilever/bend_rods_spring.mp4
```

产出（`fiction/output/videos/elastogen_cantilever/`，SSH 可看；IPC .msh 演示系列在同级 `ipc_models/`）：
- `bend_rods_nh.mp4`（10.04s）、`twist_rods_nh.mp4`（10.04s）、`bend_rods_spring.mp4`（10.04s）
- `rods_experiments_sheet.png`（三条视频全时程关键帧拼图）

参数选择依据：E 由静力悬臂挠度 δ = wL⁴/(8EI) 标定——E=5e7 时稳态挠曲约 2–3 单位，
与原图弯曲幅度同量级；t / ω 为演示观感定档（40–60 deg/s 的扭转可看清 3D）。

## 6. 与论文的边界（诚实声明）

- 材料模型与 ν 项已按 NH 解析公式实现，但 ElastoGen 是**神经网络近似的能量+局部-全局求解**，
  本文目标是「同场景、同几何、真实物理的对照动画」，不做 NN 与 fiction 的定量误差对照；
- Co-rotational / StVK 尚未实现（接口已就位，见 `NeoHookeanTetEnergy` 模板 + `ipc_scene.h` 的 material 分支）；
- PDNet 的驱动端是**软 pin**（NN 松弛近似），此处用硬 DBC 轨迹（收敛后节点精确贴轨迹）。

## 7. 后续路线（Phase 2 完整版）

1. `CorotatedTetEnergy`（SVD R 投影 + 体积项 λ(Σσ−3)²）与 `StVKTetEnergy`（Green 应变 + λ(trE)²），
   对应 Elastogen `util.get_energy` 的 cr/stvk 分支 → 完整 3 材料 × 3 ν 与 rods.pdf 逐格对齐；
2. `--force-box` 时间窗（对齐 IPC NBC 的 t_begin/t_end）；
3. 材料对比录：ν 与体积自锁（locking）对 StVK/CR 弯曲行为的影响（论文 Fig. rods-bend 的横轴内容）。
