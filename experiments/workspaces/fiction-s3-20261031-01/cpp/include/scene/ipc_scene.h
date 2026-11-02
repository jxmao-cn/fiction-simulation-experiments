#pragma once
// =============================================================================
// scene/ipc_scene.h — IPC .msh + 命令行参数的场景装配
// =============================================================================
// 职责：
//   - 解析 ipc_mesh_demo / ipc_mesh_visual 共用的命令行参数
//   - 加载 .msh → MassSpringModel → DBC（pin-box 零速固定）+ 区域外力（force-box）
//   - 汇总统计信息文本（模型规模 / 体积 / 质量 / 刚度 / 包围盒）
//
// 参数示例：
//   ipc_mesh_visual --mesh bar-186.msh
//       --pin-box 0 0 0 0.6 0.6 1.1
//       --force-box 9.4 0 0 10 0.5 1 30 0 0
//       --gravity 0 0 -9.81 --dt 0.01 --steps 300
// =============================================================================

#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "energies/RegionForce.h"
#include "energies/SpringEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"
#include "types.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace fiction {

// -----------------------------------------------------------------------------
// 命令行参数
// -----------------------------------------------------------------------------
struct IpcSceneOptions {
    std::string meshPath;                                  ///< 必填：.msh 路径
    std::string material = "spring";                       ///< 弹性模型：spring（边弹簧）| nh（NH tet FEM）
    double poisson = 0.3;                                  ///< --material nh 的泊松比
    double density = 1000.0;                               ///< 密度 kg/m^3
    double stiffness = 0.0;                                ///< 每边刚度；<=0 时 auto = youngs * 平均边长
    double youngs = 1e5;                                   ///< spring: auto 刚度用；nh: 直接作为 E
    double dt = 0.01;                                      ///< 时间步长
    int steps = 200;                                       ///< 仿真步数
    Vec3 gravity = Vec3(0.0, 0.0, -9.81);                  ///< 重力加速度
    Vec3 initVel = Vec3(0.0, 0.0, 0.0);                    ///< 全体节点初始速度
    bool enableGravity = true;                             ///< --no-gravity 关闭
    std::vector<std::array<double, 6>> pinBoxes;           ///< 零速固定区域（可重复）
    std::vector<std::array<double, 6>> driveBoxes;         ///< 扭转驱动区域（可重复，与 --twist-axis 配合）
    std::string twistAxis = "";                            ///< 扭转旋转轴：""(关) | x | y | z
    double twistRateDeg = 0.0;                             ///< 角速度（deg/s，仿真时间）
    struct ForceBox {
        std::array<double, 6> box;                         ///< [xmin,ymin,zmin,xmax,ymax,zmax]
        Vec3 accel{0.0, 0.0, 0.0};                         ///< 恒定加速度
    };
    std::vector<ForceBox> forceBoxes;                      ///< 区域外力（可重复）
    std::string exportDir;                                 ///< 非空时导出 OBJ 帧
    int exportEvery = 10;                                  ///< 每隔多少步导出一次（OBJ 与视频共用）
    std::string videoPath;                                 ///< 非空时导出 mp4（仅 ipc_mesh_demo）
    double videoFps = 25.0;                                ///< 视频帧率
    int videoWidth = 960;                                  ///< 视频宽（自动取偶）
    int videoHeight = 720;                                 ///< 视频高（自动取偶）
    bool videoOrbit = true;                                ///< 视频相机绕模型环绕（--video-static 关闭）
    bool videoKeepFrames = false;                          ///< 编码成功后保留帧目录
    double tol = 1e-2;                                     ///< Newton 容差
    int maxIter = 100;                                     ///< Newton 最大迭代数
    bool verbose = false;                                  ///< 是否打印 Newton 迭代日志
};

enum class IpcParseResult { Ok, Help, Error };

namespace detail {

inline bool parseDouble(const std::string& tok, double& out) {
    try {
        size_t pos = 0;
        out = std::stod(tok, &pos);
        return pos == tok.size();
    } catch (const std::exception&) {
        return false;
    }
}

inline bool parseInt(const std::string& tok, int& out) {
    try {
        size_t pos = 0;
        const long v = std::stol(tok, &pos);
        if (pos != tok.size() || v < 0) return false;
        out = static_cast<int>(v);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

inline bool takeDouble(const std::vector<std::string>& args, size_t& i,
                       const std::string& flag, double& out, std::string& err) {
    if (i + 1 >= args.size()) {
        err = flag + " 缺少参数值";
        return false;
    }
    if (!parseDouble(args[++i], out)) {
        err = flag + " 参数不是合法浮点数: " + args[i];
        return false;
    }
    return true;
}

inline bool takeInt(const std::vector<std::string>& args, size_t& i,
                    const std::string& flag, int& out, std::string& err) {
    if (i + 1 >= args.size()) {
        err = flag + " 缺少参数值";
        return false;
    }
    if (!parseInt(args[++i], out)) {
        err = flag + " 参数不是合法非负整数: " + args[i];
        return false;
    }
    return true;
}

inline bool takeVec3(const std::vector<std::string>& args, size_t& i,
                     const std::string& flag, Vec3& out, std::string& err) {
    if (i + 3 >= args.size()) {
        err = flag + " 需要 3 个浮点数参数";
        return false;
    }
    double v[3];
    for (int d = 0; d < 3; ++d) {
        if (!parseDouble(args[++i], v[d])) {
            err = flag + " 参数不是合法浮点数: " + args[i];
            return false;
        }
    }
    out = Vec3(v[0], v[1], v[2]);
    return true;
}

inline bool takeBox(const std::vector<std::string>& args, size_t& i,
                    const std::string& flag, std::array<double, 6>& out, std::string& err) {
    if (i + 6 >= args.size()) {
        err = flag + " 需要 6 个浮点数参数 (xmin ymin zmin xmax ymax zmax)";
        return false;
    }
    for (int d = 0; d < 6; ++d) {
        if (!parseDouble(args[++i], out[static_cast<size_t>(d)])) {
            err = flag + " 参数不是合法浮点数: " + args[i];
            return false;
        }
    }
    return true;
}

} // namespace detail

inline bool ipcOptionsWantHelp(const std::vector<std::string>& args) {
    for (const auto& a : args)
        if (a == "--help" || a == "-h") return true;
    return false;
}

inline IpcParseResult parseIpcSceneOptions(const std::vector<std::string>& args,
                                           IpcSceneOptions& opt, std::string& err) {
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "--mesh") {
            if (i + 1 >= args.size()) { err = "--mesh 缺少参数值"; return IpcParseResult::Error; }
            opt.meshPath = args[++i];
        } else if (a == "--material") {
            if (i + 1 >= args.size()) { err = "--material 缺少参数值 (spring|nh)"; return IpcParseResult::Error; }
            opt.material = args[++i];
        } else if (a == "--poisson") {
            if (!detail::takeDouble(args, i, a, opt.poisson, err)) return IpcParseResult::Error;
        } else if (a == "--density") {
            if (!detail::takeDouble(args, i, a, opt.density, err)) return IpcParseResult::Error;
        } else if (a == "--stiffness") {
            if (!detail::takeDouble(args, i, a, opt.stiffness, err)) return IpcParseResult::Error;
        } else if (a == "--youngs") {
            if (!detail::takeDouble(args, i, a, opt.youngs, err)) return IpcParseResult::Error;
        } else if (a == "--dt") {
            if (!detail::takeDouble(args, i, a, opt.dt, err)) return IpcParseResult::Error;
        } else if (a == "--steps") {
            if (!detail::takeInt(args, i, a, opt.steps, err)) return IpcParseResult::Error;
        } else if (a == "--gravity") {
            if (!detail::takeVec3(args, i, a, opt.gravity, err)) return IpcParseResult::Error;
        } else if (a == "--init-vel") {
            if (!detail::takeVec3(args, i, a, opt.initVel, err)) return IpcParseResult::Error;
        } else if (a == "--no-gravity") {
            opt.enableGravity = false;
        } else if (a == "--drive-box") {
            std::array<double, 6> box{};
            if (!detail::takeBox(args, i, a, box, err)) return IpcParseResult::Error;
            opt.driveBoxes.push_back(box);
        } else if (a == "--twist-axis") {
            if (i + 1 >= args.size()) { err = "--twist-axis 缺少参数值 (x|y|z)"; return IpcParseResult::Error; }
            opt.twistAxis = args[++i];
        } else if (a == "--twist-rate") {
            if (!detail::takeDouble(args, i, a, opt.twistRateDeg, err)) return IpcParseResult::Error;
        } else if (a == "--pin-box") {
            std::array<double, 6> box{};
            if (!detail::takeBox(args, i, a, box, err)) return IpcParseResult::Error;
            opt.pinBoxes.push_back(box);
        } else if (a == "--force-box") {
            IpcSceneOptions::ForceBox fb;
            if (!detail::takeBox(args, i, a, fb.box, err)) return IpcParseResult::Error;
            if (!detail::takeVec3(args, i, a, fb.accel, err)) return IpcParseResult::Error;
            opt.forceBoxes.push_back(fb);
        } else if (a == "--export") {
            if (i + 1 >= args.size()) { err = "--export 缺少目录参数"; return IpcParseResult::Error; }
            opt.exportDir = args[++i];
        } else if (a == "--export-every") {
            if (!detail::takeInt(args, i, a, opt.exportEvery, err)) return IpcParseResult::Error;
        } else if (a == "--video") {
            if (i + 1 >= args.size()) { err = "--video 缺少输出文件路径"; return IpcParseResult::Error; }
            opt.videoPath = args[++i];
        } else if (a == "--video-fps") {
            if (!detail::takeDouble(args, i, a, opt.videoFps, err)) return IpcParseResult::Error;
        } else if (a == "--video-size") {
            if (!detail::takeInt(args, i, a, opt.videoWidth, err)) return IpcParseResult::Error;
            if (!detail::takeInt(args, i, a, opt.videoHeight, err)) return IpcParseResult::Error;
        } else if (a == "--video-static") {
            opt.videoOrbit = false;
        } else if (a == "--video-keep-frames") {
            opt.videoKeepFrames = true;
        } else if (a == "--tol") {
            if (!detail::takeDouble(args, i, a, opt.tol, err)) return IpcParseResult::Error;
        } else if (a == "--max-iter") {
            if (!detail::takeInt(args, i, a, opt.maxIter, err)) return IpcParseResult::Error;
        } else if (a == "--verbose") {
            opt.verbose = true;
        } else if (a == "--help" || a == "-h") {
            return IpcParseResult::Help;
        } else {
            err = "未知参数: " + a;
            return IpcParseResult::Error;
        }
    }
    if (opt.meshPath.empty()) {
        err = "缺少必填参数 --mesh <path.msh>";
        return IpcParseResult::Error;
    }
    if (opt.material != "spring" && opt.material != "nh") {
        err = "--material 仅支持 spring | nh（当前: " + opt.material + "）";
        return IpcParseResult::Error;
    }
    if (opt.poisson <= -1.0 || opt.poisson >= 0.5) {
        err = "--poisson 需在 (-1, 0.5) 内";
        return IpcParseResult::Error;
    }
    if (opt.twistAxis != "" && opt.twistAxis != "x" && opt.twistAxis != "y" && opt.twistAxis != "z") {
        err = "--twist-axis 仅支持 x | y | z";
        return IpcParseResult::Error;
    }
    const bool twisting = !opt.twistAxis.empty() && opt.twistRateDeg != 0.0;
    if ((opt.twistAxis != "" || !opt.driveBoxes.empty()) && !twisting) {
        err = "启用了 --drive-box/--twist-axis 但 twist-rate 为 0 或轴为空";
        return IpcParseResult::Error;
    }
    if (twisting && opt.driveBoxes.empty()) {
        err = "--twist 驱动需要至少一个 --drive-box 指定驱动端";
        return IpcParseResult::Error;
    }
    if (opt.density <= 0.0) { err = "--density 必须为正"; return IpcParseResult::Error; }
    if (opt.dt <= 0.0) { err = "--dt 必须为正"; return IpcParseResult::Error; }
    if (opt.exportEvery <= 0) { err = "--export-every 必须为正"; return IpcParseResult::Error; }
    if (opt.videoFps <= 0.0) { err = "--video-fps 必须为正"; return IpcParseResult::Error; }
    if (opt.videoWidth <= 0 || opt.videoHeight <= 0) {
        err = "--video-size 宽高必须为正";
        return IpcParseResult::Error;
    }
    return IpcParseResult::Ok;
}

inline void printIpcSceneUsage(const char* prog, std::ostream& os) {
    os << "用法: " << prog << " --mesh <path.msh> [选项]\n"
       << "  必填:\n"
       << "    --mesh <path>            IPC .msh 四面体网格（MshIO 4.1 ASCII）\n"
       << "  物理参数:\n"
       << "    --material <spring|nh>   弹性模型：边弹簧（默认）或 Neo-Hookean tet FEM\n"
       << "    --poisson <nu>           --material nh 的泊松比（默认 0.3，(−1,0.5)）\n"
       << "    --density <rho>          密度 kg/m^3（默认 1000，IPC 默认值）\n"
       << "    --stiffness <k>          每边刚度（仅 spring）；缺省时 auto = --youngs * 平均边长\n"
       << "    --youngs <E>             spring: auto 刚度的 E（默认 1e5 Pa）；nh: 直接作为 Young's modulus\n"
       << "    --dt <h>                 时间步长（默认 0.01）\n"
       << "    --steps <n>              仿真步数（默认 200）\n"
       << "    --gravity <gx gy gz>     重力（默认 0 0 -9.81）\n"
       << "    --no-gravity             关闭重力\n"
       << "    --init-vel <vx vy vz>    全体节点初始速度（默认 0，可做纯弹簧振荡演示）\n"
       << "    --tol <tol>              Newton 容差（默认 1e-2）\n"
       << "    --max-iter <n>           Newton 最大迭代数（默认 100）\n"
       << "  边界/外力（可重复）:\n"
       << "    --pin-box xmin ymin zmin xmax ymax zmax     盒内节点零速固定（DBC）\n"
       << "    --drive-box xmin ymin zmin xmax ymax zmax   盒内节点跟随扭转轨迹（配 --twist-axis/--twist-rate）\n"
       << "    --twist-axis <x|y|z>     扭转旋转轴（默认关）\n"
       << "    --twist-rate <deg/s>     角速度（deg/s，仿真时间；应与 --drive-box 同用）\n"
       << "    --force-box xmin ymin zmin xmax ymax zmax ax ay az   盒内节点恒定加速度（NBC）\n"
       << "  输出:\n"
       << "    --export <dir>           导出表面三角面 OBJ 帧\n"
       << "    --export-every <n>       每隔 n 步导出一次（默认 10，OBJ 与视频帧共用）\n"
       << "    --video <path.mp4>       导出 mp4 视频（CPU 软光栅化 + ffmpeg，SSH 无显示环境可用，仅控制台 demo）\n"
       << "    --video-fps <fps>        视频帧率（默认 25）\n"
       << "    --video-size <W> <H>     视频分辨率（默认 960 720，自动取偶）\n"
       << "    --video-static           关闭相机环绕（默认缓慢环绕展示 3D 形变）\n"
       << "    --video-keep-frames      编码成功后保留 PNG 帧目录\n"
       << "    --verbose                打印 Newton 迭代日志\n";
}

// -----------------------------------------------------------------------------
// 场景装配
// -----------------------------------------------------------------------------
struct IpcScene {
    MshMesh mesh;
    MassSpringModel model;
    DBCInfo dbc;
    std::vector<int> pinnedNodes;                    ///< pin-box 选中的节点（去重）
    std::vector<std::vector<int>> forceNodes;        ///< 每个 force-box 选中的节点
    std::vector<int> driveNodes;                     ///< 扭转驱动的节点（不与 pin 重叠）
    Eigen::Vector3d axisPoint = Eigen::Vector3d::Zero();  ///< 旋转轴上的点（bbox 中心）
    Eigen::Vector3d axisDir = Eigen::Vector3d::UnitX();  ///< 旋转轴单位方向
};

/// 加载网格并装配 mass-spring + DBC（不做仿真）
inline IpcScene buildIpcScene(const IpcSceneOptions& opt) {
    IpcScene scene;
    scene.mesh = loadMsh(opt.meshPath);
    scene.model = buildMassSpringModel(scene.mesh, opt.density, opt.stiffness, opt.youngs);
    const int N = scene.model.numNodes();

    // pin-box → 零速 DBC（去重，避免重复罚项）
    std::unordered_set<int> pinnedSet;
    for (const auto& box : opt.pinBoxes) {
        for (int i : nodesInBox(scene.model.x0, box)) pinnedSet.insert(i);
    }
    scene.pinnedNodes.assign(pinnedSet.begin(), pinnedSet.end());
    std::sort(scene.pinnedNodes.begin(), scene.pinnedNodes.end());

    scene.dbc.indices = scene.pinnedNodes;
    scene.dbc.isDBC.assign(N, 0);
    for (int i : scene.pinnedNodes) {
        scene.dbc.targetVelocities.push_back(Vec3::Zero());
        scene.dbc.limits.push_back(scene.model.x0.row(i).transpose());
        scene.dbc.isDBC[i] = 1;
    }

    // 扭转驱动节点：drive-box 选取，剔除已被 pin 的（pin 优先）
    std::unordered_set<int> driveSet;
    for (const auto& box : opt.driveBoxes) {
        for (int i : nodesInBox(scene.model.x0, box)) {
            if (!pinnedSet.count(i)) driveSet.insert(i);
        }
    }
    scene.driveNodes.assign(driveSet.begin(), driveSet.end());
    std::sort(scene.driveNodes.begin(), scene.driveNodes.end());

    // 合并 DBC indices（pin ∪ drive），并分配绝对目标轨迹槽位
    for (int i : scene.driveNodes) {
        scene.dbc.indices.push_back(i);
        scene.dbc.targetVelocities.push_back(Vec3::Zero());
        scene.dbc.limits.push_back(scene.model.x0.row(i).transpose());
        scene.dbc.isDBC[i] = 1;
    }
    if (!scene.dbc.indices.empty()) {
        scene.dbc.targetPositions.reserve(scene.dbc.indices.size());
        for (int idx : scene.dbc.indices) {
            scene.dbc.targetPositions.push_back(scene.model.x0.row(idx).transpose());
        }
    }

    // 旋转轴：沿 twistAxis、过整模型 bbox 中心（直杆的中心轴线，与 PDNet 的绕 rest 质心一致）
    if (opt.twistAxis == "x") {
        scene.axisDir = Eigen::Vector3d::UnitX();
    } else if (opt.twistAxis == "y") {
        scene.axisDir = Eigen::Vector3d::UnitY();
    } else if (opt.twistAxis == "z") {
        scene.axisDir = Eigen::Vector3d::UnitZ();
    }
    scene.axisPoint = boxCenter(scene.model.bbox);

    // force-box → 区域节点集合
    for (const auto& fb : opt.forceBoxes)
        scene.forceNodes.push_back(nodesInBox(scene.model.x0, fb.box));

    return scene;
}

/// 按仿真时间刷新扭转 DBC 的绝对目标位置（pin → x0；drive → R(θ)·x0）
/// 在每个时间步的 stepForward 之前调用，t_end = 本步结束的仿真时刻
inline void refreshIpcDbcTargets(IpcScene& scene, const IpcSceneOptions& opt,
                                 double t_end) {
    if (scene.driveNodes.empty()) return;
    const double theta = opt.twistRateDeg * (3.14159265358979323846 / 180.0) * t_end;
    const Eigen::Matrix3d R =
        Eigen::AngleAxisd(theta, scene.axisDir).toRotationMatrix();
    std::unordered_set<int> driveLookup(scene.driveNodes.begin(), scene.driveNodes.end());
    for (size_t k = 0; k < scene.dbc.indices.size(); ++k) {
        const int idx = scene.dbc.indices[k];
        Eigen::Vector3d p = scene.model.x0.row(idx).transpose();
        if (driveLookup.count(idx)) {
            p = scene.axisPoint + R * (p - scene.axisPoint);
        }
        scene.dbc.targetPositions[k] = p;
    }
}

/// 组装能量列表（重力 + 弹性（spring | nh）+ 区域外力）
inline std::vector<EnergyPtr> makeIpcSceneEnergies(const IpcScene& scene,
                                                   const IpcSceneOptions& opt) {
    std::vector<EnergyPtr> energies;
    if (opt.enableGravity) {
        energies.push_back(std::make_shared<ExampleGravityEnergy>(scene.model.masses, opt.gravity));
    }
    if (opt.material == "nh") {
        // Neo-Hookean tet FEM（复现 Elastogen cantilever 实验，见 docs/elastogen-rods-experiments.md）
        energies.push_back(std::make_shared<NeoHookeanTetEnergy>(
            scene.mesh.tets, scene.model.x0, opt.youngs, opt.poisson));
    } else {
        energies.push_back(std::make_shared<ExampleSpringEnergy>(
            scene.model.edges, scene.model.l2, scene.model.k));
    }
    for (size_t i = 0; i < scene.forceNodes.size(); ++i) {
        energies.push_back(makeRegionForceEnergy(scene.model.masses,
                                                 scene.forceNodes[i],
                                                 opt.forceBoxes[i].accel));
    }
    return energies;
}

/// 汇总信息文本（demo 与控制台共用）
inline std::string ipcSceneSummary(const IpcScene& scene, const IpcSceneOptions& opt) {
    const MassSpringModel& m = scene.model;
    std::ostringstream os;
    os.precision(6);
    os << std::fixed;
    os << "网格文件 : " << opt.meshPath << "\n"
       << "规模     : " << m.numNodes() << " 节点, " << scene.mesh.numTets()
       << " 四面体, " << m.numEdges() << " 弹簧(唯一边), "
       << scene.mesh.numSurface() << " 表面三角面\n"
       << "体积/质量: V=" << m.volume << " m^3, 总质量=" << m.totalMass << " kg\n"
       << "包围盒   : [" << m.bbox[0] << ", " << m.bbox[1] << ", " << m.bbox[2] << "] - ["
       << m.bbox[3] << ", " << m.bbox[4] << ", " << m.bbox[5] << "]\n"
       << "弹性材料 : " << opt.material
       << (opt.material == "nh" ? " (Neo-Hookean tet FEM, nu=" + std::to_string(opt.poisson)
                                + ", E=" + std::to_string(opt.youngs) + ")"
                                : " (mass-spring)")
       << "\n"
       << "平均边长 : " << m.avgEdgeLength
       << ", 每边刚度 k=" << (m.numEdges() > 0 ? m.k[0] : 0.0)
       << (opt.stiffness > 0.0 ? " (显式)" : " (auto = youngs * l_avg)") << "\n"
       << "时间积分 : h=" << opt.dt << ", steps=" << opt.steps
       << ", tol=" << opt.tol << ", gravity="
       << (opt.enableGravity
               ? "(" + std::to_string(opt.gravity.x()) + ", " + std::to_string(opt.gravity.y())
                     + ", " + std::to_string(opt.gravity.z()) + ")"
               : std::string("off"))
       << "\n"
       << "固定节点 : " << scene.pinnedNodes.size() << " (pin-box)";
    if (!opt.pinBoxes.empty()) {
        os << " [";
        for (size_t i = 0; i < opt.pinBoxes.size(); ++i) {
            if (i) os << "; ";
            const auto& b = opt.pinBoxes[i];
            os << "(" << b[0] << "," << b[1] << "," << b[2] << ")-("
               << b[3] << "," << b[4] << "," << b[5] << ")";
        }
        os << "]";
    }
    os << "\n";
    for (size_t i = 0; i < scene.forceNodes.size(); ++i) {
        const Vec3& a = opt.forceBoxes[i].accel;
        os << "区域外力 : box" << i << " -> " << scene.forceNodes[i].size() << " 节点, accel=("
           << a.x() << ", " << a.y() << ", " << a.z() << ")\n";
    }
    if (!scene.driveNodes.empty()) {
        os << "扭转驱动 : " << scene.driveNodes.size() << " 节点, 轴=" << opt.twistAxis
           << ", rate=" << opt.twistRateDeg << " deg/s\n";
    }
    return os.str();
}

} // namespace fiction
