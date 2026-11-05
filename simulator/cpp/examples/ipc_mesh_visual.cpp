// =============================================================================
// examples/ipc_mesh_visual.cpp — IPC .msh 模型交互式可视化（GLFW + OpenGL）
// =============================================================================
// 职责：
//   - 与 ipc_mesh_demo 相同的场景装配与物理（mass-spring + 重力 + 区域外力）
//   - 用 $Surface 三角面实时渲染，pin/外力节点高亮
//   - 可选导出表面 OBJ 帧
//
// 操作：
//   左键拖拽 = 旋转视角   滚轮 = 缩放   右键(或 Shift+左键)拖拽 = 平移
//   Space = 暂停/继续     R = 复位相机   Up/Down = 每帧仿真步数 ±1   ESC = 退出
// =============================================================================

#include "io/obj_writer.h"
#include "scene/ipc_scene.h"
#include "time_integrator.h"
#include "types.h"
#include "viewer.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (fiction::ipcOptionsWantHelp(args)) {
        fiction::printIpcSceneUsage(argv[0], std::cout);
        return 0;
    }

    fiction::IpcSceneOptions opt;
    std::string err;
    if (fiction::parseIpcSceneOptions(args, opt, err) != fiction::IpcParseResult::Ok) {
        std::cerr << "[ipc_mesh_visual] 参数错误: " << err << "\n\n";
        fiction::printIpcSceneUsage(argv[0], std::cerr);
        return 1;
    }

    try {
        if (!opt.videoPath.empty()) {
            std::cerr << "[ipc_mesh_visual] 提示: 交互式版本不支持 --video，"
                         "视频导出请用 ipc_mesh_demo（CPU 软光栅化，SSH 可用）\n";
        }
        fiction::IpcScene scene = fiction::buildIpcScene(opt);
        std::cout << "===== fiction IPC 模型可视化（"
                  << (opt.material == "nh" ? "Neo-Hookean FEM" : "mass-spring") << "）=====\n"
                  << fiction::ipcSceneSummary(scene, opt) << "\n";

        fiction::TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = false; // 实时窗口不刷屏
        fiction::TimeIntegrator integrator(cfg, scene.model.masses);
        for (auto& e : fiction::makeIpcSceneEnergies(scene, opt)) {
            integrator.addEnergy(std::move(e));
        }

        // 相机按模型包围盒自动取景
        const fiction::Vec3 center = fiction::boxCenter(scene.model.bbox);
        const double diag = fiction::boxDiagonal(scene.model.bbox);
        const double radius = std::max(1.0, 1.8 * std::max(diag, 1e-6));

        fiction::Viewer::Config vcfg;
        vcfg.title = "fiction IPC mesh demo (Space=暂停 R=复位相机 Up/Down=步数 ESC=退出)";
        fiction::Viewer viewer(vcfg, center, radius);
        if (!viewer.ok()) {
            std::cerr << "[ipc_mesh_visual] Viewer 初始化失败（无显示环境？）\n";
            return 1;
        }

        // 导出目录准备
        const bool exporting = !opt.exportDir.empty();
        if (exporting) {
            std::error_code ec;
            std::filesystem::create_directories(opt.exportDir, ec);
            if (ec) {
                std::cerr << "[ipc_mesh_visual] 无法创建导出目录 " << opt.exportDir
                          << ": " << ec.message() << "\n";
                return 1;
            }
            if (!fiction::writeSurfaceObj(opt.exportDir + "/0.obj", scene.model.x0, scene.mesh.surface))
                throw std::runtime_error("OBJ 写出失败: " + opt.exportDir + "/0.obj");
        }

        // pin / 外力节点高亮掩码
        const std::vector<char> pinMask = scene.dbc.isDBC;
        std::vector<char> forceMask(scene.model.numNodes(), 0);
        for (const auto& nodes : scene.forceNodes) {
            for (int i : nodes) forceMask[static_cast<size_t>(i)] = 1;
        }

        fiction::Positions x = scene.model.x0;
        fiction::Velocities v = fiction::Velocities::Zero(scene.model.numNodes(), fiction::kDim);
        if (opt.initVel.squaredNorm() > 0.0) v.rowwise() = opt.initVel.transpose();

        const fiction::Vec3 groundOrigin(center.x(), center.y(), scene.model.bbox[2]);
        const double groundExtent = std::max(diag * 1.5, 1.0);
        int step = 0;

        while (!viewer.shouldClose()) {
            viewer.pollEvents();

            viewer.beginFrame();
            viewer.drawGround(groundOrigin, fiction::Vec3(0.0, 0.0, 1.0),
                              groundExtent, groundExtent / 10.0, 0.35f, 0.38f, 0.42f);
            viewer.drawTriangles(x, scene.mesh.surface, 0.45f, 0.62f, 0.90f);
            if (!scene.pinnedNodes.empty()) {
                viewer.drawPoints(x, pinMask, 7.0f, 1.0f, 0.25f, 0.2f);
            }
            if (!scene.forceNodes.empty()) {
                viewer.drawPoints(x, forceMask, 7.0f, 1.0f, 0.75f, 0.15f);
            }
            viewer.endFrame();

            if (!viewer.paused() && viewer.stepsPerFrame() > 0) {
                for (int k = 0; k < viewer.stepsPerFrame(); ++k) {
                    // 与控制台一致：本次求解的目标位置对应下一时刻，而非当前时刻。
                    fiction::refreshIpcDbcTargets(scene, opt, (step + 1) * opt.dt);
                    const auto state = integrator.stepForward(x, v, scene.dbc);
                    if (!state.report.converged || !state.x.allFinite() || !state.v.allFinite()) {
                        std::cerr << "[FAIL] step=" << (step + 1)
                                  << " stop_reason=" << state.report.stopReason
                                  << "; simulation stopped; completed_steps=" << step << "\n";
                        return 1;
                    }
                    x = state.x;
                    v = state.v;
                    ++step;
                    if (exporting && step % opt.exportEvery == 0) {
                        const std::string objPath = opt.exportDir + "/" + std::to_string(step) + ".obj";
                        if (!fiction::writeSurfaceObj(objPath, x, scene.mesh.surface))
                            throw std::runtime_error("OBJ 写出失败: " + objPath);
                    }
                }
            }
        }

        std::cout << std::fixed << std::setprecision(6)
                  << "退出：共仿真 " << step << " 步\n";
        if (exporting) {
            std::cout << "OBJ 导出 : " << opt.exportDir << "/（每 " << opt.exportEvery << " 步）\n";
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "[ipc_mesh_visual] 错误: " << ex.what() << "\n";
        return 1;
    }
}
